"""Custom (user-uploaded) reference database adapter.

A user can upload their own reference FASTA in the Databases stage and compare
the run's OTU/ASV representative sequences against it. Two kinds of upload are
supported and auto-detected from the FASTA headers:

* ``"sintax"`` — every header already carries a ``;tax=`` taxonomy string, so
  the file is fed to VSEARCH ``--sintax`` verbatim (same path as SILVA/PR2) and
  folded into the taxonomy consensus.
* ``"blast"`` — a plain reference FASTA (e.g. curated COI/18S sequences). We
  build a local BLAST nucleotide database with ``makeblastdb`` and align the OTU
  representatives against it with ``blastn``, keeping the single best hit per
  OTU. Every query gets a row: those with no hit are reported with status
  ``"None"`` so the UI can flag them in red.

Both paths are pure-Python around the bundled command-line tools, so the module
stays unit-testable (the BLAST path is exercised against real ``blastn`` when
available and skipped cleanly otherwise).
"""

from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from core.process import run_tool
from databases.ncbi_db import BLAST_FAILURE_SIGNATURES
from core.config import BLASTN_PATH, MAKEBLASTDB_PATH, MAKEBLASTDB_TIMEOUT_SEC, BLAST_TIMEOUT_SEC
from core.models import StepOutput


# ── FASTA inspection ─────────────────────────────────────────────────────

def _open_maybe_gzip(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _sequence_id(header: str) -> str:
    """First whitespace/`;`-delimited token of a FASTA header (the accession)."""
    token = header.split()[0] if header else ""
    return token.split(";")[0].strip()


def reference_label(header: str) -> str:
    """Human-readable taxon label for a reference header.

    Handles the three header conventions we see in practice:
    * SINTAX (``ID;tax=d:..,g:Genus,s:Genus species;``) → the species/genus field
    * NCBI/BOLD style (``ACC Genus species voucher ...``) → ``Genus species``
    * bare id → the id itself
    """
    header = header.strip()
    if not header:
        return ""
    if ";tax=" in header:
        tax = header.split(";tax=", 1)[1].rstrip(";")
        ranks = {}
        for field in tax.split(","):
            field = field.strip()
            if ":" in field:
                key, _, value = field.partition(":")
                ranks[key.strip().lower()] = value.strip()
        for key in ("s", "species"):
            if ranks.get(key):
                return ranks[key].replace("_", " ")
        for key in ("g", "genus"):
            if ranks.get(key):
                return ranks[key].replace("_", " ")
        # fall through to id when tax carries no species/genus
        return _sequence_id(header)
    parts = header.split(None, 1)
    if len(parts) == 2 and parts[1].strip():
        # Keep the first two words (Genus species) when they look like a name.
        words = parts[1].strip().split()
        if len(words) >= 2 and words[0][:1].isalpha():
            return f"{words[0]} {words[1]}"
        return parts[1].strip().split(";")[0].strip() or parts[0]
    return _sequence_id(header)


def reference_labels(reference_fasta: Path) -> Dict[str, str]:
    """Map reference sequence id → taxon label for annotating BLAST hits."""
    labels: Dict[str, str] = {}
    with _open_maybe_gzip(reference_fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                labels[_sequence_id(header)] = reference_label(header)
    return labels


# IUPAC nucleotide codes, plus the gap characters an aligned reference may carry.
# Anything else is not DNA: F, I, L, P, E, Q and friends are protein-only, and their
# presence is what tells a peptide FASTA from a nucleotide one.
_NUCLEOTIDE_ALPHABET = set("ACGTUNRYSWKMBDHV-.acgtunryswkmbdhv")
_UNAMBIGUOUS = set("ACGTUacgtu")

# A reference sequence shorter than this cannot be a useful blastn subject: megablast
# seeds on 28-base words and even blastn's default word size is 11, so a 4 bp "sequence"
# produces a database that silently matches nothing.
MIN_REFERENCE_LENGTH = 50

# A peptide can slip past the alphabet check when every one of its letters happens to be
# an IUPAC nucleotide code — K, W, V, S, Y, R, G, A, M, D, H, N all are. What it cannot do
# is look like DNA by composition. Half is deliberately loose: a reference with masked or
# padded regions is still DNA, and rejecting it would be the validator overreaching. A
# peptide made only of IUPAC letters scores near zero, so the two are not close.
_MIN_UNAMBIGUOUS_FRACTION = 0.5


def _sequence_error(seq_id: str, index: int, sequence: str) -> str:
    """Why this record's SEQUENCE is unusable, or ``""``.

    The upload used to look only at lines starting with ``>``: a file with valid headers
    and protein, empty or junk sequence lines passed validation and failed a whole stage
    later, inside ``makeblastdb``, where the user could no longer act on it.
    """
    if not sequence:
        return f"Sequence #{index} ('{seq_id}') has no bases — the header is there but the body is empty."
    invalid = {char for char in sequence if char not in _NUCLEOTIDE_ALPHABET}
    if invalid:
        sample = "".join(sorted(invalid)[:6])
        return (
            f"Sequence #{index} ('{seq_id}') contains non-nucleotide characters "
            f"('{sample}'). A custom database must be DNA — a protein FASTA "
            f"cannot be compared with blastn."
        )
    # Gaps are alignment padding, not bases: they belong in neither side of this ratio,
    # and counting them as "not A/C/G/T" would reject an aligned reference for being aligned.
    bases = [char for char in sequence if char not in "-."]
    if not bases:
        return f"Sequence #{index} ('{seq_id}') is all gaps, without a single base."
    unambiguous = sum(1 for char in bases if char in _UNAMBIGUOUS)
    if unambiguous < _MIN_UNAMBIGUOUS_FRACTION * len(bases):
        return (
            f"Sequence #{index} ('{seq_id}') is only {100 * unambiguous / len(bases):.0f}% "
            f"A/C/G/T. This does not look like DNA."
        )
    if len(bases) < MIN_REFERENCE_LENGTH:
        return (
            f"Sequence #{index} ('{seq_id}') has {len(bases)} bases; the minimum is "
            f"{MIN_REFERENCE_LENGTH}. Shorter sequences do not produce reliable alignments."
        )
    return ""


def detect_fasta_kind(path: Path) -> tuple[str, str]:
    """Return ``(kind, error)`` where kind is ``"sintax"``, ``"blast"`` or ``""``.

    ``"sintax"`` when every record header carries ``;tax=``; ``"blast"`` for any
    other well-formed FASTA. ``error`` is a user-facing message when the file is
    not a usable FASTA at all — including when its SEQUENCES are unusable, which this
    used to ignore entirely.
    """
    record_count = 0
    missing_tax = 0
    seen_ids: set[str] = set()
    duplicate = ""
    first = True
    current_id = ""
    current_seq: List[str] = []

    def _flush() -> str:
        if not record_count:
            return ""
        return _sequence_error(current_id, record_count, "".join(current_seq).strip())

    try:
        with _open_maybe_gzip(path) as handle:
            for line in handle:
                if first:
                    if not line.startswith(">"):
                        return "", "This does not look like a FASTA file (the first line does not start with '>')."
                    first = False
                if not line.startswith(">"):
                    current_seq.append(line.strip())
                    continue
                error = _flush()
                if error:
                    return "", error
                record_count += 1
                header = line[1:].strip()
                seq_id = _sequence_id(header)
                if not seq_id:
                    return "", f"Header with no identifier at sequence #{record_count}."
                if seq_id in seen_ids and not duplicate:
                    duplicate = seq_id
                seen_ids.add(seq_id)
                if ";tax=" not in header:
                    missing_tax += 1
                current_id, current_seq = seq_id, []
            error = _flush()
            if error:
                return "", error
    except (OSError, UnicodeDecodeError) as exc:
        return "", f"Could not read the file: {exc}"

    if first or record_count == 0:
        return "", "Empty FASTA file."
    if missing_tax == 0:
        # SINTAX reference: duplicated ids break VSEARCH, so reject them.
        if duplicate:
            return "", f"Duplicate ID in the custom database: '{duplicate}'."
        return "sintax", ""
    # Plain reference for BLAST. Duplicate ids are tolerable here (makeblastdb
    # renames them), so only warn via the returned kind.
    return "blast", ""


# ── Config (kept for backward compatibility with the upload callback) ────

class CustomDBConfig:
    def __init__(self, db_path: Path, sintax_cutoff: float = 0.8):
        self.db_path = Path(db_path)
        self.sintax_cutoff = sintax_cutoff
        self.kind: str = ""

    def validate(self) -> tuple[bool, str]:
        if not self.db_path.exists():
            return False, f"Custom database not found: {self.db_path}"
        suffix = self.db_path.suffix.lower()
        if suffix not in {".fasta", ".fa", ".fna", ".gz"}:
            return False, f"Unsupported format: {suffix}. Use .fasta, .fa, .fna (or .gz of any of them)"
        kind, error = detect_fasta_kind(self.db_path)
        if not kind:
            return False, error
        self.kind = kind
        return True, ""


# ── BLAST comparison engine ──────────────────────────────────────────────

_BESTHIT_COLUMNS = [
    "OTU_ID",
    "subject_id",
    "species",
    "pident",
    "qcovhsp",
    "length",
    "evalue",
    "bitscore",
    "status",
]


def blast_tools_available() -> bool:
    return bool(
        (shutil.which(BLASTN_PATH) or Path(BLASTN_PATH).is_file())
        and (shutil.which(MAKEBLASTDB_PATH) or Path(MAKEBLASTDB_PATH).is_file())
    )


def _query_ids(query_fasta: Path) -> List[str]:
    ids: List[str] = []
    with _open_maybe_gzip(query_fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(_sequence_id(line[1:].strip()))
    return ids


def build_blast_db(reference_fasta: Path, work_dir: Path) -> Path:
    """Build a nucleotide BLAST database, returning its path prefix."""
    work_dir.mkdir(parents=True, exist_ok=True)
    # makeblastdb cannot read gzip directly; decompress into the work dir first.
    if reference_fasta.suffix.lower() == ".gz":
        plain = work_dir / (reference_fasta.stem or "reference.fasta")
        with gzip.open(reference_fasta, "rt", encoding="utf-8", errors="replace") as src, plain.open("w", encoding="utf-8") as dst:
            shutil.copyfileobj(src, dst)
        source = plain
    else:
        source = reference_fasta
    db_prefix = work_dir / "custom_ref_db"
    subprocess.run(
        [MAKEBLASTDB_PATH, "-in", str(source), "-dbtype", "nucl", "-out", str(db_prefix)],
        capture_output=True,
        text=True,
        check=True,
        timeout=MAKEBLASTDB_TIMEOUT_SEC,
    )
    return db_prefix


def _parse_best_hits(blast_tsv: Path) -> Dict[str, Dict[str, str]]:
    """Best hit (max bitscore) per query from a BLAST outfmt-6 file."""
    best: Dict[str, Dict[str, str]] = {}
    if not blast_tsv.exists() or blast_tsv.stat().st_size == 0:
        return best
    with blast_tsv.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            qid = cols[0].split(";")[0].strip()
            try:
                bitscore = float(cols[11])
            except ValueError:
                bitscore = 0.0
            current = best.get(qid)
            if current is None or bitscore > current["_bitscore"]:
                best[qid] = {
                    "subject_id": cols[1].strip(),
                    "pident": cols[2].strip(),
                    "length": cols[3].strip(),
                    "evalue": cols[10].strip(),
                    "bitscore": cols[11].strip(),
                    # Older raw TSVs were written without qcovhsp; treat it as absent
                    # rather than failing to parse them.
                    "qcovhsp": cols[12].strip() if len(cols) > 12 else "",
                    "_bitscore": bitscore,
                }
    return best


def run_custom_blast(
    query_fasta: Path,
    reference_fasta: Path,
    output_dir: Path,
    *,
    perc_identity: float = 90.0,
    evalue: float = 1e-5,
    max_target_seqs: int = 5,
    qcov_hsp_perc: float = 0.0,
    threads: int = 4,
) -> StepOutput:
    """BLAST the OTU/ASV representatives against a custom reference FASTA.

    Writes ``custom_blast_besthit.csv`` with one row per query sequence; queries
    without a hit get ``status="None"`` (rendered red in the UI). Returns a
    :class:`StepOutput` whose ``metrics`` report matched/total counts.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "custom_blast.log"
    besthit_path = output_dir / "custom_blast_besthit.csv"
    raw_tsv = output_dir / "custom_blast_raw.tsv"

    if not blast_tools_available():
        message = "blastn/makeblastdb not found on PATH. Install BLAST+ to compare against a custom database."
        log_path.write_text(message, encoding="utf-8")
        return StepOutput(step_name="custom_blast", success=False, files={"log": log_path}, log=log_path, error_message=message)

    query_fasta = Path(query_fasta)
    reference_fasta = Path(reference_fasta)
    if not query_fasta.exists():
        message = f"Query FASTA not found: {query_fasta}"
        log_path.write_text(message, encoding="utf-8")
        return StepOutput(step_name="custom_blast", success=False, files={"log": log_path}, log=log_path, error_message=message)
    if not reference_fasta.exists():
        message = f"Custom database not found: {reference_fasta}"
        log_path.write_text(message, encoding="utf-8")
        return StepOutput(step_name="custom_blast", success=False, files={"log": log_path}, log=log_path, error_message=message)

    log_lines: List[str] = []
    try:
        db_prefix = build_blast_db(reference_fasta, output_dir / "db")
    except subprocess.CalledProcessError as exc:
        message = f"makeblastdb failed: {str(exc.stderr or exc.stdout or exc)[:500]}"
        log_path.write_text(message, encoding="utf-8")
        return StepOutput(step_name="custom_blast", success=False, files={"log": log_path}, log=log_path, error_message=message)
    except subprocess.TimeoutExpired:
        message = f"makeblastdb exceeded the timeout of {MAKEBLASTDB_TIMEOUT_SEC}s."
        log_path.write_text(message, encoding="utf-8")
        return StepOutput(step_name="custom_blast", success=False, files={"log": log_path}, log=log_path, error_message=message)

    command = [
        BLASTN_PATH,
        "-query",
        str(query_fasta),
        "-db",
        str(db_prefix),
        "-out",
        str(raw_tsv),
        "-outfmt",
        # qcovhsp is appended last so the existing positional column indices below
        # keep working; identity without coverage cannot distinguish a full-length
        # match from a strong hit over a short fragment of the query.
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp",
        "-evalue",
        str(evalue),
        "-perc_identity",
        str(perc_identity),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(max(1, int(threads))),
    ]
    if qcov_hsp_perc and qcov_hsp_perc > 0:
        command.extend(["-qcov_hsp_perc", str(qcov_hsp_perc)])
    # Judged by output, not exit code: this is the same blastn that exits 0 after printing
    # "(Severe Error) ... search failed" and writing an empty result. Reading only the exit
    # code here would reproduce, on the custom database, exactly the bug that made the NCBI
    # step report a refusal as `no_hits` — and a 0-byte custom_blast_raw.tsv is precisely
    # what the fix plan flagged.
    result = run_tool(command, timeout=BLAST_TIMEOUT_SEC, failure_signatures=BLAST_FAILURE_SIGNATURES)
    log_lines.append(result.output)
    if result.timed_out:
        message = f"blastn exceeded the timeout of {BLAST_TIMEOUT_SEC}s."
        log_path.write_text(message, encoding="utf-8")
        return StepOutput(step_name="custom_blast", success=False, files={"log": log_path}, log=log_path, error_message=message)
    if not result.ok:
        reason = result.failure_signature or result.output[:400] or f"exit {result.returncode}"
        message = f"blastn failed: {reason[:500]}"
        log_path.write_text("\n".join(log_lines) + f"\n{message}", encoding="utf-8")
        return StepOutput(step_name="custom_blast", success=False, files={"log": log_path}, log=log_path, error_message=message)

    best = _parse_best_hits(raw_tsv)
    labels = reference_labels(reference_fasta)
    query_ids = _query_ids(query_fasta)

    matched = 0
    rows: List[Dict[str, str]] = []
    for qid in query_ids:
        hit = best.get(qid)
        if hit:
            matched += 1
            subject_clean = _sequence_id(hit["subject_id"])
            species = labels.get(subject_clean) or labels.get(hit["subject_id"]) or subject_clean
            rows.append(
                {
                    "OTU_ID": qid,
                    "subject_id": subject_clean,
                    "species": species,
                    "pident": hit["pident"],
                    "qcovhsp": hit.get("qcovhsp", ""),
                    "length": hit["length"],
                    "evalue": hit["evalue"],
                    "bitscore": hit["bitscore"],
                    "status": "Match",
                }
            )
        else:
            rows.append(
                {
                    "OTU_ID": qid,
                    "subject_id": "",
                    "species": "None",
                    "pident": "",
                    "qcovhsp": "",
                    "length": "",
                    "evalue": "",
                    "bitscore": "",
                    "status": "None",
                }
            )

    with besthit_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_BESTHIT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    total = len(query_ids)
    log_lines.append(f"Custom BLAST: {matched}/{total} OTUs matched (perc_identity>={perc_identity}, evalue<={evalue}).")
    log_path.write_text("\n".join(l for l in log_lines if l), encoding="utf-8")

    return StepOutput(
        step_name="custom_blast",
        success=True,
        files={"custom_blast_besthit_csv": besthit_path, "custom_blast_raw_tsv": raw_tsv, "log": log_path},
        log=log_path,
        metrics={"matched": float(matched), "total": float(total), "unmatched": float(total - matched)},
    )


def get_custom_reference(marker: str) -> str:
    return f"custom:{marker.lower()}"
