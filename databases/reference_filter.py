"""Marker-filtered reference databases for VSEARCH SINTAX.

SILVA SSU (824 MB) and PR2 (334 MB) describe every small-subunit rRNA that has ever
been sequenced. A run does not need that: it needs the one marker region it actually
amplified, in the one domain it targeted. VSEARCH ``--sintax`` loads the whole
reference plus a k-mer index into RAM (~3x the FASTA), so on a small machine the full
SILVA asks for ~2.5 GB, does not get it, and the source is skipped — leaving only the
remote sources (NCBI/BOLD), which are slow and often return a species with no lineage.

Cutting the reference down to the marker fixes that at the root:

    SILVA 824 MB  --(keep Eukaryota)-->  110 MB  --(cut to the V4 amplicon)-->  ~28 MB

which fits in RAM, runs locally in seconds, and returns a full kingdom-to-species
lineage. It also makes the reference *more* correct, not less: SINTAX scores k-mer
similarity, so a 380 bp query is best matched against 380 bp references covering the
same region rather than against 1 800 bp whole-SSU sequences whose bulk is region the
query never spanned.

The amplicon is cut with the run's *own* primers (``PRIMER_PRESETS`` /
``primer_fwd``+``primer_rev``), so the reference and the query always describe the same
stretch of the gene. That is the whole point — a reference trimmed with different
primers than the reads would be worse than no trimming at all.

Everything here degrades rather than fails. If the primers are unknown, or cutadapt is
missing, or the primers barely match this reference (wrong marker for this database),
the ladder falls back a rung at a time:

    1. lineage + amplicon   (~28 MB — local, fast, full lineage)
    2. lineage only         (~110 MB — still fits, still local)
    3. the original file    (the memory gate then skips it, as it did before)

The result is cached next to the source FASTA and rebuilt only when the source or this
module's format version changes, so the ~2 minute build is paid once, not once per run.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

from core.config import CUTADAPT_PATH

# Bump to invalidate every cached reference built by an older version of this module.
FILTER_FORMAT_VERSION = "v1"

# Which slice of life each marker addresses. SILVA SSU carries all three domains, but
# an 18S primer pair only ever amplifies eukaryotes: the bacteria and archaea are 89%
# of the file and 0% of the answer. PR2 is already eukaryote-only, so this is a no-op
# there and the amplicon step does the work.
# COI is deliberately absent: BOLD is the COI reference and it is remote, and neither
# SILVA nor PR2 contains COI at all.
MARKER_LINEAGE = {
    "18S": ("d:Eukaryota",),
    "28S": ("d:Eukaryota",),
    "ITS": ("d:Eukaryota",),
}

# Bounds on the extracted amplicon. Measured on PR2's V4: median 381 bp, p99 629 bp.
# 800 keeps the divergent expansion segments (Opisthokonta, Foraminifera) that are real
# biology, while still dropping the products of a primer that matched something absurd.
AMPLICON_MIN_LENGTH = 150
AMPLICON_MAX_LENGTH = 800

# cutadapt's mismatch budget for the primer match. Primers are degenerate and reference
# sequences are imperfect; 0.15 (3 mismatches in a 20-mer) is the usual amplicon-region
# extraction setting. Too strict here silently throws away good references.
PRIMER_ERROR_RATE = 0.15

# If the primer step keeps less than this share of the lineage-filtered records, the
# primers do not belong to this reference (wrong marker, wrong database) and the result
# would be a reference that classifies nothing. Fall back to lineage-only instead.
MIN_AMPLICON_YIELD = 0.20

_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


def reverse_complement(sequence: str) -> str:
    """Reverse complement, IUPAC-aware (``R``<->``Y``, ``K``<->``M``, ...)."""
    return sequence.translate(_COMPLEMENT)[::-1]


def cutadapt_available() -> bool:
    return bool(CUTADAPT_PATH and (shutil.which(CUTADAPT_PATH) or Path(CUTADAPT_PATH).is_file()))


def normalize_marker(value) -> str:
    """Map a marker label (``Marker`` enum, ``"18S"``, ``"ssu_18s"``...) to a scope key."""
    raw = str(getattr(value, "value", value) or "").strip().upper()
    if not raw:
        return ""
    for key in MARKER_LINEAGE:
        if key in raw:
            return key
    return raw


def resolve_marker(markers: Iterable) -> str:
    """The single marker a run targets, or ``""`` if the run mixes markers.

    A mixed-marker run has no one amplicon to cut the reference to, so it gets the
    unfiltered reference rather than a reference that is wrong for half its samples.
    """
    distinct = {normalize_marker(marker) for marker in markers or ()}
    distinct.discard("")
    if len(distinct) != 1:
        return ""
    return distinct.pop()


# ``OTHER`` is this codebase's word for "we could not tell", not a marker in its own
# right: it is what ``Marker`` falls back to and what an unmatched file name produces.
# Treating it as a real scope is what let a COI/Leray run cut SILVA to "the OTHER
# marker" — a filter with no referent.
_UNKNOWN_MARKERS = {"", "OTHER"}


def marker_is_known(marker) -> bool:
    """True when ``marker`` names a real amplicon we can scope a reference to."""
    return normalize_marker(marker) not in _UNKNOWN_MARKERS


def preset_marker(params) -> str:
    """The marker implied by the run's primer preset, or ``""``.

    Imported lazily and defensively for the same reason ``cutadapt_runner`` does it:
    ``databases`` must not hard-depend on ``analysis`` to answer a question about
    configuration.
    """
    try:
        from analysis.vsearch_runner import PRESET_MARKER
    except ImportError:  # pragma: no cover - analysis package always ships with us
        return ""
    preset = str((params or {}).get("primer_preset") or "").strip()
    return normalize_marker(PRESET_MARKER.get(preset, ""))


def resolve_run_marker(markers: Iterable, params=None) -> Tuple[str, str]:
    """The marker a run targets, as ``(marker, note)``.

    The primer preset wins over the marker guessed from the sample file names. The
    preset is the one thing the user states outright about what was amplified; the
    file-name marker comes from a regex looking for ``coi``/``18s``/``its`` tokens in
    the stem, which finds nothing in a name like ``aptl64b_01`` and silently yields
    ``OTHER``. A run configured COI/Leray then had its rRNA references cut with COI
    primers — minutes spent building a reference that cannot match anything.

    ``note`` is non-empty only when the two sources disagree. The disagreement does not
    stop the run (the preset is still obeyed), but it is recorded, because a sample
    named ``..._18S_...`` amplified with COI primers is a configuration worth seeing.
    """
    from_samples = resolve_marker(markers)
    from_preset = preset_marker(params)
    if not from_preset:
        return from_samples, ""
    if from_samples and from_samples != from_preset:
        note = (
            f"marcador do preset de primers ({from_preset}) difere do inferido dos nomes "
            f"das amostras ({from_samples}); a run segue como {from_preset}"
        )
        return from_preset, note
    return from_preset, ""


def marker_filter_applies(marker, primers: Tuple[str, str] = ("", "")) -> bool:
    """Whether a run with this marker/primers gets a filtered reference.

    Lets the pre-run checks predict the *filtered* size without importing the builder
    or paying for a build — a "SILVA needs 2.5 GB and will be skipped" warning is a lie
    once the filter is going to reduce it to a tenth of that before it is ever loaded.
    """
    if not marker_is_known(marker):
        # No marker, no basis for a cut. Reaching the primer branch below with an
        # unknown marker is what let COI primers be taken to SILVA: cutadapt happily
        # cuts an rRNA reference with them and leaves almost nothing behind, and the
        # remains get loaded as if they were a reference. Let the memory gate decide
        # about the full file instead — it at least says why it skipped.
        return False
    forward, reverse = (primers or ("", ""))
    if normalize_marker(marker) in MARKER_LINEAGE:
        return True
    return bool(str(forward or "").strip() and str(reverse or "").strip() and cutadapt_available())


@dataclass(frozen=True)
class FilteredReference:
    """A reference cut down to one marker, and the story of how it got that small."""

    path: Path
    source_path: Path
    marker: str
    records: int
    source_records: int
    size_bytes: int
    source_size_bytes: int
    steps: Tuple[str, ...] = ()

    @property
    def shrink_factor(self) -> float:
        if not self.size_bytes:
            return 0.0
        return self.source_size_bytes / self.size_bytes

    def summary(self) -> str:
        """One sentence for the log and the per-source status panel."""
        return (
            f"referência filtrada para {self.marker}: {self.size_bytes / 1e6:.0f} MB "
            f"({self.records:,} sequências) de {self.source_size_bytes / 1e6:.0f} MB "
            f"({self.source_records:,}) — {'; '.join(self.steps)}"
        ).replace(",", " ")


def _filter_by_lineage(source: Path, dest: Path, lineage: Sequence[str]) -> Tuple[int, int]:
    """Keep records whose SINTAX header matches ``lineage``; return ``(kept, total)``.

    Streams a record at a time — the input is close to a gigabyte and the machine this
    exists for has less than that free. Two normalizations happen on the way through:

    * ``U`` -> ``T``. SILVA ships RNA; primers are DNA. Without this the primer match in
      the next step finds nothing. Only sequence lines are touched — taxon names in the
      header contain ``U`` too (``Unruhdinium``), and rewriting those would corrupt the
      taxonomy.
    * Wrapped sequences are unwrapped to a single line, which is what makes the length
      filter downstream mean what it says.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_dest = dest.with_suffix(dest.suffix + ".tmp")
    kept = 0
    total = 0
    with source.open("r", encoding="utf-8", errors="replace") as src, tmp_dest.open("w", encoding="utf-8") as dst:
        header: Optional[str] = None
        chunks: list[str] = []

        def flush() -> None:
            if header is not None and chunks:
                dst.write(header)
                dst.write("".join(chunks).replace("U", "T").replace("u", "t"))
                dst.write("\n")

        for line in src:
            if line.startswith(">"):
                flush()
                total += 1
                keep = not lineage or any(tag in line for tag in lineage)
                if keep:
                    kept += 1
                    header = line if line.endswith("\n") else line + "\n"
                else:
                    header = None
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
        flush()
    tmp_dest.replace(dest)
    return kept, total


def _extract_amplicon(source: Path, dest: Path, forward: str, reverse: str, threads: int = 2) -> Optional[int]:
    """Cut each reference down to the region the primers amplify. ``None`` if it failed.

    The reverse primer is given as it was ordered (5'->3' on the opposite strand), so it
    has to be reverse-complemented to be found on the same strand as the forward primer.

    Both primers are marked ``required``: cutadapt's default for a linked adapter makes
    the second one optional, which quietly keeps the whole 1 700 bp sequence whenever the
    reverse primer is missing — the exact sequences the filter exists to remove.
    """
    if not cutadapt_available() or not (forward and reverse):
        return None
    linked = f"{forward.strip()};required...{reverse_complement(reverse.strip())};required"
    command = [
        CUTADAPT_PATH,
        "-g", linked,
        "-e", str(PRIMER_ERROR_RATE),
        "--match-read-wildcards",
        "--discard-untrimmed",
        "-m", str(AMPLICON_MIN_LENGTH),
        "-M", str(AMPLICON_MAX_LENGTH),
        "-j", str(max(1, int(threads))),
        "-o", str(dest),
        str(source),
    ]
    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, OSError):
        dest.unlink(missing_ok=True)
        return None
    return _count_records(dest)


def _count_records(path: Path) -> int:
    count = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    count += 1
    except OSError:
        return 0
    return count


def _cache_paths(source: Path, marker: str, primer_tag: str) -> Tuple[Path, Path]:
    """``(fasta, sidecar)`` for this source+marker+primer combination.

    Built by string concatenation, not ``Path.with_suffix``: these names are full of
    dots (``...sintax.v2.18S.filtered.v1``) and ``with_suffix`` would treat the format
    version as the extension and replace it — dropping the very token that lets a new
    version of this module invalidate an old cache.
    """
    stem = source.name.split(".fasta")[0]
    base = f"{stem}.{marker}{primer_tag}.filtered.{FILTER_FORMAT_VERSION}"
    return source.parent / f"{base}.fasta", source.parent / f"{base}.json"


def _load_cached(fasta: Path, sidecar: Path, source: Path) -> Optional[FilteredReference]:
    """Return the cached reference, or ``None`` if it is absent or out of date."""
    if not fasta.exists() or not sidecar.exists():
        return None
    try:
        if fasta.stat().st_mtime < source.stat().st_mtime:
            return None
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    try:
        return FilteredReference(
            path=fasta,
            source_path=source,
            marker=str(payload["marker"]),
            records=int(payload["records"]),
            source_records=int(payload["source_records"]),
            size_bytes=int(payload["size_bytes"]),
            source_size_bytes=int(payload["source_size_bytes"]),
            steps=tuple(payload.get("steps") or ()),
        )
    except (KeyError, TypeError, ValueError):
        return None


def build_marker_reference(
    source_fasta: Path,
    marker,
    primers: Tuple[str, str] = ("", ""),
    threads: int = 2,
) -> Optional[FilteredReference]:
    """Build (or fetch from cache) a marker-sized copy of ``source_fasta``.

    ``None`` means "no filter applies here" — an unknown marker, or a marker this
    reference has no scope for — and the caller should use the original file.
    """
    source = Path(source_fasta)
    if not marker_is_known(marker):
        # The docstring already promised this ("an unknown marker"), but the code only
        # honoured it when there were no primers: with primers set, an OTHER run fell
        # through to the amplicon cut and spent minutes trimming a reference to a
        # marker that was never identified.
        return None
    marker_key = normalize_marker(marker)
    lineage = MARKER_LINEAGE.get(marker_key)
    forward, reverse = (primers or ("", ""))
    forward, reverse = str(forward or "").strip(), str(reverse or "").strip()
    can_cut = bool(forward and reverse and cutadapt_available())
    if not source.exists() or (not lineage and not can_cut):
        return None

    # The primers are part of the cache identity: a reference cut for V4 is not the
    # reference a V9 run needs, and handing one to the other would be silently wrong.
    primer_tag = f".{forward[:8]}_{reverse[:8]}".upper() if can_cut else ""
    fasta_path, sidecar_path = _cache_paths(source, marker_key or "ALL", primer_tag)
    cached = _load_cached(fasta_path, sidecar_path, source)
    if cached is not None:
        return cached

    source_records = 0
    steps: list[str] = []
    lineage_path = fasta_path.parent / f"{fasta_path.name}.lineage.tmp"
    try:
        kept, source_records = _filter_by_lineage(source, lineage_path, lineage or ())
        if not kept:
            return None
        # PR2 is already eukaryote-only, so on an 18S run its lineage filter removes
        # nothing. Without primers to cut the amplicon there is no second step either,
        # and carrying on would write a full-size copy of the reference to disk and call
        # it "filtered". Nothing was gained; say so by returning None.
        if kept == source_records and not can_cut:
            return None
        if lineage and kept < source_records:
            steps.append(f"{' + '.join(tag.split(':')[-1] for tag in lineage)} apenas ({kept:,} de {source_records:,})".replace(",", " "))

        records = kept
        if can_cut:
            amplicon_path = fasta_path.parent / f"{fasta_path.name}.amplicon.tmp"
            cut = _extract_amplicon(lineage_path, amplicon_path, forward, reverse, threads=threads)
            # A primer pair that matches almost nothing is a primer pair that does not
            # belong to this reference. Keeping that result would trade a database that
            # is merely too big for one that is the right size and classifies nothing.
            if cut and cut >= max(1, int(kept * MIN_AMPLICON_YIELD)):
                amplicon_path.replace(lineage_path)
                steps.append(f"recortada aos primers {forward[:12]}…/{reverse[:12]}… ({cut:,} de {kept:,})".replace(",", " "))
                records = cut
            else:
                amplicon_path.unlink(missing_ok=True)
                if cut is not None:
                    steps.append("primers não correspondem a esta base — só filtro taxonómico")
        if not steps:
            return None

        lineage_path.replace(fasta_path)
        reference = FilteredReference(
            path=fasta_path,
            source_path=source,
            marker=marker_key or "ALL",
            records=records,
            source_records=source_records,
            size_bytes=fasta_path.stat().st_size,
            source_size_bytes=source.stat().st_size,
            steps=tuple(steps),
        )
        sidecar_path.write_text(
            json.dumps(
                {
                    "marker": reference.marker,
                    "records": reference.records,
                    "source_records": reference.source_records,
                    "size_bytes": reference.size_bytes,
                    "source_size_bytes": reference.source_size_bytes,
                    "steps": list(reference.steps),
                    "source": str(source),
                    "primers": [forward, reverse] if can_cut else [],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return reference
    except OSError:
        return None
    finally:
        lineage_path.unlink(missing_ok=True)
