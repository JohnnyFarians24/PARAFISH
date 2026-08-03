"""vsearch execution helpers for OTU workflows.

The OTU workflow is split in two phases so that clustering happens ONCE across
all samples pooled together, instead of once per sample:

1. run_vsearch_filter_sample(): quality-filters one sample's reads and tags
   them with --sample <id>, which VSEARCH embeds as a `;sample=<id>;`
   annotation in each read header.
2. run_vsearch_pool_and_cluster(): concatenates every sample's filtered reads,
   dereplicates and clusters that pool ONCE, removes chimeras ONCE on the
   resulting centroids, then maps the whole pool back onto the final
   chimera-free OTU set with --usearch_global --otutabout. VSEARCH reads the
   `sample=` annotation to build a real per-sample abundance matrix natively
   (see the VSEARCH OTU-clustering tutorial), so every sample ends up counted
   against the *same* OTU numbering instead of each sample owning its own.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Dict

from core.config import MIN_UNIQUE_SIZE, VSEARCH_PATH
from core.models import RunConfig, StepOutput
from core.utils import _count_fasta_seqs


def _run(command):
    return subprocess.run(command, capture_output=True, text=True, check=True)


def _normalize_sample_headers(fasta_path: Path, sample_id: str) -> None:
    """Rewrite FASTA headers to ``>{read_id};sample={sample_id}`` (no whitespace).

    Guarantees the per-sample annotation survives VSEARCH's header truncation so
    ``--otutabout`` builds a correct read_id/OTU × sample abundance table keyed by
    the real sample id rather than the sequencer/read-name prefix.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        return
    tmp = fasta_path.with_suffix(fasta_path.suffix + ".tmp")
    with fasta_path.open("r", encoding="utf-8", errors="replace") as src, tmp.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith(">"):
                header = line[1:].strip()
                read_id = header.split()[0].split(";")[0] if header else ""
                dst.write(f">{read_id};sample={sample_id}\n")
            else:
                dst.write(line)
    tmp.replace(fasta_path)


# Primer pairs we ship by name, so a user does not have to paste sequences.
# TAReuk454FWD1 / TAReukREV3 is the standard eukaryotic 18S V4 pair — and the one
# that amplifies fish host DNA overwhelmingly in a parasite study, which is worth
# knowing before the run rather than after.
PRIMER_PRESETS = {
    "18S_V4_TAReuk": ("CCAGCASCYGCGGTAATTCC", "ACTTTCGTTCTTGATYRA"),
    "COI_Folmer": ("GGTCAACAAATCATAAAGATATTGG", "TAAACTTCAGGGTGACCAAAAAATCA"),
    "COI_Leray": ("GGWACWGGWTGAACWGTWTAYCCYCC", "TANACYTCNGGRTGNCCRAARAAYCA"),
}

# Which marker each preset amplifies. Chosen primers are the one *explicit* statement
# the user makes about what was sequenced; before this map existed the marker was
# guessed from the FASTQ file name alone, so a COI/Leray run on samples named
# `aptl64b_01` resolved to OTHER and the rRNA references were cut with COI primers.
# Kept beside PRIMER_PRESETS, not derived from the name prefix, so a preset added
# without a marker fails the paired test below instead of being silently mis-typed.
PRESET_MARKER = {
    "18S_V4_TAReuk": "18S",
    "COI_Folmer": "COI",
    "COI_Leray": "COI",
}


def _primer_strip_lengths(params: Dict) -> tuple[int, int]:
    """Bases to strip from each end of a merged amplicon.

    Accepts either explicit lengths (``primer_strip_left``/``primer_strip_right``),
    primer sequences (``primer_fwd``/``primer_rev``), or a preset name
    (``primer_preset``). Returns ``(0, 0)`` when no primers are configured — the
    pipeline never guesses, because stripping the wrong number of bases would
    quietly truncate real sequence.
    """
    # If the Preprocess stage already removed the primers with cutadapt (sequence-aware),
    # stripping again by length here would truncate real biological sequence. The
    # Preprocess stage sets this flag once cutadapt has run for the configured primers.
    if params.get("primers_removed_upstream"):
        return 0, 0

    left = params.get("primer_strip_left")
    right = params.get("primer_strip_right")
    if left is not None or right is not None:
        return int(left or 0), int(right or 0)

    forward = str(params.get("primer_fwd") or "").strip()
    reverse = str(params.get("primer_rev") or "").strip()
    preset = str(params.get("primer_preset") or "").strip()
    if not (forward or reverse) and preset in PRIMER_PRESETS:
        forward, reverse = PRIMER_PRESETS[preset]
    return len(forward), len(reverse)


def run_vsearch_filter_sample(assembled_fastq: Path, config: RunConfig, output_dir: Path, sample_id: str) -> StepOutput:
    """Quality-filter one sample's reads, tagging them for later pooling."""
    assembled = Path(assembled_fastq)
    if not assembled.exists():
        return StepOutput(step_name="vsearch_filter", success=False, files={}, error_message=f"Assembled FASTQ not found: {assembled}.")

    output_dir.mkdir(parents=True, exist_ok=True)
    filtered = output_dir / "filtered.fasta"
    log_path = output_dir / "step1_filter.log"
    params = config.vsearch_params

    command = [
        VSEARCH_PATH,
        "--fastq_filter",
        str(assembled),
        "--sample",
        str(sample_id),
        "--fastq_maxee",
        str(params.get("fastq_maxee", 1.0)),
        "--fastq_minlen",
        str(params.get("fastq_minlen", 200)),
        "--fastaout",
        str(filtered),
        "--log",
        str(log_path),
    ]

    # Strip the PCR primers. Trimmomatic's ILLUMINACLIP removes *sequencing
    # adapters*, not the amplification primers, so without this the primer — and
    # its degenerate positions — stays in every OTU and is fed to the classifier
    # and the tree as if it were biological signal. On merged amplicons the
    # forward primer occupies the first bases and the reverse primer's complement
    # the last, so a fixed-length strip is exact.
    strip_left, strip_right = _primer_strip_lengths(params)
    if strip_left:
        command.extend(["--fastq_stripleft", str(strip_left)])
    if strip_right:
        command.extend(["--fastq_stripright", str(strip_right)])
    try:
        _run(command)
    except subprocess.CalledProcessError as exc:
        return StepOutput(step_name="vsearch_filter", success=False, files={"log": log_path}, log=log_path, error_message=f"VSEARCH filter failed: {exc}")

    # Illumina FASTQ headers carry a description after a space
    # (e.g. "A01426:...:1423 1:N:0:ACGTATGG+CTCTTGTC"). VSEARCH's --sample appends
    # ";sample=<id>" AFTER that description, but VSEARCH truncates every header at
    # the first whitespace on read — so the sample tag is silently dropped and the
    # per-sample OTU table ends up grouped by the read-name prefix (the sequencer
    # id) instead of the actual sample. Rewrite each header to
    # ">{read_id};sample={sample_id}" (no whitespace) so the tag survives
    # dereplication, clustering and --otutabout.
    _normalize_sample_headers(filtered, str(sample_id))

    if _count_fasta_seqs(filtered) == 0:
        return StepOutput(
            step_name="vsearch_filter",
            success=False,
            files={"filtered_fasta": filtered, "log": log_path},
            log=log_path,
            error_message="Filtering produced 0 reads. Check fastq_maxee/fastq_minlen.",
        )

    return StepOutput(step_name="vsearch_filter", success=True, files={"filtered_fasta": filtered, "log": log_path}, log=log_path)


# The run's singleton policy is a CHOICE, not a constant: ``min_unique_size`` decides the
# smallest cluster allowed to SEED an OTU (default 1 — see core.config.MIN_UNIQUE_SIZE).
#
# Whatever the floor, dereplication only decides what may *seed* a cluster. The mapping
# step then runs ``--usearch_global`` over the *full* pooled file, so a read excluded from
# seeding is still counted against whichever OTU it matches at the identity threshold.
# That is why the two halves must stay in step: raising the floor without knowing the
# mapping counts everything anyway is what made the pipeline report 147 010 reads in some
# files and 157 527 in others (see analysis_metrics.chimera_stats).
#
# What the floor actually costs, stated precisely, because "singletons are error" is the
# folklore that got this wrong for a long time:
#
#   * floor = 1: every unique sequence may seed. Error-derived reads do NOT become their
#     own OTUs — they sit within 3% of the sequence they are errors of, so clustering
#     absorbs them. Measured on the reference run: identical output to floor = 2.
#   * floor = 2: a read that is globally unique AND >3% from every other sequence seeds
#     nothing and then matches nothing — it leaves as `unmapped_reads`. That is not the
#     error case; it is the rare-novel-taxon case, i.e. exactly what a parasite-detection
#     platform exists to find.
#
# Hence the default of 1. A floor of 2 buys clustering time (57 s -> 247 s on the
# reference run, in the other direction) and remains available; `unmapped_reads` in
# read_tracking.csv is what keeps its cost visible rather than silent.
SINGLETON_POLICY = "singletons may seed OTUs (--minuniquesize 1 by default); whatever the floor, all reads are counted in abundance via mapping"

_DEREP_DISCARDED = re.compile(r"([\d]+)\s+uniques written,\s+([\d]+)\s+clusters discarded")


def _derep_singleton_counts(derep_log: Path) -> Dict[str, float]:
    """``{uniques_written, unique_seqs_discarded}`` from VSEARCH's dereplication log."""
    try:
        match = _DEREP_DISCARDED.search(Path(derep_log).read_text(encoding="utf-8", errors="replace"))
    except OSError:
        return {}
    if not match:
        return {}
    return {"uniques_written": float(match.group(1)), "unique_seqs_discarded": float(match.group(2))}


def otu_table_totals(otu_table: Path) -> Dict[str, int]:
    """Reads per OTU from the mapping table — the run's canonical abundance.

    ``otu_table.tsv`` is written by ``--usearch_global --otutabout``, which maps every
    pooled read back onto the final OTU set. That is the one number the whole pipeline
    must agree on; anything derived from the clustering instead measures a different
    population (see :func:`_write_otu_id_map`).
    """
    totals: Dict[str, int] = {}
    otu_table = Path(otu_table)
    if not otu_table.exists():
        return totals
    import csv as _csv

    with otu_table.open("r", errors="replace") as handle:
        reader = _csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            try:
                totals[row[0].strip()] = sum(int(float(cell or 0)) for cell in row[1:])
            except ValueError:
                continue
    return totals


def _write_otu_id_map(otus_fasta: Path, otu_table: Path, out_path: Path) -> Path:
    """Map each OTU id to its abundance and the representative read it came from.

    ``size`` is the OTU's read total **from the mapping table**, not the ``;size=`` in
    the FASTA header. Those two numbers are not the same thing and used to disagree
    badly: the header carries the cluster size from ``--cluster_size --sizeout``, which
    (a) excludes every read dropped by ``--minuniquesize`` at dereplication and
    (b) counts the reads that *seeded* the cluster, not the reads that ultimately map
    to it at the identity threshold. On the reference run those two effects put OTU2 at
    14 144 by clustering and 63 260 by mapping, and made the run total 147 010 in some
    files and 157 527 in others.

    The clustering number is kept as ``cluster_size`` rather than discarded: the gap
    between the two columns is real information (it is where read redistribution shows
    up), and hiding it would only make the next person rediscover this bug.
    """
    import csv as _csv
    import re as _re

    totals = otu_table_totals(otu_table)
    rows = []
    if otus_fasta.exists():
        with otus_fasta.open("r", errors="replace") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                new_id, _, original = header.partition(" ")
                otu_id = new_id.split(";")[0]
                cluster_size = _re.search(r"size=(\d+)", new_id)
                rows.append(
                    {
                        "otu_id": otu_id,
                        # An OTU absent from the mapping table genuinely has 0 mapped
                        # reads; that is a fact worth writing, not a blank.
                        "size": totals.get(otu_id, 0),
                        "cluster_size": int(cluster_size.group(1)) if cluster_size else "",
                        "representative_read": original.strip(),
                    }
                )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = _csv.DictWriter(handle, fieldnames=["otu_id", "size", "cluster_size", "representative_read"])
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def run_vsearch_pool_and_cluster(filtered_fastas: Dict[str, Path], config: RunConfig, output_dir: Path) -> StepOutput:
    """Pool every sample's filtered reads and cluster them ONCE into a shared OTU set.

    filtered_fastas maps sample_id -> its run_vsearch_filter_sample() output
    (already tagged with --sample). Returns a StepOutput whose "otus_fasta" is
    the final, chimera-free, cross-sample-consistent OTU representative set,
    and "otu_table_tsv" is a real OTU x sample read-abundance matrix.
    """
    identity = float(config.vsearch_params.get("id", 0.97))
    if identity < 0.8 or identity > 1.0:
        return StepOutput(step_name="vsearch_pool", success=False, files={}, error_message="Parameter --id is outside the allowed range [0.8, 1.0].")

    output_dir.mkdir(parents=True, exist_ok=True)
    pooled = output_dir / "pooled_filtered.fasta"
    derep = output_dir / "derep.fasta"
    otus_raw = output_dir / "otus_raw.fasta"
    otus = output_dir / "otus.fasta"
    chimeras = output_dir / "chimeras.fasta"
    chimera_uchime = output_dir / "chimera.uchime"
    clusters_uc = output_dir / "clusters.uc"
    map_uc = output_dir / "map.uc"
    otu_table = output_dir / "otu_table.tsv"
    combined_log = output_dir / "vsearch_pool.log"

    with pooled.open("w", encoding="utf-8") as out_handle:
        for path in filtered_fastas.values():
            path = Path(path)
            if path.exists():
                out_handle.write(path.read_text(encoding="utf-8", errors="replace"))

    if _count_fasta_seqs(pooled) == 0:
        return StepOutput(step_name="vsearch_pool", success=False, files={}, error_message="No sample produced filtered reads to pool.")

    threads = str(config.vsearch_params.get("threads") or 8)
    # The smallest cluster allowed to seed an OTU; see SINGLETON_POLICY above. The
    # fallback tracks the configured default rather than repeating a literal — this line
    # said 2 while the Setup page said 1, so the "policy" depended on which entry point
    # you came in through.
    min_unique_size = int(config.vsearch_params.get("min_unique_size", MIN_UNIQUE_SIZE) or 1)

    derep_command = [VSEARCH_PATH, "--derep_fulllength", str(pooled), "--sizeout", "--output", str(derep), "--log", str(output_dir / "step2_pool_derep.log")]
    if min_unique_size > 1:
        derep_command.extend(["--minuniquesize", str(min_unique_size)])

    pieces = []
    commands = [
        derep_command,
        # --sizein/--sizeout are what make "cluster by size" and the chimera step's
        # abundance-skew test actually see abundances. Without them VSEARCH assumes
        # every sequence has abundance 1, so uchime_denovo judged parent/chimera
        # skew on dereplication counts it could not read.
        [VSEARCH_PATH, "--cluster_size", str(derep), "--id", str(identity), "--sizein", "--sizeout", "--centroids", str(otus_raw), "--uc", str(clusters_uc), "--log", str(output_dir / "step3_cluster.log")],
        # --relabel gives the surviving OTUs stable, readable ids (OTU1..OTUn).
        # Without it every downstream table is keyed by the raw Illumina header of
        # whichever read happened to be the cluster centroid
        # ("A01426:1215:HW7YLDRX5:2:2101:6289:3458"), which is not an identifier
        # anyone can cite. --relabel_keep preserves the original header after the
        # new one, so the representative read is still traceable.
        [
            VSEARCH_PATH, "--uchime_denovo", str(otus_raw), "--sizein", "--sizeout",
            "--relabel", "OTU", "--relabel_keep",
            "--nonchimeras", str(otus), "--chimeras", str(chimeras),
            "--uchimeout", str(chimera_uchime), "--log", str(output_dir / "step4_chimera.log"),
        ],
    ]
    try:
        for command in commands:
            completed = _run(command)
            pieces.append((completed.stdout or "") + "\n" + (completed.stderr or ""))
    except subprocess.CalledProcessError as exc:
        combined_log.write_text("\n".join(pieces) + "\n" + str(exc), encoding="utf-8")
        return StepOutput(step_name="vsearch_pool", success=False, files={}, log=combined_log, error_message=f"VSEARCH (pool/cluster/chimera) failed: {exc}")

    otu_map = output_dir / "otu_id_map.csv"

    if _count_fasta_seqs(otus) == 0:
        combined_log.write_text("\n".join(pieces), encoding="utf-8")
        return StepOutput(
            step_name="vsearch_pool",
            success=False,
            files={"otus_raw_fasta": otus_raw, "chimeras_fasta": chimeras, "chimera_uchime": chimera_uchime, "clusters_uc": clusters_uc, "log": combined_log},
            log=combined_log,
            error_message="Global clustering produced 0 OTUs after chimera removal. Check the --id parameter and the read quality.",
        )

    try:
        completed = _run([
            VSEARCH_PATH,
            "--usearch_global",
            str(pooled),
            "--db",
            str(otus),
            "--id",
            str(identity),
            "--uc",
            str(map_uc),
            "--otutabout",
            str(otu_table),
            "--threads",
            threads,
            "--log",
            str(output_dir / "step5_otutable.log"),
        ])
        pieces.append((completed.stdout or "") + "\n" + (completed.stderr or ""))
    except subprocess.CalledProcessError as exc:
        combined_log.write_text("\n".join(pieces) + "\n" + str(exc), encoding="utf-8")
        return StepOutput(step_name="vsearch_pool", success=False, files={"otus_fasta": otus}, log=combined_log, error_message=f"Building the OTU table failed: {exc}")

    # Written only now, and deliberately: its ``size`` is the OTU's total in the mapping
    # table above, so it cannot be produced before that table exists. It used to be
    # written right after clustering, which is precisely why it carried the clustering's
    # abundance instead of the run's canonical one.
    _write_otu_id_map(otus, otu_table, otu_map)

    # Make the singleton policy an output, not folklore. Without these numbers the gap
    # between "reads that founded OTUs" and "reads counted in the table" is invisible,
    # which is how it went unnoticed for so long.
    pooled_reads = _count_fasta_seqs(pooled)
    mapped_reads = sum(otu_table_totals(otu_table).values())
    metrics = {
        "min_unique_size": float(min_unique_size),
        "pooled_reads": float(pooled_reads),
        "mapped_reads": float(mapped_reads),
        "unmapped_reads": float(pooled_reads - mapped_reads),
        **_derep_singleton_counts(output_dir / "step2_pool_derep.log"),
    }
    pieces.append(
        f"Singleton policy: {SINGLETON_POLICY} (--minuniquesize {min_unique_size}). "
        f"{pooled_reads} reads pooled; {mapped_reads} mapped to an OTU and counted in "
        f"otu_table; {pooled_reads - mapped_reads} matched no OTU at id={identity}."
    )

    combined_log.write_text("\n".join(pieces), encoding="utf-8")
    return StepOutput(
        step_name="vsearch_pool",
        success=True,
        metrics=metrics,
        files={
            "pooled_filtered_fasta": pooled,
            "derep_fasta": derep,
            "otus_raw_fasta": otus_raw,
            "otus_fasta": otus,
            "chimeras_fasta": chimeras,
            "chimera_uchime": chimera_uchime,
            "clusters_uc": clusters_uc,
            "map_uc": map_uc,
            "otu_table_tsv": otu_table,
            "otu_id_map_csv": otu_map,
            "log": combined_log,
        },
        log=combined_log,
    )


def run_vsearch_otus(input_fasta, output_dir, vsearch_path=None, identity=0.97):
    """Compatibility wrapper returning the legacy dict-based OTU result payload."""
    vpath = vsearch_path or VSEARCH_PATH
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_path = Path(input_fasta)
    derep = out_dir / "derep.fasta"
    repseqs = out_dir / "otus.fasta"
    uc_file = out_dir / "clusters.uc"
    logs = {}

    try:
        derep_run = subprocess.run(
            [vpath, "--derep_fulllength", str(input_path), "--sizeout", "--output", str(derep)],
            capture_output=True,
            text=True,
            check=True,
        )
        logs["derep"] = (derep_run.stdout or "") + (derep_run.stderr or "")
        cluster_run = subprocess.run(
            [vpath, "--cluster_size", str(derep), "--id", str(identity), "--centroids", str(repseqs), "--uc", str(uc_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        logs["cluster"] = (cluster_run.stdout or "") + (cluster_run.stderr or "")
    except Exception as exc:
        return {"error": str(exc), "repseqs": None, "uc_file": None, "logs": logs}

    if not repseqs.exists() or repseqs.stat().st_size == 0:
        return {"error": "Clustering produced no OTUs.", "repseqs": None, "uc_file": None, "logs": logs}

    return {
        "error": None,
        "repseqs": str(repseqs),
        "uc_file": str(uc_file) if uc_file.exists() else None,
        "logs": logs,
    }

