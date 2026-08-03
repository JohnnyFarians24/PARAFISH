"""Pre-run resource checks.

The point of this module is to turn silent out-of-memory kills into sentences the
user can act on. VSEARCH's SINTAX classifier and DADA2's ``assignTaxonomy()`` both
load an entire reference database into RAM. When the reference does not fit, the
Linux OOM killer sends SIGKILL partway through: the stage dies with no error of its
own, and — because the killer picks the largest resident process — it sometimes
takes the Dash server down with it, which is what surfaces in the browser as
"Callback failed: the server did not respond".

Predicting that before the run costs nothing and is the difference between "SILVA
failed" and "SILVA needs ~2.4 GB and you have 1.0 GB".
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Optional

# VSEARCH holds the reference sequences plus a k-mer index. The index is the larger
# term, and 3x the on-disk FASTA is the conservative rule of thumb we plan against;
# under-estimating here is what gets a run killed.
SINTAX_MEMORY_MULTIPLIER = 3.0

# Rough peak RSS of the pipeline itself (Dash server, pandas frames, the pooled
# FASTA) before any reference database is loaded.
BASE_PIPELINE_MEMORY_BYTES = 600 * 1024 * 1024

# Measured floor for the DADA2 branch: an Rscript that merely loads dada2 and
# Biostrings and then exits already holds ~747 MB resident. That is the cost of the
# libraries themselves, not of the data — a two-sample, 100k-read run peaks barely
# above it. So no amount of tuning inside the R script buys headroom below this, and
# there is no point starting DADA2 without it: R gets SIGKILLed during library load.
DADA2_BASELINE_MEMORY_BYTES = 850 * 1024 * 1024

# Per-sample working memory: the trimmed FASTQ, its filtered FASTA, and the reads
# that sample contributes to the pooled dereplication.
PER_SAMPLE_MEMORY_BYTES = 120 * 1024 * 1024

# Above this, the pooled VSEARCH clustering and the DADA2 error model stop being
# the bottleneck and sheer bookkeeping does. It is a guard-rail, not a hard law of
# nature: raise it if the machine is big and the run genuinely needs to be bigger.
ABSOLUTE_MAX_SAMPLES = 96


def available_memory_bytes() -> Optional[int]:
    """Best-effort available RAM in bytes (Linux ``/proc/meminfo``)."""
    try:
        with open("/proc/meminfo", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        return None
    return None


def total_memory_bytes() -> Optional[int]:
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def _gb(value: float) -> str:
    return f"{value / 1e9:.1f} GB"


def sintax_memory_required(db_path: Path) -> Optional[int]:
    """Bytes of RAM a SINTAX run against ``db_path`` is expected to need."""
    try:
        return int(Path(db_path).stat().st_size * SINTAX_MEMORY_MULTIPLIER)
    except OSError:
        return None


def max_samples_for_memory(available: Optional[int] = None) -> int:
    """How many samples this machine can carry through one run.

    Derived from the memory left over once the pipeline's own footprint is paid
    for, then capped. Never returns less than 1 — a single sample must always be
    attemptable, even if it is going to be tight.
    """
    available = available_memory_bytes() if available is None else available
    if not available:
        return ABSOLUTE_MAX_SAMPLES
    spare = available - BASE_PIPELINE_MEMORY_BYTES
    if spare <= 0:
        return 1
    return max(1, min(ABSOLUTE_MAX_SAMPLES, int(spare // PER_SAMPLE_MEMORY_BYTES)))


def dada2_memory_shortfall(available: Optional[int] = None) -> Optional[str]:
    """Why DADA2 cannot run right now, or ``None`` if it can.

    Checked against memory free *at this moment*, not total: on a small box the thing
    that decides whether DADA2 survives is usually what else is running (an IDE and its
    extensions can easily hold 1.5 GB), not the size of the FASTQs.
    """
    available = available_memory_bytes() if available is None else available
    if available is None or available >= DADA2_BASELINE_MEMORY_BYTES:
        return None
    return (
        f"DADA2 (ASVs) needs about {_gb(DADA2_BASELINE_MEMORY_BYTES)} free just to load "
        f"the R libraries — and right now there is {_gb(available)}. This does not depend on sample "
        f"size: it is the cost of the package itself. Close whatever is using memory (VS Code "
        f"and its extensions are usually the biggest consumer) or raise the WSL2 RAM. "
        f"OTUs (VSEARCH) are unaffected and will still run."
    )


def preflight_notes(
    db_sources: Iterable[str],
    db_paths,
    n_samples: int,
    analysis_method: str = "",
    marker_filter: str = "",
    inapplicable: Iterable[str] = (),
) -> list[str]:
    """Warnings to emit before a run. Empty when all is well.

    ``marker_filter`` is the marker (e.g. ``"18S"``) that SILVA/PR2 will be cut down to
    before they are loaded — see ``databases.reference_filter``. It is passed in rather
    than looked up so this module stays free of database imports. When it is set, the
    size that matters is the filtered one, and warning about the full file would send
    the user off to fix a problem that no longer exists.

    ``inapplicable`` names sources this run's marker rules out entirely. They are
    announced and then dropped: promising a COI run that "SILVA vai ser reduzida ao teu
    marcador" describes work the run will not do, and it was the first visible symptom
    that the marker had been resolved wrongly.
    """
    notes: list[str] = []
    available = available_memory_bytes()
    total = total_memory_bytes()
    sources = {str(source).lower() for source in db_sources or ()}
    ruled_out = {str(source).lower() for source in inapplicable or ()} & sources
    for source in sorted(ruled_out):
        notes.append(
            f"{source.upper()}: does not index this run's marker — it will be marked "
            f"not applicable and skipped, without being filtered or loaded."
        )
    sources -= ruled_out

    if available is None:
        return notes

    if str(analysis_method).lower() in ("dada2", "both"):
        shortfall = dada2_memory_shortfall(available)
        if shortfall:
            notes.append(shortfall)

    # BOLD and NCBI query remote servers, one feature at a time. They are correct and
    # they are slow — a run with many OTUs spends most of its wall-clock time here, and
    # that is the etapa Databases "not finishing" when in fact it is still waiting on the
    # network. Local sources (SILVA/PR2) are fast but do not fit in this machine's RAM,
    # so on a small box the only working taxonomy is unavoidably the slow remote one.
    remote = sorted({"BOLD", "NCBI"} & {s.upper() for s in sources})
    if remote:
        notes.append(
            f"{' and '.join(remote)} query remote servers, one feature at a time — by far the slowest "
            f"stage (minutes to over an hour, depending on the number of OTUs/ASVs). The Databases "
            f"stage can look «stuck» when it is only waiting on the network. For a quick run, "
            f"deselect these sources (you get no taxonomy, but QC→Analysis finishes in seconds)."
        )
    if str(analysis_method).lower() == "both" and remote:
        notes.append(
            "Method «both»: since you generate both OTUs and ASVs, the remote sources run TWICE (once per "
            "set), which doubles the time. Choose «vsearch» or «dada2» if you only need one."
        )

    if total and total < 4 * 1024 * 1024 * 1024:
        notes.append(
            f"Total machine memory: {_gb(total)} (available now: {_gb(available)}). "
            f"That is little for large reference databases. On WSL2, raise it with a "
            f".wslconfig file on Windows: [wsl2] / memory=8GB, then 'wsl --shutdown'."
        )

    for name, attr in (("SILVA", "silva"), ("PR2", "pr2")):
        if name.lower() not in sources:
            continue
        db_path = getattr(db_paths, attr, None)
        if not db_path or not Path(db_path).exists():
            continue
        needed = sintax_memory_required(Path(db_path))
        if not needed or needed <= available:
            continue
        if marker_filter:
            notes.append(
                f"{name}: the full reference takes {_gb(Path(db_path).stat().st_size)} and would not fit "
                f"in the {_gb(available)} available, so it will be cut down to your marker "
                f"({marker_filter}) before being loaded — it becomes a fraction of the size and runs "
                f"locally. The first run with this reference takes a few minutes to prepare; "
                f"later ones reuse the filtered version."
            )
        else:
            notes.append(
                f"{name}: the reference takes {_gb(Path(db_path).stat().st_size)} and VSEARCH SINTAX "
                f"needs about {_gb(needed)} of RAM, but only {_gb(available)} is available. "
                f"The source will be skipped rather than killed mid-run by the system (OOM). "
                f"Fixes: raise the RAM, or set the run's marker and primers so the reference "
                f"is trimmed automatically."
            )

    limit = max_samples_for_memory(available)
    if n_samples > limit:
        notes.append(
            f"{n_samples} samples in this run, but with {_gb(available)} available the safe limit "
            f"is around {limit}. Split the run into smaller batches or raise the RAM."
        )

    return notes
