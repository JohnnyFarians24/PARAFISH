"""cutadapt PCR-primer removal for the Preprocess stage.

Trimmomatic's ILLUMINACLIP removes *sequencing adapters*, not the amplification
primers. Left in, the primer — degenerate bases and all — stays in every OTU/ASV and
is fed to the classifier and the tree as if it were biological signal. This step strips
the forward primer from R1 and the reverse primer from R2 on the *raw* reads, before
quality trimming and merging, using cutadapt's native IUPAC wildcard matching.

It is best-effort by design: with no primers configured, or with cutadapt missing, it
returns the *original* reads untouched (``metrics['primers_removed'] == 0.0``) so the
stage never fails on primers. In that fallback the older fixed-length strip in the
Analysis stage still applies, so primers are not silently kept when cutadapt is absent.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

from core.config import CUTADAPT_PATH, DISCARD_UNTRIMMED
from core.models import AnalysisMode, RunConfig, Sample, StepOutput

try:  # the preset table lives with the length-strip that consumes it too
    from analysis.vsearch_runner import PRIMER_PRESETS
except Exception:  # pragma: no cover - defensive
    PRIMER_PRESETS = {}


# The standard primer pair for a marker, used only when the run picks "auto".
#
# Guessing is safe *here* and nowhere else in this pipeline: cutadapt matches the primer
# as a sequence and removes it only where it actually occurs, so a wrong guess removes
# nothing. The blind length-strip in the Analysis stage must never guess, because it
# cuts N bases whether or not the primer is there.
#
# COI is deliberately absent: Folmer and Leray are both standard and amplify different
# regions, so there is no single pair to fall back on.
MARKER_DEFAULT_PRESET = {"18S": "18S_V4_TAReuk"}

# The preset value meaning "work it out from the marker". Distinct from "" — which is
# the user explicitly choosing not to strip primers, and must be obeyed.
AUTO_PRESET = "auto"


def _marker_preset(config: RunConfig) -> str:
    """The preset for the marker this run targets, or ``""`` if it is mixed/unknown."""
    markers = {str(getattr(sample.marker, "value", sample.marker) or "").upper() for sample in (getattr(config, "samples", None) or [])}
    markers.discard("")
    if len(markers) != 1:
        return ""
    return MARKER_DEFAULT_PRESET.get(markers.pop(), "")


def resolve_primers(config: RunConfig) -> Tuple[str, str]:
    """``(forward, reverse)`` primer sequences from explicit fields, a preset, or the marker.

    The preset used to default to ``""``, which is also the value of the "do not strip
    primers" option — so a run that simply never touched the control was indistinguishable
    from one that opted out, and got no primer removal. Every OTU then carried the primer
    and the heterogeneous 5' spacer in front of it, which is what split one cod across
    OTU1 and OTU2. ``"auto"`` is now the default and resolves from the marker; ``""``
    stays a real opt-out.
    """
    params = getattr(config, "vsearch_params", {}) or {}
    forward = str(params.get("primer_fwd") or "").strip()
    reverse = str(params.get("primer_rev") or "").strip()
    if forward or reverse:
        return forward, reverse

    raw_preset = params.get("primer_preset")
    # Absent (an older stored run) is treated as "auto"; an explicit "" is not.
    preset = AUTO_PRESET if raw_preset is None else str(raw_preset).strip()
    if preset == AUTO_PRESET:
        preset = _marker_preset(config)
    if preset in PRIMER_PRESETS:
        return PRIMER_PRESETS[preset]
    return "", ""


def cutadapt_available() -> bool:
    return bool(CUTADAPT_PATH and (shutil.which(CUTADAPT_PATH) or Path(CUTADAPT_PATH).is_file()))


def _counts_from_json(report_path: Path) -> Dict[str, float]:
    """Pull input/output/with-primer read counts from cutadapt's --json report."""
    try:
        data = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    read_counts = data.get("read_counts") or {}
    result: Dict[str, float] = {}
    for out_key, json_key in (("input_reads", "input"), ("output_reads", "output"), ("reads_with_primer", "read1_with_adapter")):
        value = read_counts.get(json_key)
        if isinstance(value, (int, float)):
            result[out_key] = float(value)
    return result


def _passthrough(sample: Sample, message: str = "") -> StepOutput:
    """Success carrying the ORIGINAL reads, so preprocess continues unchanged."""
    return StepOutput(
        step_name="cutadapt",
        success=True,
        files={"r1_primer": sample.paths.r1, "r2_primer": sample.paths.r2},
        metrics={"primers_removed": 0.0},
        error_message=message,
    )


def run_cutadapt(sample: Sample, config: RunConfig) -> StepOutput:
    """Remove PCR primers from one sample's raw reads before quality trimming.

    ``files['r1_primer']`` / ``files['r2_primer']`` always point at the reads to feed
    downstream — the primer-trimmed outputs when they were produced, else the originals.
    ``metrics['primers_removed']`` is ``1.0`` only when cutadapt actually ran.
    """
    forward, reverse = resolve_primers(config)
    if not (forward or reverse):
        return _passthrough(sample)
    if not cutadapt_available():
        return _passthrough(
            sample,
            "cutadapt não encontrado — primers não removidos aqui; o corte por comprimento "
            "no vsearch (Analysis) continua a aplicar-se. Instala com 'pip install cutadapt'.",
        )

    output_dir = config.output_dir / "preprocess" / sample.internal_id
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "cutadapt.log"
    report_path = output_dir / "cutadapt.json"
    threads = int((getattr(config, "vsearch_params", {}) or {}).get("threads") or 4)

    command = [CUTADAPT_PATH, "-j", str(threads), "--json", str(report_path)]
    if DISCARD_UNTRIMMED:
        # A read with no primer match cannot be shown to be an amplicon of the targeted
        # region, so it is dropped rather than carried through uncut.
        command += ["--discard-untrimmed"]
    is_pe = sample.mode == AnalysisMode.PE and sample.paths.r2 is not None
    if is_pe:
        r1_out = output_dir / f"{sample.internal_id}_R1_primer.fastq.gz"
        r2_out = output_dir / f"{sample.internal_id}_R2_primer.fastq.gz"
        if forward:
            command += ["-g", forward]
        if reverse:
            command += ["-G", reverse]
        if DISCARD_UNTRIMMED:
            # Paired-end only — cutadapt rejects --pair-filter outright in single-end
            # mode, which made the whole step fall back to untrimmed passthrough.
            # "any" drops the pair unless the amplicon is verified at BOTH ends, which
            # is the claim --discard-untrimmed exists to let the study make.
            command += ["--pair-filter=any"]
        command += ["-o", str(r1_out), "-p", str(r2_out), str(sample.paths.r1), str(sample.paths.r2)]
        files: Dict[str, Optional[Path]] = {"r1_primer": r1_out, "r2_primer": r2_out, "log": log_path, "report": report_path}
    else:
        r1_out = output_dir / f"{sample.internal_id}_primer.fastq.gz"
        if forward:
            command += ["-g", forward]
        command += ["-o", str(r1_out), str(sample.paths.r1)]
        files = {"r1_primer": r1_out, "r2_primer": None, "log": log_path, "report": report_path}

    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
        log_path.write_text((completed.stdout or "") + (completed.stderr or ""), encoding="utf-8")
    except subprocess.CalledProcessError as exc:
        log_path.write_text((exc.stdout or "") + (exc.stderr or ""), encoding="utf-8")
        # Fall back to the untrimmed reads rather than killing the whole stage.
        return _passthrough(sample, f"cutadapt falhou ({str(exc)[:200]}) — reads passam sem remoção de primers.")

    metrics: Dict[str, float] = {"primers_removed": 1.0}
    metrics.update(_counts_from_json(report_path))
    return StepOutput(step_name="cutadapt", success=True, files=files, log=log_path, metrics=metrics)
