"""Build structured run ZIP for download."""

from __future__ import annotations

import hashlib
import io
import json
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path

STAGE_ORDER = ["qc", "preprocess", "analysis", "databases", "phylo"]

# A real path is never this long. Anything longer is a message that leaked into a
# position where a path was expected, and `Path(...).exists()` raises ENAMETOOLONG on
# it rather than returning False — which is how a stage's error_message ended up
# crashing the whole ZIP download with "OSError: File name too long".
_MAX_PATH_LEN = 4096

# Keys of a serialized StepOutput that are metadata, not artifacts.
_STEP_OUTPUT_META = frozenset({"step_name", "success", "log", "metrics", "error_message"})


def _is_path_like(value) -> bool:
    return isinstance(value, str) and bool(value) and len(value) < _MAX_PATH_LEN


def _exists(path: Path) -> bool:
    """``path.exists()`` that cannot raise.

    ``Path.exists()`` only swallows a specific set of errnos; ENAMETOOLONG is not among
    them, so a stray non-path string propagates an OSError to the caller.
    """
    try:
        return path.exists()
    except OSError:
        return False


def _write(zf, seen: set, path: Path, parts: list) -> None:
    """Write one file under ``parts``, skipping anything already in the archive.

    The same artifact is reachable by several routes — ``otus.fasta`` is the pooling
    step's output, the legacy ``combined_fasta`` alias, AND the OTU feature set's
    representative FASTA — so without this the ZIP carried it three times and zipfile
    warned about duplicate names.
    """
    arcname = "/".join(parts)
    if arcname in seen:
        return
    seen.add(arcname)
    zf.write(str(path), arcname)


def _add_stage_files(zf, seen, run_id, stage, sample_id, files, subdir: str = "") -> None:
    """Add every real file in ``files`` (a name -> path map) to the archive."""
    if not isinstance(files, dict):
        return
    for key, path_str in files.items():
        if key in _STEP_OUTPUT_META or not _is_path_like(path_str):
            continue
        path = Path(path_str)
        if not _exists(path) or not path.is_file():
            continue
        parts = [run_id, stage]
        if subdir:
            parts.append(subdir)
        if sample_id:
            parts.append(sample_id)
        parts.append(path.name)
        _write(zf, seen, path, parts)


def _step_output_files(value) -> dict:
    """The ``files`` map of a serialized StepOutput, or ``{}``.

    Passing the whole StepOutput where a files map was expected is what made the
    exporter treat ``error_message`` as a filename.
    """
    if not isinstance(value, dict):
        return {}
    files = value.get("files")
    return files if isinstance(files, dict) else {}


def _add_metrics(zf, seen: set, run_id: str, stage: str, metrics_files: dict) -> None:
    """Add the metric tables, one folder per feature set.

    ``metrics_files`` is keyed ``"otu.alpha_diversity_csv"`` / ``"asv.alpha_diversity_csv"``
    with unprefixed aliases pointing at the primary set's files. Writing the aliases too
    would put two identical entries in the archive, so only the namespaced keys are
    used when any exist — that also keeps the OTU and ASV tables in separate folders
    instead of one silently overwriting the other, since they share basenames.
    """
    if not isinstance(metrics_files, dict):
        return
    namespaced = {key: value for key, value in metrics_files.items() if "." in key}
    entries = namespaced or metrics_files
    for key, path_str in entries.items():
        if not _is_path_like(path_str):
            continue
        path = Path(path_str)
        if not _exists(path) or not path.is_file():
            continue
        feature_set = key.split(".", 1)[0] if "." in key else ""
        parts = [run_id, stage, "metrics"]
        if feature_set:
            parts.append(feature_set)
        parts.append(path.name)
        _write(zf, seen, path, parts)


# ── Provenance recorded when a run starts ────────────────────────────────
#
# These feed the run_config block of progress.json (ui.callbacks_progress._run_provenance).
# None of them may raise: a run must never fail because its own description could not be
# completed. Each failure has a word of its own instead — "unknown", "unavailable".

_HASH_CHUNK = 1024 * 1024


def git_commit(repo_root=None) -> str:
    """The commit of the code that is running, or ``"unknown"`` (no git, no repository)."""
    root = Path(repo_root) if repo_root else Path(__file__).resolve().parent.parent
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(root), capture_output=True, text=True, timeout=10, check=True
        )
    except Exception:
        return "unknown"
    commit = (completed.stdout or "").strip()
    return commit or "unknown"


def _fingerprint_cache_path(path: Path) -> Path:
    return path.parent / f"{path.name}.sha256.json"


def _read_fingerprint_cache(path: Path, size: int, mtime_ns: int):
    """The cached sha256 for this exact size+mtime, or ``None``."""
    try:
        cached = json.loads(_fingerprint_cache_path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(cached, dict):
        return None
    if cached.get("size_bytes") != size or cached.get("mtime_ns") != mtime_ns:
        return None
    digest = cached.get("sha256")
    return digest if isinstance(digest, str) and digest else None


def _write_fingerprint_cache(path: Path, size: int, mtime_ns: int, digest: str) -> None:
    _fingerprint_cache_path(path).write_text(
        json.dumps({"file": path.name, "size_bytes": size, "mtime_ns": mtime_ns, "sha256": digest}, indent=2),
        encoding="utf-8",
    )


def reference_fingerprint(path) -> dict:
    """``{file, size_bytes, mtime, sha256}`` for one local reference file.

    A file NAME does not identify a reference — two downloads of "the same" release can
    differ — so the run records a sha256. Hashing a multi-hundred-MB FASTA costs seconds,
    so the digest is cached in ``<file>.sha256.json`` beside the reference and reused for
    as long as the size and the mtime are both unchanged. Where that directory is
    read-only (the Docker mounts are) the cache simply is not written and the next run
    hashes again; the digest itself is still returned.
    """
    path = Path(path)
    result = {"file": path.name, "size_bytes": None, "mtime": None, "sha256": "unavailable"}
    try:
        stat = path.stat()
        result["size_bytes"] = int(stat.st_size)
        result["mtime"] = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if not path.is_file():
            return result
        digest = _read_fingerprint_cache(path, stat.st_size, stat.st_mtime_ns)
        if digest is None:
            hasher = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(_HASH_CHUNK), b""):
                    hasher.update(chunk)
            digest = hasher.hexdigest()
            try:
                _write_fingerprint_cache(path, stat.st_size, stat.st_mtime_ns, digest)
            except OSError:
                pass
        result["sha256"] = digest
    except Exception:
        result["sha256"] = "unavailable"
    return result


def _flatten_parameters(value, prefix: str = "") -> list:
    """``{"tree": {"seed": 42}}`` -> ``[("tree.seed", 42)]`` — one line per setting."""
    if isinstance(value, dict):
        rows = []
        for key, inner in value.items():
            rows.extend(_flatten_parameters(inner, f"{prefix}.{key}" if prefix else str(key)))
        return rows
    return [(prefix, value)]


def _recorded_provenance_lines(config: dict) -> list:
    """The PROVENANCE sections a run recorded about itself; ``[]`` for an older run."""
    lines: list = []

    commit = config.get("parafish_commit")
    if commit:
        lines += ["", "PARAFISH commit", "-" * 46, f"  {commit}"]

    references = config.get("references")
    if isinstance(references, list) and references:
        lines += ["", "References", "-" * 46]
        for entry in references:
            if not isinstance(entry, dict):
                continue
            lines.append(f"  {str(entry.get('source', '')).upper():8s} {entry.get('file', '—')}")
            lines.append(f"           size   : {entry.get('size_bytes', '—')} bytes")
            lines.append(f"           mtime  : {entry.get('mtime', '—')}")
            lines.append(f"           sha256 : {entry.get('sha256', '—')}")

    parameters = config.get("parameters")
    if isinstance(parameters, dict) and parameters:
        lines += ["", "Parameters", "-" * 46]
        for section, values in parameters.items():
            lines.append(f"  [{section}]")
            for key, value in _flatten_parameters(values):
                lines.append(f"    {key} = {value}" if key else f"    {value}")
    return lines


def render_provenance(config: dict, tool_versions=None) -> str:
    """A human-readable recipe for a run, RENDERED from the run_config block progress.json
    already carries (``_run_provenance``) — a legible view of the single source, never a
    second collection that could drift from it.

    Machine-readable provenance already lives in progress.json; the value of THIS file is
    being person-readable, so a supervisor who opens the ZIP reads "Fontes: SILVA, PR2,
    NCBI · Marcador: COI" rather than a JSON dump. It ships INSIDE the deliverable ZIP so
    the recipe travels with the results instead of staying in the app.

    ``tool_versions`` is a list of ``(tool, version)`` — reproducibility is which VSEARCH /
    BLAST built the run, not only which parameters. They are resolved from the pipeline's
    configured binaries at export time; when none can be resolved the file says so plainly
    rather than looking complete without them.
    """
    config = config or {}

    if not config:
        recipe = ["  (the recipe was not recorded for this run — it predates provenance", "   recording; new runs write it into progress.json.)"]
    else:
        sources = ", ".join(str(s).upper() for s in (config.get("db_sources") or [])) or "—"
        markers = ", ".join(str(m) for m in (config.get("markers") or [])) or "—"
        primers = config.get("primers") or {}
        recipe = [
            f"  Reference sources : {sources}",
            f"  Marker(s)         : {markers}",
            f"  Analysis method   : {config.get('analysis_method', '—')}",
            f"  Branch            : {config.get('branch', '—')}",
            f"  min_unique_size   : {config.get('min_unique_size', '—')}",
            f"  Primers  fwd      : {primers.get('fwd') or '—'}",
            f"  Primers  rev      : {primers.get('rev') or '—'}",
            f"  NCBI              : {config.get('ncbi_mode', '—')}",
        ]

    lines = ["PARAFISH — run provenance", "=" * 46, "", "Recipe", "-" * 46, *recipe, "", "Tool versions", "-" * 46]
    # What the run recorded when it started wins over anything resolved now: the binary
    # on this machine at export time may no longer be the one that produced the data.
    recorded_versions = config.get("tool_versions")
    if isinstance(recorded_versions, dict) and recorded_versions:
        lines.append("  (recorded when the run started)")
        lines += [f"  {tool:12s} {version}" for tool, version in recorded_versions.items()]
    elif tool_versions:
        lines += [f"  {tool:12s} {version}" for tool, version in tool_versions]
    else:
        lines += [
            "  (not resolved in this environment — known gap)",
            "  Versions are read from the binaries configured at export",
            "  time; without them the recipe is incomplete for reproduction.",
        ]
    lines += _recorded_provenance_lines(config)
    return "\n".join(lines) + "\n"


def build_run_zip(run_id: str, stage_outputs: dict, provenance_text: str | None = None) -> io.BytesIO:
    buf = io.BytesIO()
    seen: set[str] = set()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        # The recipe travels WITH the results: a supervisor opening the ZIP can read how the
        # run was made, not just what it produced. Rendered from the run_config in
        # progress.json (see render_provenance), so it is a view of the single source.
        if provenance_text:
            zf.writestr(f"{run_id}/PROVENANCE.txt", provenance_text)
        for stage in STAGE_ORDER:
            record = stage_outputs.get(stage) or {}
            outputs = record.get("outputs") or {}

            # Per-feature-set results supersede the single-track aliases kept for older
            # readers; exporting both would duplicate every tree and taxonomy table.
            legacy_aliases = set()
            if "taxonomy_by_set" in outputs:
                legacy_aliases.add("taxonomy")
            if "phylo_by_set" in outputs:
                legacy_aliases.add("phylo")

            # Per-feature-set first, so each artifact is claimed under its OWN folder
            # (databases/otu/…) before any alias can claim it at the stage root.
            for key in ("taxonomy_by_set", "phylo_by_set"):
                for feature_set, step_output in (outputs.get(key) or {}).items():
                    _add_stage_files(zf, seen, run_id, stage, None, _step_output_files(step_output), subdir=str(feature_set))

            for key, value in outputs.items():
                if key in ("taxonomy_by_set", "phylo_by_set") or key in legacy_aliases:
                    continue

                if not isinstance(value, dict):
                    # Plain top-level path strings (combined_fasta, otu_matrix, ...).
                    if _is_path_like(value):
                        path = Path(value)
                        if _exists(path) and path.is_file():
                            _write(zf, seen, path, [run_id, stage, path.name])
                    continue

                if key == "feature_set_comparison":
                    _add_stage_files(zf, seen, run_id, stage, None, value, subdir="comparison")
                    continue

                if "files" in value:
                    # A serialized StepOutput, keyed either by sample or by tool name.
                    sample_id = None if key in ("dada2", "taxonomy", "vsearch_pool", "phylo") else key
                    _add_stage_files(zf, seen, run_id, stage, sample_id, _step_output_files(value))
                    continue

                # Preprocess-style: {sample_id: {trimmomatic: {files: ...}, pear: {...}}}
                for tool_data in value.values():
                    _add_stage_files(zf, seen, run_id, stage, key, _step_output_files(tool_data))

            # The persisted metric tables (rarefaction, alpha diversity, occupancy,
            # taxon x sample matrices) are the tables a reader actually wants, and they
            # live beside "outputs", not inside it.
            _add_metrics(zf, seen, run_id, stage, record.get("metrics_files") or {})

    buf.seek(0)
    return buf
