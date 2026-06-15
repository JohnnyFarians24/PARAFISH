"""Manifest construction, sample pairing and workflow selection helpers."""

from __future__ import annotations

import os
import re

from core.models import Sample
from core.utils import _normalize_uploaded_paths
from .sample_parsing import ParsedFastqName, parse_fastq_name, validate_manifest_rows


def _next_sample_id(existing_ids):
    i = 1
    while True:
        sid = f"S{str(i).zfill(3)}"
        if sid not in existing_ids:
            return sid
        i += 1


def _normalize_sample_name(filename):
    parsed = parse_fastq_name(filename)
    return parsed.internal_sample_id


def _detect_pe_pair_suffix(filename):
    lowered = os.path.basename(str(filename or "")).lower()
    patterns = [
        (r"_1(?:\.\w+)?$", "1"),
        (r"_2(?:\.\w+)?$", "2"),
        (r"\.r1(?:\.\w+)?$", "r1"),
        (r"\.r2(?:\.\w+)?$", "r2"),
        (r"_r1(?:\.\w+)?$", "r1"),
        (r"_r2(?:\.\w+)?$", "r2"),
        (r"\.raw_1(?:\.\w+)?$", "1"),
        (r"\.raw_2(?:\.\w+)?$", "2"),
        (r"_raw_1(?:\.\w+)?$", "1"),
        (r"_raw_2(?:\.\w+)?$", "2"),
    ]
    for pattern, suffix_type in patterns:
        if re.search(pattern, lowered):
            return suffix_type
    return None


def _is_paired_end_pair(file1, file2):
    if not file1 or not file2:
        return False
    return _normalize_sample_name(file1) == _normalize_sample_name(file2)


def _group_pe_pairs(file_list):
    if not file_list:
        return []
    processed = set()
    groups = []
    for i, file1 in enumerate(file_list):
        if i in processed:
            continue
        matched = False
        for j in range(i + 1, len(file_list)):
            if j in processed:
                continue
            file2 = file_list[j]
            if _is_paired_end_pair(file1, file2):
                groups.append(("PE", [file1, file2]))
                processed.add(i)
                processed.add(j)
                matched = True
                break
        if not matched:
            groups.append(("SE", [file1]))
            processed.add(i)
    return groups


def _ensure_pe_sample_paths(sample_info):
    if not sample_info or not isinstance(sample_info, dict):
        return {"kind": "SE", "paths": {}}
    info = sample_info.copy()
    paths = info.get("paths", {}) or {}
    kind = info.get("kind", "SE")
    if kind == "PE" and (not paths.get("r1") or not paths.get("r2")):
        if paths.get("se"):
            info["kind"] = "SE"
    return info


def _get_sample_fastq_pair(sample_id, samples_store):
    samples = samples_store or {}
    sample = samples.get(sample_id, {})
    info = _ensure_pe_sample_paths(sample)
    paths = info.get("paths", {})
    kind = info.get("kind", "SE")
    if kind == "PE" and paths.get("r1") and paths.get("r2"):
        return {"kind": "PE", "files": [paths["r1"], paths["r2"]]}
    if paths.get("se"):
        return {"kind": "SE", "files": [paths["se"]]}
    if paths.get("r1"):
        return {"kind": "SE", "files": [paths["r1"]]}
    if paths.get("r2"):
        return {"kind": "SE", "files": [paths["r2"]]}
    return {"kind": "SE", "files": []}


def _normalize_samples_store(samples_store):
    if not samples_store or not isinstance(samples_store, dict):
        return {}

    merged = {}
    seen_pair = {}

    def _pair_key_from_info(info):
        name = str((info or {}).get("pair_key") or (info or {}).get("display_name") or (info or {}).get("stem") or "")
        normalized = re.sub(r"(?i)(?:[_\.-]?(?:raw[_\.-]?)?(?:r|fwd|rev)?[ _\.-]?(?:1|2))$", "", os.path.splitext(name)[0])
        normalized = re.sub(r"[^A-Za-z0-9]+", "_", normalized).strip("_")
        return normalized

    for sid, info in list(samples_store.items()):
        info = info if isinstance(info, dict) else {}
        pair_key = info.get("pair_key") or _pair_key_from_info(info) or sid

        if pair_key in seen_pair:
            target_sid = seen_pair[pair_key]
            target = merged[target_sid]

            members = info.get("members") if isinstance(info.get("members"), list) else []
            if not members:
                members = [info.get("display_name") or info.get("display_label") or info.get("stem") or sid]
            target.setdefault("members", [])
            for member in members:
                if member and member not in target["members"]:
                    target["members"].append(member)

            for path_key, path_value in (info.get("paths") or {}).items():
                if path_value:
                    target.setdefault("paths", {})[path_key] = path_value

            for output_key, output_value in (info.get("outputs") or {}).items():
                if output_value:
                    target.setdefault("outputs", {})[output_key] = output_value

            for meta_key, meta_value in (info.get("meta") or {}).items():
                if meta_value is not None:
                    target.setdefault("meta", {})[meta_key] = meta_value

            if target.get("kind") != "PE" and info.get("kind") == "PE":
                target["kind"] = "PE"
                target["mode"] = "PE"

            if not target.get("internal_sample_id"):
                target["internal_sample_id"] = info.get("internal_sample_id") or pair_key
            if not target.get("display_label"):
                target["display_label"] = target.get("display_name") or target.get("pair_key") or target_sid
            if not target.get("display_name"):
                target["display_name"] = target.get("display_label") or target_sid
            if not target.get("marker") and info.get("marker"):
                target["marker"] = info.get("marker")
        else:
            display_name = info.get("display_name") or info.get("display_label") or info.get("stem") or sid
            merged[sid] = {
                "display_name": display_name,
                "display_label": info.get("display_label") or display_name,
                "stem": info.get("stem") or sid,
                "pair_key": pair_key,
                "internal_sample_id": info.get("internal_sample_id") or pair_key,
                "kind": info.get("kind", "SE"),
                "mode": info.get("mode", info.get("kind", "SE")),
                "marker": info.get("marker"),
                "members": list(info.get("members", [display_name])) if info.get("members") else [display_name],
                "paths": dict(info.get("paths", {}) or {}),
                "outputs": dict(info.get("outputs", {}) or {}),
                "meta": dict(info.get("meta", {}) or {}),
            }
            seen_pair[pair_key] = sid

    return merged


def _index_samples_from_dir(dataset_dir):
    candidates = []
    for root, _, files in os.walk(dataset_dir):
        for filename in files:
            if filename.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fasta", ".fa", ".fasta.gz", ".fa.gz")):
                candidates.append(os.path.join(root, filename))
    return assign_samples_from_files(sorted(candidates), current_samples={})


def _discover_dada2_preprocess_pairs(project_state):
    from core.paths import _step_base

    preprocess_root = _step_base(project_state, "preprocess")
    if not os.path.exists(preprocess_root):
        return []
    pairs = {}
    for root, _, files in os.walk(preprocess_root):
        for filename in files:
            if not filename.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                continue
            path = os.path.join(root, filename)
            sample_key = _normalize_sample_name(filename)
            pairs.setdefault(sample_key, []).append(path)
    resolved = []
    for key, items in sorted(pairs.items()):
        if len(items) >= 2:
            resolved.append((key, sorted(items)[:2]))
    return resolved


def _resolve_dada2_inputs(project_state, uploaded_files, samples_store, workflow_scope):
    selected_samples = list(_iter_workflow_samples(samples_store, workflow_scope))
    if not selected_samples:
        return {"mode": None, "fwd_files": [], "rev_files": [], "single_files": [], "source": None, "error": "No selected samples found."}

    pe_entries = []
    se_entries = []
    for sid, sample in selected_samples:
        info = _ensure_pe_sample_paths(sample)
        paths = info.get("paths", {})
        if paths.get("r1") and paths.get("r2"):
            pe_entries.append((sid, paths.get("r1"), paths.get("r2")))
        elif paths.get("se"):
            se_entries.append((sid, paths.get("se")))

    if pe_entries and se_entries:
        return {"mode": None, "fwd_files": [], "rev_files": [], "single_files": [], "source": None, "error": "Cannot mix PE and SE samples in the same DADA2 run."}
    if pe_entries:
        return {
            "mode": "PE",
            "fwd_files": [entry[1] for entry in pe_entries],
            "rev_files": [entry[2] for entry in pe_entries],
            "single_files": [],
            "source": "manifest",
            "error": None,
        }
    if se_entries:
        return {
            "mode": "SE",
            "fwd_files": [],
            "rev_files": [],
            "single_files": [entry[1] for entry in se_entries],
            "source": "manifest",
            "error": None,
        }
    return {"mode": None, "fwd_files": [], "rev_files": [], "single_files": [], "source": None, "error": "No valid DADA2 inputs found."}


def assign_samples_from_files(file_list, current_samples=None, mode_hint=None):
    samples = (current_samples or {}).copy()
    existing_ids = set(samples.keys())

    def _display_label(internal_sample_id, marker, mode):
        label = str(internal_sample_id or "sample")
        if marker:
            label = f"{label} ({marker})"
        return f"{label} ({mode})"

    for fp in file_list or []:
        if not fp:
            continue
        parsed = parse_fastq_name(fp)
        base = os.path.basename(fp)
        pair_key = parsed.internal_sample_id
        mode = str(mode_hint or parsed.mode or "SE").upper()

        sid = None
        for k, v in samples.items():
            if v.get("pair_key") == pair_key:
                sid = k
                break

        if not sid:
            sid = _next_sample_id(existing_ids)
            existing_ids.add(sid)
            samples[sid] = {
                "display_name": _display_label(pair_key, parsed.marker, mode),
                "display_label": _display_label(pair_key, parsed.marker, mode),
                "stem": os.path.splitext(base)[0],
                "pair_key": pair_key,
                "internal_sample_id": pair_key,
                "kind": mode,
                "mode": mode,
                "marker": parsed.marker,
                "members": [parsed.display_label],
                "paths": {},
                "source_files": [fp],
                "outputs": {},
                "meta": {},
            }
        else:
            if parsed.display_label not in samples[sid].get("members", []):
                samples[sid].setdefault("members", []).append(parsed.display_label)
                samples[sid]["display_name"] = _display_label(pair_key, samples[sid].get("marker") or parsed.marker, samples[sid].get("mode", mode))
                samples[sid]["display_label"] = samples[sid]["display_name"]
            samples[sid].setdefault("source_files", [])
            if fp not in samples[sid]["source_files"]:
                samples[sid]["source_files"].append(fp)
            if not samples[sid].get("internal_sample_id"):
                samples[sid]["internal_sample_id"] = pair_key
            if not samples[sid].get("marker") and parsed.marker:
                samples[sid]["marker"] = parsed.marker

        if parsed.mode == "PE" and parsed.mate == 1:
            samples[sid]["paths"]["r1"] = fp
            samples[sid]["kind"] = "PE"
            samples[sid]["mode"] = "PE"
        elif parsed.mode == "PE" and parsed.mate == 2:
            samples[sid]["paths"]["r2"] = fp
            samples[sid]["kind"] = "PE"
            samples[sid]["mode"] = "PE"
        else:
            samples[sid]["paths"]["se"] = fp
            if not samples[sid]["paths"].get("r1") and not samples[sid]["paths"].get("r2"):
                samples[sid]["kind"] = "SE"
                samples[sid]["mode"] = "SE"

        if "display_name" not in samples[sid] or not samples[sid]["display_name"]:
            marker = samples[sid].get("marker")
            samples[sid]["display_name"] = _display_label(samples[sid].get("internal_sample_id") or pair_key, marker, samples[sid].get("mode", "SE"))
            samples[sid]["display_label"] = samples[sid]["display_name"]

    return samples


def _iter_workflow_samples(samples_store, workflow_scope):
    scope = workflow_scope or {}
    selected_ids = scope.get("samples_selected", [])
    samples = samples_store or {}
    if not selected_ids or "__ALL__" in selected_ids:
        return list(samples.items())
    result = []
    for sid in selected_ids:
        if sid in samples:
            result.append((sid, samples[sid]))
    return result


def _get_workflow_scope_samples(workflow_scope):
    scope = workflow_scope or {}
    samples = scope.get("samples_selected", [])
    if not samples or "__ALL__" in samples:
        return []
    return samples


def _normalize_samples_for_manifest(samples_store, workflow_scope=None):
    rows = []
    for sid, sample in _iter_workflow_samples(samples_store, workflow_scope):
        sample = sample or {}
        paths = sample.get("paths", {}) if isinstance(sample, dict) else {}
        rows.append(
            {
                "internal_sample_id": sample.get("internal_sample_id") or sample.get("pair_key") or sid,
                "display_label": sample.get("display_label") or sample.get("display_name") or sid,
                "marker": sample.get("marker"),
                "mode": sample.get("mode") or sample.get("kind") or "SE",
                "forward": paths.get("r1"),
                "reverse": paths.get("r2"),
                "single": paths.get("se"),
            }
        )
    return rows


def build_manifest(samples: list[ParsedFastqName]) -> list[dict[str, str]]:
    return [
        {
            "display_label": sample.display_label,
            "internal_sample_id": sample.internal_sample_id,
            "mode": sample.mode,
            "marker": sample.marker,
        }
        for sample in samples
    ]


def validate_manifest(manifest: list[dict[str, str]]) -> bool:
    required = {"display_label", "internal_sample_id", "mode", "marker"}
    return all(isinstance(record, dict) and required.issubset(record.keys()) for record in manifest)
