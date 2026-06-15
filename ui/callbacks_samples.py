"""Callbacks related to the sample upload stage."""

from __future__ import annotations

import base64
import io
from datetime import datetime, timezone
import re
from collections import OrderedDict
import tempfile
import os
import shutil
import subprocess
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

from dash import ALL, Input, Output, State, callback_context, html, no_update
from dash.exceptions import PreventUpdate

from core.models import AnalysisBranch
from ui.runtime_state import (
	build_run_config,
	build_sample_store_entry,
	decode_upload_content,
	project_run_dir,
	resolve_default_dada2_reference,
	serialize_run_config,
	timestamp_slug,
)
from ui.components.sample_table import render_sample_table


PAIR_PATTERNS = (
	re.compile(r"(?i)(?P<base>.*?)(?:[_\-.]?r(?P<mate>[12]))$"),
	re.compile(r"(?i)(?P<base>.*?)(?:[_\-.]?read(?P<mate>[12]))$"),
)


def _normalize_uploaded_files(contents, filenames):
	if not contents or not filenames:
		return []
	if not isinstance(contents, list):
		contents = [contents]
	if not isinstance(filenames, list):
		filenames = [filenames]
	return [{"content": content, "filename": filename} for content, filename in zip(contents, filenames)]


def _sample_key(filename: str) -> str:
	base, _ = _split_filename(filename)
	return _slug(base)


def _split_filename(filename: str) -> tuple[str, str | None]:
	name = PurePosixPath(filename).name
	stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name, flags=re.IGNORECASE)
	for pattern in PAIR_PATTERNS:
		match = pattern.match(stem)
		if match:
			base = match.group("base").rstrip("_-. ")
			return base or stem, str(match.group("mate"))
	return stem, None


def _slug(value: str) -> str:
	value = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
	return value.lower() or "sample"


def parse_sample_names(filenames: list[str]):
    from collections import defaultdict

    from analysis.sample_parsing import parse_fastq_name

    parsed = [parse_fastq_name(f) for f in filenames]

    # Agrupa por (base normalizada sem sufixos raw/R1/R2, marker)
    groups: dict[tuple, dict] = defaultdict(lambda: {"files": [], "parsed": []})

    for filename, p in zip(filenames, parsed):
        # Remove "raw" residual do base
        base = re.sub(r"(?i)[_\.-]?raw$", "", p.internal_sample_id).strip("_.-") or p.internal_sample_id
        key = (base, p.marker or "OTHER")
        groups[key]["files"].append(filename)
        groups[key]["parsed"].append(p)

    samples = []
    idx = 0
    for (base, marker), group in groups.items():
        parsed_list = group["parsed"]
        files_list = group["files"]

        # Separa SE (extendedFrags/assembled) de PE
        se_files = [(f, p) for f, p in zip(files_list, parsed_list) if p.mode == "SE"]
        pe_files = [(f, p) for f, p in zip(files_list, parsed_list) if p.mode == "PE"]

        # Cada SE é uma amostra individual
        for f, p in se_files:
            idx += 1
            internal = f"{_slug(base)}_{idx:02d}"
            samples.append({
                "internal_sample_id": internal,
                "internal_id": internal,
                "display_label": p.display_label,
                "mode": "SE",
                "kind": "SE",
                "marker": marker,
                "status": "ready",
                "selected": True,
                "paths": {"se": f},
                "r1": f,
                "r2": None,
            })

        # PE: agrupa por pares (mate 1 + mate 2)
        r1_files = [(f, p) for f, p in pe_files if p.mate == 1]
        r2_files = [(f, p) for f, p in pe_files if p.mate == 2]
        # Emparelha por ordem se contagens iguais, fallback para alfabético
        pairs = list(zip(
            sorted(r1_files, key=lambda x: x[0]),
            sorted(r2_files, key=lambda x: x[0]),
        ))
        for (f1, p1), (f2, p2) in pairs:
            idx += 1
            internal = f"{_slug(base)}_{idx:02d}"
            samples.append({
                "internal_sample_id": internal,
                "internal_id": internal,
                "display_label": p1.display_label,
                "mode": "PE",
                "kind": "PE",
                "marker": marker,
                "status": "ready",
                "selected": True,
                "paths": {"r1": f1, "r2": f2},
                "r1": f1,
                "r2": f2,
            })
        # PE sem par (ímpar): SE fallback
        if len(r1_files) != len(r2_files):
            unpaired = r1_files[len(pairs):] + r2_files[len(pairs):]
            for f, p in unpaired:
                idx += 1
                internal = f"{_slug(base)}_{idx:02d}"
                samples.append({
                    "internal_sample_id": internal,
                    "internal_id": internal,
                    "display_label": p.display_label,
                    "mode": "SE",
                    "kind": "SE",
                    "marker": marker,
                    "status": "ready",
                    "selected": True,
                    "paths": {"se": f},
                    "r1": f,
                    "r2": None,
                })

    return samples


def _write_uploaded_files(entries: list[dict], project_data: dict | None) -> tuple[Path, list[dict]]:
	run_root = project_run_dir(project_data)
	raw_dir = run_root / "input"
	raw_dir.mkdir(parents=True, exist_ok=True)
	written = []
	for entry in entries:
		filename = PurePosixPath(entry["filename"]).name
		file_path = raw_dir / filename
		content = entry["content"]
		file_bytes = decode_upload_content(content) if isinstance(content, str) else bytes(content)
		file_path.write_bytes(file_bytes)
		written.append({"filename": filename, "path": str(file_path), "size": file_path.stat().st_size})
	return run_root, written


def _extract_and_collect_fastqs(contents_list, filenames_list, temp_dir):
	"""
	Recebe ficheiros uploaded (base64). Extrai ZIPs.
	Retorna lista de (nome_ficheiro, bytes) de todos os FASTQs encontrados.
	"""
	del temp_dir
	fastq_files = []
	fastq_extensions = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

	for content, filename in zip(contents_list or [], filenames_list or []):
		_, content_b64 = content.split(",", 1)
		raw_bytes = base64.b64decode(content_b64)

		if filename.lower().endswith(".zip"):
			with zipfile.ZipFile(io.BytesIO(raw_bytes)) as zf:
				for member in zf.namelist():
					member_path = Path(member)
					if member_path.name.startswith(".") or member.endswith("/"):
						continue
					suffix = "".join(member_path.suffixes[-2:]).lower()
					if suffix in fastq_extensions or member_path.suffix.lower() in fastq_extensions:
						fastq_files.append((member_path.name, zf.read(member)))
		elif filename.lower().endswith(".rar"):
			fastq_files.extend(_extract_rar(raw_bytes))
		elif any(filename.lower().endswith(ext) for ext in fastq_extensions):
			fastq_files.append((filename, raw_bytes))

	return fastq_files


def _extract_rar(content_bytes: bytes) -> list[tuple[str, bytes]]:
	"""Extract FASTQ-like files from RAR archive bytes using system tools.

	This attempts to use `unrar` first, falling back to `bsdtar` if available.
	Returns a list of (filename, bytes) for extracted FASTQ files.
	"""
	fastq_extensions = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}
	with tempfile.TemporaryDirectory(prefix="parafish_rar_") as tmp_dir:
		# write the uploaded bytes to a temp rar file
		temp_rar = Path(tmp_dir) / "upload.rar"
		temp_rar.write_bytes(content_bytes)
		# try unrar then bsdtar
		unrar_bin = shutil.which("unrar")
		bsdtar_bin = shutil.which("bsdtar")
		extracted = False
		try:
			if unrar_bin:
				# extract quietly to tmp_dir
				subprocess.run([unrar_bin, "x", "-y", str(temp_rar), str(tmp_dir)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
				extracted = True
			elif bsdtar_bin:
				subprocess.run([bsdtar_bin, "-xvf", str(temp_rar), "-C", str(tmp_dir)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
				extracted = True
		except Exception:
			# extraction failed
			extracted = False
		if not extracted:
			return []
		# walk tmp_dir and collect fastq files
		results: list[tuple[str, bytes]] = []
		for root, _, files in os.walk(tmp_dir):
			for name in files:
				path = Path(root) / name
				# skip the temporary rar itself
				if path == temp_rar:
					continue
				suffix = "".join(path.suffixes[-2:]).lower()
				if suffix in fastq_extensions or path.suffix.lower() in fastq_extensions:
					try:
						results.append((name, path.read_bytes()))
					except Exception:
						# skip unreadable
						continue
		return results


def _decode_preview(content: str | None) -> str:
	if not content:
		return ""
	try:
		if isinstance(content, str):
			_, encoded = content.split(",", 1)
			decoded = base64.b64decode(encoded)
		else:
			decoded = bytes(content)
		return decoded[:120].decode("utf-8", errors="ignore")
	except Exception:
		return ""


def register_callbacks(app: Any) -> None:
	@app.callback(
		Output("store-samples", "data"),
		Output("store-run-config", "data", allow_duplicate=True),
		Output("samples-upload-status", "children"),
		Input("upload-fastq", "contents"),
		State("upload-fastq", "filename"),
		State("store-project", "data"),
		State("store-samples", "data"),
		prevent_initial_call=True,
	)
	def _ingest_upload(contents, filenames, project_data, existing_samples_data):
		entries = _normalize_uploaded_files(contents, filenames)
		if not entries:
			return no_update, no_update, no_update
		has_rar = any(str(entry["filename"]).lower().endswith(".rar") for entry in entries)
		fastq_entries = _extract_and_collect_fastqs(
			[entry["content"] for entry in entries],
			[entry["filename"] for entry in entries],
			temp_dir=None,
		)
		if not fastq_entries:
			message = "RAR not supported, please use ZIP" if has_rar else "No FASTQ files found in the upload."
			warning_class = "status-pill badge badge-red" if has_rar else "status-pill"
			return no_update, no_update, html.Div(message, className=warning_class)
		normalized_entries = [{"filename": filename, "content": raw_bytes} for filename, raw_bytes in fastq_entries]
		run_root, written_files = _write_uploaded_files(normalized_entries, project_data)
		filenames_only = [filename for filename, _ in fastq_entries]
		samples = parse_sample_names(filenames_only)
		for sample in samples:
			base = sample["display_label"]
			for item in written_files:
				item_base, mate = _split_filename(item["filename"])
				if item_base == base and mate == "1":
					sample["r1_path"] = item["path"]
				if item_base == base and mate == "2":
					sample["r2_path"] = item["path"]
			if "r1_path" not in sample:
				match = next((item for item in written_files if _split_filename(item["filename"])[0] == base), None)
				if match:
					sample["r1_path"] = match["path"]
					sample["r2_path"] = None
			if not sample.get("r2_path"):
				sample["r2_path"] = None
		run_config = build_run_config(project_data, {"samples": samples})
		run_config.dada2_params["ref_fasta"] = resolve_default_dada2_reference()
		if getattr(run_config, "analysis_method", "vsearch") in {"dada2", "both"} and not run_config.dada2_params.get("ref_fasta"):
			run_config.dada2_params["ref_fasta"] = None
		preview = _decode_preview(normalized_entries[0]["content"])
		status_children = [
			html.Div("Upload complete", className="status-pill"),
			html.Div(f"{len(samples)} sample group(s) detected in {run_root}", className="mini-note"),
		]
		if preview:
			status_children.append(html.Div(preview[:80], className="mini-note"))
		existing_samples = (existing_samples_data or {}).get("samples") or []
		existing_ids = {s["internal_id"] for s in existing_samples if s.get("internal_id")}
		new_samples = [s for s in samples if s.get("internal_id") not in existing_ids]
		merged_samples = existing_samples + new_samples

		existing_files = (existing_samples_data or {}).get("files") or []
		merged_files = existing_files + written_files
		for sample in merged_samples:
			r1_name = PurePosixPath(sample["r1"]).name if sample.get("r1") else None
			r2_name = PurePosixPath(sample["r2"]).name if sample.get("r2") else None
			for item in written_files:
				if r1_name and item["filename"] == r1_name:
					sample["r1_path"] = item["path"]
					sample["paths"]["r1"] = item["path"]
				if r2_name and item["filename"] == r2_name:
					sample["r2_path"] = item["path"]
					sample["paths"]["r2"] = item["path"]
			# SE fallback
			if not sample.get("r1_path"):
				se_name = PurePosixPath(sample["r1"]).name if sample.get("r1") else None
				match = next((item for item in written_files if item["filename"] == se_name), None)
				if match:
					sample["r1_path"] = match["path"]
					sample["paths"]["r1"] = match["path"]
			sample["r2_path"] = sample.get("r2_path") or None
		return (
			{
				"files": merged_files,
				"samples": merged_samples,
				"run_root": str(run_root),
				"updated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
			},
			serialize_run_config(run_config),
			status_children,
		)

	@app.callback(
		Output("samples-table", "children"),
		Input("store-samples", "data"),
	)
	def _render_samples_table(samples_data):
		samples = (samples_data or {}).get("samples") if isinstance(samples_data, dict) else samples_data
		if not samples:
			return html.Div(render_sample_table([], selectable=True), className="card")
		return html.Div(render_sample_table(samples, selectable=True), className="card")

	@app.callback(
		Output("store-samples", "data", allow_duplicate=True),
		Output("samples-table", "children", allow_duplicate=True),
		Output("samples-upload-status", "children", allow_duplicate=True),
		Input("btn-clear-samples", "n_clicks"),
		prevent_initial_call=True,
	)
	def _clear_samples(n_clicks):
		if not n_clicks:
			raise PreventUpdate
		return {}, render_sample_table([], selectable=True), html.Div("Samples cleared.", className="status-pill")

	@app.callback(
		Output("store-samples", "data", allow_duplicate=True),
		Input({"type": "sample-checkbox", "index": ALL}, "value"),
		State("store-samples", "data"),
		prevent_initial_call=True,
	)
	def _update_sample_selection(checkbox_values, samples_data):
		if not callback_context.triggered or samples_data is None:
			raise PreventUpdate
		samples = (samples_data or {}).get("samples") or []
		triggered_ids = [t["prop_id"].split(".")[0] for t in callback_context.triggered]
		for raw_id, values in zip(triggered_ids, checkbox_values):
			try:
				import json
				id_dict = json.loads(raw_id)
				sample_index = id_dict.get("index")
			except Exception:
				continue
			for sample in samples:
				if sample.get("internal_id") == sample_index:
					sample["selected"] = bool(values)
		return {**samples_data, "samples": samples}
