"""Reusable generic helpers with no stage-specific logic."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile


def _path_is_valid(path):
	return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _read_json(path, default=None):
	if not os.path.exists(path):
		return default
	try:
		with open(path, "r", encoding="utf-8") as f:
			return json.load(f)
	except Exception:
		return default


def _write_json(path, payload):
	os.makedirs(os.path.dirname(path), exist_ok=True)
	fd, tmp_path = tempfile.mkstemp(prefix="tmp_session_", suffix=".json", dir=os.path.dirname(path))
	try:
		with os.fdopen(fd, "w", encoding="utf-8") as f:
			json.dump(payload, f, indent=2)
		os.replace(tmp_path, path)
	finally:
		if os.path.exists(tmp_path):
			os.remove(tmp_path)


def _write_text_file(path, content):
	if not path:
		return None
	os.makedirs(os.path.dirname(path), exist_ok=True)
	with open(path, "w", encoding="utf-8") as handle:
		handle.write(content if content is not None else "")
	return path


def _is_within_directory(root_dir, target_path):
	root = os.path.abspath(root_dir)
	target = os.path.abspath(target_path)
	return target == root or target.startswith(root + os.sep)


def _safe_extract_zip(zip_path, dest_dir):
	os.makedirs(dest_dir, exist_ok=True)
	with zipfile.ZipFile(zip_path, "r") as zf:
		for member in zf.infolist():
			member_name = member.filename
			if not member_name or member_name.endswith("/"):
				continue
			final_path = os.path.abspath(os.path.join(dest_dir, member_name))
			if not _is_within_directory(dest_dir, final_path):
				raise ValueError(f"Unsafe ZIP path detected: {member_name}")
		zf.extractall(dest_dir)


def _safe_extract_archive(archive_path, dest_dir):
	lower = archive_path.lower()
	if lower.endswith(".zip"):
		_safe_extract_zip(archive_path, dest_dir)
		return
	if lower.endswith(".rar"):
		os.makedirs(dest_dir, exist_ok=True)
		unrar_bin = shutil.which("unrar")
		if unrar_bin:
			res = subprocess.run(
				[unrar_bin, "x", "-y", archive_path, dest_dir],
				capture_output=True,
				text=True,
			)
			if res.returncode == 0:
				return
		bsdtar_bin = shutil.which("bsdtar")
		if bsdtar_bin:
			res = subprocess.run(
				[bsdtar_bin, "-xf", archive_path, "-C", dest_dir],
				capture_output=True,
				text=True,
			)
			if res.returncode == 0:
				return
		raise RuntimeError("RAR extraction failed. Install 'unrar' or 'bsdtar'.")
	raise ValueError("Unsupported archive format. Please upload ZIP or RAR.")


def _download_href(file_path):
	file_path = os.path.abspath(file_path)
	upload_dir = os.path.abspath(os.path.join(os.getcwd(), "uploads"))
	if file_path.startswith(upload_dir + os.sep):
		rel = os.path.relpath(file_path, upload_dir).replace(os.sep, "/")
		return f"/uploads/{rel}"
	cwd = os.path.abspath(os.getcwd()) + os.sep
	if file_path.startswith(cwd):
		rel = os.path.relpath(file_path, os.getcwd()).replace(os.sep, "/")
		return f"/files/{rel}"
	return "#"


def _normalize_uploaded_paths(filenames):
	if not filenames:
		return []
	if isinstance(filenames, str):
		filenames = [filenames]
	if not isinstance(filenames, list):
		return []
	normalized = []
	for item in filenames:
		if item:
			normalized.append(str(item))
	return normalized


def _as_selection_list(value):
	if value in (None, "", []):
		return []
	if isinstance(value, list):
		return [str(v) for v in value if v]
	return [str(value)]


def _normalize_active_sample_value(value):
	selected = _as_selection_list(value)
	if not selected or "__ALL__" in selected:
		return "__ALL__"
	if len(selected) == 1:
		return selected[0]
	return selected


def _ctx_uploaded_files_store():
	try:
		import dash

		current = (dash.callback_context.states or {}).get("uploaded-files.data") or {}
	except Exception:
		current = {}
	if not isinstance(current, dict):
		return {}
	return current.copy()


def _merge_uploaded_files_from_ctx(patch):
	merged = _ctx_uploaded_files_store()
	if not isinstance(patch, dict):
		return merged
	merged.update(patch)
	return merged


def _resolve_local_source(candidates):
	for path in candidates:
		if path and os.path.exists(path) and os.path.getsize(path) > 0:
			return path
	return None


def _resolve_step_input(candidate_paths):
	for source_label, path in candidate_paths or []:
		if path and os.path.exists(path) and os.path.getsize(path) > 0:
			return source_label, path
	return None


def _count_fasta_seqs(path):
	if not path or not os.path.exists(path):
		return 0
	count = 0
	try:
		with open(path, "r", encoding="utf-8") as fh:
			for line in fh:
				if line.startswith(">"):
					count += 1
	except Exception:
		return 0
	return count


def _get_directory_size(path):
	if not path or not os.path.exists(path):
		return 0
	total = 0
	for root, _, files in os.walk(path):
		for name in files:
			file_path = os.path.join(root, name)
			if os.path.isfile(file_path):
				try:
					total += os.path.getsize(file_path)
				except OSError:
					continue
	return total


def _download_file(url, dest_path):
	os.makedirs(os.path.dirname(dest_path), exist_ok=True)
	urllib.request.urlretrieve(url, dest_path)
	return dest_path


def _gunzip_file(gz_path, out_path):
	with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
		shutil.copyfileobj(src, dst)
	return out_path


def _file_sha256(path):
	h = hashlib.sha256()
	with open(path, "rb") as f:
		while True:
			chunk = f.read(1024 * 1024)
			if not chunk:
				break
			h.update(chunk)
	return h.hexdigest()


def _mirror_output_to_results(source_path, results_path):
	if not source_path or not results_path:
		return results_path
	if not os.path.exists(source_path) or os.path.getsize(source_path) == 0:
		return results_path
	os.makedirs(os.path.dirname(results_path), exist_ok=True)
	try:
		if os.path.abspath(source_path) != os.path.abspath(results_path):
			shutil.copy2(source_path, results_path)
	except Exception:
		return results_path
	return results_path


def prefix_fasta_ids(fasta_text: str, prefix: str) -> str:
	"""Prefix all FASTA headers with a sample identifier."""
	lines = []
	for line in fasta_text.splitlines():
		if line.startswith(">"):
			lines.append(f">{prefix}_{line[1:]}")
		else:
			lines.append(line)
	return "\n".join(lines) + "\n"
