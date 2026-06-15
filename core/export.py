"""Build structured run ZIP for download."""

from __future__ import annotations

import io
import zipfile
from pathlib import Path


STAGE_ORDER = ["qc", "preprocess", "analysis", "databases", "phylo"]


def build_run_zip(run_id: str, stage_outputs: dict) -> io.BytesIO:
	buf = io.BytesIO()
	with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
		for stage in STAGE_ORDER:
			record = stage_outputs.get(stage) or {}
			outputs = record.get("outputs") or {}
			for sample_id, sample_data in outputs.items():
				if sample_id in ("combined_fasta", "dada2", "taxonomy"):
					# top-level stage outputs, not per sample
					_add_stage_files(zf, run_id, stage, None, sample_data if isinstance(sample_data, dict) else {})
					continue
				if not isinstance(sample_data, dict):
					continue
				if "files" in sample_data:
					# QC-style: {sample_id: {files: {...}, metrics: {...}}}
					_add_stage_files(zf, run_id, stage, sample_id, sample_data.get("files") or {})
				else:
					# Preprocess-style: {sample_id: {trimmomatic: {files: {...}}, pear: {files: {...}}}}
					for tool, tool_data in sample_data.items():
						if isinstance(tool_data, dict):
							_add_stage_files(zf, run_id, stage, sample_id, tool_data.get("files") or {})
		# combined_fasta and taxonomy as top-level files
		for stage in ("analysis", "databases"):
			record = stage_outputs.get(stage) or {}
			outputs = record.get("outputs") or {}
			cf = outputs.get("combined_fasta")
			if cf and Path(cf).exists():
				zf.write(cf, f"{run_id}/{stage}/{Path(cf).name}")
			tax = outputs.get("taxonomy")
			if isinstance(tax, dict):
				_add_stage_files(zf, run_id, stage, None, tax.get("files") or {})
	buf.seek(0)
	return buf


def _add_stage_files(zf, run_id, stage, sample_id, files: dict):
	for key, path_str in files.items():
		if not path_str or not isinstance(path_str, str):
			continue
		p = Path(path_str)
		if not p.exists():
			continue
		if sample_id:
			arc_path = f"{run_id}/{stage}/{sample_id}/{p.name}"
		else:
			arc_path = f"{run_id}/{stage}/{p.name}"
		zf.write(str(p), arc_path)