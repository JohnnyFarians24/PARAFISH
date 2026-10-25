"""What a Databases table download hands over, and the per-database assignment counts.

Pure functions, no Dash: the callbacks pass in what the DataTable reports
(``selected_rows``, ``derived_virtual_indices``) and get back row positions.
"""

from __future__ import annotations

from typing import Iterable, Optional


def table_view_state(data, selected_rows, virtual_indices) -> dict:
	"""The part of a DataTable's state a download needs, as stored between callbacks."""
	return {
		"n_rows": len(data or []),
		"selected_rows": list(selected_rows or []),
		"virtual_indices": None if virtual_indices is None else list(virtual_indices),
	}


def rows_to_export(n_rows: int, selected_rows: Optional[Iterable] = None, virtual_indices: Optional[Iterable] = None) -> Optional[list]:
	"""Row positions to export, or ``None`` for "the whole table, exactly as before".

	* Rows ticked in the table: those rows, in the order the table shows them (a ticked row
	  that the table filter now hides is still exported, after the visible ones: it was
	  chosen explicitly).
	* Nothing ticked, table filtered: the rows the filter leaves, in the order shown.
	* Nothing ticked, nothing filtered (sorting alone does not count): ``None``, and the
	  caller keeps its previous behaviour byte for byte.

	Positions outside ``range(n_rows)`` are ignored, so a stale selection from a previous
	render can never pick the wrong row; the caller also checks ``n_rows``.
	"""
	try:
		total = int(n_rows)
	except (TypeError, ValueError):
		return None

	def _valid(values) -> list:
		out = []
		for value in values or []:
			if isinstance(value, bool):
				continue
			if isinstance(value, int) and 0 <= value < total:
				out.append(value)
		return out

	selected = _valid(selected_rows)
	view = None if virtual_indices is None else _valid(virtual_indices)
	if selected:
		chosen = set(selected)
		ordered = [index for index in (view or []) if index in chosen]
		seen = set(ordered)
		ordered += [index for index in sorted(chosen) if index not in seen]
		return ordered
	if view is not None and len(set(view)) < total:
		return view
	return None


def export_positions(state, n_rows: int) -> Optional[list]:
	"""``rows_to_export`` from a stored ``table_view_state``, or ``None`` if it is stale.

	The stored state belongs to the table as it was rendered. If the table the download
	rebuilds has a different number of rows (the filters changed and the table has not
	re-reported), the state is not trusted and the whole table is exported, as before.
	"""
	if not isinstance(state, dict):
		return None
	if state.get("n_rows") != n_rows:
		return None
	return rows_to_export(n_rows, state.get("selected_rows"), state.get("virtual_indices"))


def source_assignment_summary(frame):
	"""One row per database: how many features it assigned and how many reads they hold.

	``frame`` is the ranked hits table (one row per feature x hit, with ``source``,
	``feature_id`` and ``reads``). A feature counts for a source when the source has at
	least one hit for it. Reads are those of the feature (the same for all its rows), so a
	feature's reads are counted once per source, and ``% of reads`` is relative to all the
	features present in ``frame``. Sorted by features, then reads, descending.
	"""
	import pandas as pd

	columns = ["Database", "Features assigned", "Reads in those features", "% of reads"]
	if frame is None or frame.empty or "source" not in frame.columns:
		return pd.DataFrame(columns=columns)
	id_col = "feature_id" if "feature_id" in frame.columns else ("ASV_ID" if "ASV_ID" in frame.columns else None)
	if id_col is None:
		return pd.DataFrame(columns=columns)

	work = frame[[id_col, "source"] + (["reads"] if "reads" in frame.columns else [])].copy()
	work["source"] = work["source"].fillna("").astype(str).str.strip()
	work = work[work["source"] != ""]
	if "reads" in work.columns:
		work["reads"] = pd.to_numeric(work["reads"], errors="coerce")
	else:
		work["reads"] = float("nan")

	per_feature_reads = work.groupby(id_col)["reads"].first()
	total_reads = per_feature_reads.sum(min_count=1)

	rows = []
	for source, group in work.groupby("source"):
		features = group[id_col].drop_duplicates()
		reads = per_feature_reads.reindex(features).sum(min_count=1)
		pct = (100.0 * reads / total_reads) if pd.notna(reads) and pd.notna(total_reads) and total_reads else float("nan")
		rows.append(
			{
				"Database": source,
				"Features assigned": int(len(features)),
				"Reads in those features": "" if pd.isna(reads) else int(reads),
				"% of reads": "" if pd.isna(pct) else round(float(pct), 2),
			}
		)
	result = pd.DataFrame(rows, columns=columns)
	sort_reads = pd.to_numeric(result["Reads in those features"], errors="coerce").fillna(-1)
	order = sorted(range(len(result)), key=lambda i: (-result["Features assigned"].iloc[i], -sort_reads.iloc[i], result["Database"].iloc[i]))
	return result.iloc[order].reset_index(drop=True)
