"""Reusable generic helpers with no stage-specific logic."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re as _re
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


# Bumped whenever the SINTAX conversion changes meaning, so cached reference
# databases built by an older mapping are rebuilt instead of silently reused.
# v3: PR2 no longer forces its Division/Class onto the Linnean p:/c: ranks, and its
# unresolved-lineage placeholders are blanked. A cache built by v2 still carries
# "Phylum=Opisthokonta", so the name must change for it to be rebuilt.
# v4: SILVA stops claiming p/c/o/f/g positionally on lineages that declare no ranks
# (it was emitting Phylum=Amorphea, Genus=Teleostei). Same reason: v3 caches are wrong.
SINTAX_FORMAT_VERSION = "v4"

_SINTAX_BASE_RANKS = ("d", "p", "c", "o", "f", "g", "s")

# PR2 v5 lineages carry nine fixed levels:
#   Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species
#
# Only three of those are safely equivalent to the Linnean ranks the consensus table
# is keyed on. PR2's backbone is its own, and it is not Linnean:
#
#   * Supergroup, Division and Subdivision have no Linnean counterpart at all.
#   * PR2's "Class"/"Order"/"Family" are not Linnean classes/orders/families. On the
#     metazoan branch they hold clades and placeholders — Craniata is a subphylum,
#     and Craniata_X is PR2's marker for "unresolved child of Craniata".
#
# The old mapping sent Division onto p: and Class onto c:, which is what produced
# Phylum=Opisthokonta and Class=Craniata for a fish, made taxa_table_phylum report
# "Opisthokonta" for the whole run, and made conflict_report flag a kingdom/phylum/
# class conflict on *every* OTU — an artefact of mixing two backbones against a
# Linnean source, not a disagreement about biology.
#
# So the non-equivalent ranks are dropped (``None``) rather than invented onto a slot
# they do not occupy. This is a real trade: for protists — PR2's actual domain — its
# Class (e.g. Dinophyceae) IS a genuine class, and it is lost here. The Linnean middle
# ranks must come from a Linnean source (NCBI/SILVA/BOLD). Preserving them properly
# needs PR2's own ranks as first-class columns in the consensus schema, which is a
# larger change than this one.
PR2_RANK_LABELS = ("d", None, None, None, None, None, None, "g", "s")

# PR2 names an unresolved child by suffixing its parent: Craniata -> Craniata_X ->
# Craniata_XX -> Craniata_XXX, and "<parent>_sp." for an unnamed species. These are
# not taxa. Carried through, they make a placeholder read as an identification —
# "Genus=Craniata_XXX, Species=Craniata_XXX_sp." says nothing except "a craniate".
_PR2_PLACEHOLDER = _re.compile(r"(_X+|_sp\.?)$", _re.IGNORECASE)


def is_pr2_placeholder(value) -> bool:
	"""True for PR2's ``*_X``/``*_XX``/``*_sp.`` unresolved-lineage markers."""
	return bool(_PR2_PLACEHOLDER.search(str(value or "").strip()))


def _clean_taxon(value) -> str:
	"""Strip the delimiters SINTAX's own tax= field uses."""
	value = str(value or "").strip()
	for char in (",", ";", ":"):
		value = value.replace(char, "_")
	return "_".join(value.split())


def silva_rank_labels(depth: int) -> list:
	"""Rank letters for a SILVA lineage of ``depth`` levels. Only the ends are claimed.

	SILVA's ``tax_silva`` lineages vary in depth and **declare no ranks at all**: the
	file gives an ordered path, and only two positions have a guaranteed meaning — the
	first field is the domain and the last is the organism name.

	The intermediate levels are therefore left unranked (``None``). Claiming them
	positionally, as this used to (``p, c, o, f`` for the first four, ``g`` for the
	last-but-one), does not recover SILVA's ranks — it fabricates them. On the
	reference run that produced Phylum=Amorphea (a supergroup), Class=Obazoa (a clade)
	and Genus=Teleostei (an infraclass) for a fish, which is the same error PR2's
	Division→Phylum mapping made, from the other direction. Two sources inventing
	different middle ranks is also what made conflict_report flag a kingdom/phylum/
	class conflict on every OTU.

	The cost is that SILVA now contributes only domain and organism. That is what SILVA
	actually knows without its separate rank file (``tax_slv_*.txt``, which maps each
	taxon path to a declared rank and is not shipped here); parsing that file is what
	would let the middle ranks come back honestly.

	The returned list always has exactly ``depth`` entries, so no field is truncated.
	"""
	if depth <= 0:
		return []
	if depth == 1:
		return ["d"]
	if depth == 2:
		return ["d", "s"]
	return ["d", *([None] * (depth - 2)), "s"]


# The ranks SILVA's rank file can declare that have a Linnean column in the consensus.
# Everything else it declares (major_clade, superphylum, subclass, kingdom, ...) has no
# slot and is left out rather than squeezed into a neighbouring one. "kingdom" is left
# out on purpose: the SINTAX d: field already fills the kingdom column with the domain.
SILVA_DECLARED_RANK_LETTERS = {
	"domain": "d",
	"phylum": "p",
	"class": "c",
	"order": "o",
	"family": "f",
	"genus": "g",
}


def load_silva_rank_map(path) -> dict:
	"""Read SILVA's rank file (``tax_slv_ssu_<release>.txt``, plain or gzipped).

	Each line is ``path<TAB>taxid<TAB>rank<TAB>remark<TAB>release``, where ``path`` is the
	lineage up to and including that taxon, ``;``-terminated (``Eukaryota;Amorphea;``).
	Returns ``{path: rank}`` with the rank lower-cased. Lines without a path or a rank are
	skipped. An empty map means the file is not a SILVA rank file; the caller refuses it.
	"""
	path = str(path)
	opener = gzip.open if path.lower().endswith(".gz") else open
	rank_map = {}
	with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
		for line in handle:
			parts = line.rstrip("\r\n").split("\t")
			if len(parts) < 3:
				continue
			lineage = parts[0].strip()
			rank = parts[2].strip().lower()
			if not lineage or not rank:
				continue
			if not lineage.endswith(";"):
				lineage += ";"
			rank_map[lineage] = rank
	return rank_map


def silva_declared_labels(fields: list, rank_map: dict) -> list:
	"""Rank letters for a SILVA lineage, taken from the ranks SILVA itself declares.

	The counterpart of :func:`silva_rank_labels` when the rank file is available. Each
	intermediate field is looked up by its accumulated path (``A;B;C;``), which is how the
	rank file keys taxa, so a name reused on different branches resolves on its own
	branch. Only the ranks in ``SILVA_DECLARED_RANK_LETTERS`` are kept; each letter is
	used once and only in descending order (d, p, c, o, f, g), so a malformed path cannot
	put a family above an order. The two ends keep their fixed meaning: the first field is
	the domain and the last is the organism name (it is not in the rank file).

	Always returns exactly ``len(fields)`` entries, like ``silva_rank_labels``.
	"""
	depth = len(fields)
	if depth <= 0:
		return []
	labels = [None] * depth
	labels[0] = "d"
	if depth == 1:
		return labels
	order = ("d", "p", "c", "o", "f", "g")
	last_index = 0
	for position in range(1, depth - 1):
		lineage = ";".join(fields[: position + 1]) + ";"
		letter = SILVA_DECLARED_RANK_LETTERS.get(rank_map.get(lineage, ""))
		if not letter:
			continue
		index = order.index(letter)
		if index <= last_index:
			continue
		labels[position] = letter
		last_index = index
	labels[-1] = "s"
	return labels


def taxonomy_fields_to_sintax(fields: list, rank_labels=None) -> str:
	"""Convert an ordered list of taxonomy levels into a VSEARCH SINTAX tax= value.

	``rank_labels`` gives the SINTAX letter for each incoming field; a ``None``
	entry drops it. Pass the database's own mapping (``PR2_RANK_LABELS``,
	``silva_rank_labels``) whenever the source's levels are not plain Linnean —
	positional assignment is only correct for an ordered kingdom..species list.
	"""
	if rank_labels is None:
		labels = list(_SINTAX_BASE_RANKS)
		if len(fields) > len(labels):
			labels += [f"r{i}" for i in range(len(labels) + 1, len(fields) + 1)]
	else:
		labels = list(rank_labels)
	parts = []
	for label, value in zip(labels, fields):
		if not label:
			continue
		value = _clean_taxon(value)
		if value:
			parts.append(f"{label}:{value}")
	return ",".join(parts)


def prefix_fasta_ids(fasta_text: str, prefix: str) -> str:
	"""Prefix all FASTA headers with a sample identifier."""
	lines = []
	for line in fasta_text.splitlines():
		if line.startswith(">"):
			lines.append(f">{prefix}_{line[1:]}")
		else:
			lines.append(line)
	return "\n".join(lines) + "\n"
