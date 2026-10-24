"""Resolve NCBI taxids into a full ranked lineage via Entrez.

BLAST reports a scientific name and a taxid per hit, but no lineage — so a
species-level call arrives with Kingdom..Family blank, which empties the sunburst,
the per-rank breakdown, and the community-composition plot. This module fills those
ranks by asking NCBI's taxonomy database for the lineage of each unique taxid.

It is deliberately best-effort. With no ``ENTREZ_EMAIL`` set, or with no network, it
returns whatever is cached plus a plain-language warning, and the upper ranks stay
blank rather than the whole Databases stage failing. Results are cached on disk keyed
by taxid, so re-runs and later runs never re-fetch a lineage already seen.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from core.config import get_entrez_config

# NCBI rank names -> the pipeline's column keys. NCBI carries "superkingdom"
# (Eukaryota/Bacteria/Archaea); the classic "kingdom" (Metazoa, Fungi, Viridiplantae)
# is the more informative top level, so when both are present "kingdom" wins by being
# assigned later in lineage order.
_RANK_MAP = {
    "superkingdom": "kingdom",
    "kingdom": "kingdom",
    "phylum": "phylum",
    "class": "class_",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
}

# The ranks this module fills. Genus/species are usually already recovered from the
# BLAST binomial; they are included so a taxid whose binomial was unparseable still
# gets them from NCBI.
LINEAGE_KEYS = ("kingdom", "phylum", "class_", "order", "family", "genus", "species")

_TAXID_RE = re.compile(r"^\d+$")

# Words that mark a name as a placeholder rather than a taxon. NCBI Taxonomy files
# environmental sequences under names such as "uncultured eukaryote" with rank
# "species", and SILVA writes the same kind of label ("uncultured_eukaryote"); taken at
# face value they count an unidentified sequence as resolved to species.
_PLACEHOLDER_MARKERS = frozenset(
    {"uncultured", "unidentified", "unclassified", "environmental", "sample", "clone", "isolate", "metagenome", "symbiont"}
)


def is_placeholder_name(name) -> bool:
    """True for "uncultured eukaryote", "uncultured_marine_eukaryote", "environmental sample"..."""
    tokens = str(name or "").replace("_", " ").lower().split()
    return bool(_PLACEHOLDER_MARKERS & {token.strip(".,;()") for token in tokens})


def _clean_taxids(taxids: Iterable[str]) -> List[str]:
    """Unique, order-preserving numeric taxids. ``"12345;67890"`` -> two ids."""
    seen: Dict[str, None] = {}
    for raw in taxids:
        for token in str(raw or "").replace(";", " ").replace(",", " ").split():
            token = token.strip()
            if _TAXID_RE.match(token):
                seen.setdefault(token, None)
    return list(seen)


def _load_cache(cache_path: Optional[Path]) -> dict:
    if not cache_path:
        return {}
    try:
        data = json.loads(Path(cache_path).read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _save_cache(cache_path: Optional[Path], cache: dict) -> None:
    if not cache_path:
        return
    try:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        Path(cache_path).write_text(json.dumps(cache), encoding="utf-8")
    except OSError:
        pass


def _lineage_from_record(record: dict) -> Dict[str, str]:
    ranks: Dict[str, str] = {}
    for node in record.get("LineageEx", []) or []:
        rank = str(node.get("Rank") or "").lower()
        key = _RANK_MAP.get(rank)
        name = str(node.get("ScientificName") or "").strip()
        if key and name:
            ranks[key] = name
    # The record's own terminal taxon (usually the species itself).
    own_key = _RANK_MAP.get(str(record.get("Rank") or "").lower())
    own_name = str(record.get("ScientificName") or "").strip()
    if own_key and own_name:
        ranks.setdefault(own_key, own_name)
    return ranks


def resolve_lineages(
    taxids: Iterable[str],
    *,
    cache_path: Optional[Path] = None,
    batch_size: int = 200,
) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """Map each taxid -> ``{kingdom, phylum, ..., species}``. Best-effort.

    Returns ``(lineages, warnings)``. Taxids that could not be resolved simply do not
    appear in ``lineages``; ``warnings`` says why (no email, offline, ...) so the UI
    can explain that the upper ranks are blank on purpose rather than by a bug.
    """
    wanted = _clean_taxids(taxids)
    warnings: List[str] = []
    cache = _load_cache(cache_path)
    lineages: Dict[str, Dict[str, str]] = {tid: cache[tid] for tid in wanted if tid in cache}
    missing = [tid for tid in wanted if tid not in cache]
    if not missing:
        return lineages, warnings

    cfg = get_entrez_config()
    if not cfg["email"]:
        warnings.append(
            "Taxonomic lineage (Kingdom→Family) not resolved: ENTREZ_EMAIL is missing. "
            "Set ENTREZ_EMAIL in .env so NCBI fills the full lineage; "
            "without it, the ranks above genus stay blank."
        )
        return lineages, warnings

    try:
        from Bio import Entrez
    except ImportError:
        warnings.append("Biopython unavailable: NCBI lineage not resolved.")
        return lineages, warnings

    Entrez.email = cfg["email"]
    # NCBI asks automated callers to name their tool alongside the email.
    Entrez.tool = "PARAFISH"
    if cfg["api_key"]:
        Entrez.api_key = cfg["api_key"]

    fetched_any = False
    for start in range(0, len(missing), batch_size):
        chunk = missing[start : start + batch_size]
        try:
            handle = Entrez.efetch(db="taxonomy", id=",".join(chunk), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
        except Exception as exc:  # network, HTTP, parse — all best-effort
            warnings.append(
                f"NCBI lineage not resolved for {len(missing) - start} taxid(s) "
                f"(offline or network error?): {str(exc)[:200]}"
            )
            break
        for record in records:
            ranks = _lineage_from_record(record)
            if not ranks:
                continue
            ids = [str(record.get("TaxId") or "").strip()]
            # A taxid merged into another comes back under the new record with the old
            # id listed in AkaTaxIds; register both so the queried id resolves.
            ids.extend(str(aka).strip() for aka in (record.get("AkaTaxIds") or []))
            for tid in ids:
                if tid:
                    lineages[tid] = ranks
                    cache[tid] = ranks
                    fetched_any = True

    if fetched_any:
        _save_cache(cache_path, cache)
    return lineages, warnings


def fill_records(records: Iterable[dict], lineages: Dict[str, Dict[str, str]]) -> None:
    """Fill blank rank fields of each record from its taxid's lineage, in place.

    Only blank fields are written, so a name the source already asserted (e.g. the
    species from the BLAST binomial) is never overwritten by the lineage lookup.
    """
    for record in records:
        ranks = lineages.get(str(record.get("taxid") or "").strip())
        if not ranks:
            continue
        for key in LINEAGE_KEYS:
            value = ranks.get(key)
            if key in ("genus", "species") and is_placeholder_name(value):
                continue
            if not str(record.get(key) or "").strip() and value:
                record[key] = value
