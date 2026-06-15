"""Plot and tree conversion helpers for phylogenetic outputs."""

from __future__ import annotations

from pathlib import Path


def biopython_to_cytoscape_elements(tree) -> list[dict]:
    """
    Convert a Bio.Phylo tree object to Cytoscape.js elements for Dash Cytoscape.

    Returns a list of node and edge dicts compatible with dash-cytoscape.
    """
    try:
        from Bio import Phylo
    except Exception:
        return []

    elements = []
    node_id = [0]

    def traverse(clade, parent_id=None):
        current_id = f"node_{node_id[0]}"
        node_id[0] += 1
        label = clade.name or current_id
        elements.append({"data": {"id": current_id, "label": label}})
        if parent_id is not None:
            branch_length = clade.branch_length or 0.0
            elements.append({
                "data": {
                    "source": parent_id,
                    "target": current_id,
                    "weight": round(branch_length, 6),
                }
            })
        for child in clade.clades:
            traverse(child, current_id)

    traverse(tree.root)
    return elements


def build_phylo_plot_data(tree_path: str | Path) -> dict:
    """
    Build a plotting payload from a Newick tree file.

    Returns:
        {
            "newick": str,           # raw newick string for JS renderers
            "node_count": int,
            "leaf_count": int,
            "cytoscape_elements": list,
        }
    """
    tree_path = Path(tree_path)
    if not tree_path.exists() or tree_path.stat().st_size == 0:
        return {"error": "Tree file not found or empty", "newick": "", "node_count": 0, "leaf_count": 0, "cytoscape_elements": []}

    newick = tree_path.read_text(encoding="utf-8", errors="replace").strip()

    try:
        from io import StringIO
        from Bio import Phylo
        tree = Phylo.read(StringIO(newick), "newick")
        node_count = len(list(tree.find_clades()))
        leaf_count = len(tree.get_terminals())
        cytoscape_elements = biopython_to_cytoscape_elements(tree)
    except Exception:
        node_count = newick.count("(")
        leaf_count = 0
        cytoscape_elements = []

    return {
        "newick": newick,
        "node_count": node_count,
        "leaf_count": leaf_count,
        "cytoscape_elements": cytoscape_elements,
    }
