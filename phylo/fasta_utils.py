"""FASTA sanitization utilities shared across phylo runners."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

VALID_IUPAC = set("ACGTRYSWKMBDHVN")


def sanitize_fasta_for_phylo(fasta_path: Path, output_path: Path) -> Dict[str, str]:
    """Rename headers to ASV IDs, remove short sequences, normalize bases."""
    mapping: Dict[str, str] = {}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    current_header = None
    current_sequence: list[str] = []
    index = 1

    def flush(handle):
        nonlocal index, current_header, current_sequence
        if current_header is None:
            return
        sequence = "".join(current_sequence).upper()
        sequence = "".join(b if b in VALID_IUPAC else "N" for b in sequence)
        if len(sequence) < 100:
            current_header = None
            current_sequence = []
            return
        new_name = f"ASV{index}"
        mapping[new_name] = current_header
        handle.write(f">{new_name}\n{sequence}\n")
        index += 1
        current_header = None
        current_sequence = []

    with fasta_path.open("r", encoding="utf-8", errors="replace") as src, \
         output_path.open("w", encoding="utf-8") as dst:
        for line in src:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush(dst)
                current_header = line[1:].strip()
                continue
            current_sequence.append(line)
        flush(dst)

    return mapping
