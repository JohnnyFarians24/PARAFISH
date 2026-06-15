"""Custom reference database adapter."""

from __future__ import annotations
from pathlib import Path
from typing import Optional


class CustomDBConfig:
    def __init__(self, db_path: Path, sintax_cutoff: float = 0.8):
        self.db_path = Path(db_path)
        self.sintax_cutoff = sintax_cutoff

    def validate(self) -> tuple[bool, str]:
        if not self.db_path.exists():
            return False, f"Custom DB not found: {self.db_path}"
        suffix = self.db_path.suffix.lower()
        if suffix not in {".fasta", ".fa", ".gz"}:
            return False, f"Unsupported format: {suffix}"
        return True, ""


def get_custom_reference(marker: str) -> str:
    return f"custom:{marker.lower()}"
