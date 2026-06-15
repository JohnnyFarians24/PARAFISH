"""Legacy compatibility entrypoint for the modular PARAFISH platform."""

from __future__ import annotations

from analysis.manifest import assign_samples_from_files
from analysis.sample_parsing import ParsedFastqName, parse_fastq_name, validate_manifest_rows
from app import create_app, main

__all__ = [
    "ParsedFastqName",
    "assign_samples_from_files",
    "create_app",
    "main",
    "parse_fastq_name",
    "validate_manifest_rows",
]


if __name__ == "__main__":
    main()