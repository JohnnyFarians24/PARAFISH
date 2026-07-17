"""Storage retention and cleanup helpers.

Only the run-cleanup surface is wired into the runtime (see
``_apply_cleanup_post_run``, called from ``ui/callbacks_progress.py``). The
legacy session-metadata and per-step manifest helpers were removed as dead code.
"""

from __future__ import annotations

import gzip
import os
import shutil
import time

from .paths import _project_exports_dir


class StorageConfig:
    """Configuration for storage management policies."""

    MAX_RUNS_TO_KEEP = int(os.getenv("MAX_RUNS_TO_KEEP", "1"))
    MAX_RUN_AGE_DAYS = int(os.getenv("MAX_RUN_AGE_DAYS", "30"))
    AUTO_DELETE_TEMP = str(os.getenv("AUTO_DELETE_TEMP", "1")).strip().lower() in {"1", "true", "yes", "on"}
    TEMP_FILE_MAX_AGE_DAYS = int(os.getenv("TEMP_FILE_MAX_AGE_DAYS", "3"))
    DELETE_INTERMEDIATE_STEPS = [
        step.strip()
        for step in os.getenv("DELETE_INTERMEDIATE_STEPS", "qc,preprocess").split(",")
        if step.strip()
    ]
    MAX_EXPORTS_TO_KEEP = int(os.getenv("MAX_EXPORTS_TO_KEEP", "10"))


def _apply_cleanup_post_run(project_dir, run_dir):
    """Apply cleanup after a run completes."""
    if not project_dir or not isinstance(project_dir, str):
        return False
    if not run_dir or not isinstance(run_dir, str):
        return False
    try:
        StorageManager.cleanup_old_runs(project_dir)
        if StorageConfig.AUTO_DELETE_TEMP:
            StorageManager.cleanup_temp_files(project_dir)
        StorageManager.compact_run_outputs(run_dir)
        StorageManager.cleanup_export_archives(project_dir, run_dir)
        return True
    except Exception:
        return False


class StorageManager:
    """Manages storage cleanup and retention policies."""

    @staticmethod
    def cleanup_old_runs(project_dir):
        runs_dir = os.path.join(project_dir, "runs")
        if not os.path.exists(runs_dir):
            return []
        deleted = []
        try:
            entries = []
            for entry in os.listdir(runs_dir):
                entry_path = os.path.join(runs_dir, entry)
                if os.path.isdir(entry_path):
                    entries.append((os.path.getmtime(entry_path), entry_path, entry))
            entries.sort(reverse=True)
            for _, entry_path, entry in entries[StorageConfig.MAX_RUNS_TO_KEEP:]:
                try:
                    shutil.rmtree(entry_path)
                    deleted.append(entry)
                except Exception:
                    pass
        except Exception:
            pass
        return deleted

    @staticmethod
    def cleanup_temp_files(project_dir):
        deleted = []
        temp_patterns = [
            ("steps/*/temp_*", "step temp files"),
            ("runs/*/steps/*/tmp_*", "run temp files"),
            ("uploads/*_temp_*", "upload temp files"),
        ]
        now = time.time()
        max_age_seconds = StorageConfig.TEMP_FILE_MAX_AGE_DAYS * 86400
        try:
            import glob

            for pattern, _ in temp_patterns:
                full_pattern = os.path.join(project_dir, pattern)
                for filepath in glob.glob(full_pattern):
                    try:
                        if now - os.path.getmtime(filepath) > max_age_seconds:
                            if os.path.isfile(filepath):
                                os.remove(filepath)
                            elif os.path.isdir(filepath):
                                shutil.rmtree(filepath)
                            deleted.append(os.path.basename(filepath))
                    except Exception:
                        pass
        except Exception:
            pass
        return deleted

    @staticmethod
    def compact_run_outputs(run_dir):
        compressed_files = []
        try:
            import glob

            for pattern in ("steps/*/intermediate/*.fastq", "steps/*/intermediate/*.fq"):
                full_pattern = os.path.join(run_dir, pattern)
                for filepath in glob.glob(full_pattern):
                    try:
                        if os.path.getsize(filepath) > 10 * 1024 * 1024:
                            gz_path = filepath + ".gz"
                            if not os.path.exists(gz_path):
                                with open(filepath, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                os.remove(filepath)
                                compressed_files.append(os.path.basename(filepath))
                    except Exception:
                        pass
        except Exception:
            pass
        return compressed_files

    @staticmethod
    def cleanup_export_archives(project_dir, run_dir=None):
        deleted = []
        try:
            if run_dir and os.path.isdir(run_dir):
                export_staging_dir = os.path.join(run_dir, "export_staging")
                if os.path.isdir(export_staging_dir):
                    shutil.rmtree(export_staging_dir, ignore_errors=True)
                for name in os.listdir(run_dir):
                    file_path = os.path.join(run_dir, name)
                    if os.path.isfile(file_path) and name.lower().endswith(".zip"):
                        try:
                            os.remove(file_path)
                            deleted.append(file_path)
                        except OSError:
                            pass

            exports_dir = _project_exports_dir(project_dir) if project_dir else None
            if exports_dir and os.path.isdir(exports_dir):
                zip_entries = []
                for name in os.listdir(exports_dir):
                    file_path = os.path.join(exports_dir, name)
                    if os.path.isfile(file_path) and name.lower().endswith(".zip"):
                        try:
                            zip_entries.append((os.path.getmtime(file_path), file_path))
                        except OSError:
                            continue

                zip_entries.sort(key=lambda item: item[0], reverse=True)
                for _, file_path in zip_entries[StorageConfig.MAX_EXPORTS_TO_KEEP :]:
                    try:
                        os.remove(file_path)
                        deleted.append(file_path)
                    except OSError:
                        pass
        except Exception:
            pass
        return deleted
