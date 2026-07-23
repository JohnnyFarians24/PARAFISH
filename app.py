"""Future entrypoint for the fish parasite analysis platform."""

from __future__ import annotations

import json as _json
import threading
from importlib import import_module
from pathlib import Path
from typing import Any


_pipeline_threads: dict[str, threading.Thread] = {}
_progress_lock = threading.Lock()


def create_app() -> Any:
	"""Create and return the Dash application instance."""
	try:
		dash_module = import_module("dash")
	except ImportError as exc:  # pragma: no cover
		raise RuntimeError("Dash is required to create the web application.") from exc

	from ui.layout import get_layout
	from ui.callbacks_project import register_callbacks as register_project_callbacks
	from ui.callbacks_samples import register_callbacks as register_samples_callbacks
	from ui.callbacks_analysis import register_callbacks as register_analysis_callbacks
	from ui.callbacks_databases import register_callbacks as register_database_callbacks
	from ui.callbacks_phylo import register_callbacks as register_phylo_callbacks
	from ui.callbacks_preprocess import register_callbacks as register_preprocess_callbacks
	from ui.callbacks_setup import register_callbacks as register_setup_callbacks
	from ui.callbacks_progress import register_callbacks as register_progress_callbacks
	from ui.callbacks_results import register_callbacks as register_results_callbacks
	from ui.callbacks_results_summary import register_callbacks as register_results_summary_callbacks
	from ui.callbacks_qc import register_callbacks as register_qc_callbacks

	flask_module = import_module("flask")
	server = flask_module.Flask(__name__)

	def _resolve_within_runs(candidate: str) -> Path | None:
		"""Resolve ``candidate`` and confirm it stays inside the runs root.

		Returns the resolved path when it is contained in ``get_runs_dir()``,
		otherwise ``None``. Guards the download routes against directory
		traversal (e.g. ``run_dir=/etc`` or ``../../secrets``).
		"""
		from core.config import get_runs_dir

		if not candidate:
			return None
		try:
			runs_root = get_runs_dir().resolve()
			resolved = Path(candidate).resolve()
		except (OSError, ValueError):
			return None
		if resolved != runs_root and runs_root not in resolved.parents:
			return None
		return resolved

	@server.route("/download/run_zip")
	def _download_run_zip():
		from flask import abort, request, send_file
		from core.export import build_run_zip

		run_dir = request.args.get("run_dir", "")
		if not run_dir:
			abort(400)
		resolved_run_dir = _resolve_within_runs(run_dir)
		if resolved_run_dir is None:
			abort(403)
		run_dir = str(resolved_run_dir)
		progress_path = resolved_run_dir / "progress.json"
		if not progress_path.exists():
			abort(404)
		try:
			state = _json.loads(progress_path.read_text(encoding="utf-8"))
		except Exception:
			abort(500)
		run_id = Path(run_dir).name
		stage_outputs = {
			"qc": state.get("qc_record"),
			"preprocess": state.get("preprocess_record"),
			"analysis": state.get("analysis_record"),
			"databases": state.get("db_record"),
			"phylo": state.get("phylo_record"),
		}
		buf = build_run_zip(run_id, stage_outputs)
		return send_file(buf, mimetype="application/zip", as_attachment=True, download_name=f"{run_id}.zip")

	@server.route("/health")
	def _health():
		"""Cheap liveness probe for the connection banner.

		The banner needs to distinguish "the server is thinking" from "the server is
		gone". It has to be gone-detectable from the browser, because when the OOM
		killer takes the process down there is nobody left to report it — which is
		exactly the state the old wall of ten identical "server did not respond"
		toasts was failing to communicate.
		"""
		from flask import jsonify

		return jsonify({"ok": True})

	# Serve run artifacts under GET /runs/<path:relpath>
	@server.route("/runs/<path:relpath>")
	def _serve_runs(relpath: str):
		from flask import abort, send_from_directory
		from core.config import get_runs_dir
		root = get_runs_dir()
		# Defense in depth: reject any path that escapes the runs root before
		# delegating to send_from_directory (which also enforces safe joins).
		if _resolve_within_runs(str(root / relpath)) is None:
			abort(403)
		return send_from_directory(str(root), relpath, as_attachment=True)
	server.config["MAX_CONTENT_LENGTH"] = 4 * 1024 * 1024 * 1024
	app = dash_module.Dash(__name__, server=server, suppress_callback_exceptions=True)
	app.layout = get_layout()
	app.validation_layout = dash_module.html.Div(
		[
			get_layout(),
		],
	)
	from ui.callbacks_ordination import register as register_ordination
	from ui.components.connection_banner import register_clientside

	register_clientside(app)
	# One ordination implementation, three mounts: full controls on Analysis, a
	# read-only compact plot on Results.
	register_ordination(app, prefix="ord-analysis")
	register_ordination(app, prefix="ord-results", compact=True)
	register_ordination(app, prefix="ord-phylo")
	register_project_callbacks(app)
	register_samples_callbacks(app)
	register_setup_callbacks(app)
	register_progress_callbacks(app)
	register_results_callbacks(app)
	register_results_summary_callbacks(app)
	register_qc_callbacks(app)
	register_preprocess_callbacks(app)
	register_analysis_callbacks(app)
	register_database_callbacks(app)
	register_phylo_callbacks(app)

	from ui.callbacks_phylo import register_phylo_kpis

	register_phylo_kpis(app)
	return app


def main() -> None:
	"""Run the development server for manual local testing."""
	app = create_app()
	# Startup reconciliation: a run whose progress.json still says "running" has no
	# owning thread now (the server is only just starting), so it is an OOM-killed or
	# abandoned orphan. Mark those interrupted before serving — otherwise the newest one
	# masquerades as the active run and the disk-restore rebuilds BLANK Eixo B pages from
	# it, indistinguishable from that bug regressing. Wired here, not in create_app, so
	# the test suite (which builds the app but does not start the server) never mutates
	# real runs. The predicate keys on the live-thread map, empty at startup.
	from core.run_registry import reconcile_interrupted_runs

	with _progress_lock:
		_live_run_ids = {rid for rid, thread in _pipeline_threads.items() if thread.is_alive()}
	reconcile_interrupted_runs(is_live=lambda run_id: run_id in _live_run_ids)
	# dev_tools_ui stays ON. Turning it off to suppress the wall of identical
	# "server did not respond" toasts also suppressed *genuine* callback exceptions, and a
	# real bug then failed silently — the panel simply never updated and the app looked
	# dead with nothing to go on. That is a far worse failure mode than a noisy overlay.
	#
	# The duplicate-toast problem is handled where it belongs: the connection banner
	# reports a dead server once, with a retry count. A crashed callback is a different
	# thing and must remain visible.
	app.run(
		debug=True,
		threaded=True,
		use_reloader=False,
		dev_tools_props_check=False,
	)


if __name__ == "__main__":
	main()


