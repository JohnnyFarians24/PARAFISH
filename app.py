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
	from ui.callbacks_qc import register_callbacks as register_qc_callbacks

	flask_module = import_module("flask")
	server = flask_module.Flask(__name__)

	@server.route("/download/run_zip")
	def _download_run_zip():
		from flask import abort, request, send_file
		from core.export import build_run_zip

		run_dir = request.args.get("run_dir", "")
		if not run_dir:
			abort(400)
		progress_path = Path(run_dir) / "progress.json"
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

	# Serve run artifacts under GET /runs/<path:relpath>
	@server.route("/runs/<path:relpath>")
	def _serve_runs(relpath: str):
		from flask import send_from_directory
		from core.config import get_runs_dir
		root = str(get_runs_dir())
		# Prevent directory traversal by resolving path inside root
		return send_from_directory(root, relpath, as_attachment=True)
	server.config["MAX_CONTENT_LENGTH"] = 4 * 1024 * 1024 * 1024
	app = dash_module.Dash(__name__, server=server, suppress_callback_exceptions=True)
	app.layout = get_layout()
	app.validation_layout = dash_module.html.Div(
		[
			get_layout(),
		],
	)
	register_project_callbacks(app)
	register_samples_callbacks(app)
	register_setup_callbacks(app)
	register_progress_callbacks(app)
	register_results_callbacks(app)
	register_qc_callbacks(app)
	register_preprocess_callbacks(app)
	register_analysis_callbacks(app)
	register_database_callbacks(app)
	register_phylo_callbacks(app)
	return app


def main() -> None:
	"""Run the development server for manual local testing."""
	app = create_app()
	app.run(debug=True, threaded=True, use_reloader=False)


if __name__ == "__main__":
	main()


