"""Project and page routing callbacks."""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from typing import Any

from dash import Input, Output, State, ctx, html, no_update
from dash.exceptions import PreventUpdate
from dash import ALL

from core.project_store import load_all_projects, load_project_bundle, save_project_bundle
from ui.pipeline_state import get_stage_record, record_belongs_to, stage_status
from ui.components.sidebar import render_sidebar
from ui.components.topbar import render_topbar
from ui.runtime_state import resolve_project_data, restore_run_state

from ui.pages.page_analysis import layout as analysis_layout
from ui.pages.page_databases import layout as databases_layout
from ui.pages.page_phylo import layout as phylo_layout
from ui.pages.page_preprocess import layout as preprocess_layout
from ui.pages.page_project import layout as project_layout
from ui.pages.page_qc import layout as qc_layout
from ui.pages.page_results import layout as results_layout
from ui.pages.page_samples import layout as samples_layout
from ui.pages.page_setup import layout as setup_layout
from ui.pages.page_progress import layout as progress_layout


ROUTES = {
	"/": project_layout,
	"/project": project_layout,
	"/setup": setup_layout,
	"/samples": samples_layout,
	"/qc": qc_layout,
	"/preprocess": preprocess_layout,
	"/analysis": analysis_layout,
	"/databases": databases_layout,
	"/phylo": phylo_layout,
	"/progress": progress_layout,
	"/results": results_layout,
}


EXAMPLE_PROJECT = {
	"project_name": "Trout Gill Parasites",
	"mode": "Complete",
	"notes": "Example project loaded from the dashboard.",
	"source": "example",
	"created_at": "2026-05-22T00:00:00Z",
}


def _stage_states(pathname: str | None, project_data: dict | None, samples_data: dict | None, step_outputs: dict | None) -> list[str]:
	pathname = pathname or "/project"
	project_ready = bool(resolve_project_data(project_data).get("project_name"))
	samples_ready = bool(samples_data)
	qc_record = get_stage_record(step_outputs, "qc")
	preprocess_record = get_stage_record(step_outputs, "preprocess")
	analysis_record = get_stage_record(step_outputs, "analysis")
	databases_record = get_stage_record(step_outputs, "databases")
	phylo_record = get_stage_record(step_outputs, "phylo")
	enabled_records = [
		qc_record if True else None,
    	preprocess_record if True else None,
    	analysis_record if True else None,
    	databases_record if True else None,
    	phylo_record if True else None,
	]
	active_records = [r for r in enabled_records if r is not None]
	results_ready = bool(active_records) and all(stage_status(r) == "success" for r in active_records)
	
	stage_rows = [
		("/project", "done" if project_ready else "active"),
		("/samples", "done" if samples_ready else ("active" if project_ready else "pending")),
		("/setup", "active" if pathname == "/setup" else ("done" if samples_ready else "pending")),
		("/qc", stage_status(qc_record)),
		("/preprocess", stage_status(preprocess_record)),
		("/analysis", stage_status(analysis_record)),
		("/databases", stage_status(databases_record)),
		("/phylo", stage_status(phylo_record)),
		("/progress", "active" if pathname == "/progress" else "pending"),
		("/results", "done" if results_ready else "pending"),
	]

	states: list[str] = []
	failure_seen = False
	for route, status in stage_rows:
		if status == "failure":
			failure_seen = True
		if pathname == route and status not in {"success", "failure"}:
			states.append("active")
			continue
		if status == "success":
			states.append("done")
			continue
		if status == "failure":
			states.append("error")
			continue
		if failure_seen:
			states.append("pending")
			continue
		states.append("pending" if status in {"pending", "blocked"} else "done")
	return states


def _results_status(step_outputs: dict | None, enabled_stages: dict | None) -> str:
	enabled_stages = enabled_stages or {}
	qc_record = get_stage_record(step_outputs, "qc")
	preprocess_record = get_stage_record(step_outputs, "preprocess")
	analysis_record = get_stage_record(step_outputs, "analysis")
	databases_record = get_stage_record(step_outputs, "databases")
	phylo_record = get_stage_record(step_outputs, "phylo")
	relevant = [
		qc_record if enabled_stages.get("qc", True) else None,
		preprocess_record if enabled_stages.get("preprocess", True) else None,
		analysis_record if enabled_stages.get("analysis", True) else None,
		databases_record if enabled_stages.get("databases", True) else None,
		phylo_record if enabled_stages.get("phylo", True) else None,
	]
	active_records = [record for record in relevant if record is not None]
	if not active_records:
		return "pending"
	return "done" if active_records and all(stage_status(record) == "success" for record in active_records) else "pending"


def _selected_sample_count(samples_data) -> int:
	samples = (samples_data or {}).get("samples") if isinstance(samples_data, dict) else None
	return sum(1 for sample in (samples or []) if (sample or {}).get("selected", True))


def _active_project_banner(project_store, samples_data):
	"""A prominent 'which project am I working on' strip for Samples and Setup.

	The active project used to show only in the corner badge, so it was easy to create
	a fresh project, upload into a different one that was still active, and run the
	wrong data without noticing. This states it in-flow, at the top of the page.
	"""
	active = resolve_project_data(project_store)
	name = active.get("project_name") or active.get("name")
	if not name:
		return html.Div(
			[
				html.Span("⚠ No active project", className="badge badge-amber"),
				html.A("Create or select a project →", href="/project", style={"color": "var(--blue)", "marginLeft": "10px"}),
			],
			className="control-row",
			style={"marginBottom": "12px", "alignItems": "center"},
		)
	return html.Div(
		[
			html.Span("Active project", className="logo-subtext"),
			html.Span(str(name), className="badge badge-green"),
			html.Span(f"{_selected_sample_count(samples_data)} sample(s)", className="status-pill"),
		],
		className="control-row",
		style={"marginBottom": "12px", "alignItems": "center", "gap": "10px"},
	)


def _run_target_note(project_store, samples_data):
	"""One line, next to Run, naming exactly which project and how many samples run."""
	active = resolve_project_data(project_store)
	name = active.get("project_name") or active.get("name")
	if not name:
		return html.Div(html.Span("No active project — create one before running.", className="badge badge-amber"), className="mini-note")
	count = _selected_sample_count(samples_data)
	if count == 0:
		return html.Div(
			[html.Span("You are about to run project "), html.B(str(name)), html.Span(" — "), html.Span("no samples selected", className="badge badge-amber")],
			className="mini-note",
		)
	return html.Div([html.Span("You are about to run project "), html.B(str(name)), html.Span(f" · {count} sample(s)")], className="mini-note")


def _was_clicked(triggered) -> bool:
	"""True only for a real click: a re-rendered button reports n_clicks None (or 0)."""
	return any(bool(item.get("value")) for item in (triggered or []))


def register_callbacks(app: Any) -> None:

	@app.callback(
		Output("store-project", "data", allow_duplicate=True),
		Output("store-samples", "data", allow_duplicate=True),
		Output("url", "pathname", allow_duplicate=True),
		Input("btn-create-project", "n_clicks"),
		State("input-project-name", "value"),
		State("textarea-project-notes", "value"),
		State("store-project", "data"),
		prevent_initial_call=True,
	)
	def _create_project(create_clicks, project_name, notes, current_store):
		del create_clicks
		project_name = (project_name or "").strip()
		if not project_name:
			return no_update, no_update, no_update
		store = current_store or {"active_id": None, "projects": {}}
		projects = dict(store.get("projects", {}))
		proj_id = f"proj_{uuid.uuid4().hex[:6]}"
		projects[proj_id] = {
			"id": proj_id,
			"name": project_name,
			"project_name": project_name,
			"mode": "Complete",
			"notes": notes or "",
			"created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
			"samples": [],
			"status": "active",
			"run_id": f"run_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}",
		}
		store = {"active_id": proj_id, "projects": projects}
		save_project_bundle(proj_id, projects[proj_id], {})
		# Move the user to the next step. Creating the project used to leave them on this
		# page with no navigation, no confirmation and a session panel still reading
		# "empty until a project is selected" — so the button looked broken even though it
		# had worked. A new project has no samples, and samples are what it needs next.
		return store, {}, "/samples"

	@app.callback(
		Output("store-project", "data", allow_duplicate=True),
		Output("store-samples", "data", allow_duplicate=True),
		Output("url", "pathname", allow_duplicate=True),
		Input("btn-load-example", "n_clicks"),
		State("store-project", "data"),
		prevent_initial_call=True,
	)
	def _load_example(n_clicks, current_store):
		if not n_clicks:
			return no_update, no_update, no_update
		store = current_store or {"active_id": None, "projects": {}}
		projects = dict(store.get("projects", {}))
		proj_id = "proj_example"
		projects[proj_id] = {
			"id": proj_id,
			"name": "Example - Trout COI",
			"project_name": "Example - Trout COI",
			"mode": "Complete",
			"notes": "Demo project with synthetic COI reads.",
			"created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
			"samples": [],
			"status": "active",
			"run_id": "example_run",
		}
		store = {"active_id": proj_id, "projects": projects}
		return store, {}, "/samples"

	@app.callback(
		Output("project-session-state-info", "children"),
		Input("store-project", "data"),
		Input("store-samples", "data"),
		Input("store-run-config", "data"),
	)
	def _render_session_state(project_store, samples_data, run_config_data):
		"""What is loaded right now, and what to do next.

		This panel used to be a hard-coded placeholder — no callback ever wrote to it, so
		it read "Session state is empty until a project is selected" forever, including
		immediately after you had just selected a project. That is what made a working
		"Create project" look like a dead button.
		"""
		from ui.components.cards import status_badge

		store = project_store or {}
		active_id = store.get("active_id")
		project = (store.get("projects") or {}).get(active_id) if active_id else None

		if not project:
			return html.Div(
				[
					html.Div("No active project.", className="empty-state-sub"),
					html.Div("Create a project above, or load one of the recent ones.", className="mini-note"),
				]
			)

		samples = [s for s in (samples_data or {}).get("samples", []) if s.get("selected", True)]
		run_id = (run_config_data or {}).get("run_id")

		# The next action, stated plainly. A project with no samples cannot run, and the
		# UI never said so.
		if not samples:
			next_step = html.Div(
				[
					html.Span("Next step: "),
					html.A("add samples →", href="/samples", style={"color": "var(--blue)"}),
				],
				className="mini-note",
			)
		else:
			next_step = html.Div(
				[
					html.Span("Next step: "),
					html.A("set up and run →", href="/setup", style={"color": "var(--blue)"}),
				],
				className="mini-note",
			)

		return html.Div(
			[
				html.Div(
					[
						html.Div(project.get("project_name") or project.get("name") or active_id, style={"fontWeight": "600"}),
						status_badge("success" if samples else "pending", label="ready" if samples else "no samples"),
					],
					style={"display": "flex", "alignItems": "center", "gap": "10px"},
				),
				html.Div(
					[
						html.Div(f"{len(samples)} sample(s) selected", className="status-pill"),
						html.Div(f"Modo: {project.get('mode', 'Complete')}", className="status-pill"),
						*([html.Div(f"Run: {run_id}", className="status-pill")] if run_id else []),
					],
					className="chip-list",
					style={"marginTop": "8px"},
				),
				next_step,
			]
		)

	@app.callback(
		Output("recent-projects-list", "children"),
		Input("store-project", "data"),
	)
	def _render_recent(store):
		return _render_recent_projects(store)

	@app.callback(
		Output("stat-projects-count", "children"),
		Output("stat-samples-count", "children"),
		Output("stat-modes-count", "children"),
		Input("store-project", "data"),
		Input("store-samples", "data"),
	)
	def _update_project_stats(store, samples_data):
		projects = (store or {}).get("projects") or {}
		n_projects = len(projects)
		n_samples = len(samples_data or {})
		modes = {str(project.get("mode", "")).strip() for project in projects.values()}
		modes.discard("")
		n_modes = len(modes)

		def _card(value, label):
			return [html.Div(str(value), className="stat-value"), html.Div(label, className="stat-label")]

		return _card(n_projects, "Projects"), _card(n_samples, "Samples"), _card(n_modes, "Active modes")

	@app.callback(
		Output("store-run-config", "data", allow_duplicate=True),
		Output("store-qc-outputs", "data", allow_duplicate=True),
		Output("store-preprocess-outputs", "data", allow_duplicate=True),
		Output("store-analysis-outputs", "data", allow_duplicate=True),
		Output("store-databases-outputs", "data", allow_duplicate=True),
		Output("store-phylo-outputs", "data", allow_duplicate=True),
		Input("store-project", "data"),
		State("store-run-config", "data"),
		State("store-qc-outputs", "data"),
		State("store-preprocess-outputs", "data"),
		State("store-analysis-outputs", "data"),
		State("store-databases-outputs", "data"),
		State("store-phylo-outputs", "data"),
		prevent_initial_call="initial_duplicate",
	)
	def _restore_last_run(project_data, run_config, qc, preprocess, analysis, databases, phylo):
		"""Reload the project's newest run from disk when the session stores are empty.

		The stage stores are per-session, so restarting the app or opening a new tab
		wiped every result even though the run was still on disk — the Results page
		then had no stages to offer. Restoring is skipped whenever any store already
		holds something, so a run in progress is never clobbered.
		"""
		if any(bool(store) for store in (run_config, qc, preprocess, analysis, databases, phylo)):
			raise PreventUpdate
		restored = restore_run_state(project_data)
		if not restored:
			raise PreventUpdate
		return (
			restored["run_config"],
			restored["qc"],
			restored["preprocess"],
			restored["analysis"],
			restored["databases"],
			restored["phylo"],
		)

	@app.callback(
		Output("store-qc-outputs", "data", allow_duplicate=True),
		Output("store-preprocess-outputs", "data", allow_duplicate=True),
		Output("store-analysis-outputs", "data", allow_duplicate=True),
		Output("store-databases-outputs", "data", allow_duplicate=True),
		Output("store-phylo-outputs", "data", allow_duplicate=True),
		Input("store-run-config", "data"),
		State("store-qc-outputs", "data"),
		State("store-preprocess-outputs", "data"),
		State("store-analysis-outputs", "data"),
		State("store-databases-outputs", "data"),
		State("store-phylo-outputs", "data"),
		prevent_initial_call=True,
	)
	def _drop_records_from_another_run(run_config, qc, preprocess, analysis, databases, phylo):
		"""Discard stage records that a DIFFERENT run produced.

		The stage stores are per-session and the run-start callback writes only
		``store-run-config`` — nothing ever cleared these five, so starting a new run
		left the previous run's records in place with ``status: "success"``. Every stage
		page reads the store before it reads the disk, so the tabs stayed unlocked and
		showed the old run's tables while the new run was still executing.

		Checked here, once, rather than at the ~50 places that read a record: this
		fires whenever the active run changes, so every reader downstream sees either
		this run's record or nothing, without any of them learning about run identity.

		Not a clear-on-start: a record that matches is left alone (``no_update``), which
		is what lets a finished run's pages survive a re-render. Only a mismatch — or a
		record too old to say which run it came from — is dropped.
		"""
		run_id = (run_config or {}).get("run_id")
		if not run_id:
			# No active run to compare against. Dropping here would clobber the restore
			# during initial load, and there is no wrong data to protect against yet.
			raise PreventUpdate
		results = [
			{} if store and not record_belongs_to(store, run_id) else no_update
			for store in (qc, preprocess, analysis, databases, phylo)
		]
		if all(value is no_update for value in results):
			raise PreventUpdate
		return tuple(results)

	@app.callback(
		Output("qc-run-banner", "children"),
		Output("preprocess-run-banner", "children"),
		Output("analysis-run-banner", "children"),
		Output("databases-run-banner", "children"),
		Output("phylo-run-banner", "children"),
		Input("store-run-config", "data"),
		Input("store-qc-outputs", "data"),
		Input("store-preprocess-outputs", "data"),
		Input("store-analysis-outputs", "data"),
		Input("store-databases-outputs", "data"),
		Input("store-phylo-outputs", "data"),
	)
	def _stage_run_banners(run_config, qc, preprocess, analysis, databases, phylo):
		"""Say outright that a stage has not run in the ACTIVE run.

		Without this, dropping another run's records only turns the old tables into a
		blank page — which reads as the restore having broken, the exact confusion the
		Eixo B work exists to prevent. A stage with no record while a run is live gets a
		sentence saying so; a stage with no record and no live run says nothing, and the
		page's own empty state stands.
		"""
		from core.run_registry import load_run_state

		run_id = (run_config or {}).get("run_id")
		output_dir = (run_config or {}).get("output_dir")
		if not run_id or not output_dir:
			return "", "", "", "", ""
		state = load_run_state(output_dir)
		if state.get("run_id") != run_id or state.get("done") or state.get("interrupted"):
			# Not a live run: nothing is pending, so there is nothing to wait for.
			return "", "", "", "", ""

		stages = state.get("stages") or {}

		def _banner(store, stage):
			if record_belongs_to(store, run_id):
				return ""
			status = str(stages.get(stage) or "pending")
			detail = "running now" if status == "running" else "has not started yet"
			return html.Div(
				f"This stage has not run in this run yet ({detail}). "
				f"Results appear once the run finishes.",
				className="badge badge-amber",
				style={"marginBottom": "8px"},
			)

		return (
			_banner(qc, "qc"),
			_banner(preprocess, "preprocess"),
			_banner(analysis, "analysis"),
			_banner(databases, "databases"),
			_banner(phylo, "phylo"),
		)

	@app.callback(
		Output("page-project", "style"),
		Output("page-samples", "style"),
		Output("page-setup", "style"),
		Output("page-qc", "style"),
		Output("page-preprocess", "style"),
		Output("page-analysis", "style"),
		Output("page-databases", "style"),
		Output("page-phylo", "style"),
		Output("page-progress", "style"),
		Output("page-results", "style"),
		Input("url", "pathname"),
	)
	def _route(pathname):
		pathname = pathname or "/project"
		visible = {"display": "block"}
		hidden = {"display": "none"}
		# Exactly one page is visible at a time. The Results page used to also be
		# shown whenever the run had finished, so it appended itself to the bottom
		# of every other stage page.
		return (
			visible if pathname in {"/", "/project"} else hidden,
			visible if pathname == "/samples" else hidden,
			visible if pathname == "/setup" else hidden,
			visible if pathname == "/qc" else hidden,
			visible if pathname == "/preprocess" else hidden,
			visible if pathname == "/analysis" else hidden,
			visible if pathname == "/databases" else hidden,
			visible if pathname == "/phylo" else hidden,
			visible if pathname == "/progress" else hidden,
			visible if pathname == "/results" else hidden,
		)


	def _render_recent_projects(store):
		store = store or {}
		# Union the browser store with projects persisted on disk, so previously
		# created projects show up even in a fresh browser / after a restart.
		projects = dict(load_all_projects())
		projects.update(store.get("projects") or {})
		active_id = store.get("active_id")
		if not projects:
			return html.Div("No projects yet.", className="empty-state-sub")
		rows = []
		for proj in projects.values():
			proj_id = proj.get("id") or proj.get("project_name") or "project"
			rows.append(
				html.Tr(
					[
						html.Td(proj_id, className="mono-cell"),
						html.Td(proj.get("project_name") or proj.get("name") or ""),
						html.Td(html.Span(proj.get("mode", "Complete"), className="badge badge-blue")),
						html.Td(html.Span(proj.get("status", "active"), className="badge badge-green")),
						html.Td(html.Button("Load", id={"type": "btn-load-recent", "index": proj_id}, className="btn-parafish")),
					],
					className="is-done" if proj_id == active_id else "",
				)
			)
		return html.Table(
			[
				html.Thead(html.Tr([html.Th("ID"), html.Th("Name"), html.Th("Mode"), html.Th("Status"), html.Th("Action")])),
				html.Tbody(rows),
			],
			className="parafish-table",
		)

	@app.callback(
		Output("store-project", "data", allow_duplicate=True),
		Output("store-samples", "data", allow_duplicate=True),
		Output("url", "pathname", allow_duplicate=True),
		Input({"type": "btn-load-recent", "index": ALL}, "n_clicks"),
		State("store-project", "data"),
		prevent_initial_call=True,
	)
	def _load_recent(n_clicks_list, current_store):
		triggered = ctx.triggered_id
		if not isinstance(triggered, dict) or triggered.get("type") != "btn-load-recent":
			return no_update, no_update, no_update
		# The table is re-rendered on every change of store-project, and Dash fires this
		# callback for the new buttons with n_clicks None. Acting on that switched the
		# active project to whichever button came first (and replaced the samples with
		# that project's), so a project just created or selected was silently undone.
		if not _was_clicked(ctx.triggered):
			return no_update, no_update, no_update
		proj_id = triggered.get("index")
		store = current_store or {"active_id": None, "projects": {}}
		projects = dict(load_all_projects())
		projects.update(store.get("projects") or {})
		if proj_id not in projects:
			return no_update, no_update, no_update
		# Restore the project's persisted samples so a reopened project is usable,
		# not wiped (the store is per-browser and does not keep per-project samples).
		disk_project, disk_samples = load_project_bundle(proj_id)
		if disk_project:
			projects[proj_id] = disk_project
		store = {"active_id": proj_id, "projects": projects}
		return store, (disk_samples or {}), no_update

	@app.callback(
		Output("topbar", "children"),
		Input("store-project", "data"),
		Input("store-run-config", "data"),
	)
	def sync_topbar(store, run_config_data):
		active_project = resolve_project_data(store)
		project_name = active_project.get("project_name", "")
		# The run-config is the source of truth for the mode the pipeline actually runs
		# under: the Setup page writes ``workflow_mode`` there when "Run Pipeline" is
		# pressed. The project record's ``mode`` is only the creation-time default
		# ("Complete") and would leave the badge stuck on it even after a Special run.
		mode = ""
		if isinstance(run_config_data, dict):
			mode = run_config_data.get("workflow_mode") or ""
		mode = mode or active_project.get("mode") or "Complete"
		return render_topbar(project_name=project_name, mode=str(mode).capitalize(), docker_ok=False)

	@app.callback(
		Output("samples-active-project", "children"),
		Input("store-project", "data"),
		Input("store-samples", "data"),
	)
	def _samples_active_project(project_store, samples_data):
		return _active_project_banner(project_store, samples_data)

	@app.callback(
		Output("setup-active-project", "children"),
		Input("store-project", "data"),
		Input("store-samples", "data"),
	)
	def _setup_active_project(project_store, samples_data):
		return _active_project_banner(project_store, samples_data)

	@app.callback(
		Output("setup-run-target", "children"),
		Input("store-project", "data"),
		Input("store-samples", "data"),
	)
	def _setup_run_target(project_store, samples_data):
		return _run_target_note(project_store, samples_data)

	@app.callback(
		Output("sidebar", "children"),
		Input("url", "pathname"),
		Input("store-project", "data"),
		Input("store-samples", "data"),
		Input("store-qc-outputs", "data"),
		Input("store-preprocess-outputs", "data"),
		Input("store-analysis-outputs", "data"),
		Input("store-databases-outputs", "data"),
		Input("store-phylo-outputs", "data"),
		State("store-run-config", "data"),
	)
	def _sync_shell(pathname, project_data, samples_data, qc_outputs, preprocess_outputs, analysis_outputs, databases_outputs, phylo_outputs, run_config_data):
		pathname = pathname or "/project"
		step_outputs = {
			"qc": qc_outputs or {},
			"preprocess": preprocess_outputs or {},
			"analysis": analysis_outputs or {},
			"databases": databases_outputs or {},
			"phylo": phylo_outputs or {},
		}
		enabled_stages = (run_config_data or {}).get("enabled_stages", {}) if isinstance(run_config_data, dict) else {}
		step_states = _stage_states(pathname, project_data, samples_data, step_outputs)
		if enabled_stages:
			results_ready = _results_status(step_outputs, enabled_stages)
			step_states[-1] = "done" if results_ready == "done" else step_states[-1]
		return render_sidebar(active_path=pathname, step_states=step_states)