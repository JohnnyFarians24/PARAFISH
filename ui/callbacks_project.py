"""Project and page routing callbacks."""

from __future__ import annotations

import uuid
from datetime import datetime, timezone
from typing import Any

from dash import Input, Output, State, ctx, html, no_update
from dash import ALL

from ui.pipeline_state import get_stage_record, stage_status
from ui.components.sidebar import render_sidebar
from ui.components.topbar import render_topbar
from ui.runtime_state import resolve_project_data

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
		return store, {}, no_update

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
		return store, {}, no_update

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
		State("store-run-config", "data"),
		State("store-qc-outputs", "data"),
		State("store-preprocess-outputs", "data"),
		State("store-analysis-outputs", "data"),
		State("store-databases-outputs", "data"),
		State("store-phylo-outputs", "data"),
	)
	def _route(pathname, run_config_data, qc_outputs, preprocess_outputs, analysis_outputs, databases_outputs, phylo_outputs):
		print("ROUTE PATHNAME:", pathname, flush=True)
		pathname = pathname or "/project"
		visible = {"display": "block"}
		hidden = {"display": "none"}
		step_outputs = {
			"qc": qc_outputs or {},
			"preprocess": preprocess_outputs or {},
			"analysis": analysis_outputs or {},
			"databases": databases_outputs or {},
			"phylo": phylo_outputs or {},
		}
		enabled_stages = (run_config_data or {}).get("enabled_stages", {}) if isinstance(run_config_data, dict) else {}
		results_state = _results_status(step_outputs, enabled_stages)
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
			visible if pathname == "/results" or results_state == "done" else hidden,
		)


	def _render_recent_projects(store):
		store = store or {}
		projects = store.get("projects") or {}
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
		proj_id = triggered.get("index")
		store = current_store or {"active_id": None, "projects": {}}
		if proj_id in (store.get("projects") or {}):
			store = {"active_id": proj_id, "projects": dict(store.get("projects", {}))}
			return store, {}, no_update
		return no_update, no_update, no_update

	@app.callback(
		Output("topbar", "children"),
		Input("store-project", "data"),
	)
	def sync_topbar(store):
		active_project = resolve_project_data(store)
		project_name = active_project.get("project_name", "")
		mode = active_project.get("mode", "Complete")
		return render_topbar(project_name=project_name, mode=mode, docker_ok=False)

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