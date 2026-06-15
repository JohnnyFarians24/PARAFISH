"""Project page."""

from __future__ import annotations

from dash import dcc, html

def layout():
	return html.Div(
		[
			html.Div(
				[
					html.Div([html.Div("Project workspace", className="page-title"), html.Div("Create a run, load examples, and keep session metadata aligned.", className="page-sub")], className="page-header"),
					html.Div(
						[
							html.Div(id="stat-projects-count", className="stat-card"),
							html.Div(id="stat-samples-count", className="stat-card"),
							html.Div(id="stat-modes-count", className="stat-card"),
						],
						className="grid-3",
					),
				],
				className="card card-body fade-in",
			),
			html.Div(
				[
					html.Div(
						[
							html.Div([html.Div("Project setup", className="card-title")], className="card-header"),
							html.Div(
								[
										html.Div([html.Label("Project name"), dcc.Input(id="input-project-name", type="text", placeholder="Trout Gill Parasites")]),
										html.Div([html.Label("Notes"), dcc.Textarea(id="textarea-project-notes", placeholder="Capture protocol and run context here.")]),
										html.Div([html.Button("Create project", id="btn-create-project", className="btn-parafish btn-run"), html.Button("Load example", id="btn-load-example", className="btn-parafish")], className="control-row"),
								],
								className="card-body",
							),
						],
						className="card fade-in",
					),
					html.Div(
						[
							html.Div([html.Div("Recent projects", className="card-title")], className="card-header"),
									html.Div(id="recent-projects-list", className="card-body"),
							html.Div(
								id="project-session-state-info",
								children=[
									html.Div("Session state is empty until a project is selected.", className="empty-state-title"),
									html.Div("This area will track the active run, open logs, and the current execution branch.", className="empty-state-sub"),
								],
								className="empty-state",
							),
						],
						className="card fade-in",
					),
				],
				className="grid-2",
			),
		],
		className="page-stack",
	)