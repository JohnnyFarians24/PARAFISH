"""Topbar component."""

from __future__ import annotations

from dash import html


def _tools_badge():
	"""Live badge reflecting which external tools are actually on PATH."""
	try:
		from core.config import check_core_tools

		health = check_core_tools()
	except Exception:
		return html.Span("Tools —", className="badge badge-sand", title="Could not check tools")
	available, total = health["available"], health["total"]
	missing = [name for name, ok in health["tools"].items() if not ok]
	if health["ok"]:
		return html.Span("Tools ready", className="badge badge-green", title="All bioinformatics tools found on PATH")
	label = f"Tools {available}/{total}"
	klass = "badge-amber" if available > 0 else "badge-red"
	return html.Span(label, className=f"badge {klass}", title="Missing: " + ", ".join(missing))


def render_topbar(project_name: str = "", mode: str = "Complete", docker_ok: bool = False):
	return html.Div(
		[
			html.Div(
				[
					html.Img(src="/assets/logo.png", alt="Parafish", className="app-logo"),
					html.Div(
						[
							html.Div("Parafish", className="logo-text"),
							html.Div("Modular sequencing workflow", className="logo-subtext"),
						],
						className="brand-stack",
					),
				],
				className="brand",
			),
			html.Div(
				[
					html.Span(f"{mode} mode", className="badge badge-amber"),
					_tools_badge(),
					html.Span("Project", className="logo-subtext"),
					html.Span(project_name or "No project active", className="badge badge-sand"),
				],
				className="topbar-right",
			),
		],
		className="topbar-inner",
	)