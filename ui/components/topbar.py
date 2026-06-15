"""Topbar component."""

from __future__ import annotations

from dash import html


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
					html.Span("UI shell", className="badge badge-blue"),
					html.Span("Docker ready" if docker_ok else "Docker offline", className=f"badge {'badge-green' if docker_ok else 'badge-red'}"),
					html.Span("Project", className="logo-subtext"),
					html.Span(project_name or "No project active", className="badge badge-sand"),
				],
				className="topbar-right",
			),
		],
		className="topbar-inner",
	)