"""Databases page."""

from __future__ import annotations

from dash import dcc, html

from ui.components.collapsible import details
from ui.components.db_toggle import db_toggle


def feature_set_bar():
	"""Which feature set every table on this page refers to.

	OTUs and ASVs are shown one at a time, never side by side in the same table: a
	97% cluster and an exact variant are different units, and a row that silently
	mixes them is wrong in a way nobody catches downstream.
	"""
	return html.Div(
		[
			html.Div("Feature set", className="card-title"),
			html.Div(
				"Every table below refers to the chosen set. OTUs and ASVs are never mixed in "
				"the same table — to compare them, use the “OTU vs ASV” panel.",
				className="mini-note",
			),
			dcc.RadioItems(
				id="databases-feature-set",
				options=[{"label": " OTUs (VSEARCH)", "value": "otu"}, {"label": " ASVs (DADA2)", "value": "asv"}],
				value="otu",
				inline=True,
				className="chip-list",
				inputStyle={"marginRight": "6px"},
				labelStyle={"marginRight": "18px"},
			),
			html.Div(id="databases-feature-set-note", className="mini-note"),
		],
		className="card card-body",
	)


def layout():
	return html.Div(
		[
			html.Div([html.Div("Databases", className="page-title"), html.Div("Enable and inspect SILVA, PR2, NCBI and BOLD sources.", className="page-sub")], className="page-header"),
			html.Div(id="databases-run-banner"),
			html.Div([html.Button("Run Databases", id="btn-run-databases", className="btn-run")], className="control-row"),
			feature_set_bar(),
			html.Div(
				[
					html.Div([
						html.Div("Active databases", className="card-title"),
						html.Div(id="databases-summary", className="empty-state card"),
						db_toggle("SILVA 138.2 · SSU", True, "db-silva"),
						db_toggle("PR2 5.1.1", True, "db-pr2"),
						db_toggle("BOLD", True, "db-bold"),
						db_toggle("NCBI nt", False, "db-ncbi"),
						html.Hr(),
						html.Div("Custom Database", className="card-title"),
						dcc.Upload(
							id="upload-custom-db",
							children=html.Div(["Drag and drop or ", html.A("select a FASTA file"), " (.fasta, .fa, .fasta.gz, max 500 MB)"]),
							style={"borderStyle": "dashed", "borderWidth": "1px", "borderRadius": "4px", "padding": "12px", "textAlign": "center"},
							max_size=500 * 1024 * 1024,
							multiple=False,
						),
						html.Div(id="custom-db-upload-status"),
					], className="card card-body"),
					html.Div(
						[
							html.Div("Taxonomy consensus", className="card-title"),
							dcc.Loading(html.Div(id="databases-consensus"), type="dot"),
						],
						className="card card-body",
					),
				],
				className="grid-2",
			),
			html.Div(
				[
					html.Div("Comparison with custom database (BLAST)", className="card-title"),
					html.Div(
						"After running Databases with a custom reference uploaded, each OTU/ASV is BLASTed against it. "
						"Matches show the reference species; OTUs with no hit are flagged None (red). Filter by minimum "
						"identity or species to compare dozens to hundreds of samples.",
						className="mini-note",
					),
					html.Div(
						[
							html.Div(
								[
									html.Label("Minimum identity (%)", className="mini-note"),
									dcc.Slider(
										id="custom-cmp-minid",
										min=0,
										max=100,
										step=1,
										value=0,
										marks={0: "0", 90: "90", 97: "97", 100: "100"},
										tooltip={"placement": "bottom", "always_visible": False},
									),
								],
								style={"flex": "1", "minWidth": "220px"},
							),
							html.Div(
								[
									html.Label("Filter by species", className="mini-note"),
									dcc.Input(
										id="custom-cmp-species",
										type="text",
										placeholder="e.g. Diplostomum",
										debounce=True,
										style={"width": "100%"},
									),
								],
								style={"flex": "1", "minWidth": "220px"},
							),
						],
						style={"display": "flex", "gap": "16px", "flexWrap": "wrap", "marginBottom": "12px"},
					),
					dcc.Loading(html.Div(id="custom-comparison"), type="dot"),
				],
				className="card card-body",
			),
			html.Div(
				[
					html.Div("Best hits per feature", className="card-title"),
					html.Div(
						"Every candidate from every source, ranked by % identity then by % coverage. "
						"hit_rank = 1 is the best. The runners-up are kept on purpose: the top BLAST hit is "
						"often an unnamed environmental sequence, while the 2nd or 3rd is the species that "
						"actually matters — that choice is yours.",
						className="mini-note",
					),
					html.Div(
						[
							html.Div(
								[
									html.Label("% min identity", className="mini-note"),
									dcc.Slider(
										id="hits-min-identity",
										min=0,
										max=100,
										step=1,
										value=0,
										marks={0: "0", 90: "90", 97: "97", 100: "100"},
										tooltip={"placement": "bottom", "always_visible": False},
									),
								],
								style={"flex": "1", "minWidth": "220px"},
							),
							html.Div(
								[
									html.Label("% min coverage", className="mini-note"),
									dcc.Slider(
										id="hits-min-coverage",
										min=0,
										max=100,
										step=1,
										value=0,
										marks={0: "0", 80: "80", 100: "100"},
										tooltip={"placement": "bottom", "always_visible": False},
									),
								],
								style={"flex": "1", "minWidth": "220px"},
							),
							html.Div(
								[
									html.Label("How many hits per feature", className="mini-note"),
									dcc.Dropdown(
										id="hits-top-n",
										options=[
											{"label": "Best only (1)", "value": 1},
											{"label": "Top 3", "value": 3},
											{"label": "Top 5", "value": 5},
											{"label": "All", "value": 0},
										],
										value=3,
										clearable=False,
									),
								],
								style={"flex": "1", "minWidth": "200px"},
							),
						],
						style={"display": "flex", "gap": "16px", "flexWrap": "wrap", "marginBottom": "12px"},
					),
					html.Div(
						"Tick rows to download only those; with none ticked, the download holds every "
						"row the table shows (after its filters).",
						className="mini-note",
					),
					html.Div(
						[
							html.Button("Download table (CSV)", id="btn-download-hits", className="btn-parafish"),
							dcc.Download(id="download-hits-csv"),
						],
						className="control-row",
					),
					dcc.Store(id="hits-table-state", storage_type="memory"),
					dcc.Loading(html.Div(id="databases-hits"), type="dot"),
				],
				className="card card-body",
			),
			# The deliverable: one row per feature, the consensus call beside the best hit
			# of EACH evidence family. Rendered straight from taxonomy_summary.csv (written
			# by the pipeline) — the same file the download hands over, so the screen and
			# the exported table can never diverge.
			html.Div(
				[
					html.Div("Per-feature summary", className="card-title"),
					html.Div(
						"One row per feature: the consensus call (LCA), and the best hit of each evidence "
						"family — alignment (identity) and SINTAX (bootstrap) — in SEPARATE columns. There is "
						"no single winning score: when the two families disagree, you see both and decide.",
						className="mini-note",
					),
					html.Div(
						"Sorted by reads. Tick rows to download only those; with none ticked and no "
						"filter, the download is the pipeline's own taxonomy_summary.csv.",
						className="mini-note",
					),
					html.Div(
						[
							html.Button("Download summary (CSV)", id="btn-download-summary", className="btn-parafish"),
							dcc.Download(id="download-summary-csv"),
						],
						className="control-row",
					),
					dcc.Store(id="feature-summary-table-state", storage_type="memory"),
					dcc.Loading(html.Div(id="databases-feature-summary"), type="dot"),
				],
				className="card card-body",
			),
			html.Div(
				[
					html.Div("OTU vs ASV", className="card-title"),
					html.Div(
						"The two sets side by side, in separate tables. The bridge between them is the TAXON, "
						"never the feature id: a 97% OTU typically contains several ASVs, so merging them into "
						"one table would produce meaningless counts.",
						className="mini-note",
					),
					dcc.Loading(html.Div(id="databases-feature-comparison"), type="dot"),
				],
				className="card card-body",
			),
			# Layer 1 — the hero. "How deep did the classification get?" is the question
			# this whole page exists to answer, so it is never behind a click.
			html.Div(
				[
					html.Div("Taxonomic resolution", className="card-title"),
					html.Div("How deep the consensus reached, counted per feature and per read. They rarely coincide.", className="mini-note"),
					dcc.Loading(html.Div(id="databases-resolution"), type="dot"),
				],
				className="card card-body",
			),
			# Layer 2 — audit detail. Open it when you go looking for it.
			details(
				"Assignment confidence and reliability",
				dcc.Loading(html.Div(id="databases-confidence"), type="dot"),
				subtitle="identity and coverage come from the sources that align; SINTAX gives a bootstrap",
			),
			details("Conflicts between sources", html.Div(id="databases-conflicts")),
			details("Negative controls and contamination", dcc.Loading(html.Div(id="databases-contaminants"), type="dot")),
			details("Taxon prevalence", dcc.Loading(html.Div(id="databases-prevalence"), type="dot")),
			details("Taxonomic composition and sample comparison", html.Div(id="databases-visuals")),
		],
		className="page-stack",
	)