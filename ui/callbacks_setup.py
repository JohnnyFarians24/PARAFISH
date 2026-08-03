"""Callbacks for the unified pipeline setup page."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from core.config import MIN_UNIQUE_SIZE
from core.models import PhyloMode
from ui.callbacks_databases import SUPPORTED_CUSTOM_DB_SUFFIXES, _is_supported_custom_db_filename
from ui.components.sample_table import render_sample_table
from ui.runtime_state import build_run_config, decode_upload_content, enum_from_value, resolve_project_data, serialize_run_config, slugify


def _samples_list(samples_store):
    if isinstance(samples_store, list):
        return samples_store
    if isinstance(samples_store, dict):
        if isinstance(samples_store.get("samples"), list):
            return samples_store["samples"]
        return list(samples_store.values()) if samples_store else []
    return []


def register_callbacks(app: Any) -> None:
    """Register setup page callbacks."""

    @app.callback(
        Output("setup-enable-qc", "value"),
        Output("setup-enable-preprocess", "value"),
        Output("setup-enable-analysis", "value"),
        Output("setup-enable-databases", "value"),
        Output("setup-enable-phylo", "value"),
        Output("setup-enable-qc", "disabled"),
        Output("setup-enable-preprocess", "disabled"),
        Output("setup-enable-analysis", "disabled"),
        Output("setup-enable-databases", "disabled"),
        Output("setup-enable-phylo", "disabled"),
        Output("setup-run-status", "children"),
        Input("setup-workflow-mode", "value"),
        State("setup-enable-qc", "value"),
        State("setup-enable-preprocess", "value"),
        State("setup-enable-analysis", "value"),
        State("setup-enable-databases", "value"),
        State("setup-enable-phylo", "value"),
    )
    def _validate_workflow_mode(mode, qc, preprocess, analysis, databases, phylo):
        values = [bool(qc), bool(preprocess), bool(analysis), bool(databases), bool(phylo)]

        if mode == "complete":
            values = [True, True, True, True, True]
            disabled = [True, True, True, True, True]
            return (*values, *disabled, "")

        if mode == "one":
            enabled_indexes = [idx for idx, val in enumerate(values) if val]
            if len(enabled_indexes) == 0:
                values[0] = True
            elif len(enabled_indexes) > 1:
                keep = enabled_indexes[0]
                for idx in range(len(values)):
                    values[idx] = idx == keep
            disabled = [False, False, False, False, False]
            return (*values, *disabled, "")

        # special mode
        enabled_count = sum(1 for v in values if v)
        warning = ""
        if enabled_count < 2:
            warning = html.Div("Special mode requires at least 2 enabled stages.", className="status-pill badge badge-amber")
        disabled = [False, False, False, False, False]
        return (*values, *disabled, warning)

    @app.callback(
        Output("setup-samples-table", "children"),
        Output("setup-sample-count", "children"),
        Input("store-samples", "data"),
        Input("setup-primer-preset", "value"),
    )
    def _render_samples_table(samples_store, primer_preset):
        from databases.reference_filter import resolve_run_marker

        samples = _samples_list(samples_store)
        if not samples:
            return html.Div("No samples uploaded.", className="empty-state-sub"), "0"
        selected = sum(1 for s in samples if s.get("selected", True))
        # The column showed the file-name guess, which the run stopped acting on once
        # the preset decides the marker. On the COI/Leray run it read "NONE" while the
        # Databases stage went on to treat the run as OTHER and cut SILVA with COI
        # primers — a readout that agreed with nothing the pipeline did.
        run_marker, _ = resolve_run_marker(
            (s.get("marker") for s in samples), {"primer_preset": primer_preset}
        )
        return (
            render_sample_table(samples, selectable=True, run_marker=run_marker),
            f"{selected} selected / {len(samples)} total",
        )

    @app.callback(
        Output("setup-db-recommendation", "children"),
        Input("store-samples", "data"),
        Input("setup-primer-preset", "value"),
    )
    def _db_recommendation(samples_data, primer_preset):
        from databases.reference_filter import resolve_run_marker

        samples = (samples_data or {}).get("samples") or []
        if not samples:
            return "Upload samples to see database recommendations."
        run_marker, note = resolve_run_marker(
            (s.get("marker") for s in samples), {"primer_preset": primer_preset}
        )
        recommendation = {
            "COI": "COI — BOLD (primary), NCBI (secondary). SILVA/PR2 do not index COI.",
            "18S": "18S — PR2 (primary), SILVA (secondary).",
            "ITS": "ITS — SILVA (primary), PR2 (secondary).",
            "28S": "28S — SILVA (primary).",
        }.get(run_marker)
        if not recommendation:
            # Said "Upload samples..." here before, with the samples sitting on screen:
            # every marker resolved to OTHER, OTHER was discarded, and an empty set was
            # reported as an empty upload.
            return (
                "Marker not identified from the primer preset or the file names — "
                "pick a primer preset, or enable all sources."
            )
        return f"{recommendation} | {note}" if note else recommendation


    @app.callback(
        Output("setup-custom-db-status", "children"),
        Output("store-custom-db", "data"),
        Input("setup-upload-custom-db", "contents"),
        State("setup-upload-custom-db", "filename"),
        State("store-project", "data"),
        prevent_initial_call=True,
    )
    def _handle_setup_custom_db(contents, filename, project_data):
        from databases.custom_db import CustomDBConfig

        if not contents or not filename:
            raise PreventUpdate
        if not _is_supported_custom_db_filename(filename):
            return html.Div(
                f"Unsupported file. Use {', '.join(SUPPORTED_CUSTOM_DB_SUFFIXES)}.",
                className="badge badge-amber",
            ), no_update
        project = resolve_project_data(project_data)
        project_name = project.get("project_name") or "project"
        db_dir = Path("projects") / slugify(project_name) / "custom_dbs"
        db_dir.mkdir(parents=True, exist_ok=True)
        dest = db_dir / Path(filename).name
        try:
            dest.write_bytes(decode_upload_content(contents))
        except Exception as exc:
            return html.Div(f"Could not read the upload: {exc}", className="badge badge-amber"), no_update
        config = CustomDBConfig(dest)
        valid, error = config.validate()
        if not valid:
            dest.unlink(missing_ok=True)
            return html.Div(f"Invalid file: {error}", className="badge badge-amber"), no_update
        if config.kind == "blast":
            note = "Sequences will be compared against it with BLAST (blastn)."
        else:
            note = "Headers carry taxonomy (;tax=); it will be used for SINTAX classification."
        status = html.Div(
            [
                html.Div(f"✓ {filename} loaded", style={"color": "#3fb950", "fontWeight": "700"}),
                html.Div(note, className="empty-state-sub"),
                html.Div("Run the pipeline (or the Databases stage) to see the comparison.", className="empty-state-sub"),
            ],
            style={"padding": "8px 10px", "borderRadius": "8px", "background": "rgba(63,185,80,0.10)", "border": "1px solid rgba(63,185,80,0.35)"},
        )
        return status, {"path": str(dest), "kind": config.kind, "filename": filename}

    @app.callback(
        Output("setup-preflight", "children"),
        Input("store-samples", "data"),
        Input("setup-db-sources", "value"),
        Input("setup-analysis-method", "value"),
        Input("setup-primer-preset", "value"),
        Input("setup-primer-fwd", "value"),
        Input("setup-primer-rev", "value"),
    )
    def _render_preflight(samples_data, db_sources, analysis_method, primer_preset, primer_fwd, primer_rev):
        """What this machine can and cannot do with the run as configured.

        The OOM killer gives no warning and leaves no error of its own: it SIGKILLs
        whichever process is biggest, which on this pipeline is either a VSEARCH
        SINTAX search or the Dash server itself (that is what surfaces in the browser
        as "Callback failed: the server did not respond"). Saying so before the run
        is the only way the user can act on it.

        The primers are inputs because they change the verdict: with them, SILVA/PR2 are
        cut to the amplicon and fit; without them, they are skipped.
        """
        from analysis.vsearch_runner import PRIMER_PRESETS
        from core.config import get_db_paths
        from core.preflight import available_memory_bytes, max_samples_for_memory, preflight_notes, total_memory_bytes
        from databases.reference_filter import marker_filter_applies, resolve_run_marker

        samples = [s for s in (samples_data or {}).get("samples", []) if s.get("selected", True)]
        n_samples = len(samples)
        available = available_memory_bytes()
        total = total_memory_bytes()
        limit = max_samples_for_memory(available)

        # Same resolution the run uses, preset first: a preflight that predicted OTHER
        # while the run resolved COI would advise on a filter that never happens.
        marker, _ = resolve_run_marker(
            (s.get("marker") for s in samples), {"primer_preset": primer_preset}
        )
        forward = str(primer_fwd or "").strip()
        reverse = str(primer_rev or "").strip()
        if not (forward or reverse):
            # Mirror what the run will actually do, including "auto" resolving from the
            # marker — otherwise the preflight advises on primers the run will not use.
            from preprocess.cutadapt_runner import AUTO_PRESET, MARKER_DEFAULT_PRESET

            preset = AUTO_PRESET if primer_preset is None else str(primer_preset).strip()
            if preset == AUTO_PRESET:
                preset = MARKER_DEFAULT_PRESET.get(marker, "")
            if preset in PRIMER_PRESETS:
                forward, reverse = PRIMER_PRESETS[preset]
        marker_filter = marker if marker_filter_applies(marker, (forward, reverse)) else ""

        from databases.taxonomy_runner import SOURCE_SUPPORTED_MARKERS, _marker_incompatible

        inapplicable = [key for key in SOURCE_SUPPORTED_MARKERS if _marker_incompatible(key, marker)]
        notes = preflight_notes(
            db_sources or [],
            get_db_paths(),
            n_samples,
            analysis_method or "",
            marker_filter=marker_filter,
            inapplicable=inapplicable,
        )

        pills = [
            html.Div(f"{n_samples} sample(s) selected", className="status-pill"),
            html.Div(f"Recommended maximum: {limit}", className="status-pill"),
        ]
        if total and available:
            pills.append(html.Div(f"RAM: {available / 1e9:.1f} GB livre de {total / 1e9:.1f} GB", className="status-pill"))

        children = [
            html.Div("Machine resources", className="card-title"),
            html.Div(pills, className="chip-list"),
        ]
        for note in notes:
            children.append(
                html.Div(
                    note,
                    className="mini-note",
                    style={
                        "padding": "8px 10px",
                        "borderRadius": "8px",
                        "background": "rgba(210,153,34,0.10)",
                        "border": "1px solid rgba(210,153,34,0.35)",
                        "marginTop": "8px",
                    },
                )
            )
        if not notes:
            children.append(html.Div("No memory problems expected for this configuration.", className="mini-note"))
        return html.Div(children, className="card card-body", style={"marginBottom": "12px"})

    @app.callback(
        Output("vsearch-params-panel", "style"),
        Output("dada2-params-panel", "style"),
        Input("setup-analysis-method", "value"),
    )
    def _toggle_analysis_panels(method):
        if not method:
            return no_update, no_update
        if method == "vsearch":
            return {"display": "block"}, {"display": "none"}
        if method == "dada2":
            return {"display": "none"}, {"display": "block"}
        return {"display": "block"}, {"display": "block"}


    @app.callback(
        Output("store-run-config", "data", allow_duplicate=True),
        Output("store-project", "data", allow_duplicate=True),
        Output("url", "pathname", allow_duplicate=True),
        Output("setup-run-status", "children", allow_duplicate=True),
        Output("store-samples", "data", allow_duplicate=True),
        Output("progress-interval", "disabled", allow_duplicate=True),
        Input("btn-run-pipeline", "n_clicks"),
        State("store-project", "data"),
        State("store-samples", "data"),
        State("setup-workflow-mode", "value"),
        State("setup-enable-qc", "value"),
        State("setup-enable-preprocess", "value"),
        State("setup-enable-analysis", "value"),
        State("setup-enable-databases", "value"),
        State("setup-enable-phylo", "value"),
        State("setup-analysis-method", "value"),
        State("setup-vsearch-id", "value"),
        State("setup-vsearch-maxee", "value"),
        State("setup-vsearch-minlen", "value"),
        State("setup-vsearch-sintax", "value"),
        State("setup-vsearch-minsize", "value"),
        State("setup-primer-preset", "value"),
        State("setup-primer-fwd", "value"),
        State("setup-primer-rev", "value"),
        State("setup-dada2-truncF", "value"),
        State("setup-dada2-truncR", "value"),
        State("setup-dada2-maxeeF", "value"),
        State("setup-dada2-maxeeR", "value"),
        State("setup-pear-min-overlap", "value"),
        State("setup-pear-min-length", "value"),
        State("setup-pear-pvalue", "value"),
        State("setup-trim-leading", "value"),
        State("setup-trim-trailing", "value"),
        State("setup-trim-slidingwindow", "value"),
        State("setup-trim-minlen", "value"),
        State("setup-db-sources", "value"),
        State("setup-db-consensus", "value"),
        State("setup-species-threshold", "value"),
        State("setup-sintax-cutoff", "value"),
        State("setup-ncbi-evalue", "value"),
        State("setup-ncbi-maxhits", "value"),
        State("setup-bold-db", "value"),
        State("setup-bold-mode", "value"),
        State("setup-phylo-mode", "value"),
        State("setup-raxml-bootstraps", "value"),
        State("setup-threads", "value"),
        State("setup-workers", "value"),
        State("store-custom-db", "data"),
        prevent_initial_call=True,
    )
    def _run_pipeline(n_clicks,
                      project_data,
                      samples_data,
                      workflow_mode,
                      enable_qc,
                      enable_preprocess,
                      enable_analysis,
                      enable_databases,
                      enable_phylo,
                      analysis_method,
                      v_id,
                      v_maxee,
                      v_minlen,
                      v_sintax,
                      v_minsize,
                      primer_preset,
                      primer_fwd,
                      primer_rev,
                      d_truncF,
                      d_truncR,
                      d_maxeeF,
                      d_maxeeR,
                      pear_min_overlap,
                      pear_min_length,
                      pear_pvalue,
                      trim_leading,
                      trim_trailing,
                      trim_slidingwindow,
                      trim_minlen,
                      db_sources,
                      db_consensus,
                      species_threshold,
                        sintax_cutoff,
                        ncbi_evalue,
                        ncbi_maxhits,
                        bold_db,
                        bold_mode,
                      phylo_mode,
                      raxml_bootstraps,
                      threads,
                      workers,
                      custom_db_store,
    ):
        import sys
        import traceback

        try:
            return _run_pipeline_inner(
                n_clicks,
                project_data,
                samples_data,
                workflow_mode,
                enable_qc,
                enable_preprocess,
                enable_analysis,
                enable_databases,
                enable_phylo,
                analysis_method,
                v_id,
                v_maxee,
                v_minlen,
                v_sintax,
                v_minsize,
                primer_preset,
                primer_fwd,
                primer_rev,
                d_truncF,
                d_truncR,
                d_maxeeF,
                d_maxeeR,
                pear_min_overlap,
                pear_min_length,
                pear_pvalue,
                trim_leading,
                trim_trailing,
                trim_slidingwindow,
                trim_minlen,
                db_sources,
                db_consensus,
                species_threshold,
                sintax_cutoff,
                ncbi_evalue,
                ncbi_maxhits,
                bold_db,
                bold_mode,
                phylo_mode,
                raxml_bootstraps,
                threads,
                workers,
                custom_db_store,
            )
        except PreventUpdate:
            raise
        except Exception as exc:
            print(f"CALLBACK ERROR: {exc}", file=sys.stderr, flush=True)
            traceback.print_exc(file=sys.stderr)
            return no_update, no_update, no_update, str(exc), no_update, no_update

    def _run_pipeline_inner(n_clicks,
                            project_data,
                            samples_data,
                            workflow_mode,
                            enable_qc,
                            enable_preprocess,
                            enable_analysis,
                            enable_databases,
                            enable_phylo,
                            analysis_method,
                            v_id,
                            v_maxee,
                            v_minlen,
                            v_sintax,
                            v_minsize,
                            primer_preset,
                            primer_fwd,
                            primer_rev,
                            d_truncF,
                            d_truncR,
                            d_maxeeF,
                            d_maxeeR,
                            pear_min_overlap,
                            pear_min_length,
                            pear_pvalue,
                            trim_leading,
                            trim_trailing,
                            trim_slidingwindow,
                            trim_minlen,
                            db_sources,
                            db_consensus,
                            species_threshold,
                            sintax_cutoff,
                            ncbi_evalue,
                            ncbi_maxhits,
                            bold_db,
                            bold_mode,
                            phylo_mode,
                            raxml_bootstraps,
                            threads,
                            workers,
                            custom_db_store=None,
    ):
        if not n_clicks:
            raise PreventUpdate

        def _is_enabled(value):
            if isinstance(value, (list, tuple, set)):
                return bool(value)
            return bool(value)

        enabled_stages = {
            "qc": _is_enabled(enable_qc),
            "preprocess": _is_enabled(enable_preprocess),
            "analysis": _is_enabled(enable_analysis),
            "databases": _is_enabled(enable_databases),
            "phylo": _is_enabled(enable_phylo),
        }
        

        if workflow_mode == "special" and sum(1 for v in enabled_stages.values() if v) < 2:
            warning = html.Div("Special mode requires at least 2 enabled stages.", className="status-pill badge badge-amber")
            return no_update, no_update, no_update, warning, no_update, no_update

        all_samples = _samples_list(samples_data)
        selected_samples = [s for s in all_samples if s.get("selected", True)]
        if not selected_samples:
            warning = html.Div(
                "Select at least one sample before running the pipeline.",
                className="status-pill badge badge-amber",
                )
            return no_update, no_update, no_update, warning, no_update, no_update
        samples_data_filtered = {**(samples_data or {}), "samples": selected_samples}
        run_config = build_run_config(project_data, samples_data_filtered)
        run_config.enabled_stages = enabled_stages

        # Preprocess params
        run_config.pear_params.update({
            "min_overlap": pear_min_overlap,
            "min_length": pear_min_length,
            "p_value": pear_pvalue,
        })
        window_size, window_quality = run_config.trimmomatic_params.get("window_size", 4), run_config.trimmomatic_params.get("window_quality", 20)
        if trim_slidingwindow and ":" in str(trim_slidingwindow):
            raw_size, raw_quality = str(trim_slidingwindow).split(":", 1)
            try:
                window_size, window_quality = int(raw_size), int(raw_quality)
            except ValueError:
                pass
        run_config.trimmomatic_params.update({
            "leading": trim_leading,
            "trailing": trim_trailing,
            "window_size": window_size,
            "window_quality": window_quality,
            "minlen": trim_minlen,
        })

        # VSEARCH params
        run_config.vsearch_params.update({
            "id": v_id,
            "fastq_maxee": v_maxee,
            "fastq_minlen": v_minlen,
            # Empty control -> the configured default, not a literal typed here. This
            # line's hardcoded 1 disagreed with the runner's hardcoded 2, so the policy
            # depended on the entry point.
            "min_unique_size": int(v_minsize) if v_minsize else MIN_UNIQUE_SIZE,
            # PCR primers: Trimmomatic strips sequencing adapters, not these. When set,
            # the Preprocess stage removes them with cutadapt (sequence-aware, before
            # trimming); if cutadapt is absent the Analysis stage falls back to a
            # fixed-length strip. Left in entirely, the primer's degenerate positions
            # would be clustered as if they were biology.
            "primer_preset": primer_preset or "",
            "primer_fwd": (primer_fwd or "").strip(),
            "primer_rev": (primer_rev or "").strip(),
            "threads": int(threads) if threads else run_config.vsearch_params.get("threads", 8),
        })

        # DADA2 params
        run_config.dada2_params.update({
            "truncLen_F": d_truncF,
            "truncLen_R": d_truncR,
            "maxEE_F": d_maxeeF,
            "maxEE_R": d_maxeeR,
        })

        # Databases
        run_config.vsearch_params["sintax_cutoff"] = float(sintax_cutoff) if sintax_cutoff else 0.8
        run_config.bold_params["db"] = int(bold_db) if bold_db else 1
        run_config.bold_params["mode"] = int(bold_mode) if bold_mode else 1
        run_config.vsearch_params["ncbi_evalue"] = ncbi_evalue or "1e-10"
        run_config.vsearch_params["ncbi_maxhits"] = int(ncbi_maxhits) if ncbi_maxhits else 5
        run_config.phylo_mode = enum_from_value(PhyloMode, phylo_mode, run_config.phylo_mode)
        run_config.phylo_params["bs_trees"] = int(raxml_bootstraps) if raxml_bootstraps else 10

        # Enabled steps as simple flags
        run_config.pear_params["enabled"] = bool(enable_preprocess)
        run_config.trimmomatic_params["enabled"] = bool(enable_preprocess)
        run_config.vsearch_params["enabled"] = bool(enable_analysis)
        run_config.dada2_params["enabled"] = bool(enable_analysis)
        run_config.bold_params["enabled"] = bool(enable_databases)
        run_config.phylo_params["enabled"] = bool(enable_phylo)
        
        
        # Additional DB and worker settings
        run_config.vsearch_params["db_consensus"] = db_consensus
        run_config.vsearch_params["species_threshold"] = species_threshold
        run_config.analysis_method = analysis_method or run_config.analysis_method
        run_config.vsearch_params["workers"] = int(workers) if workers else None

        # Custom (user-uploaded) reference database — carried from the upload store
        # so "Run Pipeline" compares each OTU/ASV against it in the Databases stage.
        custom_db_store = custom_db_store or {}
        custom_db_path = custom_db_store.get("path")
        if custom_db_path and Path(custom_db_path).exists():
            run_config.vsearch_params["custom_db_path"] = custom_db_path
            run_config.vsearch_params["custom_db_kind"] = custom_db_store.get("kind") or ""
        

        project_store = dict(project_data or {})
        projects = dict(project_store.get("projects", {}))
        active_id = project_store.get("active_id")
        if active_id and active_id in projects:
            projects[active_id] = dict(projects[active_id])
            projects[active_id]["mode"] = workflow_mode
            projects[active_id]["samples"] = selected_samples
            project_store = {"active_id": active_id, "projects": projects}
        
        config_data = serialize_run_config(run_config)
        config_data["db_sources"] = list(db_sources or [])
        config_data["workflow_mode"] = workflow_mode
        # "Run Pipeline" configures the run and navigates to Progress, but does NOT
        # launch it. Starting a pipeline consumes hours of remote NCBI/BOLD queries and,
        # on this box (3 GB, active OOM killer), a second concurrent run is a guaranteed
        # OOM — so a run must begin only on an explicit click of the Progress page's
        # "Start Pipeline" button, the single deliberate trigger. run_id and output_dir
        # travel in config_data, so that click launches exactly this configuration. The
        # poller stays disabled until then; "Start Pipeline" re-arms it.
        return config_data, project_store, "/progress", "", samples_data, True