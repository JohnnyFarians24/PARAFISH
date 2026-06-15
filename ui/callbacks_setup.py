"""Callbacks for the unified pipeline setup page."""

from __future__ import annotations

from typing import Any

from dash import Input, Output, State, html, no_update
from dash.exceptions import PreventUpdate

from core.models import PhyloMode
from ui.components.sample_table import render_sample_table
from ui.runtime_state import build_run_config, enum_from_value, serialize_run_config


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
    )
    def _render_samples_table(samples_store):
        samples = _samples_list(samples_store)
        if not samples:
            return html.Div("No samples uploaded.", className="empty-state-sub"), "0"
        selected = sum(1 for s in samples if s.get("selected", True))
        return render_sample_table(samples, selectable=True), f"{selected} selected / {len(samples)} total"

    @app.callback(
        Output("setup-db-recommendation", "children"),
        Input("store-samples", "data"),
    )
    def _db_recommendation(samples_data):
        samples = (samples_data or {}).get("samples") or []
        markers = {str(s.get("marker") or "OTHER").upper() for s in samples}
        markers.discard("OTHER")
        if not markers:
            return "Upload samples to see database recommendations."
        lines = []
        if "COI" in markers:
            lines.append("COI detected — BOLD (primary), NCBI (secondary).")
        if "18S" in markers:
            lines.append("18S detected — PR2 (primary), SILVA (secondary).")
        if "ITS" in markers:
            lines.append("ITS detected — SILVA (primary), PR2 (secondary).")
        if "28S" in markers:
            lines.append("28S detected — SILVA (primary).")
        unknown = markers - {"COI", "18S", "ITS", "28S"}
        if unknown:
            lines.append(f"Unknown markers ({', '.join(unknown)}) — enable all sources.")
        return " | ".join(lines)


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
    ):
        import sys
        import traceback

        print(
            "SETUP RUN CLICKED:",
            {
                "n_clicks": n_clicks,
                "workflow_mode": workflow_mode,
                "analysis_method": analysis_method,
                "enable_qc": enable_qc,
                "enable_preprocess": enable_preprocess,
                "enable_analysis": enable_analysis,
                "enable_databases": enable_databases,
                "enable_phylo": enable_phylo,
            },
            flush=True,
        )

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
            )
        except PreventUpdate:
            raise
        except Exception as exc:
            print(f"CALLBACK ERROR: {exc}", file=sys.stderr, flush=True)
            traceback.print_exc(file=sys.stderr)
            return no_update, no_update, no_update, str(exc), no_update

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
            return no_update, no_update, no_update, warning

        all_samples = _samples_list(samples_data)
        selected_samples = [s for s in all_samples if s.get("selected", True)]
        if not selected_samples:
            warning = html.Div(
                "Select at least one sample before running the pipeline.",
                className="status-pill badge badge-amber",
                )
            return no_update, no_update, no_update, warning
        samples_data_filtered = {**(samples_data or {}), "samples": selected_samples}
        run_config = build_run_config(project_data, samples_data_filtered)
        run_config.enabled_stages = enabled_stages

        print(
            "RUN CONFIG BUILT:",
            {
                "run_id": run_config.run_id,
                "n_samples": len(run_config.samples),
                "branch": getattr(run_config.branch, "value", str(run_config.branch)),
                "analysis_method": run_config.analysis_method,
                "enabled_stages": run_config.enabled_stages,
            },
            flush=True,
        )

        # Preprocess params
        run_config.pear_params.update({
            "min_overlap": pear_min_overlap,
            "min_length": pear_min_length,
            "pvalue": pear_pvalue,
        })
        run_config.trimmomatic_params.update({
            "leading": trim_leading,
            "trailing": trim_trailing,
            "slidingwindow": trim_slidingwindow,
            "minlen": trim_minlen,
        })

        # VSEARCH params
        run_config.vsearch_params.update({
            "id": v_id,
            "maxee": v_maxee,
            "minlen": v_minlen,
            "sintax": v_sintax,
            "threads": int(threads) if threads else None,
        })

        # DADA2 params
        run_config.dada2_params.update({
            "truncF": d_truncF,
            "truncR": d_truncR,
            "maxeeF": d_maxeeF,
            "maxeeR": d_maxeeR,
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
        

        project_store = dict(project_data or {})
        projects = dict(project_store.get("projects", {}))
        active_id = project_store.get("active_id")
        if active_id and active_id in projects:
            projects[active_id] = dict(projects[active_id])
            projects[active_id]["mode"] = workflow_mode
            projects[active_id]["samples"] = selected_samples
            project_store = {"active_id": active_id, "projects": projects}
        
        print("DB SOURCES RAW:", db_sources, flush=True)
        print("RUN CONFIG VSEARCH PARAMS:", run_config.vsearch_params, flush=True)

        config_data = serialize_run_config(run_config)
        config_data["db_sources"] = list(db_sources or [])
        config_data["workflow_mode"] = workflow_mode
        print(
            "SETUP RETURN:",
            {
                "config_has_db_sources": "db_sources" in config_data,
                "config_keys": list(config_data.keys()),
                "project_active_id": project_store.get("active_id"),
                "project_keys": list(project_store.keys()),
                "pathname": "/progress",
            },
            flush=True,
        )
        return config_data, project_store, "/progress", "", samples_data