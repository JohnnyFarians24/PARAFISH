"""vsearch execution helpers for OTU workflows."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Dict, Optional

from core.config import VSEARCH_PATH
from core.models import RunConfig, StepOutput
from core.utils import _count_fasta_seqs


def _parse_vsearch_uc_to_matrix(uc_file, sample_counts):
    if not uc_file or not os.path.exists(uc_file) or not sample_counts:
        return None

    otu_counts = {}
    otu_order = []
    with open(uc_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            otu_key = parts[1]
            query_label = parts[8]
            if "|" not in query_label:
                continue
            sample_id = query_label.split("|", 1)[0]
            if sample_id not in sample_counts:
                continue
            if otu_key not in otu_counts:
                otu_counts[otu_key] = {}
                otu_order.append(otu_key)
            otu_counts[otu_key][sample_id] = otu_counts[otu_key].get(sample_id, 0) + 1

    if not otu_counts:
        return None

    rows = []
    for otu_key in otu_order:
        row = {sid: otu_counts[otu_key].get(sid, 0) for sid in sample_counts}
        row["OTU"] = f"OTU_{len(rows) + 1}"
        rows.append(row)

    import pandas as pd

    matrix = pd.DataFrame(rows)
    if matrix.empty:
        return None
    return matrix.set_index("OTU")


def _run(command):
    return subprocess.run(command, capture_output=True, text=True, check=True)


def run_vsearch(assembled_fastq: Path, config: RunConfig, output_dir: Path, db_sintax: Optional[Path] = None) -> StepOutput:
    """Run the VSEARCH OTU workflow from a prepared assembled FASTQ file."""
    identity = float(config.vsearch_params.get("id", 0.97))
    if identity < 0.8 or identity > 1.0:
        return StepOutput(step_name="vsearch", success=False, files={}, error_message="Parâmetro --id fora do intervalo permitido [0.8, 1.0].")

    assembled = Path(assembled_fastq)
    if not assembled.exists():
        return StepOutput(step_name="vsearch", success=False, files={}, error_message=f"FASTQ assembled não encontrado: {assembled}.")

    output_dir.mkdir(parents=True, exist_ok=True)
    filtered = output_dir / "filtered.fasta"
    derep = output_dir / "derep.fasta"
    otus = output_dir / "otus.fasta"
    nonchimeras = output_dir / "nonchimeras.fasta"
    chimeras = output_dir / "chimeras.fasta"
    chimera_uchime = output_dir / "chimera.uchime"
    taxonomy = output_dir / "taxonomy.sintax"
    otu_table = output_dir / "otu_table.biom"
    clusters_uc = output_dir / "clusters.uc"
    map_uc = output_dir / "map.uc"
    combined_log = output_dir / "vsearch.log"

    params = config.vsearch_params
    pieces = []
    sintax_warning = None

    commands = [
        [
            VSEARCH_PATH,
            "--fastq_filter",
            str(assembled),
            "--fastq_maxee",
            str(params.get("fastq_maxee", 1.0)),
            "--fastq_minlen",
            str(params.get("fastq_minlen", 200)),
            "--fastaout",
            str(filtered),
            "--log",
            str(output_dir / "step1_filter.log"),
        ],
        [
            VSEARCH_PATH,
            "--derep_fulllength",
            str(filtered),
            "--sizeout",
            "--output",
            str(derep),
            "--log",
            str(output_dir / "step2_derep.log"),
        ],
        [
            VSEARCH_PATH,
            "--cluster_size",
            str(derep),
            "--id",
            str(identity),
            "--centroids",
            str(otus),
            "--uc",
            str(clusters_uc),
            "--log",
            str(output_dir / "step3_cluster.log"),
        ],
        [
            VSEARCH_PATH,
            "--uchime_denovo",
            str(otus),
            "--nonchimeras",
            str(nonchimeras),
            "--chimeras",
            str(chimeras),
            "--uchimeout",
            str(chimera_uchime),
            "--log",
            str(output_dir / "step4_chimera.log"),
        ],
    ]

    try:
        for command in commands:
            completed = _run(command)
            pieces.append((completed.stdout or "") + "\n" + (completed.stderr or ""))
    except subprocess.CalledProcessError as exc:
        combined_log.write_text(str(exc), encoding="utf-8")
        return StepOutput(step_name="vsearch", success=False, files={}, log=combined_log, error_message=f"VSEARCH falhou: {exc}")

    if _count_fasta_seqs(otus) == 0:
        combined_log.write_text("\n".join(pieces), encoding="utf-8")
        return StepOutput(
            step_name="vsearch",
            success=False,
            files={
                "filtered_fasta": filtered,
                "derep_fasta": derep,
                "otus_fasta": otus,
                "nonchimeras_fasta": nonchimeras,
                "chimeras_fasta": chimeras,
                "chimera_uchime": chimera_uchime,
                "taxonomy_sintax": taxonomy,
                "otu_table_biom": otu_table,
                "clusters_uc": clusters_uc,
                "log": combined_log,
            },
            log=combined_log,
            error_message="Clustering produziu 0 OTUs. Verifica o parâmetro --id e a qualidade das reads.",
        )

    if db_sintax is not None and db_sintax.exists():
        try:
            completed = _run([
                VSEARCH_PATH,
                "--sintax",
                str(nonchimeras),
                "--db",
                str(db_sintax),
                "--sintax_cutoff",
                str(params.get("sintax_cutoff", 0.8)),
                "--tabbedout",
                str(taxonomy),
                "--threads",
                str(params.get("threads", 8)),
                "--log",
                str(output_dir / "step5_sintax.log"),
            ])
            pieces.append((completed.stdout or "") + "\n" + (completed.stderr or ""))
        except subprocess.CalledProcessError as exc:
            sintax_warning = f"SINTAX falhou (ignorado): {exc}"
            pieces.append(sintax_warning)
    else:
        sintax_warning = "db_sintax não encontrado; passo SINTAX ignorado."

    try:
        completed = _run([
            VSEARCH_PATH,
            "--usearch_global",
            str(filtered),
            "--db",
            str(nonchimeras),
            "--id",
            str(identity),
            "--uc",
            str(map_uc),
            "--otutabout",
            str(otu_table),
            "--threads",
            str(params.get("threads", 8)),
            "--log",
            str(output_dir / "step6_otutable.log"),
        ])
        pieces.append((completed.stdout or "") + "\n" + (completed.stderr or ""))
    except subprocess.CalledProcessError as exc:
        combined_log.write_text("\n".join(pieces) + "\n" + str(exc), encoding="utf-8")
        return StepOutput(step_name="vsearch", success=False, files={}, log=combined_log, error_message=f"OTU table falhou: {exc}")

    combined_log.write_text("\n".join(pieces), encoding="utf-8")
    return StepOutput(
        step_name="vsearch",
        success=True,
        files={
            "filtered_fasta": filtered,
            "derep_fasta": derep,
            "otus_fasta": otus,
            "nonchimeras_fasta": nonchimeras,
            "chimeras_fasta": chimeras,
            "chimera_uchime": chimera_uchime,
            "taxonomy_sintax": taxonomy,
            "otu_table_biom": otu_table,
            "clusters_uc": clusters_uc,
            "log": combined_log,
        },
        log=combined_log,
        error_message=sintax_warning,
    )


def run_vsearch_otus(input_fasta, output_dir, vsearch_path=None, identity=0.97):
    """Compatibility wrapper returning the legacy dict-based OTU result payload."""
    vpath = vsearch_path or VSEARCH_PATH
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    input_path = Path(input_fasta)
    derep = out_dir / "derep.fasta"
    repseqs = out_dir / "otus.fasta"
    uc_file = out_dir / "clusters.uc"
    logs = {}

    try:
        derep_run = subprocess.run(
            [vpath, "--derep_fulllength", str(input_path), "--sizeout", "--output", str(derep)],
            capture_output=True,
            text=True,
            check=True,
        )
        logs["derep"] = (derep_run.stdout or "") + (derep_run.stderr or "")
        cluster_run = subprocess.run(
            [vpath, "--cluster_size", str(derep), "--id", str(identity), "--centroids", str(repseqs), "--uc", str(uc_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        logs["cluster"] = (cluster_run.stdout or "") + (cluster_run.stderr or "")
    except Exception as exc:
        return {"error": str(exc), "repseqs": None, "uc_file": None, "logs": logs}

    if not repseqs.exists() or repseqs.stat().st_size == 0:
        return {"error": "Clustering produced no OTUs.", "repseqs": None, "uc_file": None, "logs": logs}

    return {
        "error": None,
        "repseqs": str(repseqs),
        "uc_file": str(uc_file) if uc_file.exists() else None,
        "logs": logs,
    }

