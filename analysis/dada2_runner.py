"""dada2 execution helpers for ASV workflows."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from core.config import RSCRIPT_PATH
from core.models import AnalysisMode, Marker, RunConfig, Sample, StepOutput
from core.utils import _count_fasta_seqs


def _gather_trimmed_files(samples: List[Sample]):
    fwd = []
    rev = []
    single = []
    for sample in samples:
        if sample.mode == AnalysisMode.PE and sample.paths.r1 and sample.paths.r2:
            fwd.append(sample.paths.r1)
            rev.append(sample.paths.r2)
        elif sample.paths.r1:
            single.append(sample.paths.r1)
    return fwd, rev, single


def run_dada2(samples: List[Sample], config: RunConfig, output_dir: Optional[Path] = None) -> StepOutput:
    ref_fasta = config.dada2_params.get("ref_fasta") or ""
    if not ref_fasta:
        raise RuntimeError("DADA2 requires a reference FASTA database (ref_fasta).")

    output_dir = output_dir or (config.output_dir / "analysis" / config.run_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    script_path = output_dir / "dada2_script.R"
    log_path = output_dir / "dada2_run.log"
    asv_table = output_dir / "asv_table.rds"
    abundance_csv = output_dir / "asv_abundance.csv"
    taxonomy_csv = output_dir / "taxonomy.csv"
    asv_fasta = output_dir / "asv_sequences.fasta"
    tracking_csv = output_dir / "read_tracking.csv"

    fwd_files, rev_files, single_files = _gather_trimmed_files(samples)
    path_filtered = output_dir / "filtered"
    path_filtered.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        if sample.paths.r1 and sample.paths.r1.exists():
            shutil.copy2(sample.paths.r1, path_filtered / sample.paths.r1.name)
        if sample.mode == AnalysisMode.PE and sample.paths.r2 and sample.paths.r2.exists():
            shutil.copy2(sample.paths.r2, path_filtered / sample.paths.r2.name)

    trunc_len_f = int(config.dada2_params.get("truncLen_F", 0))
    trunc_len_r = int(config.dada2_params.get("truncLen_R", 0))
    has_its = any(sample.marker == Marker.ITS for sample in samples)
    has_real_pe = any(
        sample.mode == AnalysisMode.PE
        and sample.paths.r2 is not None
        and sample.paths.r2.exists()
        for sample in samples
    )
    if not has_real_pe:
        trunc_len_f = 0
        trunc_len_r = 0
    warning = None
    r_script = f"""
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
ref_fasta <- args[3]
mode <- toupper(args[4])
truncLen_F <- as.integer(args[5])
truncLen_R <- as.integer(args[6])
maxEE_F <- as.numeric(args[7])
maxEE_R <- as.numeric(args[8])
truncQ <- as.integer(args[9])
multithread <- as.logical(args[10])

if (!dir.exists(input_dir)) {{
    stop("Input directory not found: ", input_dir)
}}
if (!dir.exists(output_dir)) {{
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}}

normalize_sample_name <- function(path) {{
    name <- basename(path)
    name <- sub("[.]fastq[.]gz$", "", name, ignore.case = TRUE)
    name <- sub("[.]fq[.]gz$", "", name, ignore.case = TRUE)
    name <- sub("[.]fastq$", "", name, ignore.case = TRUE)
    name <- sub("[.]fq$", "", name, ignore.case = TRUE)
    name <- sub("[.]gz$", "", name, ignore.case = TRUE)
    name <- sub("([._-]trimmed)$", "", name, ignore.case = TRUE)
    name <- sub("([._-]R?[12])$", "", name, ignore.case = TRUE)
    name <- sub("([._-](forward|reverse))$", "", name, ignore.case = TRUE)
    name <- sub("([._-](raw)?[._-]?[12])$", "", name, ignore.case = TRUE)
    name
}}

list_input_files <- function(dir_path) {{
    files <- list.files(dir_path, full.names = TRUE, recursive = FALSE)
    files[grepl("[.](fastq|fq)([.]gz)?$", files, ignore.case = TRUE)]
}}

build_taxonomy_string <- function(tax_row) {{
    values <- as.character(tax_row)
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 0) {{
        return("")
    }}
    paste(values, collapse = ";")
}}

write_outputs <- function(seqtab_nochim, taxa, track, output_dir) {{
    saveRDS(seqtab_nochim, file.path(output_dir, "seqtab_nochim.rds"))
    write.csv(as.data.frame(t(seqtab_nochim)), file.path(output_dir, "asv_abundance.csv"), quote = FALSE)
    if (exists("taxa") && !is.null(taxa) && nrow(taxa) > 0) {{
        taxonomy_df <- data.frame(
            ASV_ID = rownames(taxa),
            Taxonomy = apply(taxa, 1, build_taxonomy_string),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        write.csv(taxonomy_df, file.path(output_dir, "taxonomy.csv"), row.names = FALSE, quote = FALSE)
    }} else {{
        write.csv(data.frame(ASV_ID = character(), Taxonomy = character()), file.path(output_dir, "taxonomy.csv"), row.names = FALSE, quote = FALSE)
    }}
    if (!is.null(track) && nrow(track) > 0) {{
        track_df <- as.data.frame(track)
        track_df <- cbind(Sample = rownames(track_df), track_df)
        write.csv(track_df, file.path(output_dir, "read_tracking.csv"), row.names = FALSE, quote = FALSE)
    }}
    if (ncol(seqtab_nochim) > 0) {{
        asv_seqs <- DNAStringSet(colnames(seqtab_nochim))
        names(asv_seqs) <- paste0("ASV", seq_len(length(asv_seqs)))
        writeXStringSet(asv_seqs, file.path(output_dir, "asv_sequences.fasta"), format = "fasta")
    }} else {{
        writeLines(character(), file.path(output_dir, "asv_sequences.fasta"))
    }}
}}

input_files <- list_input_files(input_dir)
if (length(input_files) == 0) {{
    stop("No FASTQ inputs found in: ", input_dir)
}}

if (mode == "PE") {{
    fnFs <- input_files[grepl("_R1_trimmed[^/]*[.]fastq([.]gz)?$|_R1_trimmed[.]fastq([.]gz)?$", input_files, ignore.case = TRUE)]
    fnRs <- input_files[grepl("_R2_trimmed[^/]*[.]fastq([.]gz)?$|_R2_trimmed[.]fastq([.]gz)?$", input_files, ignore.case = TRUE)]
    if (length(fnFs) == 0 || length(fnRs) == 0 || length(fnFs) != length(fnRs)) {{
        stop("PE mode requires matching forward and reverse FASTQ files in the input directory.")
    }}
    fwd_names <- vapply(fnFs, normalize_sample_name, character(1))
    rev_names <- vapply(fnRs, normalize_sample_name, character(1))
    paired <- intersect(fwd_names, rev_names)
    if (length(paired) == 0) {{
        stop("PE mode could not match any forward/reverse pairs.")
    }}
    fnFs <- fnFs[match(paired, fwd_names)]
    fnRs <- fnRs[match(paired, rev_names)]
    sample_names <- paired

    filt_path <- file.path(output_dir, "filtered")
    if (!dir.exists(filt_path)) {{
        dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)
    }}
    filtFs <- file.path(filt_path, paste0(sample_names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample_names, "_R_filt.fastq.gz"))

    ft_multithread <- if (truncLen_F == 0 && truncLen_R == 0) FALSE else multithread
    out <- filterAndTrim(
        fnFs, filtFs,
        fnRs, filtRs,
        truncLen = c(truncLen_F, truncLen_R),
        maxEE = c(maxEE_F, maxEE_R),
        truncQ = truncQ,
        multithread = ft_multithread,
        compress = TRUE,
        verbose = TRUE
    )
    exists_filt <- file.exists(filtFs) & file.exists(filtRs)
    if (sum(exists_filt) == 0) {{
        stop("All samples were removed by filterAndTrim. Try reducing truncLen or increasing maxEE.")
    }}
    out_all <- out
    filtFs <- filtFs[exists_filt]
    filtRs <- filtRs[exists_filt]
    sample_names <- sample_names[exists_filt]

    errF <- learnErrors(filtFs, multithread = FALSE)
    errR <- learnErrors(filtRs, multithread = FALSE)

    names(filtFs) <- sample_names
    names(filtRs) <- sample_names
    derepFs <- derepFastq(filtFs)
    derepRs <- derepFastq(filtRs)
    names(derepFs) <- sample_names
    names(derepRs) <- sample_names

    dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
    dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
    seqtab <- makeSequenceTable(mergers)
    seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
    if (nchar(ref_fasta) > 0 && file.exists(ref_fasta)) {{
        taxa <- assignTaxonomy(seqtab_nochim, ref_fasta, multithread = FALSE)
    }} else {{
        taxa <- matrix(NA, nrow = ncol(seqtab_nochim), ncol = 1, dimnames = list(colnames(seqtab_nochim), "Kingdom"))
    }}

    out_filt <- out_all[exists_filt, , drop = FALSE]
    track <- cbind(
        out_filt,
        sapply(dadaFs, getN),
        sapply(mergers, getN),
        rowSums(seqtab_nochim)
    )
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample_names

    write_outputs(seqtab_nochim, taxa, track, output_dir)
}} else {{
    fnFs <- input_files
    sample_names <- vapply(fnFs, normalize_sample_name, character(1))

    filt_path <- file.path(output_dir, "filtered")
    if (!dir.exists(filt_path)) {{
        dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)
    }}
    filtFs <- file.path(filt_path, paste0(sample_names, "_filt.fastq.gz"))

    ft_multithread <- if (truncLen_F == 0) FALSE else multithread
    out <- filterAndTrim(
        fnFs, filtFs,
        truncLen = truncLen_F,
        maxEE = maxEE_F,
        truncQ = truncQ,
        multithread = ft_multithread,
        compress = TRUE,
        verbose = TRUE
    )
    exists_filt <- file.exists(filtFs)
    if (sum(exists_filt) == 0) {{
        stop("All samples were removed by filterAndTrim. Try reducing truncLen or increasing maxEE.")
    }}
    out_all <- out
    filtFs <- filtFs[exists_filt]
    sample_names <- sample_names[exists_filt]

    errF <- learnErrors(filtFs, multithread = FALSE)
    names(filtFs) <- sample_names
    derepFs <- derepFastq(filtFs)
    names(derepFs) <- sample_names

    dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
    seqtab <- makeSequenceTable(dadaFs)
    seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
    if (nchar(ref_fasta) > 0 && file.exists(ref_fasta)) {{
        taxa <- assignTaxonomy(seqtab_nochim, ref_fasta, multithread = FALSE)
    }} else {{
        taxa <- matrix(NA, nrow = ncol(seqtab_nochim), ncol = 1, dimnames = list(colnames(seqtab_nochim), "Kingdom"))
    }}

    out_filt <- out_all[exists_filt, , drop = FALSE]
    track <- cbind(
        out_filt,
        sapply(dadaFs, getN),
        rowSums(seqtab_nochim)
    )
    colnames(track) <- c("input", "filtered", "denoised", "nonchim")
    rownames(track) <- sample_names

    write_outputs(seqtab_nochim, taxa, track, output_dir)
}}
"""
    # If there is no real paired-end input, keep truncation at 0 so DADA2 uses the natural read length.
    has_real_pe = any(
        sample.mode == AnalysisMode.PE and sample.paths.r2 and sample.paths.r2.exists()
        for sample in samples
    )
    if not has_real_pe:
        trunc_len_f = 0
        trunc_len_r = 0
    if has_its:
        trunc_len_f = 0
        trunc_len_r = 0
        warning = "ITS detectado: truncLen configurado para 0,0."
    script_path.write_text(r_script, encoding="utf-8")

    command = [
        RSCRIPT_PATH,
        "--vanilla",
        str(script_path),
        str(path_filtered),
        str(output_dir),
        str(ref_fasta),
                str("PE" if any(sample.mode == AnalysisMode.PE for sample in samples) else "SE"),
                str(trunc_len_f),
                str(trunc_len_r),
                str(config.dada2_params.get("maxEE_F", 2)),
                str(config.dada2_params.get("maxEE_R", 2)),
                str(config.dada2_params.get("truncQ", 2)),
                "false",
    ]

    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
        log_path.write_text((completed.stdout or "") + "\n" + (completed.stderr or ""), encoding="utf-8")
    except subprocess.CalledProcessError as exc:
        log_path.write_text((exc.stdout or "") + "\n" + (exc.stderr or "") + "\n" + str(exc), encoding="utf-8")
        return StepOutput(
            step_name="dada2",
            success=False,
            files={
                "asv_table_rds": asv_table,
                "asv_abundance_csv": abundance_csv,
                "taxonomy_csv": taxonomy_csv,
                "asv_sequences_fasta": asv_fasta,
                "read_tracking_csv": tracking_csv,
                "log": log_path,
            },
            log=log_path,
            error_message=f"DADA2 falhou: {exc}",
        )

    if not tracking_csv.exists():
        return StepOutput(
            step_name="dada2",
            success=False,
            files={},
            log=log_path,
            error_message="DADA2 não gerou tracking de reads. Verifica o log R.",
        )

    if not asv_fasta.exists() or _count_fasta_seqs(asv_fasta) == 0:
        return StepOutput(
            step_name="dada2",
            success=False,
            files={},
            log=log_path,
            error_message="DADA2 não gerou ASVs válidas.",
        )

    return StepOutput(
        step_name="dada2",
        success=True,
        files={
            "asv_table_rds": asv_table,
            "asv_abundance_csv": abundance_csv,
            "taxonomy_csv": taxonomy_csv,
            "asv_sequences_fasta": asv_fasta,
            "read_tracking_csv": tracking_csv,
            "log": log_path,
        },
        log=log_path,
        error_message=warning,
    )

