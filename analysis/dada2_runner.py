"""dada2 execution helpers for ASV workflows.

DADA2 produces the ASV feature set: denoised, chimera-free amplicon sequence
variants plus an ASV x sample abundance matrix. It deliberately does NOT assign
taxonomy by default.

Taxonomy for ASVs is produced by the shared taxonomy stage (SINTAX / BLAST /
BOLD), exactly as it is for VSEARCH OTUs. Two reasons:

1. Parity. ``assignTaxonomy()`` is a naive-Bayes k-mer classifier: it reports a
   bootstrap, never a % identity or a query coverage. Routing ASVs through the
   same aligners as OTUs is the only way the two feature sets can carry the same
   evidence columns and be compared on equal terms.
2. Memory. ``assignTaxonomy()`` loads the whole reference into the R heap
   (several GB for a PR2/SILVA-sized reference), which is far more than most
   workstations have spare. The aligner path streams instead.

``ref_fasta`` therefore remains supported but strictly optional: set it only to
get DADA2's own native bootstrap alongside the shared taxonomy.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from core.config import RSCRIPT_PATH
from core.models import AnalysisMode, Marker, RunConfig, Sample, StepOutput
from core.utils import _count_fasta_seqs


# assignTaxonomy() holds the reference k-mer table in the R heap. The blow-up
# over the on-disk FASTA is large; this multiplier is the conservative rule of
# thumb used to refuse the call before the OS kills the whole process.
ASSIGN_TAXONOMY_MEMORY_MULTIPLIER = 6.0

# DADA2's own default for learnErrors(). The error model converges long before it
# needs this much data, so on a small machine it is mostly wasted resident memory.
DEFAULT_ERR_NBASES = 100_000_000
MIN_ERR_NBASES = 10_000_000

# Roughly how much of the R heap a base of error-model training data costs.
_ERR_BYTES_PER_BASE = 12
# Never hand the error model more than this share of free RAM: denoising still has to
# happen afterwards, in the same process.
_ERR_MEMORY_SHARE = 0.25


def learn_errors_nbases() -> int:
    """How many bases ``learnErrors()`` may read, given the RAM actually free.

    The full 1e8 default is a big resident block on a small box, and it is the easiest
    thing to give back: the model is fitting a handful of transition-rate parameters
    and converges on far less data. Denoising accuracy is unaffected — what changes is
    only how much of the run is sampled to fit the model.
    """
    available = _available_memory_bytes()
    if not available:
        return DEFAULT_ERR_NBASES
    budget = int(available * _ERR_MEMORY_SHARE / _ERR_BYTES_PER_BASE)
    return max(MIN_ERR_NBASES, min(DEFAULT_ERR_NBASES, budget))


def _available_memory_bytes() -> Optional[int]:
    """Best-effort available RAM in bytes (Linux /proc/meminfo)."""
    try:
        with open("/proc/meminfo", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        return None
    return None


def assign_taxonomy_memory_shortfall(ref_fasta: Path) -> Optional[str]:
    """Explain why ``ref_fasta`` cannot be loaded by assignTaxonomy(), or ``None``.

    Running it anyway means R is SIGKILLed partway through and the whole DADA2
    stage is lost, so a reference that cannot fit is skipped up front and the ASVs
    are still produced.
    """
    try:
        size = Path(ref_fasta).stat().st_size
    except OSError:
        return None
    available = _available_memory_bytes()
    if available is None:
        return None
    needed = size * ASSIGN_TAXONOMY_MEMORY_MULTIPLIER
    if needed <= available:
        return None
    return (
        f"assignTaxonomy() skipped: the reference ({size / 1e9:.1f} GB) needs about "
        f"{needed / 1e9:.1f} GB of RAM and only {available / 1e9:.1f} GB is available. "
        f"The ASVs are still produced; taxonomy comes from the Databases stage "
        f"(SINTAX/BLAST/BOLD), which does not load the whole reference into memory."
    )


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


R_SCRIPT = r"""
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
err_nbases <- as.numeric(args[11])

if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir)
}
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

getN <- function(x) sum(getUniques(x))

normalize_sample_name <- function(path) {
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
}

list_input_files <- function(dir_path) {
    files <- list.files(dir_path, full.names = TRUE, recursive = FALSE)
    files[grepl("[.](fastq|fq)([.]gz)?$", files, ignore.case = TRUE)]
}

# Every downstream table (taxonomy, abundance, tree) joins on the ASV id, so the
# id must be minted ONCE here and reused everywhere. Writing the raw DNA sequence
# as the abundance-matrix row key instead — which is what t(seqtab) gives you for
# free — makes the matrix impossible to join against a taxonomy keyed by "ASV1".
write_outputs <- function(seqtab_nochim, taxa, track, output_dir) {
    saveRDS(seqtab_nochim, file.path(output_dir, "asv_table.rds"))

    asv_seqs <- colnames(seqtab_nochim)
    n_asv <- length(asv_seqs)
    asv_ids <- if (n_asv > 0) paste0("ASV", seq_len(n_asv)) else character(0)
    sizes <- if (n_asv > 0) colSums(seqtab_nochim) else numeric(0)

    # ASV x sample abundance, keyed by ASV id. Mirrors the OTU branch's
    # otu_table.tsv so both feature sets feed the same metrics code.
    abundance <- as.data.frame(t(seqtab_nochim), check.names = FALSE)
    rownames(abundance) <- NULL
    abundance <- cbind(ASV_ID = asv_ids, abundance)
    write.csv(abundance, file.path(output_dir, "asv_abundance.csv"), row.names = FALSE, quote = FALSE)
    write.table(abundance, file.path(output_dir, "asv_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

    # Keeps "ASV7" traceable back to the exact sequence it stands for.
    write.csv(
        data.frame(asv_id = asv_ids, size = as.integer(sizes), sequence = asv_seqs, stringsAsFactors = FALSE),
        file.path(output_dir, "asv_id_map.csv"),
        row.names = FALSE,
        quote = FALSE
    )

    if (!is.null(taxa) && nrow(taxa) > 0) {
        # assignTaxonomy() keys its result by sequence; re-key it to the ASV id.
        idx <- match(asv_seqs, rownames(taxa))
        tax_strings <- vapply(seq_along(asv_seqs), function(i) {
            if (is.na(idx[i])) return("")
            values <- as.character(taxa[idx[i], ])
            values <- values[!is.na(values) & nzchar(values)]
            if (length(values) == 0) "" else paste(values, collapse = ";")
        }, character(1))
        write.csv(
            data.frame(ASV_ID = asv_ids, Taxonomy = tax_strings, stringsAsFactors = FALSE),
            file.path(output_dir, "taxonomy.csv"),
            row.names = FALSE,
            quote = FALSE
        )
    } else {
        write.csv(
            data.frame(ASV_ID = character(), Taxonomy = character()),
            file.path(output_dir, "taxonomy.csv"),
            row.names = FALSE,
            quote = FALSE
        )
    }

    if (!is.null(track) && nrow(track) > 0) {
        track_df <- as.data.frame(track)
        track_df <- cbind(Sample = rownames(track_df), track_df)
        write.csv(track_df, file.path(output_dir, "read_tracking.csv"), row.names = FALSE, quote = FALSE)
    }

    if (n_asv > 0) {
        seqs <- DNAStringSet(asv_seqs)
        # ";size=" mirrors the VSEARCH OTU header so the shared taxonomy and phylo
        # code can read an abundance off either feature set's FASTA identically.
        names(seqs) <- paste0(asv_ids, ";size=", as.integer(sizes))
        writeXStringSet(seqs, file.path(output_dir, "asv_sequences.fasta"), format = "fasta")
    } else {
        writeLines(character(), file.path(output_dir, "asv_sequences.fasta"))
    }
}

input_files <- list_input_files(input_dir)
if (length(input_files) == 0) {
    stop("No FASTQ inputs found in: ", input_dir)
}

assign_taxa <- function(seqtab_nochim) {
    if (nchar(ref_fasta) > 0 && file.exists(ref_fasta)) {
        return(assignTaxonomy(seqtab_nochim, ref_fasta, multithread = FALSE))
    }
    # No reference: ASVs are still fully valid. Taxonomy is assigned downstream by
    # the shared SINTAX/BLAST/BOLD stage, which is also what gives them the
    # identity/coverage columns the OTUs get.
    NULL
}

if (mode == "PE") {
    fnFs <- input_files[grepl("_R1_trimmed[^/]*[.]fastq([.]gz)?$|_R1_trimmed[.]fastq([.]gz)?$", input_files, ignore.case = TRUE)]
    fnRs <- input_files[grepl("_R2_trimmed[^/]*[.]fastq([.]gz)?$|_R2_trimmed[.]fastq([.]gz)?$", input_files, ignore.case = TRUE)]
    if (length(fnFs) == 0 || length(fnRs) == 0 || length(fnFs) != length(fnRs)) {
        stop("PE mode requires matching forward and reverse FASTQ files in the input directory.")
    }
    fwd_names <- vapply(fnFs, normalize_sample_name, character(1))
    rev_names <- vapply(fnRs, normalize_sample_name, character(1))
    paired <- intersect(fwd_names, rev_names)
    if (length(paired) == 0) {
        stop("PE mode could not match any forward/reverse pairs.")
    }
    fnFs <- fnFs[match(paired, fwd_names)]
    fnRs <- fnRs[match(paired, rev_names)]
    sample_names <- paired

    filt_path <- file.path(output_dir, "filtered")
    if (!dir.exists(filt_path)) {
        dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)
    }
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
    if (sum(exists_filt) == 0) {
        stop("All samples were removed by filterAndTrim. Try reducing truncLen or increasing maxEE.")
    }
    out_all <- out
    filtFs <- filtFs[exists_filt]
    filtRs <- filtRs[exists_filt]
    sample_names <- sample_names[exists_filt]

    names(filtFs) <- sample_names
    names(filtRs) <- sample_names

    # nbases caps how much data the error model reads. The default (1e8) is far more
    # than it needs to converge and is a large slice of the peak memory on a small box.
    errF <- learnErrors(filtFs, nbases = err_nbases, multithread = FALSE, randomize = TRUE)
    errR <- learnErrors(filtRs, nbases = err_nbases, multithread = FALSE, randomize = TRUE)

    # Sample-by-sample, DADA2's own "big data" workflow. Dereplicating every sample at
    # once (derepFastq on the whole vector) holds every unique sequence of every sample
    # in memory simultaneously, which is what the OOM killer was reaping. Here only one
    # sample's dereps are live at a time and the merged result — much smaller — is all
    # that is kept.
    mergers <- vector("list", length(sample_names))
    names(mergers) <- sample_names
    denoised_f <- integer(length(sample_names))
    denoised_r <- integer(length(sample_names))
    merged_n <- integer(length(sample_names))

    for (i in seq_along(sample_names)) {
        sam <- sample_names[i]
        cat("Denoising sample", sam, "\n")
        derep_f <- derepFastq(filtFs[[i]])
        dada_f <- dada(derep_f, err = errF, multithread = FALSE)
        derep_r <- derepFastq(filtRs[[i]])
        dada_r <- dada(derep_r, err = errR, multithread = FALSE)
        merger <- mergePairs(dada_f, derep_f, dada_r, derep_r, verbose = TRUE)

        mergers[[sam]] <- merger
        denoised_f[i] <- getN(dada_f)
        denoised_r[i] <- getN(dada_r)
        merged_n[i] <- getN(merger)

        rm(derep_f, derep_r, dada_f, dada_r, merger)
        gc()
    }

    seqtab <- makeSequenceTable(mergers)
    rm(mergers)
    gc()
    seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
    taxa <- assign_taxa(seqtab_nochim)

    out_filt <- out_all[exists_filt, , drop = FALSE]
    track <- cbind(out_filt, denoised_f, denoised_r, merged_n, rowSums(seqtab_nochim))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample_names

    write_outputs(seqtab_nochim, taxa, track, output_dir)
} else {
    fnFs <- input_files
    sample_names <- vapply(fnFs, normalize_sample_name, character(1))

    filt_path <- file.path(output_dir, "filtered")
    if (!dir.exists(filt_path)) {
        dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)
    }
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
    if (sum(exists_filt) == 0) {
        stop("All samples were removed by filterAndTrim. Try reducing truncLen or increasing maxEE.")
    }
    out_all <- out
    filtFs <- filtFs[exists_filt]
    sample_names <- sample_names[exists_filt]

    names(filtFs) <- sample_names
    errF <- learnErrors(filtFs, nbases = err_nbases, multithread = FALSE, randomize = TRUE)

    # Same one-sample-at-a-time loop as the PE branch, for the same reason.
    dada_list <- vector("list", length(sample_names))
    names(dada_list) <- sample_names
    denoised_n <- integer(length(sample_names))

    for (i in seq_along(sample_names)) {
        sam <- sample_names[i]
        cat("Denoising sample", sam, "\n")
        derep_f <- derepFastq(filtFs[[i]])
        dada_f <- dada(derep_f, err = errF, multithread = FALSE)
        dada_list[[sam]] <- dada_f
        denoised_n[i] <- getN(dada_f)
        rm(derep_f)
        gc()
    }

    seqtab <- makeSequenceTable(dada_list)
    rm(dada_list)
    gc()
    seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
    taxa <- assign_taxa(seqtab_nochim)

    out_filt <- out_all[exists_filt, , drop = FALSE]
    track <- cbind(out_filt, denoised_n, rowSums(seqtab_nochim))
    colnames(track) <- c("input", "filtered", "denoised", "nonchim")
    rownames(track) <- sample_names

    write_outputs(seqtab_nochim, taxa, track, output_dir)
}
"""


def run_dada2(samples: List[Sample], config: RunConfig, output_dir: Optional[Path] = None) -> StepOutput:
    """Denoise reads into ASVs. Taxonomy is left to the shared taxonomy stage."""
    from core.preflight import dada2_memory_shortfall

    output_dir = output_dir or (config.output_dir / "analysis" / config.run_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Refuse before spending a minute loading R only to be SIGKILLed by the OOM killer
    # partway through. The message names the real culprit; "signal 9" does not.
    shortfall = dada2_memory_shortfall()
    if shortfall:
        log_path = output_dir / "dada2_run.log"
        log_path.write_text(shortfall, encoding="utf-8")
        return StepOutput(step_name="dada2", success=False, files={"log": log_path}, log=log_path, error_message=shortfall)

    script_path = output_dir / "dada2_script.R"
    log_path = output_dir / "dada2_run.log"
    asv_table_rds = output_dir / "asv_table.rds"
    abundance_csv = output_dir / "asv_abundance.csv"
    asv_table_tsv = output_dir / "asv_table.tsv"
    asv_id_map = output_dir / "asv_id_map.csv"
    taxonomy_csv = output_dir / "taxonomy.csv"
    asv_fasta = output_dir / "asv_sequences.fasta"
    tracking_csv = output_dir / "read_tracking.csv"

    warnings: List[str] = []

    # A reference is optional. When one is configured, only use it if it actually
    # fits in RAM — otherwise R gets SIGKILLed and the ASVs are lost along with it.
    ref_fasta = str(config.dada2_params.get("ref_fasta") or "")
    if ref_fasta and not Path(ref_fasta).exists():
        warnings.append(f"DADA2 reference not found ({ref_fasta}); assignTaxonomy() skipped.")
        ref_fasta = ""
    if ref_fasta:
        shortfall = assign_taxonomy_memory_shortfall(Path(ref_fasta))
        if shortfall:
            warnings.append(shortfall)
            ref_fasta = ""

    path_filtered = output_dir / "filtered_input"
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
        sample.mode == AnalysisMode.PE and sample.paths.r2 is not None and sample.paths.r2.exists()
        for sample in samples
    )
    # Truncating to a fixed length only makes sense when there is a real read pair to
    # overlap; ITS amplicons vary in length by design and must never be truncated.
    if not has_real_pe or has_its:
        trunc_len_f = 0
        trunc_len_r = 0
    if has_its:
        warnings.append("ITS detected: truncLen set to 0,0.")

    err_nbases = int(config.dada2_params.get("err_nbases") or learn_errors_nbases())
    if err_nbases < DEFAULT_ERR_NBASES:
        warnings.append(
            f"Error model limited to {err_nbases / 1e6:.0f} Mbases (normally "
            f"{DEFAULT_ERR_NBASES / 1e6:.0f}) to fit in the available memory. "
            f"That is enough sampling for the model to converge; it does not affect the ASVs."
        )

    script_path.write_text(R_SCRIPT, encoding="utf-8")

    files = {
        "asv_table_rds": asv_table_rds,
        "asv_abundance_csv": abundance_csv,
        "asv_table_tsv": asv_table_tsv,
        "asv_id_map_csv": asv_id_map,
        "taxonomy_csv": taxonomy_csv,
        "asv_sequences_fasta": asv_fasta,
        "read_tracking_csv": tracking_csv,
        "log": log_path,
    }

    command = [
        RSCRIPT_PATH,
        "--vanilla",
        str(script_path),
        str(path_filtered),
        str(output_dir),
        ref_fasta,
        "PE" if has_real_pe else "SE",
        str(trunc_len_f),
        str(trunc_len_r),
        str(config.dada2_params.get("maxEE_F", 2)),
        str(config.dada2_params.get("maxEE_R", 2)),
        str(config.dada2_params.get("truncQ", 2)),
        "false",
        str(int(err_nbases)),
    ]

    preamble = ("\n".join(warnings) + "\n\n") if warnings else ""
    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
        log_path.write_text(preamble + (completed.stdout or "") + "\n" + (completed.stderr or ""), encoding="utf-8")
    except subprocess.CalledProcessError as exc:
        log_path.write_text(preamble + (exc.stdout or "") + "\n" + (exc.stderr or "") + "\n" + str(exc), encoding="utf-8")
        detail = str(exc.stderr or exc.stdout or exc)[:400]
        # A negative return code is the OS killing R, which on a denoising run is
        # almost always the out-of-memory killer.
        if exc.returncode is not None and exc.returncode < 0:
            detail = (
                f"R was killed by the system (signal {-exc.returncode}), almost certainly for "
                f"lack of memory (OOM). Reduce the number of samples per run, or raise the RAM."
            )
        return StepOutput(step_name="dada2", success=False, files={"log": log_path}, log=log_path, error_message=f"DADA2 failed: {detail}")

    if not tracking_csv.exists():
        return StepOutput(step_name="dada2", success=False, files={"log": log_path}, log=log_path, error_message="DADA2 produced no read tracking. Check the R log.")

    if not asv_fasta.exists() or _count_fasta_seqs(asv_fasta) == 0:
        return StepOutput(step_name="dada2", success=False, files={"log": log_path}, log=log_path, error_message="DADA2 produced no valid ASVs.")

    return StepOutput(
        step_name="dada2",
        success=True,
        files={key: path for key, path in files.items() if path.exists()},
        log=log_path,
        error_message="\n".join(warnings) if warnings else None,
    )
