#!/usr/bin/env python
"""Cross-file consistency checks for a PARAFISH run. Prints PASS/FAIL per check.

This is the non-regression net for the phased fix plan. Every check compares files
written by *different* functions that must agree; none of them re-implements the
pipeline's logic, they only assert that two independent writers tell the same story.
That is deliberate — the whole class of bug being hunted here is "what one function
writes is not what the next one reads".

The five checks:

  (a) read totals agree across otu_table / otu_matrix / chimera_stats /
      taxonomic_resolution / every taxa_table_*
  (b) otu_id_map.size == that OTU's otu_table row sum
  (c) taxa_table_class agrees with the Class column of tax_table, collapsed over
      otu_table abundances
  (d) no invented middle ranks (Phylum=Opisthokonta, Class=Craniata)
  (e) source_status matches the assignments each source actually made

Every check runs PER FEATURE SET actually present in the run (OTU and/or ASV). This
matters more than it looks: a verifier that only reads ``databases/otu/`` does not FAIL
on an ASV run — it SKIPs, silently, and the run exits 0. A blind instrument goes quiet,
it does not scream, so an ASV run's "verify passed" would be validation against nothing.
The checks resolve the set's own paths (``databases/<set>/``, the set's abundance table)
and label their name with the set, so ``(d/asv)`` is a real, biting check.

Checks whose inputs are absent report SKIP, not PASS: a run without the databases
stage must not look like a clean bill of health. Exit code is non-zero only on FAIL.

Usage:
    python scripts/verify.py <run_dir> [--quiet]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Optional

# Ranks that must never be invented into a Linnean slot. Opisthokonta is a PR2
# supergroup and Craniata a subphylum-ish clade; both appeared because a source's
# own backbone was force-fitted into Phylum/Class.
BANNED = {"phylum": {"opisthokonta", "amorphea", "obazoa"}, "class": {"craniata"}}

UNCLASSIFIED = {"", "na", "nan", "none", "unclassified", "unassigned", "unknown"}

# The raw output each source leaves behind, relative to databases/<set>/. Presence of
# any of these means the source RAN, whatever source_status does or does not say about
# it — which is the point: a source that ran and was then dropped on the floor is
# invisible to a check that only walks source_status's own rows.
SOURCE_RAW_OUTPUTS = {
    "SILVA": ["silva/silva_sintax_result.csv"],
    "PR2": ["pr2/pr2_sintax_result.csv"],
    "NCBI": ["ncbi/blast_result.tsv"],
    "BOLD": [
        "bold/boldigger3_data/sanitized_input_identification_result.parquet.snappy",
        "bold/boldigger3_data/sanitized_input_identification_result.xlsx",
    ],
    "CUSTOM": ["custom_blast/custom_blast_raw.tsv"],
}


class Report:
    def __init__(self, quiet: bool = False):
        self.rows: list[tuple[str, str, str]] = []
        self.quiet = quiet

    def add(self, status: str, name: str, detail: str = "") -> None:
        self.rows.append((status, name, detail))
        if not self.quiet:
            print(f"  {status:4}  {name}" + (f"\n        {detail}" if detail else ""))

    @property
    def failed(self) -> int:
        return sum(1 for status, _, _ in self.rows if status == "FAIL")

    @property
    def skipped(self) -> int:
        return sum(1 for status, _, _ in self.rows if status == "SKIP")


def _read_csv(path: Path, delimiter: str = ",") -> list[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    return list(csv.DictReader(text.splitlines(), delimiter=delimiter))


def _norm(value) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in UNCLASSIFIED else text


def _matrix(path: Path, delimiter: str) -> tuple[list[str], list[str], dict[str, list[float]]]:
    """(feature_ids, sample_names, {feature: counts}) from a feature x sample table."""
    if not path.exists() or path.stat().st_size == 0:
        return [], [], {}
    rows = list(csv.reader(path.read_text(encoding="utf-8", errors="replace").splitlines(), delimiter=delimiter))
    rows = [r for r in rows if r and any(c.strip() for c in r)]
    if not rows:
        return [], [], {}
    samples = [c.strip() for c in rows[0][1:]]
    features, values = [], {}
    for row in rows[1:]:
        feature = row[0].strip()
        features.append(feature)
        values[feature] = [float(c or 0) for c in row[1:]]
    return features, samples, values


def _total(values: dict[str, list[float]]) -> float:
    return sum(sum(counts) for counts in values.values())


# ── feature-set layout ───────────────────────────────────────────────────
#
# Each stage writes one directory per feature set, and the two sets name their files
# differently on the analysis side (VSEARCH pools OTUs; DADA2 denoises ASVs). These
# resolvers are the whole ASV-awareness: a check reads THIS set's files, so it verifies
# an ASV run instead of skipping it. The databases side is a clean ``databases/<set>/``
# substitution; only the analysis-side names differ.

FEATURE_SETS = ("otu", "asv")


def _abundance_table(run: Path, fset: str) -> tuple[Path, str]:
    """(path, delimiter) of the feature × sample table for a set."""
    if fset == "asv":
        return run / "analysis/dada2/asv_table.tsv", "\t"
    return run / "analysis/pooled/otu_table.tsv", "\t"


def _matrix_table(run: Path, fset: str) -> tuple[Path, str]:
    """(path, delimiter) of the set's other abundance table (the .csv sibling)."""
    if fset == "asv":
        return run / "analysis/dada2/asv_abundance.csv", ","
    return run / "analysis/otu_matrix.csv", ","


def _id_map(run: Path, fset: str) -> tuple[Path, str]:
    """(path, id_column) of the per-feature size map. Both sets carry one; the id column differs."""
    if fset == "asv":
        return run / "analysis/dada2/asv_id_map.csv", "asv_id"
    return run / "analysis/pooled/otu_id_map.csv", "otu_id"


def _db_dir(run: Path, fset: str) -> Path:
    return run / "databases" / fset


def _analysis_metric(run: Path, fset: str, name: str) -> Path:
    return run / "analysis" / "metrics" / fset / name


def _present_sets(run: Path) -> list[str]:
    """The feature sets this run actually produced, by the presence of the abundance table.

    Empty only when neither exists; callers default to a single OTU pass so an empty run
    still SKIPs loudly rather than falling silent.
    """
    return [fset for fset in FEATURE_SETS if _abundance_table(run, fset)[0].exists()]


def _requested_sources(run: Path) -> Optional[set[str]]:
    """The db_sources the run was ASKED for, from progress.json's run_config (upper-cased).

    This is the record of what was REQUESTED, written at init — distinct from what arrived,
    which is exactly the distinction check (f) needs to tell a config choice from a death.
    Returns ``None`` when provenance does not record the request (older runs), so the caller
    can fall back rather than mistake "unrecorded" for "nothing requested".
    """
    try:
        progress = json.loads((run / "progress.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    db_sources = (progress.get("run_config") or {}).get("db_sources")
    if db_sources is None:
        return None
    return {str(source).strip().upper() for source in db_sources if str(source).strip()}


# ── (a) read totals ──────────────────────────────────────────────────────


def check_read_totals(run: Path, report: Report, fset: str) -> None:
    name = f"(a/{fset}) read totals agree"
    table_path, table_delim = _abundance_table(run, fset)
    if not table_path.exists():
        report.add("SKIP", name, f"no {table_path.relative_to(run)}")
        return
    _, _, values = _matrix(table_path, table_delim)
    canonical = _total(values)
    totals = {table_path.name: canonical}

    matrix_path, matrix_delim = _matrix_table(run, fset)
    if matrix_path.exists():
        _, _, mvals = _matrix(matrix_path, matrix_delim)
        totals[matrix_path.name] = _total(mvals)

    # DADA2 removes chimeras internally, so there is no ASV chimera_stats — the leg is OTU-only.
    if fset == "otu":
        chimera = _read_csv(_analysis_metric(run, fset, "chimera_stats.csv"))
        if chimera:
            totals["chimera_stats.total_reads"] = float(chimera[0].get("total_reads") or 0)

    resolution = _read_csv(_db_dir(run, fset) / "metrics/taxonomic_resolution.csv")
    if resolution:
        totals["taxonomic_resolution.n_reads"] = sum(float(r.get("n_reads") or 0) for r in resolution)

    taxa_dir = _db_dir(run, fset) / "taxa_tables"
    if taxa_dir.exists():
        for path in sorted(taxa_dir.glob("taxa_table_*.csv")):
            if path.name.endswith("_relative.csv"):
                continue
            _, _, tvals = _matrix(path, ",")
            totals[path.name] = _total(tvals)

    bad = {name_: value for name_, value in totals.items() if abs(value - canonical) > 0.5}
    detail = " · ".join(f"{name_}={int(value)}" for name_, value in totals.items())
    if bad:
        report.add("FAIL", name, f"canonical({table_path.name})={int(canonical)} · {detail}")
    else:
        report.add("PASS", name, f"all {len(totals)} sources = {int(canonical)}")


# ── (b) otu_id_map.size ──────────────────────────────────────────────────


def check_otu_id_map(run: Path, report: Report, fset: str) -> None:
    name = f"(b/{fset}) id_map.size == table row sum"
    id_map, id_col = _id_map(run, fset)
    table_path, table_delim = _abundance_table(run, fset)
    if not id_map.exists() or not table_path.exists():
        report.add("SKIP", name, "missing input")
        return
    _, _, values = _matrix(table_path, table_delim)
    rows = _read_csv(id_map)
    if not rows:
        report.add("SKIP", name, f"empty {id_map.name}")
        return
    mismatches = []
    for row in rows:
        feature = (row.get(id_col) or "").strip()
        size = float(row.get("size") or 0)
        expected = sum(values.get(feature, []))
        if abs(size - expected) > 0.5:
            mismatches.append(f"{feature}: size={int(size)} row_sum={int(expected)}")
    if mismatches:
        report.add("FAIL", name, " · ".join(mismatches[:6]))
    else:
        report.add("PASS", name, f"{len(rows)}/{len(rows)} features match")


# ── (c) taxa_table_class vs tax_table ────────────────────────────────────


def check_class_table(run: Path, report: Report, fset: str) -> None:
    name = f"(c/{fset}) taxa_table_class agrees with tax_table"
    taxa_class = _db_dir(run, fset) / "taxa_tables/taxa_table_class.csv"
    tax_table = _db_dir(run, fset) / "phyloseq/tax_table.tsv"
    table_path, table_delim = _abundance_table(run, fset)
    if not (taxa_class.exists() and tax_table.exists() and table_path.exists()):
        report.add("SKIP", name, "missing input")
        return

    _, samples, values = _matrix(table_path, table_delim)
    classes = {
        (row.get("feature_id") or "").strip(): _norm(row.get("Class"))
        for row in _read_csv(tax_table, delimiter="\t")
    }

    # Collapse otu_table by tax_table's Class exactly as the exporter should have.
    expected: dict[str, list[float]] = {}
    for feature, counts in values.items():
        label = classes.get(feature, "") or "Unclassified"
        acc = expected.setdefault(label, [0.0] * len(counts))
        for j, count in enumerate(counts):
            acc[j] += count

    _, actual_samples, actual = _matrix(taxa_class, ",")

    # Compare per sample NAME, not per column position: the two tables are written by
    # different functions and nothing guarantees they order samples the same way. A
    # positional compare would pass on transposed columns, which is a real bug class.
    if set(actual_samples) != set(samples):
        report.add(
            "FAIL",
            name,
            f"sample sets differ: abundance={samples} taxa_table_class={actual_samples}",
        )
        return
    index = [actual_samples.index(sample) for sample in samples]

    problems = []
    for label in set(expected) | set(actual):
        exp = expected.get(label)
        act = actual.get(label)
        if exp is None:
            problems.append(f"{label}: in taxa_table_class only ({int(sum(act))})")
        elif act is None:
            problems.append(f"{label}: in tax_table only ({int(sum(exp))})")
        else:
            aligned = [act[i] for i in index]
            if any(abs(a - b) > 0.5 for a, b in zip(exp, aligned)):
                problems.append(f"{label}: tax_table={[int(v) for v in exp]} taxa_table={[int(v) for v in aligned]}")
    if problems:
        report.add("FAIL", name, " · ".join(problems[:6]))
    else:
        report.add("PASS", name, f"{len(expected)} class(es) match")


# ── (d) no invented middle ranks ─────────────────────────────────────────


def check_banned_ranks(run: Path, report: Report, fset: str) -> None:
    name = f"(d/{fset}) no invented middle ranks"
    tax_table = _db_dir(run, fset) / "phyloseq/tax_table.tsv"
    consensus = _db_dir(run, fset) / "taxonomy_consensus.csv"
    sources = []
    if tax_table.exists():
        sources.append((tax_table, _read_csv(tax_table, delimiter="\t"), {"phylum": "Phylum", "class": "Class"}))
    if consensus.exists():
        sources.append((consensus, _read_csv(consensus), {"phylum": "Consensus_Phylum", "class": "Consensus_Class"}))
    if not sources:
        report.add("SKIP", name, "no taxonomy output")
        return
    hits = []
    for path, rows, columns in sources:
        for row in rows:
            for rank, column in columns.items():
                value = _norm(row.get(column)).lower()
                if value and value in BANNED[rank]:
                    feature = row.get("feature_id") or row.get("ASV_ID") or "?"
                    hits.append(f"{path.name}:{feature} {rank}={value}")
    if hits:
        report.add("FAIL", name, " · ".join(hits[:6]))
    else:
        report.add("PASS", name, "no Opisthokonta/Amorphea/Obazoa/Craniata")


# ── (e) ran-on-disk vs source_status vs consumed-by-consensus ────────────


def check_source_status(run: Path, report: Report, fset: str) -> None:
    """Three-way: what RAN, what source_status ADMITS, what the consensus USED.

    Two legs are not enough. A source can produce output, be listed with an honest
    status, make real assignments — and still never reach the consensus; and the
    consensus can credit a source that assigned nothing. Both are invisible to any
    check that only compares status against assignments, and both are live bugs:
    NCBI ran with the right answer and was never listed, while BOLD was credited as
    Best_Source on all 7 features having assigned none of them.
    """
    check_name = f"(e/{fset}) source_status matches assignments"
    db = _db_dir(run, fset)
    status_path = db / "source_status.csv"
    assign_path = db / "source_assignments.csv"
    consensus_path = db / "taxonomy_consensus.csv"
    if not status_path.exists() or not assign_path.exists():
        report.add("SKIP", check_name, "missing input")
        return

    # source -> {feature: True} for REAL assignments (a named taxon at a real rank).
    assigned_by_source: dict[str, set[str]] = {}
    for row in _read_csv(assign_path):
        source = (row.get("source") or "").strip().upper()
        if not source:
            continue
        assigned_by_source.setdefault(source, set())
        if _norm(row.get("taxon")) and _norm(row.get("rank")):
            assigned_by_source[source].add((row.get("ASV_ID") or "").strip())
    counts = {source: len(features) for source, features in assigned_by_source.items()}

    problems = []

    # Leg 1: status vs the assignments that source actually made, and vs its raw output.
    reported = set()
    for row in _read_csv(status_path):
        source = (row.get("source") or "").strip().upper()
        status = (row.get("status") or "").strip().lower()
        reported.add(source)
        assigned = counts.get(source, 0)
        if status in ("skipped", "not_applicable"):
            # Never ran, or ran against the wrong marker: assignments say nothing here.
            continue
        if status == "ok" and assigned == 0:
            problems.append(f"{source}: status=ok but 0 assignments")
        elif status in ("no_hits", "failed") and assigned > 0:
            problems.append(f"{source}: status={status} but {assigned} assignments")
        elif status == "no_hits":
            # "no_hits" asserts the source ran and genuinely found nothing. With an
            # empty or missing raw output that assertion is unfounded — it is a failure,
            # and the two must not be spelled the same. A remote BLAST killed server-side
            # reported no_hits off a 0-byte result and an error-filled log.
            candidates = SOURCE_RAW_OUTPUTS.get(source, [])
            raw = next((db / name for name in candidates if (db / name).exists()), None)
            if candidates and (raw is None or raw.stat().st_size == 0):
                where = raw.name if raw else "/".join(candidates[0].split("/")[-1:])
                problems.append(f"{source}: status=no_hits but raw output is empty/absent ({where}) — that is a failure")

    # Leg 2: raw output present on disk vs what source_status admits to.
    #
    # This deliberately does NOT claim the source "ran and was dropped" — it cannot see
    # that. The Databases stage can be re-run into an existing run directory with a
    # different db_sources selection, which leaves the previous selection's raw output
    # behind; that is how a stale blast_result.tsv sat next to a source_status that
    # never mentions NCBI and read exactly like a discarded result. Both causes are
    # worth failing on (either way the directory contradicts the report), but the
    # message must state what is observed, not why.
    for source, candidates in SOURCE_RAW_OUTPUTS.items():
        raw = next((db / name for name in candidates if (db / name).exists()), None)
        if raw is None or source in reported:
            continue
        problems.append(
            f"{source}: raw output present ({raw.name}, {raw.stat().st_size} bytes) but absent from source_status"
            f" — this source was not run in the invocation that wrote source_status, or its output is stale"
            + (f"; it made {counts.get(source, 0)} assignments" if counts.get(source) else "")
        )

    # Leg 3: what the consensus CREDITS and USES vs what each source assigned, per feature.
    #
    # Two distinct questions, and conflating them produced a false positive worth
    # recording: Best_Source names the strongest *alignment*, and SINTAX sources (SILVA,
    # PR2) carry no identity, so they can never be Best_Source when a BLAST source is
    # present. "SILVA is never Best_Source" is therefore normal, not a bug — SILVA feeds
    # the LCA regardless. Participation is read from Consensus_Sources; attribution from
    # Best_Source. Per feature, not per run: a source that assigned something *somewhere*
    # must still not be credited for a feature it said nothing about.
    if consensus_path.exists():
        credited_nothing = []
        unused = {source for source, count in counts.items() if count > 0}
        for row in _read_csv(consensus_path):
            feature = (row.get("ASV_ID") or "").strip()
            for source in (s.strip().upper() for s in (row.get("Consensus_Sources") or "").split("|") if s.strip()):
                unused.discard(source)
            best = _norm(row.get("Best_Source")).upper()
            if not best:
                continue
            if feature not in assigned_by_source.get(best, set()):
                credited_nothing.append(f"{feature}<-{best}")
        if credited_nothing:
            problems.append(
                f"consensus credits a source with no assignment for that feature: {', '.join(credited_nothing[:5])}"
                + (f" (+{len(credited_nothing) - 5} more)" if len(credited_nothing) > 5 else "")
            )
        for source in sorted(unused):
            problems.append(f"{source}: made {counts[source]} assignments but fed no feature's consensus")

    if problems:
        report.add("FAIL", check_name, " · ".join(problems))
    else:
        summary = " · ".join(f"{source}={count}" for source, count in sorted(counts.items()))
        report.add("PASS", check_name, summary or "no sources")


# ── (f) a Linnean backbone actually reached the consensus ────────────────


# Sources that carry a full Linnean lineage (Phylum..Family). SILVA and PR2 are excluded
# deliberately: their own backbones are not Linnean (PR2 has Supergroup/Division/
# Subdivision), which is why forcing them into those columns invented Opisthokonta as a
# "Phylum" — and why, with them as the only sources, Phylum..Family come out empty.
LINNEAN_SOURCES = {"NCBI", "BOLD", "CUSTOM"}


def check_linnean_backbone(run: Path, report: Report, fset: str) -> None:
    """A run whose only Linnean source died is coherent — and useless. Fail it.

    The other five checks test whether the files agree with each other, and they all pass
    on a run where NCBI was refused: nothing contradicts anything, every rank is honestly
    empty, and the report reads 5/5. That is the one way this suite could mislead, so the
    question is asked directly: did anything actually deliver a Linnean lineage?
    """
    name = f"(f/{fset}) a Linnean backbone reached the consensus"
    consensus_path = _db_dir(run, fset) / "taxonomy_consensus.csv"
    if not consensus_path.exists():
        report.add("SKIP", name, "no taxonomy output")
        return

    rows = _read_csv(consensus_path)
    if not rows:
        report.add("FAIL", name, "taxonomy_consensus.csv is empty")
        return

    contributing = set()
    for row in rows:
        for source in (s.strip().upper() for s in (row.get("Consensus_Sources") or "").split("|") if s.strip()):
            contributing.add(source)
    linnean = contributing & LINNEAN_SOURCES

    # The rank columns those sources exist to fill.
    filled = sum(
        1 for row in rows if any(_norm(row.get(f"Consensus_{rank}")) for rank in ("Phylum", "Class", "Order", "Family"))
    )
    if linnean:
        report.add(
            "PASS",
            name,
            f"{'|'.join(sorted(linnean))} · {filled}/{len(rows)} features with Phylum..Family",
        )
        return

    # No Linnean source reached the consensus. Whether that is a FAILURE or simply
    # not-applicable depends on what was REQUESTED, read from provenance — not from what
    # arrived. An alarm that fires on a source nobody asked for (SILVA+PR2-only runs, where
    # empty Phylum..Family is a marker choice) trains the user to ignore it, and then it stops
    # guarding the case it exists for: a Linnean source that WAS requested and died. Same
    # distinction check (e) draws between skipped/not_applicable and failed, fed by run_config.
    requested = _requested_sources(run)
    requested_linnean = (requested & LINNEAN_SOURCES) if requested is not None else None
    if requested is not None and not requested_linnean:
        report.add(
            "SKIP",
            name,
            f"no Linnean source was requested (db_sources={'|'.join(sorted(requested)) or 'none'}); "
            "empty Phylum..Family is a config choice, not a failure",
        )
        return

    status_note = ""
    for row in _read_csv(_db_dir(run, fset) / "source_status.csv"):
        if (row.get("source") or "").strip().upper() in LINNEAN_SOURCES:
            status_note += f" · {row.get('source')}={row.get('status')}"
    asked = f"requested {'|'.join(sorted(requested_linnean))} but " if requested_linnean else ""
    report.add(
        "FAIL",
        name,
        f"{asked}no Linnean source fed any feature's consensus (saw: {'|'.join(sorted(contributing)) or 'none'})"
        f"{status_note} — a requested Linnean source did not reach the consensus",
    )


def check_results_reader(run: Path, report: Report) -> None:
    """The UI's reader must see the numbers this run actually wrote.

    Every other check compares two files. This one compares a file to the CODE THAT
    READS IT, because that is where the whole Phase 3/4 taxonomy can be correct on disk
    and absent from the user's screen — the reader looking at a path the pipeline does
    not write returns "nothing", which is indistinguishable from an honest empty run.

    It is deliberately anchored on the run's own frozen contract files rather than on a
    fixture: a fixture is written by whoever wrote the reader and can share its wrong
    belief about the layout, which is exactly how 401 tests passed while the Results page
    was blank.
    """
    name = "(g) the Results page reads this run"

    chimera = _read_csv(run / "analysis/metrics/otu/chimera_stats.csv")
    resolution = _read_csv(run / "databases/otu/metrics/taxonomic_resolution.csv")
    alpha = _read_csv(run / "analysis/metrics/otu/alpha_diversity.csv")
    if not chimera or not resolution or not alpha:
        report.add("SKIP", name, "run lacks the analysis/taxonomy metrics this reads")
        return

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from analysis.results_summary import clustered_read_total, sample_dashboard, species_classified_total

    # What the run's own tables say, computed here from the files directly.
    total = float(chimera[0].get("total_reads") or 0)
    chimeric = float(chimera[0].get("chimeric_reads") or 0)
    want_clustered = int(max(total - chimeric, 0))
    want_classified = next((int(float(r.get("n_reads") or 0)) for r in resolution if (r.get("rank") or "") == "species"), 0)
    want_richness = {r["sample"]: int(float(r.get("observed") or 0)) for r in alpha}

    problems = []

    got_clustered = clustered_read_total(run)
    if got_clustered != want_clustered:
        problems.append(f"clustered_read_total -> {got_clustered}, chimera_stats says {want_clustered}")

    got_classified = species_classified_total(run)
    if got_classified != want_classified:
        problems.append(f"species_classified_total -> {got_classified}, taxonomic_resolution says {want_classified}")

    rows = {row["sample"]: row for row in sample_dashboard(run)}
    for sample, expected in want_richness.items():
        got = (rows.get(sample) or {}).get("richness")
        if got != expected:
            problems.append(f"dashboard[{sample}].richness -> {got}, alpha_diversity says {expected}")
    # A run that classified anything must not report every sample at 0% to species.
    if want_classified and rows and not any((row.get("pct_species") or 0) > 0 for row in rows.values()):
        problems.append(f"dashboard pct_species is 0 for every sample, yet {want_classified} reads reached species")

    if problems:
        report.add(
            "FAIL",
            name,
            "the reader does not see this run's output:\n        - " + "\n        - ".join(problems),
        )
        return
    report.add(
        "PASS",
        name,
        f"clustered={got_clustered} · classified={got_classified} · richness matches alpha_diversity on {len(want_richness)} sample(s)",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    run = args.run_dir
    if not run.exists():
        print(f"FAIL  run dir does not exist: {run}")
        return 2

    print(f"verify: {run}")
    report = Report(quiet=args.quiet)
    # (a)–(f) run PER FEATURE SET actually present, so an ASV run is verified rather than
    # skipped into a false green. Default to a single OTU pass when neither set exists, so an
    # empty run still SKIPs loudly instead of producing no checks at all.
    for fset in _present_sets(run) or ["otu"]:
        check_read_totals(run, report, fset)
        check_otu_id_map(run, report, fset)
        check_class_table(run, report, fset)
        check_banned_ranks(run, report, fset)
        check_source_status(run, report, fset)
        check_linnean_backbone(run, report, fset)
    # (g) tests the Results reader, which is chimera-based (OTU-only) — it resolves its own
    # paths, so it runs once rather than per set.
    check_results_reader(run, report)

    passed = sum(1 for s, _, _ in report.rows if s == "PASS")
    print(f"\n{passed} passed · {report.failed} failed · {report.skipped} skipped")
    return 1 if report.failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
