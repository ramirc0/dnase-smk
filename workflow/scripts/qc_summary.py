"""Gather a run's QC into one long table and one figure.

Rows are ``sample, section, metric, value``. Every input file is named
``<sample>.<rest>``, and sample names never contain dots.
"""

import argparse
import gzip
import json
from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl

# Panels of the summary figure: section, metric, title.
PANELS = [
    ("flagstat.nuclear", "passed.total", "Nuclear reads"),
    ("duplicates", "PERCENT_DUPLICATION", "Duplicate fraction"),
    ("preseq", "preseq-est-max", "Preseq distinct reads"),
    ("hotspot3", "SPOT", "SPOT"),
    ("hotspot3", "peaks.fdr0.001", "Peaks, FDR 0.001"),
    ("footprints", "fdr0.001", "Footprints, FDR 0.001"),
]


def build_parser():
    """Build the qc_summary.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc_summary.py flags.
    """
    p = argparse.ArgumentParser(description="Gather a run's QC.")
    for flag, what in (
        ("--cutadapt", "cutadapt JSON, <sample>.<layout>.cutadapt.json"),
        ("--flagstat", "samtools flagstat JSON, <sample>.<bam>.flagstat.json"),
        ("--stats", "samtools stats, <sample>.<bam>.stats.txt"),
        ("--bamcounts", "bamcounts.py counts"),
        ("--markdup", "samtools markdup JSON"),
        ("--preseq", "preseq targets"),
        ("--insert-sizes", "insert_sizes.py TSV"),
        ("--fit-stats", "hotspot3 fit_stats.bed.gz"),
        ("--total-cutcounts", "hotspot3 total_cutcounts"),
        ("--calls", "hotspot3 hotspots and peaks, <sample>.<kind>.fdr<F>.bed.gz"),
        ("--footprints", "footprints, <sample>.footprints.fdr<F>.bed.gz"),
        ("--dispersion-model", "ftd learn_dm JSON"),
    ):
        p.add_argument(flag, nargs="*", default=[], help=f"{what} files.")
    p.add_argument("--table", required=True, help="TSV to write.")
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def name_parts(path):
    """Split a QC file name at its dots.

    Parameters
    ----------
    path : str
        File named ``<sample>.<rest>``.

    Returns
    -------
    list of str
        The name's dot-separated parts, the sample first.
    """
    return Path(path).name.split(".")


def tsv_pairs(path):
    """Read a two-column ``key<TAB>value`` file.

    Parameters
    ----------
    path : str
        File to read. A ``metric<TAB>value`` header is skipped.

    Returns
    -------
    list of tuple of str
        One (key, value) per line.
    """
    with open(path) as f:
        pairs = [tuple(line.rstrip("\n").split("\t")[:2]) for line in f if line.strip()]
    return [p for p in pairs if p != ("metric", "value")]


def cutadapt_rows(path):
    """Read a cutadapt JSON report as trimming rows.

    Parameters
    ----------
    path : str
        ``<sample>.<layout>.cutadapt.json``.

    Returns
    -------
    list of tuple
        ``(sample, section, metric, value)`` rows under ``trimming.<layout>``.
    """
    sample, layout = name_parts(path)[:2]
    with open(path) as f:
        report = json.load(f)
    counts = report["read_counts"]
    metrics = {
        "input_reads": counts["input"],
        "output_reads": counts["output"],
        "casava_filtered": counts["filtered"]["casava_filtered"],
        "too_short": counts["filtered"]["too_short"],
        "read1_with_adapter": counts["read1_with_adapter"],
        "read2_with_adapter": counts["read2_with_adapter"],
        "input_bp": report["basepair_counts"]["input"],
        "output_bp": report["basepair_counts"]["output"],
    }
    return [(sample, f"trimming.{layout}", k, v) for k, v in metrics.items()]


def flagstat_rows(path):
    """Read a samtools flagstat JSON.

    Parameters
    ----------
    path : str
        ``<sample>.<bam>.flagstat.json``.

    Returns
    -------
    list of tuple
        ``(sample, section, metric, value)`` rows under ``flagstat.<bam>``,
        with ``passed.`` and ``failed.`` metrics for QC-passed and QC-failed
        reads.
    """
    sample, bam = name_parts(path)[:2]
    with open(path) as f:
        stats = json.load(f)
    return [
        (sample, f"flagstat.{bam}", f"{label}.{k}", v)
        for label, key in (("passed", "QC-passed reads"), ("failed", "QC-failed reads"))
        for k, v in stats[key].items()
    ]


def stats_rows(path):
    """Read the summary numbers of a samtools stats report.

    Parameters
    ----------
    path : str
        ``<sample>.<bam>.stats.txt``.

    Returns
    -------
    list of tuple
        ``(sample, section, metric, value)`` rows under ``stats.<bam>``, one
        per ``SN`` line.
    """
    sample, bam = name_parts(path)[:2]
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("SN\t"):
                _, key, value = line.rstrip("\n").split("\t")[:3]
                rows.append((sample, f"stats.{bam}", key.rstrip(":"), value))
    return rows


def markdup_rows(path):
    """Read a samtools markdup JSON, with Picard's PERCENT_DUPLICATION.

    Parameters
    ----------
    path : str
        ``<sample>.markdup.json``.

    Returns
    -------
    list of tuple
        ``(sample, section, metric, value)`` rows under ``duplicates``.
        markdup counts pairs in reads, so Picard's formula reduces to
        duplicate reads over examined reads.
    """
    sample = name_parts(path)[0]
    with open(path) as f:
        stats = json.load(f)
    stats.pop("COMMAND")
    examined = stats["PAIRED"] + stats["SINGLE"]
    duplicates = stats["DUPLICATE PAIR"] + stats["DUPLICATE SINGLE"]
    stats["PERCENT_DUPLICATION"] = duplicates / examined if examined else 0.0
    return [(sample, "duplicates", k, v) for k, v in stats.items()]


def pair_rows(section, path):
    """Read a two-column file as rows under one section.

    Parameters
    ----------
    section : str
        Section of the rows.
    path : str
        ``<sample>.<rest>`` file of ``key<TAB>value`` lines.

    Returns
    -------
    list of tuple
        ``(sample, section, metric, value)`` rows.
    """
    sample = name_parts(path)[0]
    return [(sample, section, k, v) for k, v in tsv_pairs(path)]


def fit_stats_rows(path):
    """Read hotspot3's genome-wide SPOT score from its fit stats.

    Parameters
    ----------
    path : str
        ``<sample>.fit_stats.bed.gz``. Every row repeats the genome-wide SPOT.

    Returns
    -------
    list of tuple
        ``SPOT`` and ``SPOT_std`` rows under ``hotspot3``.
    """
    sample = name_parts(path)[0]
    with gzip.open(path, "rt") as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        row = dict(zip(header, f.readline().rstrip("\n").split("\t"), strict=True))
    return [(sample, "hotspot3", k, row[k]) for k in ("SPOT", "SPOT_std")]


def total_rows(path):
    """Read hotspot3's total cutcounts.

    Parameters
    ----------
    path : str
        ``<sample>.total_cutcounts``.

    Returns
    -------
    list of tuple
        One ``total_cutcounts`` row under ``hotspot3``.
    """
    return [
        (
            name_parts(path)[0],
            "hotspot3",
            "total_cutcounts",
            Path(path).read_text().strip(),
        )
    ]


def count_rows(section, name, path):
    """Count the intervals of a BED.gz.

    Parameters
    ----------
    section : str
        Section of the row.
    name : callable
        Takes the file name's parts and returns the metric name.
    path : str
        ``<sample>.<rest>.bed.gz``. ``#`` lines are not counted.

    Returns
    -------
    list of tuple
        One ``(sample, section, metric, count)`` row.
    """
    with gzip.open(path, "rt") as f:
        n = sum(1 for line in f if not line.startswith("#"))
    parts = name_parts(path)
    return [(parts[0], section, name(parts), n)]


def dispersion_rows(path):
    """Record where a sample's dispersion model is.

    Parameters
    ----------
    path : str
        ``<sample>.dm.json``.

    Returns
    -------
    list of tuple
        One ``dispersion_model`` row under ``footprints``.
    """
    return [(name_parts(path)[0], "footprints", "dispersion_model", path)]


def main(argv=None):
    """Gather every QC file into the table, then draw the summary panels.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    readers = (
        (args.cutadapt, cutadapt_rows),
        (args.flagstat, flagstat_rows),
        (args.stats, stats_rows),
        (args.bamcounts, partial(pair_rows, "bamcounts")),
        (args.markdup, markdup_rows),
        (args.preseq, partial(pair_rows, "preseq")),
        (args.insert_sizes, partial(pair_rows, "insert_sizes")),
        (args.fit_stats, fit_stats_rows),
        (args.total_cutcounts, total_rows),
        # <sample>.<kind>.fdr0.<digits>.bed.gz
        (
            args.calls,
            partial(count_rows, "hotspot3", lambda p: f"{p[1]}.{p[2]}.{p[3]}"),
        ),
        (
            args.footprints,
            partial(count_rows, "footprints", lambda p: f"{p[2]}.{p[3]}"),
        ),
        (args.dispersion_model, dispersion_rows),
    )
    rows = [row for paths, read in readers for path in paths for row in read(path)]
    table = pl.DataFrame(
        [(s, sec, m, "" if v is None else str(v)) for s, sec, m, v in rows],
        schema=["sample", "section", "metric", "value"],
        orient="row",
    )
    table.write_csv(args.table, separator="\t")
    plot(table, args.plot)


def plot(table, path):
    """Draw one bar panel per summary metric, a bar per sample.

    Parameters
    ----------
    table : polars.DataFrame
        The QC table.
    path : str
        SVG to write. The PNG goes beside it.
    """
    pk.use("plain")
    fig, axes = plt.subplots(2, len(PANELS) // 2, figsize=(4 * len(PANELS) // 2, 6))
    for ax, (section, metric, title) in zip(axes.flat, PANELS, strict=True):
        values = table.filter(
            (pl.col("section") == section) & (pl.col("metric") == metric)
        ).with_columns(pl.col("value").cast(pl.Float64, strict=False))
        ax.barh(values["sample"].to_list(), values["value"].fill_null(0).to_list())
        ax.set(title=title)
        ax.invert_yaxis()
        pk.despine(ax, categorical_y=True)
    pk.save_figure(fig, path)
    plt.close(fig)


if __name__ == "__main__":
    main()
