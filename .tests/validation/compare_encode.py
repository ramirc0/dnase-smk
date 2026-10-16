"""Print the validation run's QC and calls beside ENCODE's.

Two references:

- ENCODE's v3 QC, saved in ``encode_qc.tsv``;
- ENCODE pipeline ENCPL336PGV's 2025 files in the Weng Lab mirror: hotspot,
  peak and footprint counts and overlaps, and the correlation of the
  normalized density bigWigs.

Run from the repository root::

    pixi run --manifest-path workflow/envs/pixi.toml -e test \
        python .tests/validation/compare_encode.py results/encode-validation
"""

import argparse
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
import pyBigWig

HERE = Path(__file__).parent
MIRROR = Path("/data/projects/encode/data")

# Our section and metric, ENCODE's section and metric, and a factor on ours.
QC = [
    ("flagstat.unfiltered", "passed.total", "flagstat.unfiltered", "total", 1),
    ("flagstat.unfiltered", "passed.mapped", "flagstat.unfiltered", "mapped", 1),
    (
        "flagstat.unfiltered",
        "passed.properly paired",
        "flagstat.unfiltered",
        "paired_properly",
        1,
    ),
    ("flagstat.nuclear", "passed.total", "flagstat.nuclear", "total", 1),
    ("flagstat.nuclear", "passed.duplicates", "flagstat.nuclear", "duplicates", 1),
    ("duplicates", "PERCENT_DUPLICATION", "duplicates", "Percent Duplication", 1),
    ("duplicates", "DUPLICATE PAIR", "duplicates", "Read Pair Duplicates", 0.5),
    ("duplicates", "DUPLICATE SINGLE", "duplicates", "Unpaired Read Duplicates", 1),
    ("duplicates", "ESTIMATED_LIBRARY_SIZE", "duplicates", "Estimated Library Size", 1),
    ("insert_sizes", "MEDIAN_INSERT_SIZE", "insert_sizes", "MEDIAN_INSERT_SIZE", 1),
    ("insert_sizes", "MEAN_INSERT_SIZE", "insert_sizes", "MEAN_INSERT_SIZE", 1),
    ("insert_sizes", "STANDARD_DEVIATION", "insert_sizes", "STANDARD_DEVIATION", 1),
    ("insert_sizes", "READ_PAIRS", "insert_sizes", "READ_PAIRS", 1),
    ("hotspot3", "SPOT", "hotspot1", "SPOT1", 1),
]

# ENCPL336PGV's files of each validation replicate, by experiment and kind.
ENCODE_2025 = {
    "UHT_rep1": (
        "ENCSR693UHT",
        {
            "hotspots.fdr0.05": "ENCFF553ZVW",
            "peaks.fdr0.01": "ENCFF961ETW",
            "peaks.fdr0.001": "ENCFF116FRB",
            "footprints.fdr0.01": "ENCFF845KSE",
            "footprints.fdr0.001": "ENCFF478LTM",
            "density": "ENCFF189QHJ",
        },
    ),
    "IEA_rep1": (
        "ENCSR426IEA",
        {
            "hotspots.fdr0.05": "ENCFF856ZDV",
            "peaks.fdr0.01": "ENCFF121CQK",
            "peaks.fdr0.001": "ENCFF176PEH",
            "footprints.fdr0.01": "ENCFF568GRM",
            "footprints.fdr0.001": "ENCFF233PFY",
            "density": "ENCFF233QNW",
        },
    ),
    "IEA_rep2": (
        "ENCSR426IEA",
        {
            "hotspots.fdr0.05": "ENCFF799DGR",
            "peaks.fdr0.01": "ENCFF809TPU",
            "peaks.fdr0.001": "ENCFF732GEQ",
            "density": "ENCFF663CEC",
        },
    ),
    "RWU_rep1": (
        "ENCSR420RWU",
        {
            "hotspots.fdr0.05": "ENCFF390JSJ",
            "peaks.fdr0.01": "ENCFF660HML",
            "peaks.fdr0.001": "ENCFF921RUU",
            "density": "ENCFF255OAP",
        },
    ),
}


def build_parser():
    """Build the compare_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the compare_encode.py arguments.
    """
    p = argparse.ArgumentParser(description="Compare a validation run with ENCODE.")
    p.add_argument("run", type=Path, help="Run directory, results/<run_id>.")
    p.add_argument(
        "--bin", type=int, default=1000, help="Bin size for the density correlation."
    )
    return p


def read_table(path):
    """Read a ``sample, section, metric, value`` TSV.

    Parameters
    ----------
    path : pathlib.Path
        QC table, ``#`` lines skipped.

    Returns
    -------
    dict of tuple to str
        Value by (sample, section, metric).
    """
    table = pl.read_csv(path, separator="\t", comment_prefix="#", infer_schema_length=0)
    return {(r[0], r[1], r[2]): r[3] for r in table.iter_rows()}


def qc_rows(ours, theirs):
    """Pair our QC with ENCODE's v3 QC.

    Parameters
    ----------
    ours, theirs : dict of tuple to str
        Values by (sample, section, metric).

    Returns
    -------
    polars.DataFrame
        One row per sample and metric both have, with their ratio.
    """
    rows = []
    for sample in sorted({k[0] for k in theirs}):
        for section, metric, e_section, e_metric, factor in QC:
            a, b = (
                ours.get((sample, section, metric)),
                theirs.get((sample, e_section, e_metric)),
            )
            if a is None or b is None:
                continue
            a, b = float(a) * factor, float(b)
            rows.append((sample, f"{section}.{metric}", a, b, a / b if b else None))
    for (sample, section, metric), b in sorted(theirs.items()):
        if section == "preseq" and (sample, section, metric) in ours:
            a = float(ours[sample, section, metric])
            rows.append((sample, f"preseq.{metric}", a, float(b), a / float(b)))
    return pl.DataFrame(
        rows, schema=["sample", "metric", "ours", "encode_v3", "ratio"], orient="row"
    )


def bed_count(path):
    """Count the intervals of a BED, gzipped or not.

    Parameters
    ----------
    path : pathlib.Path
        BED file. ``#`` lines are not counted.

    Returns
    -------
    int
        Number of intervals.
    """
    out = subprocess.run(
        f"zcat -f {path} | grep -vc '^#'",
        shell=True,
        capture_output=True,
        text=True,
        check=False,
    )
    return int(out.stdout)


def overlapping(a, b):
    """Count the intervals of one BED that overlap another.

    Parameters
    ----------
    a, b : pathlib.Path
        BED files.

    Returns
    -------
    int
        Intervals of `a` with at least one base in `b`.
    """
    out = subprocess.run(
        f"bedtools intersect -u -a <(zcat -f {a} | grep -v '^#') -b <(zcat -f {b} | grep -v '^#') | wc -l",
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
        check=True,
    )
    return int(out.stdout)


def ours_path(run, sample, kind):
    """Return our file of a kind from `ENCODE_2025`.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    sample : str
        Sample ID.
    kind : str
        ``hotspots.fdrF``, ``peaks.fdrF``, ``footprints.fdrF`` or ``density``.

    Returns
    -------
    pathlib.Path
        Our file.
    """
    if kind == "density":
        return run / f"{sample}/hotspot3/{sample}.normalized_density.bw"
    what, fdr = kind.split(".", 1)
    if what == "footprints":
        return run / f"{sample}/footprints/{sample}.footprints.{fdr}.bed.gz"
    return run / f"{sample}/hotspot3/{fdr}/{sample}.{kind}.bed.gz"


def binned(path, bin_size, chroms):
    """Average a bigWig in fixed bins.

    Parameters
    ----------
    path : pathlib.Path
        Signal track.
    bin_size : int
        Bin width in bp.
    chroms : dict of str to int
        Contigs and sizes to bin.

    Returns
    -------
    numpy.ndarray
        Mean per bin, contigs in order, empty bins as 0.
    """
    with pyBigWig.open(str(path)) as bw:
        parts = [
            bw.stats(c, 0, size, nBins=size // bin_size, type="mean")
            for c, size in chroms.items()
        ]
    return np.nan_to_num(np.array([v for p in parts for v in p], dtype=float))


def density_rows(run, bin_size):
    """Correlate our normalized density with ENCODE's 2025 one.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    bin_size : int
        Bin width in bp.

    Returns
    -------
    polars.DataFrame
        Pearson of log1p means and Spearman of means, per sample.
    """
    rows = []
    for sample, (experiment, files) in ENCODE_2025.items():
        ours = ours_path(run, sample, "density")
        theirs = MIRROR / experiment / f"{files['density']}.bigWig"
        if not ours.exists() or not theirs.exists():
            continue
        with pyBigWig.open(str(theirs)) as bw:
            chroms = {c: s for c, s in bw.chroms().items() if c != "chrY"}
        a, b = binned(ours, bin_size, chroms), binned(theirs, bin_size, chroms)
        pearson = np.corrcoef(np.log1p(a), np.log1p(b))[0, 1]
        ranks = pl.DataFrame({"a": a, "b": b}).select(pl.all().rank())
        spearman = np.corrcoef(ranks["a"], ranks["b"])[0, 1]
        rows.append((sample, len(a), pearson, spearman))
    return pl.DataFrame(
        rows, schema=["sample", "bins", "pearson_log1p", "spearman"], orient="row"
    )


def call_rows(run):
    """Count and overlap our calls with ENCODE's 2025 ones.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.

    Returns
    -------
    polars.DataFrame
        Counts of both, and the fraction of each set that overlaps the other.
    """
    rows = []
    for sample, (experiment, files) in ENCODE_2025.items():
        for kind, accession in files.items():
            ours, theirs = (
                ours_path(run, sample, kind),
                MIRROR / experiment / f"{accession}.bed.gz",
            )
            if kind == "density" or not ours.exists() or not theirs.exists():
                continue
            n_ours, n_theirs = bed_count(ours), bed_count(theirs)
            rows.append(
                (
                    sample,
                    kind,
                    n_ours,
                    n_theirs,
                    overlapping(theirs, ours) / n_theirs,
                    overlapping(ours, theirs) / n_ours,
                )
            )
    return pl.DataFrame(
        rows,
        schema=[
            "sample",
            "kind",
            "ours",
            "encode_2025",
            "encode_in_ours",
            "ours_in_encode",
        ],
        orient="row",
    )


def main(argv=None):
    """Print the QC, call and density comparisons.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    ours = read_table(args.run / "qc.tsv")
    theirs = read_table(HERE / "encode_qc.tsv")
    with pl.Config(tbl_rows=-1, tbl_width_chars=200, fmt_str_lengths=60):
        print("QC against ENCODE v3\n", qc_rows(ours, theirs))
        print("Calls against ENCODE 2025\n", call_rows(args.run))
        print(
            f"Normalized density against ENCODE 2025, {args.bin} bp bins\n",
            density_rows(args.run, args.bin),
        )


if __name__ == "__main__":
    main()
