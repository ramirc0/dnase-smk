"""Compute Picard's insert size metrics and stampipes' two derived values from a histogram.

The metrics follow Picard CollectInsertSizeMetrics 2.8.1: median, MAD, the
widths, minimum and maximum come from the whole histogram. The histogram is
then trimmed to median + 10 MAD, and the mean, standard deviation and
stampipes' ``insert-ls-ratio`` and ``insert-ft-eleven`` come from the rest.
"""

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk

DEVIATIONS = 10
WIDTHS = (10, 20, 30, 40, 50, 60, 70, 80, 90, 99)


def build_parser():
    """Build the insert_sizes.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the insert_sizes.py arguments.
    """
    p = argparse.ArgumentParser(description="Compute insert size metrics.")
    p.add_argument("histogram", help="Two columns: insert size, then pairs.")
    p.add_argument("--output", required=True, help="Metrics TSV to write.")
    p.add_argument("--plot", required=True, help="Histogram SVG; a PNG goes beside it.")
    return p


def read_histogram(path):
    """Read an insert size histogram.

    Parameters
    ----------
    path : str
        Two whitespace-separated columns: insert size, then pairs.

    Returns
    -------
    dict of int to int
        Pairs by insert size, in increasing size.
    """
    with open(path) as f:
        rows = (line.split() for line in f)
        return dict(sorted((int(k), int(v)) for k, v in rows))


def median(hist):
    """Return the median of a histogram, as htsjdk's Histogram.getMedian.

    Parameters
    ----------
    hist : dict of int to int
        Counts by value, in increasing value.

    Returns
    -------
    float
        Mean of the two middle values when the count is even.
    """
    count = sum(hist.values())
    if count == 0:
        return 0.0
    low = count / 2 if count % 2 == 0 else math.ceil(count / 2)
    high = low + 1 if count % 2 == 0 else low
    total = 0
    low_value = high_value = None
    for value, n in hist.items():
        total += n
        if low_value is None and total >= low:
            low_value = value
        if high_value is None and total >= high:
            high_value = value
            break
    return (low_value + high_value) / 2


def picard_metrics(hist):
    """Compute Picard's InsertSizeMetrics fields, plus the mode.

    Parameters
    ----------
    hist : dict of int to int
        Pairs by insert size, in increasing size.

    Returns
    -------
    dict of str to float
        Metric by Picard column name, and ``MODE_INSERT_SIZE``.
    trimmed : dict of int to int
        `hist` without the sizes above median + 10 MAD.
    """
    total = sum(hist.values())
    mid = median(hist)
    deviations = {}
    for value, n in hist.items():
        d = abs(value - mid)
        deviations[d] = deviations.get(d, 0) + n
    mad = median(dict(sorted(deviations.items())))
    metrics = {
        "MEDIAN_INSERT_SIZE": mid,
        "MODE_INSERT_SIZE": max(hist, key=hist.get),
        "MEDIAN_ABSOLUTE_DEVIATION": mad,
        "MIN_INSERT_SIZE": min(hist),
        "MAX_INSERT_SIZE": max(hist),
    }
    low = high = mid
    covered = 0
    widths = {}
    while low >= min(hist) or high <= max(hist):
        covered += hist.get(int(low), 0)
        if low != high:
            covered += hist.get(int(high), 0)
        for p in WIDTHS:
            if p not in widths and covered / total >= p / 100:
                widths[p] = int(high - low) + 1
        low -= 1
        high += 1
    trimmed = {k: v for k, v in hist.items() if k <= int(mid + DEVIATIONS * mad)}
    n = sum(trimmed.values())
    mean = sum(k * v for k, v in trimmed.items()) / n
    var = sum(v * (k - mean) ** 2 for k, v in trimmed.items()) / (n - 1)
    metrics |= {
        "MEAN_INSERT_SIZE": mean,
        "STANDARD_DEVIATION": math.sqrt(var),
        "READ_PAIRS": total,
    }
    metrics |= {f"WIDTH_OF_{p}_PERCENT": widths.get(p, 0) for p in WIDTHS}
    return metrics, trimmed


def stampipes_metrics(hist):
    """Compute stampipes' insert-ls-ratio and insert-ft-eleven.

    Parameters
    ----------
    hist : dict of int to int
        Trimmed insert size histogram.

    Returns
    -------
    dict of str to float
        ``insert-ls-ratio``, pairs at 120 over pairs at 60, and
        ``insert-ft-eleven``, the real part of the 11th Fourier coefficient,
        less the mean coefficient, of the differenced max-normalized histogram
        over the 105 sizes from its mode. NaN where undefined.
    """
    counts = np.zeros(max(hist) + 1)
    for k, v in hist.items():
        counts[k] = v
    ratio = counts[120] / counts[60] if len(counts) > 120 and counts[60] else math.nan
    norm = counts / counts.max()
    window = norm[norm.argmax() : norm.argmax() + 105]
    ft = np.fft.rfft(np.diff(window))
    eleven = (ft[11] - ft.mean()).real if len(ft) > 11 else math.nan
    return {"insert-ls-ratio": ratio, "insert-ft-eleven": eleven}


def plot(hist, path):
    """Plot an insert size histogram.

    Parameters
    ----------
    hist : dict of int to int
        Pairs by insert size.
    path : str
        SVG to write; a PNG goes beside it.
    """
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(list(hist), list(hist.values()))
    ax.set(xlabel="insert size (bp)", ylabel="read pairs")
    pk.despine(ax)
    pk.save_figure(fig, path)


def main(argv=None):
    """Write the insert size metrics and plot of a histogram.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    hist = read_histogram(args.histogram)
    metrics, trimmed = picard_metrics(hist)
    metrics |= stampipes_metrics(trimmed)
    with open(args.output, "w") as out:
        out.write("metric\tvalue\n")
        out.writelines(f"{name}\t{value}\n" for name, value in metrics.items())
    plot(trimmed, args.plot)


if __name__ == "__main__":
    main()
