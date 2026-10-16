"""Write the duplicate histogram that preseq reads, as stampipes' mark_dups.py --hist.

Records MUST be in coordinate order. A template is counted once, at its first
record. Templates whose first records share a contig, a start and a TLEN form
one set, whatever their strand. The histogram counts sets by size.
"""

import argparse
import itertools
from collections import Counter

import pysam


def build_parser():
    """Build the dup_hist.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the dup_hist.py arguments.
    """
    p = argparse.ArgumentParser(description="Write preseq's duplicate histogram.")
    p.add_argument("bam", help="Coordinate-sorted BAM.")
    p.add_argument("output", help="Histogram to write: set size, then sets.")
    return p


def histogram(bam):
    """Count duplicate sets by size.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Open coordinate-sorted BAM.

    Returns
    -------
    collections.Counter
        Number of sets by set size.
    """
    hist = Counter()
    # Paired templates whose first record was counted and whose mate is still ahead.
    waiting = set()
    starts = itertools.groupby(
        bam.fetch(until_eof=True), lambda r: (r.reference_id, r.reference_start)
    )
    for _, records in starts:
        sets = Counter()
        mates = []
        for r in records:
            if r.query_name in waiting:
                mates.append(r.query_name)
                continue
            if r.is_paired:
                waiting.add(r.query_name)
            sets[r.template_length] += 1
        hist.update(sets.values())
        waiting.difference_update(mates)
    return hist


def main(argv=None):
    """Write the duplicate histogram of a BAM.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        hist = histogram(bam)
    with open(args.output, "w") as out:
        out.writelines(f"{size}\t{hist[size]}\n" for size in sorted(hist))


if __name__ == "__main__":
    main()
