"""Tally read counts of a BAM under stampipes' bamcounts.py labels.

The output holds one ``label<TAB>count`` line per count, in five blocks, each
sorted by label: the filter counts, reads per contig, then MAPQ, flag and read
length tallies.
"""

import argparse
from collections import Counter

import pysam

LABELS = (
    "u",
    "u-pf",
    "u-pf-n",
    "u-pf-n-mm{mm}",
    "u-pf-n-mm{mm}-mito",
    "mm",
    "nm",
    "qc-flagged",
    "duplicate",
    "duplicate-nuclear",
    "nuclear-align",
    "autosomal-align",
    "paired-aligned",
    "paired-nuclear-align",
    "paired-autosomal-align",
    "all-aligned",
    "all-mapq-filter",
)


def build_parser():
    """Build the bamcounts.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the bamcounts.py arguments.
    """
    p = argparse.ArgumentParser(description="Tally stampipes' BAM read counts.")
    p.add_argument("bam", help="BAM to count.")
    p.add_argument("output", help="Counts file to write.")
    p.add_argument("--min_mapping_quality", type=int, default=10, help="Lowest MAPQ.")
    p.add_argument("--max_mismatches", type=int, default=2, help="Highest NM.")
    return p


def paired_passes(read, counts):
    """Return whether a paired record counts as a properly paired alignment.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Mapped paired record that passed the MAPQ filter.
    counts : collections.Counter
        Filter counts, updated in place.

    Returns
    -------
    bool
        True when both mates are mapped, 0x2 is set and 0x200 is not.
    """
    if read.flag & 0xC or not read.flag & 0x2:
        return False
    if read.is_qcfail:
        counts["paired-aligned-qcfail"] += 1
        return False
    counts["paired-aligned"] += 1
    return True


def tally(bam, min_mapq, max_mismatches):
    """Count a BAM's records as stampipes' bamcounts.py does.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Open BAM.
    min_mapq : int
        Lowest passing MAPQ.
    max_mismatches : int
        Highest passing NM.

    Returns
    -------
    list of collections.Counter
        Filter counts, reads per contig, and the MAPQ, flag and read length
        tallies.
    """
    mm = f"u-pf-n-mm{max_mismatches}"
    counts = Counter({label.format(mm=max_mismatches): 0 for label in LABELS})
    contigs, mapqs, flags, lengths = Counter(), Counter(), Counter(), Counter()
    for read in bam.fetch(until_eof=True):
        mapqs[f"mapq-{read.mapping_quality}"] += 1
        flags[f"samflag-{read.flag}"] += 1
        lengths[f"readlength-{read.query_length}"] += 1
        if read.has_tag("XT") and read.get_tag("XT") == "R":
            counts["mm"] += 1
        if read.is_qcfail:
            counts["qc-flagged"] += 1
        if read.is_unmapped:
            counts["nm"] += 1
            continue
        counts["all-aligned"] += 1
        if read.mapping_quality < min_mapq:
            counts["all-mapq-filter"] += 1
            continue
        if read.is_paired and not paired_passes(read, counts):
            continue
        contig = read.reference_name
        nuclear = contig not in ("chrM", "chrC")
        autosomal = nuclear and contig not in ("chrX", "chrY", "chrZ", "chrW")
        if nuclear:
            counts["nuclear-align"] += 1
            counts["autosomal-align"] += autosomal
            if read.is_paired:
                counts["paired-nuclear-align"] += 1
                counts["paired-autosomal-align"] += autosomal
        if read.is_duplicate:
            counts["duplicate"] += 1
            counts["duplicate-nuclear"] += nuclear
        counts["u"] += 1
        lengths[f"aligned-readlength-{read.query_length}"] += 1
        if read.is_qcfail:
            continue
        counts["u-pf"] += 1
        if "N" in (read.query_sequence or ""):
            continue
        counts["u-pf-n"] += 1
        if read.has_tag("NM") and read.get_tag("NM") > max_mismatches:
            continue
        counts[mm] += 1
        contigs[contig] += 1
        if contig != "chrM":
            counts[f"{mm}-mito"] += 1
    return [counts, contigs, mapqs, flags, lengths]


def main(argv=None):
    """Count a BAM and write the tallies.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam:
        blocks = tally(bam, args.min_mapping_quality, args.max_mismatches)
    with open(args.output, "w") as out:
        for block in blocks:
            out.writelines(f"{label}\t{block[label]}\n" for label in sorted(block))


if __name__ == "__main__":
    main()
