"""Flag reads that fail QC (0x200) and proper pairs (0x2), as stampipes' filter_reads.py.

Records MUST arrive grouped by name, as an aligner writes them. Each group is
judged on its primary records:

- a record fails when unmapped, with MAPQ below the minimum, with more than the
  maximum edit distance, or with its 5' end soft-clipped;
- a pair fails unless both mates pass and face each other on one contig, with
  nonzero insert sizes that sum to zero and are at most the maximum;
- a lone mate of a paired read fails.

A passing pair gets 0x2 set and 0x200 cleared. Every other primary record gets
0x2 cleared, and 0x200 set when it fails. Secondary and supplementary records
get 0x200 set and 0x2 cleared.
"""

import argparse
import itertools

import pysam

PROPER = 0x2
QCFAIL = 0x200
NOT_PRIMARY = 0x900


def build_parser():
    """Build the flag_reads.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the flag_reads.py arguments.
    """
    p = argparse.ArgumentParser(description="Set QC-fail and proper-pair flags.")
    p.add_argument("input", help="Name-grouped BAM.")
    p.add_argument("output", help="BAM to write, or - for uncompressed BAM on stdout.")
    p.add_argument("--min_mapq", type=int, default=10, help="Lowest passing MAPQ.")
    p.add_argument("--max_mismatches", type=int, default=2, help="Highest passing NM.")
    p.add_argument(
        "--max_insert_size", type=int, default=750, help="Highest passing |TLEN|."
    )
    return p


def read_passes(read, min_mapq, max_mismatches):
    """Return whether a primary record passes the single-record rules.

    Parameters
    ----------
    read : pysam.AlignedSegment
        Primary record.
    min_mapq : int
        Lowest passing MAPQ.
    max_mismatches : int
        Highest passing NM.

    Returns
    -------
    bool
        False when the record is unmapped, below `min_mapq`, above
        `max_mismatches`, or soft-clipped at its 5' end, whose cut site would
        then be off.
    """
    if (
        read.is_unmapped
        or read.mapping_quality < min_mapq
        or read.get_tag("NM") > max_mismatches
    ):
        return False
    op, _ = read.cigartuples[-1 if read.is_reverse else 0]
    return op != pysam.CSOFT_CLIP


def pair_passes(read1, read2, max_insert_size):
    """Return whether two passing mates form a proper pair.

    Parameters
    ----------
    read1, read2 : pysam.AlignedSegment
        Primary mates, each passing `read_passes`.
    max_insert_size : int
        Highest passing |TLEN|.

    Returns
    -------
    bool
        True when the mates lie on opposite strands of one contig, with nonzero
        TLENs that sum to zero and are at most `max_insert_size`.
    """
    return (
        read1.is_reverse != read2.is_reverse
        and read1.reference_id == read2.reference_id
        and read1.template_length != 0
        and read1.template_length + read2.template_length == 0
        and abs(read1.template_length) <= max_insert_size
    )


def flag_group(records, args):
    """Set the QC-fail and proper-pair flags of one read's records in place.

    Parameters
    ----------
    records : list of pysam.AlignedSegment
        Every record of one read name.
    args : argparse.Namespace
        Parsed thresholds.
    """
    primary = [r for r in records if not r.flag & NOT_PRIMARY]
    passes = [read_passes(r, args.min_mapq, args.max_mismatches) for r in primary]
    paired = len(primary) == 2 and all(r.is_paired for r in primary)
    if paired:
        proper = all(passes) and pair_passes(*primary, args.max_insert_size)
        ok = [proper, proper]
    else:
        proper = False
        ok = [p and not r.is_paired for r, p in zip(primary, passes, strict=True)]
    for r in records:
        r.flag &= ~(PROPER | QCFAIL)
        if r.flag & NOT_PRIMARY:
            r.flag |= QCFAIL
    for r, good in zip(primary, ok, strict=True):
        if proper:
            r.flag |= PROPER
        if not good:
            r.flag |= QCFAIL


def main(argv=None):
    """Flag every record of a name-grouped BAM and write it in the same order.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with pysam.AlignmentFile(args.input, "rb", check_sq=False) as bam:
        mode = "wbu" if args.output == "-" else "wb"
        with pysam.AlignmentFile(args.output, mode, template=bam) as out:
            for _, group in itertools.groupby(
                bam.fetch(until_eof=True), lambda r: r.query_name
            ):
                records = list(group)
                flag_group(records, args)
                for r in records:
                    out.write(r)


if __name__ == "__main__":
    main()
