"""Compare each step with ENCODE's DNase-seq pipeline v3.0.0, run in its image."""

import subprocess

import pysam
import pytest
from conftest import GENOME, LAYOUTS, encode, lanes

# ENCODE's crop: cut every sequence and quality line to the read length.
CROP = "awk 'NR % 2 == 0 {{ print substr($0, 1, {n}) }} NR % 2 != 0'"


@pytest.fixture(scope="session")
def encode_index():
    """Build ENCODE's bwa 0.7.12 index of the fixture genome.

    Returns
    -------
    list of pathlib.Path
        The FASTA link and its index files.
    """
    out = encode("bwa_index", [GENOME], f"bwa index {GENOME.name}")
    return [out / GENOME.name] + [
        out / f"{GENOME.name}.{ext}" for ext in ("amb", "ann", "bwt", "pac", "sa")
    ]


def encode_trim(sample, layout):
    """Trim a sample's lanes of one layout as ENCODE does.

    Paired-end lanes go through cutadapt 2.10 with ENCODE's settings. Then every
    read is cropped to the read length.

    Parameters
    ----------
    sample : str
        Fixture sample.
    layout : {"pe", "se"}
        Layout.

    Returns
    -------
    list of pathlib.Path
        Cropped FASTQs, R1 then R2 for paired end.
    """
    rows = lanes(sample, layout)
    crop = CROP.format(n=rows[0]["read_length"])
    r1 = " ".join(r["r1"].name for r in rows)
    command = f"cat {r1} > R1.fastq.gz\n"
    if layout == "pe":
        r2 = " ".join(r["r2"].name for r in rows)
        command += (
            f"cat {r2} > R2.fastq.gz\n"
            f"cutadapt -a {rows[0]['adapter_r1']} -A {rows[0]['adapter_r2']} --cores=8 "
            "--pair-adapters --minimum-length 1 --error-rate 0.1 "
            "--output trim.R1.fastq.gz --paired-output trim.R2.fastq.gz "
            "R1.fastq.gz R2.fastq.gz > trimstats.txt\n"
        )
        mates = {"R1": "trim.R1", "R2": "trim.R2"}
    else:
        mates = {"R1": "R1"}
    for mate, source in mates.items():
        command += (
            f"pigz -dc {source}.fastq.gz | {crop} | pigz -n > crop.{mate}.fastq.gz\n"
        )
    inputs = [r["r1"] for r in rows] + [r["r2"] for r in rows if r["r2"]]
    out = encode(f"trim/{sample}.{layout}", inputs, command)
    return [out / f"crop.{mate}.fastq.gz" for mate in mates]


def fastq_reads(path):
    """Read a FASTQ, gzipped or not, as a set of reads.

    Parameters
    ----------
    path : pathlib.Path
        FASTQ file.

    Returns
    -------
    set of tuple of str
        One (name, sequence, quality) per read, the name cut at its first space.
    """
    with pysam.FastxFile(str(path)) as f:
        return {(r.name, r.sequence, r.quality) for r in f}


def bam_fastqs(bam, directory):
    """Write the primary reads of a BAM as FASTQs, in sequencing orientation.

    Parameters
    ----------
    bam : pathlib.Path
        Aligned BAM.
    directory : pathlib.Path
        Where to write the FASTQs.

    Returns
    -------
    dict of str to pathlib.Path
        FASTQ by layout and mate: ``pe.R1`` and ``pe.R2`` hold paired reads,
        ``se.R1`` the single-end ones.
    """
    paths = {k: directory / f"{k}.fastq" for k in ("pe.R1", "pe.R2", "se.R1")}
    subprocess.run(
        ["samtools", "fastq", "-n", "-F", "0x900"]
        + [
            "-1",
            str(paths["pe.R1"]),
            "-2",
            str(paths["pe.R2"]),
            "-0",
            str(paths["se.R1"]),
        ]
        + ["-s", "/dev/null", str(bam)],
        check=True,
        capture_output=True,
    )
    return paths


def five_prime_ends(bam):
    """Return the 5' end of each read that passes filter_reads.py's record rules.

    A record passes when it is primary and mapped, with MAPQ >= 10 and NM <= 2.

    Parameters
    ----------
    bam : pathlib.Path
        Aligned BAM.

    Returns
    -------
    dict of tuple to tuple
        (contig, 5' position, reverse) by (name, is_read2).
    """
    ends = {}
    with pysam.AlignmentFile(str(bam), check_sq=False) as f:
        for r in f:
            if r.flag & 0x904 or r.mapping_quality < 10 or r.get_tag("NM") > 2:
                continue
            end = r.reference_end if r.is_reverse else r.reference_start
            ends[r.query_name, r.is_read2] = (r.reference_name, end, r.is_reverse)
    return ends


@pytest.mark.parametrize(("sample", "layout"), LAYOUTS)
def test_trimming(run, tmp_path, sample, layout):
    """Our trimming keeps >= 99.9% of reads identical to ENCODE's cutadapt 2.10 and awk."""
    ours = bam_fastqs(run / f"{sample}/align/{sample}.unfiltered.bam", tmp_path)
    for theirs in encode_trim(sample, layout):
        mate = theirs.name.split(".")[1]
        ours_reads = fastq_reads(ours[f"{layout}.{mate}"])
        theirs_reads = fastq_reads(theirs)
        shared = len(ours_reads & theirs_reads)
        print(f"{sample}.{layout} {mate}: {shared} of {len(theirs_reads)} identical")
        assert shared >= 0.999 * max(len(ours_reads), len(theirs_reads))


@pytest.mark.parametrize(("sample", "layout"), LAYOUTS)
def test_alignment(run, encode_index, sample, layout):
    """>= 95% of reads that pass in both share their 5' end and strand with bwa aln."""
    crops = encode_trim(sample, layout)
    names = " ".join(c.name for c in crops)
    aln = "bwa aln -Y -n 0.04 -l 32 -t 16 {fa} {fq} > {fq}.sai\n"
    command = "".join(aln.format(fa=GENOME.name, fq=c.name) for c in crops)
    pair = "sampe -a 750" if layout == "pe" else "samse"
    sais = " ".join(f"{c.name}.sai" for c in crops)
    command += f"bwa {pair} -n 10 {GENOME.name} {sais} {names} | samtools view -b -o aln.bam -\n"
    out = encode(f"align/{sample}.{layout}", encode_index + crops, command)
    theirs = five_prime_ends(out / "aln.bam")
    ours = five_prime_ends(run / f"{sample}/align/{sample}.unfiltered.bam")
    both = ours.keys() & theirs.keys()
    same = sum(ours[k] == theirs[k] for k in both)
    print(
        f"{sample}.{layout}: {len(both)} pass in both, {same} share their 5' end; "
        f"{len(ours)} pass ours, {len(theirs)} pass ENCODE's"
    )
    assert same >= 0.95 * len(both)
