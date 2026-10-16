"""Compare each step with ENCODE's DNase-seq pipeline v3.0.0, run in its image."""

import itertools
import json
import subprocess

import pysam
import pytest
from conftest import (
    GENOME,
    LAYOUTS,
    REFS,
    ROOT,
    RUN,
    SAMPLES,
    encode,
    first_difference,
    lanes,
)
from filelock import FileLock

FLAG_READS = ROOT / "workflow/scripts/flag_reads.py"
OURS = RUN.parent / "tests"

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


def flagged(sample):
    """Flag a sample's unfiltered BAM with flag_reads.py and sort it by coordinate.

    The BAM is rebuilt when it is older than the unfiltered BAM.

    Parameters
    ----------
    sample : str
        Fixture sample.

    Returns
    -------
    pathlib.Path
        Flagged, coordinate-sorted BAM, before duplicate marking.
    """
    unfiltered = RUN / f"{sample}/align/{sample}.unfiltered.bam"
    out = OURS / f"flagged/{sample}.bam"
    out.parent.mkdir(parents=True, exist_ok=True)
    with FileLock(f"{out}.lock"):
        if not out.exists() or out.stat().st_mtime < unfiltered.stat().st_mtime:
            subprocess.run(
                f"python {FLAG_READS} {unfiltered} - | samtools sort -@ 4 -o {out}.part - "
                f"&& mv {out}.part {out}",
                shell=True,
                check=True,
            )
    return out


def clean_groups(sample):
    """Write the name groups of an unfiltered BAM that stampipes can judge.

    stampipes pairs adjacent records, so a group passes when it holds only
    primary records and none is soft-clipped at its 5' end.

    Parameters
    ----------
    sample : str
        Fixture sample.

    Returns
    -------
    pathlib.Path
        BAM of the passing groups, in aligner order.
    """
    unfiltered = RUN / f"{sample}/align/{sample}.unfiltered.bam"
    out = OURS / f"clean/{sample}.bam"
    out.parent.mkdir(parents=True, exist_ok=True)
    with FileLock(f"{out}.lock"):
        if out.exists() and out.stat().st_mtime > unfiltered.stat().st_mtime:
            return out
        with (
            pysam.AlignmentFile(str(unfiltered), check_sq=False) as bam,
            pysam.AlignmentFile(f"{out}.part", "wb", template=bam) as dest,
        ):
            for _, group in itertools.groupby(bam, lambda r: r.query_name):
                records = list(group)
                if all(
                    not r.flag & 0x900
                    and (
                        r.is_unmapped
                        or r.cigartuples[-1 if r.is_reverse else 0][0]
                        != pysam.CSOFT_CLIP
                    )
                    for r in records
                ):
                    for r in records:
                        dest.write(r)
        (OURS / f"clean/{sample}.bam.part").rename(out)
    return out


def nuclear_list():
    """Write the nuclear contig list that stampipes' filter_reads.py reads.

    Returns
    -------
    pathlib.Path
        One contig name per line, from the workflow's nuclear chrom sizes.
    """
    out = OURS / "nuclear_chroms.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    names = [
        line.split("\t")[0]
        for line in (REFS / "nuclear.chrom.sizes").read_text().splitlines()
    ]
    if not out.exists() or out.read_text().split() != names:
        out.write_text("".join(f"{n}\n" for n in names))
    return out


def flags(bam, drop=0):
    """Read the name and flag of every record of a BAM, in file order.

    Parameters
    ----------
    bam : pathlib.Path
        BAM to read.
    drop : int, default 0
        Flag bits to clear before comparing.

    Returns
    -------
    list of str
        ``name<TAB>flag`` per record.
    """
    with pysam.AlignmentFile(str(bam), check_sq=False) as f:
        return [f"{r.query_name}\t{r.flag & ~drop}" for r in f]


def picard_metrics(path):
    """Read the metrics row of a Picard metrics file.

    Parameters
    ----------
    path : pathlib.Path
        Picard METRICS_FILE.

    Returns
    -------
    dict of str to str
        Value by column, from the row under the ``## METRICS CLASS`` line.
    """
    rows = path.read_text().splitlines()
    i = next(i for i, line in enumerate(rows) if line.startswith("## METRICS CLASS"))
    return dict(zip(rows[i + 1].split("\t"), rows[i + 2].split("\t"), strict=False))


def count(*args):
    """Count records with ``samtools view -c``.

    Parameters
    ----------
    *args : str or pathlib.Path
        samtools view arguments, ending with the BAM.

    Returns
    -------
    int
        Number of records selected.
    """
    done = subprocess.run(
        ["samtools", "view", "-c", *map(str, args)],
        check=True,
        capture_output=True,
        text=True,
    )
    return int(done.stdout)


@pytest.mark.parametrize("sample", SAMPLES)
def test_flag_reads(run, tmp_path, sample):
    """Our flag_reads.py sets the same flags as stampipes' filter_reads.py, less 0x1000."""
    clean = clean_groups(sample)
    contigs = nuclear_list()
    theirs = encode(
        f"flag/{sample}",
        [clean, contigs],
        f"python3 $(which filter_reads.py) {clean.name} flagged.bam {contigs.name}",
    )
    ours = tmp_path / "flagged.bam"
    subprocess.run(["python", str(FLAG_READS), str(clean), str(ours)], check=True)
    assert (
        first_difference(flags(ours), flags(theirs / "flagged.bam", drop=0x1000))
        is None
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_markdup(run, sample):
    """Our duplicate marking matches Picard's duplicates within 0.1% and its PERCENT_DUPLICATION within 0.001."""
    bam = flagged(sample)
    regex = "'[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'"
    picard = "java -Xmx16G -XX:ParallelGCThreads=1 -jar $(which picard.jar)"
    out = encode(
        f"markdup/{sample}",
        [bam],
        f"{picard} RevertOriginalBaseQualitiesAndAddMateCigar INPUT={bam.name} OUTPUT=mc.bam "
        "RESTORE_ORIGINAL_QUALITIES=false VALIDATION_STRINGENCY=SILENT MAX_RECORDS_TO_EXAMINE=0\n"
        f"{picard} MarkDuplicatesWithMateCigar INPUT=mc.bam OUTPUT=marked.bam METRICS_FILE=metrics.txt "
        f"ASSUME_SORTED=true MINIMUM_DISTANCE=300 READ_NAME_REGEX={regex} VALIDATION_STRINGENCY=SILENT\n",
    )
    stats = json.loads((run / f"{sample}/filter/{sample}.markdup.json").read_bytes())
    ours = stats["DUPLICATE TOTAL"]
    theirs = count("-f", "0x400", out / "marked.bam")
    percent = (stats["DUPLICATE PAIR"] + stats["DUPLICATE SINGLE"]) / (
        stats["PAIRED"] + stats["SINGLE"]
    )
    picard_percent = float(picard_metrics(out / "metrics.txt")["PERCENT_DUPLICATION"])
    # Picard keeps the copy with the longest mapped length, samtools the one with the best
    # base qualities. Sets mix QC-fail reads, so the nuclear BAM's duplicate count moves more.
    nuclear = count("-f", "0x400", run / f"{sample}/filter/{sample}.nuclear.bam")
    picard_nuclear = count(
        "-f", "0x400", "-F", "0x200", "-L", REFS / "nuclear.bed", out / "marked.bam"
    )
    print(
        f"{sample}: duplicates {ours} vs {theirs}, nuclear {nuclear} vs {picard_nuclear}; "
        f"PERCENT_DUPLICATION {percent:.5f} vs {picard_percent:.5f}"
    )
    assert abs(ours - theirs) <= 0.001 * theirs
    assert abs(percent - picard_percent) <= 0.001


@pytest.mark.parametrize("sample", SAMPLES)
def test_nuclear(run, sample):
    """The nuclear BAM holds the reads ENCODE's -F 512 and -F 4096 filters keep."""
    bam = flagged(sample)
    nuclear = set((REFS / "nuclear.chrom.sizes").read_text().split()[::2])
    marked = OURS / f"nonnuclear/{sample}.bam"
    marked.parent.mkdir(parents=True, exist_ok=True)
    with (
        pysam.AlignmentFile(str(bam)) as f,
        pysam.AlignmentFile(str(marked), "wb", template=f) as dest,
    ):
        for r in f:
            if not r.is_unmapped and r.reference_name not in nuclear:
                r.flag |= 0x1000
            dest.write(r)
    out = encode(
        f"nuclear/{sample}",
        [marked],
        f"samtools view -F 512 -b {marked.name} > filtered.bam\n"
        "samtools view -F 4096 -b filtered.bam > nuclear.bam\n",
    )

    def reads(path):
        with pysam.AlignmentFile(str(path)) as f:
            return sorted(
                f"{r.query_name}\t{r.flag & 0xC1}\t{r.reference_name}\t{r.reference_start}"
                for r in f
            )

    assert (
        first_difference(
            reads(run / f"{sample}/filter/{sample}.nuclear.bam"),
            reads(out / "nuclear.bam"),
        )
        is None
    )
