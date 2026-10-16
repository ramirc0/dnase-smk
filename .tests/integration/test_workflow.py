"""Check sheet validation and the shape of the job graph."""

import re

import pytest
from conftest import CONFIG, ROOT, snakemake

HEADER = "sample_id\tr1\tr2\tadapter_r1\tadapter_r2\tread_length\tpreseq_defects\tfootprints\n"
PE = "S\tr1.fastq.gz\tr2.fastq.gz\tACGT\tTTGA\t36\t\t\n"
SE = "S\tse.fastq.gz\t\t\t\t36\t\t\n"


def dry_run(tmp_path, sheet, *args):
    """Dry-run the workflow on a sheet in a scratch directory.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Scratch directory.
    sheet : str
        Sample sheet text.
    *args : str
        Extra Snakemake arguments.

    Returns
    -------
    subprocess.CompletedProcess
        The finished Snakemake process.
    """
    (tmp_path / "samples.tsv").write_text(sheet)
    (tmp_path / "extra.yaml").write_text("samples: samples.tsv\n")
    (tmp_path / "resources").symlink_to(ROOT / "resources")
    (tmp_path / ".tests").symlink_to(ROOT / ".tests")
    for fastq in re.findall(r"\S+\.fastq\.gz", sheet):
        (tmp_path / fastq).touch()
    return snakemake(
        "--snakefile",
        ROOT / "workflow/Snakefile",
        "--directory",
        tmp_path,
        "--configfile",
        CONFIG,
        "extra.yaml",
        "-n",
        *args,
        cwd=tmp_path,
    )


@pytest.mark.parametrize(
    ("sheet", "message"),
    [
        ("sample_id\tr2\n", r"missing columns \['r1'\]"),
        (HEADER, "no rows"),
        (HEADER + PE.replace("S\t", "S.1\t"), "must match"),
        (HEADER + PE.replace("r1.fastq.gz", ""), "empty r1"),
        (HEADER + PE.replace("r1.fastq.gz", "r1.fastq"), "must be gzipped"),
        (
            HEADER + PE.replace("\t36\t", "\tlong\t"),
            "read_length must be a positive integer",
        ),
        (
            HEADER + PE.replace("\t\t\n", "\tyes\t\n"),
            "preseq_defects must be true or false",
        ),
        (HEADER + PE.replace("\t\t\n", "\t\tno\n"), "footprints must be true or false"),
        (HEADER + PE.replace("ACGT", ""), "need adapter_r1 and adapter_r2"),
        (HEADER + PE.replace("ACGT", "acgt"), "need adapter_r1 and adapter_r2"),
        (HEADER + PE + PE.replace("TTGA", "TTGC"), "disagree on adapters"),
        (HEADER + PE + SE.replace("\t36\t", "\t76\t"), "disagree on read_length"),
    ],
)
def test_bad_sheet(tmp_path, sheet, message):
    """A malformed sheet fails the dry run with a message naming the problem."""
    done = dry_run(tmp_path, sheet)
    assert done.returncode != 0
    assert re.search(message, done.stdout + done.stderr), done.stderr[-2000:]


def job_counts(done):
    """Read the job count of each rule from a dry run's job stats.

    Parameters
    ----------
    done : subprocess.CompletedProcess
        A finished dry run.

    Returns
    -------
    dict of str to int
        Jobs per rule.
    """
    stats = (done.stdout + done.stderr).split("Job stats:")[1]
    return {
        m[1]: int(m[2]) for m in re.finditer(r"^(\w+)\s+(\d+)$", stats, re.MULTILINE)
    }


def test_mixed_layouts(tmp_path):
    """A mixed sample aligns each layout once and skips insert sizes; PE-only gets them."""
    sheet = HEADER + PE + SE + PE.replace("S\t", "P\t") + SE.replace("S\t", "U\t")
    done = dry_run(tmp_path, sheet)
    assert done.returncode == 0, done.stderr[-2000:]
    jobs = job_counts(done)
    assert jobs["align"] == 4
    assert jobs["merge"] == 3
    assert jobs["insert_sizes"] == 1
    assert jobs["bam_qc"] == 6
    assert jobs["footprints"] == 3


def test_footprints_off(tmp_path):
    """A footprints value of false drops the sample's dispersion model and footprints."""
    sheet = HEADER + SE.replace("\t\t\n", "\t\tfalse\n") + SE.replace("S\t", "T\t")
    done = dry_run(tmp_path, sheet)
    assert done.returncode == 0, done.stderr[-2000:]
    jobs = job_counts(done)
    assert jobs["dispersion_model"] == 1
    assert jobs["footprints"] == 1
