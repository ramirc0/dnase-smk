"""Helpers for the integration tests: the workflow run, single rules, and ENCODE's tasks.

pytest-xdist runs tests in parallel, so the workflow run and each ENCODE task hold
a file lock: the first worker does the work, the others wait and reuse it.
"""

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import polars as pl
import pytest
from filelock import FileLock

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ".tests/integration/config.yaml"
RUN = ROOT / "results-fixture/fixture"
REFS = ROOT / "results-fixture/refs"
FIXTURE = ROOT / "resources/fixture"
ENCODE_RUNS = ROOT / "results-fixture/encode"
SIF = Path(
    os.environ.get(
        "ENCODE_SIF",
        "/data/zusers/ramirezc/dnase-smk-refs/containers/dnase-seq-pipeline_v3.0.0.sif",
    )
)
GENOME = FIXTURE / "genome.fa"
SHEET = pl.read_csv(
    ROOT / ".tests/integration/samples.tsv",
    separator="\t",
    infer_schema_length=0,
    comment_prefix="#",
).fill_null("")
LAYOUTS = sorted(
    {(r["sample_id"], "pe" if r["r2"] else "se") for r in SHEET.iter_rows(named=True)}
)


def lanes(sample, layout):
    """Return a sample's sheet rows of one layout.

    Parameters
    ----------
    sample : str
        Sample ID from the fixture sheet.
    layout : {"pe", "se"}
        Layout.

    Returns
    -------
    list of dict
        Sheet rows in sheet order, with FASTQ paths made absolute.
    """
    rows = SHEET.filter(
        (pl.col("sample_id") == sample) & ((pl.col("r2") != "") == (layout == "pe"))
    ).iter_rows(named=True)
    return [
        {**r, "r1": ROOT / r["r1"], "r2": ROOT / r["r2"] if r["r2"] else None}
        for r in rows
    ]


def snakemake(*args, cwd=ROOT):
    """Run the launcher's Snakemake.

    Parameters
    ----------
    *args : str or pathlib.Path
        Snakemake arguments.
    cwd : pathlib.Path, default ROOT
        Working directory.

    Returns
    -------
    subprocess.CompletedProcess
        The finished process, with stdout and stderr captured as text.
    """
    return subprocess.run(
        ["pixi", "run", "--manifest-path", str(ROOT / "pixi.toml"), "snakemake"]
        + [str(a) for a in args],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.fixture(scope="session")
def run():
    """Run the workflow on the fixture, or skip when the fixture is not built.

    Returns
    -------
    pathlib.Path
        Run directory, ``results-fixture/fixture``.
    """
    if not GENOME.exists():
        pytest.skip("fixture missing: run .tests/integration/make_fixture.sh")
    RUN.parent.mkdir(parents=True, exist_ok=True)
    with FileLock(RUN.parent / "run.lock"):
        # Some tests read temp files, such as the per-layout BAMs.
        done = snakemake(
            "--profile", "profiles/local", "--configfile", CONFIG, "--notemp"
        )
    assert done.returncode == 0, done.stderr[-5000:]
    return RUN


def encode(task, inputs, command):
    """Run commands in ENCODE's v3.0.0 image, unless its inputs are older.

    Parameters
    ----------
    task : str
        Output directory under ``results-fixture/encode``.
    inputs : list of pathlib.Path
        Files the task reads. Each is linked into the output directory under
        its own name, so ENCODE never writes next to the workflow's outputs.
    command : str
        Shell command run in the output directory inside the container.

    Returns
    -------
    pathlib.Path
        The output directory.
    """
    out = ENCODE_RUNS / task
    stamp = out / ".done"
    out.parent.mkdir(parents=True, exist_ok=True)
    with FileLock(f"{out}.lock"):
        if stamp.exists() and all(
            p.stat().st_mtime < stamp.stat().st_mtime for p in inputs
        ):
            return out
        if not SIF.exists():
            pytest.skip(f"ENCODE image missing: {SIF}")
        shutil.rmtree(out, ignore_errors=True)
        out.mkdir(parents=True)
        for p in inputs:
            (out / p.name).symlink_to(p)
        subprocess.run(
            [
                "apptainer",
                "exec",
                "--bind",
                "/data",
                str(SIF),
                "bash",
                "-euo",
                "pipefail",
                "-c",
                command,
            ],
            cwd=out,
            check=True,
        )
        stamp.touch()
    return out


def lines(path):
    """Read a text file, gzipped or not, as a list of lines.

    Parameters
    ----------
    path : pathlib.Path
        File to read.

    Returns
    -------
    list of str
        Lines without their newlines.
    """
    with path.open("rb") as f:
        gzipped = f.read(2) == b"\x1f\x8b"
    opener = gzip.open if gzipped else open
    with opener(path, "rt") as f:
        return f.read().splitlines()


def first_difference(ours, theirs):
    """Describe where two long line lists first differ.

    pytest's own diff of million-line lists takes minutes, so tests compare
    with this instead.

    Parameters
    ----------
    ours, theirs : list of str
        Lines to compare.

    Returns
    -------
    str or None
        None when the lists are equal, else the first differing line of each
        and both lengths.
    """
    if ours == theirs:
        return None
    i = next(
        (i for i, (a, b) in enumerate(zip(ours, theirs, strict=False)) if a != b),
        min(len(ours), len(theirs)),
    )
    a = ours[i] if i < len(ours) else "<end>"
    b = theirs[i] if i < len(theirs) else "<end>"
    return f"line {i}: ours {a!r}, theirs {b!r}; lengths {len(ours)} vs {len(theirs)}"


def run_rule(workdir, rule, target, files, config=None):
    """Run one rule alone in a scratch working directory.

    Every other rule is disallowed, so the files given stand in for upstream
    outputs. The scratch directory links to the repository's fixture and
    `.tests`, so the fixture config resolves unchanged.

    Parameters
    ----------
    workdir : pathlib.Path
        Scratch directory, e.g. pytest's ``tmp_path``.
    rule : str or list of str
        Rule or rules to run.
    target : str
        Output to request, relative to ``workdir``.
    files : dict of str to bytes or pathlib.Path
        Inputs to place, by path relative to ``workdir``: bytes are written,
        paths are copied.
    config : str, optional
        YAML merged over the fixture config.

    Returns
    -------
    subprocess.CompletedProcess
        The finished Snakemake process.
    """
    for name in ("resources", ".tests"):
        (workdir / name).symlink_to(ROOT / name)
    for name, content in files.items():
        dest = workdir / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(content, Path):
            shutil.copy(content, dest)
        else:
            dest.write_bytes(content)
    rules = [rule] if isinstance(rule, str) else rule
    configs = [CONFIG]
    if config:
        (workdir / "extra.yaml").write_text(config)
        configs.append("extra.yaml")
    return snakemake(
        "--snakefile",
        ROOT / "workflow/Snakefile",
        "--directory",
        workdir,
        "--profile",
        ROOT / "profiles/local",
        "--configfile",
        *configs,
        "--allowed-rules",
        *rules,
        "--forcerun",
        *rules,
        "--retries",
        "0",
        target,
        cwd=workdir,
    )
