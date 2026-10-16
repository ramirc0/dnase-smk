# Config, sample sheet, and shared helpers.

import re
from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run
REFDIR = f"{RESULTS}/refs"              # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
TMP_DIR = config["tmp_dir"]

REFERENCES = config["references"]
BWA_INDEX = f"{REFDIR}/bwa_mem3/reference"

# ENCODE pipeline ENCPL336PGV's hotspot3 and footprint thresholds.
HOTSPOT_FDRS = ("0.05", "0.01", "0.001")
FOOTPRINT_FDRS = ("0.01", "0.001")

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)
HOTSPOT3_ENV = pixi(workspace="../envs", env="hotspot3", locked=True)
FOOTPRINTS_ENV = pixi(workspace="../envs", env="footprints", locked=True)
PRESEQ_ENV = pixi(workspace="../envs", env="preseq", locked=True)


# The pixi plugin installs the env inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install it once before any job starts.
onstart:
    shell("pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked")


# samtools sort -T and the per-job scratch directories live here.
Path(TMP_DIR).mkdir(parents=True, exist_ok=True)

# Dots separate name parts in file names, so names never contain one.
NAME = r"[A-Za-z0-9][A-Za-z0-9_-]*"

_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: sample_id, r1, [r2], [adapter_r1], [adapter_r2], [read_length], "
        "[preseq_defects], [footprints]; see config/samples.example.tsv."
    )


def _fail(message):
    """Raise a sample sheet error.

    Parameters
    ----------
    message : str
        What is wrong.

    Raises
    ------
    WorkflowError
        Always, prefixed with the sheet path.
    """
    raise WorkflowError(f"{_SHEET}: {message}")


# All-string read; optional columns default to "".
_OPTIONAL = ("r2", "adapter_r1", "adapter_r2", "read_length", "preseq_defects", "footprints")
_sheet = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")
_missing = {"sample_id", "r1"} - set(_sheet.columns)
if _missing:
    _fail(f"missing columns {sorted(_missing)}")
_sheet = _sheet.with_columns(
    pl.lit("").alias(c) for c in _OPTIONAL if c not in _sheet.columns
).with_columns(pl.all().fill_null("").str.strip_chars())
if _sheet.is_empty():
    _fail("no rows")

SAMPLES = {}  # sample -> {"read_length", "defects", "footprints", "adapters"}
LANES = {}  # (sample, layout) -> {"r1", "r2"}
for row in _sheet.iter_rows(named=True):
    sample = row["sample_id"]
    if not re.fullmatch(NAME, sample):
        _fail(f"sample_id {sample!r} must match {NAME}")
    if not row["r1"]:
        _fail(f"{sample}: empty r1")
    for fastq in (row["r1"], row["r2"]):
        if fastq and not fastq.endswith(".gz"):
            _fail(f"{sample}: {fastq} must be gzipped (.gz)")
    if not re.fullmatch(r"([1-9][0-9]*)?", row["read_length"]):
        _fail(f"{sample}: read_length must be a positive integer")
    for column in ("preseq_defects", "footprints"):
        if row[column] not in ("", "true", "false"):
            _fail(f"{sample}: {column} must be true or false")
    layout = "pe" if row["r2"] else "se"
    attributes = {
        "read_length": int(row["read_length"] or 0),
        "defects": row["preseq_defects"] == "true",
        "footprints": row["footprints"] != "false",
    }
    s = SAMPLES.setdefault(sample, {**attributes, "adapters": None})
    if {k: s[k] for k in attributes} != attributes:
        _fail(f"{sample}: rows disagree on read_length, preseq_defects or footprints")
    if layout == "pe":
        adapters = (row["adapter_r1"], row["adapter_r2"])
        if not all(re.fullmatch(r"[ACGTN]+", a) for a in adapters):
            _fail(f"{sample}: paired-end rows need adapter_r1 and adapter_r2 (ACGTN)")
        if s["adapters"] not in (None, adapters):
            _fail(f"{sample}: paired-end rows disagree on adapters")
        s["adapters"] = adapters
    lanes = LANES.setdefault((sample, layout), {"r1": [], "r2": []})
    lanes["r1"].append(row["r1"])
    if row["r2"]:
        lanes["r2"].append(row["r2"])

LAYOUTS = {s: [layout for layout in ("pe", "se") if (s, layout) in LANES] for s in SAMPLES}
PAIRED_ONLY = [s for s in SAMPLES if LAYOUTS[s] == ["pe"]]


wildcard_constraints:
    sample=NAME,
    layout=r"pe|se",
    bam=r"unfiltered|nuclear",


def prefix(sample, step):
    """Return the output prefix of a sample's step.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    step : str
        Step directory, such as `align` or `hotspot3`.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<sample>/<step>/<sample>`.
    """
    return f"{OUTDIR}/{sample}/{step}/{sample}"
