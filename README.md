# dnase-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for DNase-seq. It
ports the [ENCODE DNase-seq pipeline](https://github.com/ENCODE-DCC/dnase-seq-pipeline)
v3.0.0 to maintained tools. Read handling follows v3.0.0. Peaks and footprints
follow ENCODE pipeline ENCPL336PGV (2025), which calls them with hotspot3 and
footprint-tools. Per replicate it runs FASTQs to an unfiltered BAM, a nuclear
BAM, hotspots and peaks, footprints, and a QC table.

## Layout

```
config/
  config.yaml.template   paths and references; copy to config.yaml (gitignored)
  samples.example.tsv    one row per FASTQ (pair): sample_id, r1, r2, adapters, ...
workflow/
  Snakefile              includes the rule files, defines rule all
  rules/
    common.smk           config, sample sheet, envs, constants
    reference.smk        FASTA copy and index, nuclear contigs, mappable bases, bwa-mem3 index
    align.smk            trim and crop, align per layout, merge
    filter.smk           QC flags, duplicate marking, nuclear BAM
    qc.smk               samtools stats, bamcounts, duplicate histogram, preseq, insert sizes
    peaks.smk            hotspot3
    footprints.smk       dispersion model, footprints
    report.smk           run-level QC table and figure
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         rule envs (pixi workspace): default, hotspot3, footprints, preseq, test
.tests/
  download.sh            raw FASTQs, mappable bases, bias model and ENCODE's image
  integration/           chr22 fixture run, step by step against ENCODE's v3.0.0 image
  validation/            four replicates on GRCh38, compared with ENCODE's QC and 2025 calls
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env, lint and test tasks
pyproject.toml           ruff + numpydoc config
```

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/local           # dry run
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
```

bwa-mem3 needs AVX2. The SLURM profile holds Weng Lab cluster settings, such as
partitions and the `cascadelake` constraint. Adapt them for another cluster with
the [SLURM executor docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Tune `set-resources` from the `benchmarks/<run_id>/` TSVs after a first run.
Indexing GRCh38 needs about 70 GB of memory. hotspot3 needs about 80 GB at 6
CPUs.

`tmp_dir` holds per-job scratch and defaults to `.tmp` in the workdir.

## Testing

Both tests use data on the Weng Lab cluster, so they run only there. Download
the inputs once, about 28 GB, then build the fixture:

```bash
bash .tests/download.sh
pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/integration/make_fixture.sh
pixi run test
```

The fixture holds the chr22 reads, and every tenth chrM read, of three ENCODE
replicates. ENCODE's v3 unfiltered BAMs name the reads. They come from the raw
FASTQs.

| Sample | Replicate | Reads |
|---|---|---|
| UHT | ENCSR693UHT rep1 | PE 76, 2 lanes |
| IEA | ENCSR426IEA rep1 | PE 36, 4 lanes, and SE 36, 1 lane |
| RWU | ENCSR420RWU rep1 | SE 36, 2 lanes, preseq in defects mode |

The tests run the workflow on the fixture. Then they run ENCODE's tasks from the
`encodedcc/dnase-seq-pipeline:v3.0.0` image on the same inputs, and compare. Set
`ENCODE_SIF` to use another copy of the image. Trimming, read flags, bamcounts,
the duplicate histogram, the preseq targets, the nuclear read set and footprint
thresholding match ENCODE's exactly. The other steps agree within the bounds in
the tests' docstrings. `test_workflow.py` checks sheet validation and the job
graph.

The validation runs four replicates on full GRCh38 under SLURM: ENCSR693UHT
rep1, both ENCSR426IEA replicates and ENCSR420RWU rep1.
`compare_encode.py` prints their QC beside ENCODE's v3 QC in `encode_qc.tsv`. It
also compares the calls and density with ENCODE's 2025 files.

```bash
snakemake --profile profiles/slurm --configfile .tests/validation/config.yaml
pixi run --manifest-path workflow/envs/pixi.toml -e test \
    python .tests/validation/compare_encode.py results/encode-validation
```

Against ENCODE's 2025 files, the share of ENCODE's calls that overlap ours, and the
Pearson correlation of log1p normalized density in 1 kb bins:

| Replicate | Hotspots, FDR 0.05 | Peaks, FDR 0.001 | Footprints, FDR 0.001 | Density |
|---|---|---|---|---|
| ENCSR693UHT rep1 | 97.3% | 98.0% | 92.8% | 0.996 |
| ENCSR426IEA rep1 | 95.9% | 97.8% | 93.3% | 0.998 |
| ENCSR426IEA rep2 | 94.7% | 97.8% | no ENCODE file | 0.997 |
| ENCSR420RWU rep1 | 92.2% | 94.9% | no ENCODE file | 0.997 |

Call counts lie within 4% of ENCODE's. Against ENCODE's v3 QC, nuclear read
counts lie within 3%, and preseq targets within 2% up to 150M reads. UHT's
median insert size is equal, and its mean within 0.2%. The duplicate fraction of
UHT and RWU is within 0.001 of ENCODE's. IEA's is lower: 0.161 against 0.179 in
rep1, and 0.083 against 0.098 in rep2. On the fixture, samtools and Picard mark
the same duplicates in one BAM, so the gap most likely comes from where
bwa-mem3 and `aln` place reads. The cause is not yet shown.

## Inputs

The sample sheet is a TSV with one row per FASTQ or FASTQ pair. Rows that share
a `sample_id` are lanes of one replicate. One replicate MAY mix single-end and
paired-end lanes. Lines starting with `#` are comments.

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Letters, digits, `_` and `-`. No dots. |
| `r1` | yes | Gzipped FASTQ. |
| `r2` | no | Gzipped mate FASTQ. Empty means single-end. |
| `adapter_r1`, `adapter_r2` | paired-end rows | The library's adapters, trimmed as pairs. |
| `read_length` | no | Crop reads to this length. Empty means no crop. |
| `preseq_defects` | no | `true` runs `preseq lc_extrap -defects`. |
| `footprints` | no | `false` skips the dispersion model and footprints. |

`read_length`, `preseq_defects` and `footprints` MUST agree on every row of a
sample. The adapters MUST agree on its paired-end rows. Take adapters from
ENCODE's trimming reports. The v3.0.0 README gives the P5 sequence as
ENCSR693UHT's R2 adapter, while ENCODE's run trimmed with the one in
`config/samples.example.tsv`.

`ftd learn_dm` fails on a shallow sample, as its piecewise fit has too few
points. Set `footprints` to `false` for such a sample.

## Config

`config/config.yaml.template` lists every key. `samples`, `outdir`, `run_id` and
`tmp_dir` are paths. `references` names:

- `fasta`, the GRCh38 no-alt analysis set FASTA
  (`GRCh38_no_alt_analysis_set_GCA_000001405.15`);
- `nuclear_chroms`, the contigs of the nuclear BAM and of hotspot3. Names
  missing from the FASTA are skipped;
- `mappable_bases`, Altius' K36 mappable bases less the extended blacklist,
  which hotspot3's README names;
- `bias_model`, footprint-tools' `vierstra_et_al.6mer-model.txt`.

## Steps

References, in `reference.smk`: `reference` copies and indexes the FASTA and
writes the nuclear chrom sizes and BED. `mappable_bases` re-compresses the
mappable bases as BGZF and indexes them. `bwa_mem3_index` indexes the FASTA.

Alignment, in `align.smk`: `align` trims each layout of a sample. Paired-end
reads lose adapters with `cutadapt --pair-adapters -e 0.1 -m 1`. All reads are
cropped to `read_length`. bwa-mem3 then aligns them. `merge` joins the layouts
into the unfiltered BAM.

Filtering, in `filter.smk`: `filter` flags reads that fail QC with
`flag_reads.py`, a port of stampipes' `filter_reads.py`. It then marks
duplicates with `samtools markdup --include-fails`. It keeps reads that pass QC
on nuclear contigs. Nothing removes duplicates. hotspot3 and ftd count them, as
hotspot2 did.

QC, in `qc.smk`: `bam_qc` runs `samtools stats` and `flagstat` on both BAMs.
`bamcounts` ports stampipes' counts. `dup_hist` ports the duplicate histogram of
stampipes' `mark_dups.py`, and `preseq` extrapolates it with ENCODE's targets.
`insert_sizes` computes Picard's insert size metrics and stampipes' two derived
values for paired-end-only samples.

Peaks, in `peaks.smk`: `hotspot3` calls hotspots and peaks at FDR 0.05, 0.01 and
0.001, and writes the normalized density.

Footprints, in `footprints.smk`: `dispersion_model` learns ftd's model on the
FDR 0.05 hotspots. `footprints` runs `ftd detect`, then thresholds at FDR 0.01
and 0.001, pads by 3 bp and merges.

Report, in `report.smk`: `qc_summary` gathers every metric into `qc.tsv` and
plots the key ones in `qc.svg`.

## Outputs

Each sample writes to `results/<run_id>/<sample>/<step>/`, with files named
`<sample>...`. References go to `results/refs/`. Logs and benchmarks mirror the
run under `logs/<run_id>/` and `benchmarks/<run_id>/`.

| Step directory | Key outputs |
|---|---|
| `align/` | Unfiltered BAM, in aligner order; cutadapt JSON per layout |
| `filter/` | Nuclear BAM and index; markdup JSON |
| `qc/` | samtools stats and flagstat of both BAMs, bamcounts, preseq curve and targets, insert sizes and plot |
| `hotspot3/` | Cutcounts, fit stats with SPOT, normalized density and background bigWigs; `fdr<F>/` hotspots and peaks as bed.gz and bigBed |
| `footprints/` | Dispersion model; footprints at FDR 0.01 and 0.001 as bed.gz and bigBed |

`results/<run_id>/qc.tsv` holds every metric of the run, one per row:
`sample, section, metric, value`. `qc.svg` plots the key ones.

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule envs live in a second pixi workspace,
`workflow/envs/`. Its `pixi.lock` pins every package.

- `default` runs most rules.
- `hotspot3` holds hotspot3 1.0.1. It pins NumPy 1.26.0, pandas 2.3 and
  scipy 1.13.
- `footprints` builds footprint-tools from git, commit `e839b4e`, with
  genome-tools 1.0.5.
- `preseq` holds preseq 3.2.0. Its builds need an old zlib.
- `test` adds pytest, pytest-xdist and seqkit.

pixi MUST be on `PATH`. An `onstart` hook installs the envs before the first
job.

Add a dependency with `pixi add <pkg>` in `workflow/envs/`. After a hand edit,
run `pixi lock` there (rules use `locked=True`). Any change to the manifest or
lock reruns every rule. So does moving the workdir, because the env hash includes
the workspace's absolute path.

Three dependencies fail the 18-month maintenance test. footprint-tools last
changed in June 2024. No other tool implements its method, and ENCODE still runs
it. preseq 3.2.0 dates from 2022, though its repository still gets commits.
scipy 1.13 is pinned for hotspot3, as noted below.

## Conventions

The workflow follows the
[Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
and [Snakemake's best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html):

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation uses `:q` quoting.
- Config reaches `shell` through `params:` (a `lambda`), never directly.
- All paths are relative.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Shell blocks end with `sync {output:q}`, so the submit host sees the outputs
  on the shared filesystem without waiting.
- Plots use [poikilos](https://github.com/ramirc0/poikilos) with the `plain`
  theme. `pk.save_figure()` writes the SVG and PNG pair. Rules declare both.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and
  `numpydoc lint`. Neither tool parses `.smk` files.

## Differences from ENCODE

Read handling follows v3.0.0, and the output set follows ENCPL336PGV, except
where listed here. Figures are from the chr22 fixture.

- bwa-mem3 `mem -K 100000000` replaces bwa 0.7.12 `aln`, `sampe` and `samse`.
  `-K` makes the output independent of the thread count. Of the reads that pass
  QC in both, 99.8% or more share their 5' end and strand.
- `flag_reads.py` also fails secondary and supplementary records, and primary
  records soft-clipped at their 5' end, whose cut site would be off. `aln` never
  clipped. It does not set the nonstandard 0x1000 flag. The nuclear BAM keeps
  reads on `nuclear_chroms` by position instead.
- `cutadapt --discard-casava` replaces `bwa aln -Y`. cutadapt 5.2 replaces 2.10.
  The trimmed reads are identical.
- `samtools markdup --include-fails -d 100` replaces Picard 2.8.1
  MarkDuplicatesWithMateCigar. Duplicate totals and PERCENT_DUPLICATION match
  Picard's to within 2 reads. Picard keeps the copy with the longest mapped
  length and samtools the one with the best base qualities. Sets mix QC-fail
  reads, so the nuclear BAM's duplicate count is about 1% lower.
- `samtools fixmate` recomputes TLEN after filtering, so insert sizes can move
  by a few bp.
- Insert sizes come from the |TLEN| of first mates, with Picard's formulas.
  Median and MAD match. Mean and SD agree within 0.01%. There is no PDF.
- preseq 3.2.0 replaces 2.0.3. Its curve lies within 0.2% of 2.0.3's up to 10
  times the observed depth, and the targets match.
- The preseq targets skip `nan` rows, which `-defects` can write past some depth.
  stampipes' awk took them as larger than every target. In defects mode, the
  250M target and the maximum still move with the curve's tail.
- hotspot3 1.0.1 replaces hotspot1 and hotspot2, as in ENCPL336PGV. There is no
  SPOT1 or SPOT2. `qc.tsv` reports hotspot3's SPOT and its SD, which do not
  compare with either. There are no narrowPeak, allcalls or starch files.
  Of ENCODE's 2025 FDR 0.001 peaks on chr22, 96–98% overlap ours.
- hotspot3 runs its cutcount script through `/bin/sh` with `set -o pipefail`,
  which dash rejects. The rule runs that script with bash and passes the
  cutcounts in.
- One mappable bases file serves all read lengths. v3.0.0 used separate K36 and
  K76 files. Reads are still cropped to `read_length`.
- footprint-tools runs from master, commit `e839b4e`, past 1.3.7. The
  thresholds skip the bedgraph's `#` lines, which ENCODE's awk passes as a
  bogus interval, and clamp padded starts at 0. ENCODE's 2025 footprints use a
  100M-read subsample. These use every read. Of ENCODE's chr22
  footprints, 93–95% overlap ours.
- `ftd detect` on more than one thread does not reproduce exactly, even with
  `--seed`. Its workers share one random generator and one batch queue.
- The `bamcounts` `mm` count is 0, because `mem` writes no `XT:A:R` tag.
- hotspot3 1.0.1 needs pandas below 3, whose copy-on-write rejects its in-place
  updates. It needs scipy below 1.14, whose Boost `betainc` returns NaN instead
  of 0 on float32 underflow.
- Not ported: the QC-fail-filtered BAM with every contig, the flagged and
  marked BAM, UMI parsing of read names, and genomes other than GRCh38.

## References

Daley, Timothy, and Andrew D. Smith. 2013. “Predicting the Molecular Complexity of Sequencing Libraries.” Nature Methods 10 (4): 325–27. https://doi.org/10.1038/nmeth.2375.

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2): giab008. https://doi.org/10.1093/gigascience/giab008.

John, Sam, Peter J. Sabo, Robert E. Thurman, et al. 2011. “Chromatin Accessibility Pre-Determines Glucocorticoid Receptor Binding Patterns.” Nature Genetics 43 (3): 264–68. https://doi.org/10.1038/ng.759.

Li, Heng. 2013. “Aligning Sequence Reads, Clone Sequences and Assembly Contigs with BWA-MEM.” arXiv:1303.3997. https://arxiv.org/abs/1303.3997.

Martin, Marcel. 2011. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” EMBnet.Journal 17 (1): 10–12. https://doi.org/10.14806/ej.17.1.200.

Meuleman, Wouter, Alexander Muratov, Eric Rynes, et al. 2020. “Index and Biological Spectrum of Human DNase I Hypersensitive Sites.” Nature 584 (7820): 244–51. https://doi.org/10.1038/s41586-020-2559-3.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10: 33. https://doi.org/10.12688/f1000research.29032.3.

Vierstra, Jeff, John Lazar, Richard Sandstrom, et al. 2020. “Global Reference Mapping of Human Transcription Factor Footprints.” Nature 583 (7818): 729–36. https://doi.org/10.1038/s41586-020-2528-x.

Software without a paper: [bwa-mem3](https://github.com/fg-labs/bwa-mem3),
[hotspot3](https://github.com/vierstralab/hotspot3),
[footprint-tools](https://github.com/vierstralab/footprint-tools) and
[stampipes](https://github.com/StamLab/stampipes).
