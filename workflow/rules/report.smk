# Run-level QC table and figure.


def qc_inputs():
    """Return every QC file of the run, by `qc_summary.py` flag.

    Returns
    -------
    dict of str to list of str
        Paths per flag, keyed with underscores. File names start with the
        sample.
    """
    files = {}
    for sample in SAMPLES:
        qc = prefix(sample, "qc")
        hotspot3 = prefix(sample, "hotspot3")
        add = {
            "cutadapt": [f"{prefix(sample, 'align')}.{layout}.cutadapt.json" for layout in LAYOUTS[sample]],
            "flagstat": [f"{qc}.{bam}.flagstat.json" for bam in BAMS],
            "stats": [f"{qc}.{bam}.stats.txt" for bam in BAMS],
            "bamcounts": [f"{qc}.bamcounts.txt"],
            "markdup": [f"{prefix(sample, 'filter')}.markdup.json"],
            "preseq": [f"{qc}.preseq_targets.txt"],
            "insert_sizes": [f"{qc}.insert_sizes.tsv"] * (sample in PAIRED_ONLY),
            "fit_stats": [f"{hotspot3}.fit_stats.bed.gz"],
            "total_cutcounts": [f"{hotspot3}.total_cutcounts"],
            "calls": [f for f in hotspot3_calls(sample) if f.endswith(".bed.gz")],
            "footprints": [],
            "dispersion_model": [],
        }
        if SAMPLES[sample]["footprints"]:
            add["footprints"] = [f for f in footprint_files(sample) if f.endswith(".bed.gz")]
            add["dispersion_model"] = [f"{prefix(sample, 'footprints')}.dm.json"]
        for flag, paths in add.items():
            files.setdefault(flag, []).extend(paths)
    return files


rule qc_summary:
    input:
        **qc_inputs(),
    output:
        table=f"{OUTDIR}/qc.tsv",
        plot=f"{OUTDIR}/qc.svg",
        plot_png=f"{OUTDIR}/qc.png",
    log:
        f"{LOGDIR}/qc_summary.txt",
    benchmark:
        f"{BENCHDIR}/qc_summary.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/qc_summary.py \
            --cutadapt {input.cutadapt:q} \
            --flagstat {input.flagstat:q} \
            --stats {input.stats:q} \
            --bamcounts {input.bamcounts:q} \
            --markdup {input.markdup:q} \
            --preseq {input.preseq:q} \
            --insert-sizes {input.insert_sizes:q} \
            --fit-stats {input.fit_stats:q} \
            --total-cutcounts {input.total_cutcounts:q} \
            --calls {input.calls:q} \
            --footprints {input.footprints:q} \
            --dispersion-model {input.dispersion_model:q} \
            --table {output.table:q} \
            --plot {output.plot:q}

        sync {output:q}
        """
