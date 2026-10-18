# Read-level QC, as ENCODE's: stats of both BAMs, bamcounts, preseq and insert sizes.

BAMS = {
    "unfiltered": f"{OUTDIR}/{{sample}}/align/{{sample}}.unfiltered.bam",
    "nuclear": f"{OUTDIR}/{{sample}}/filter/{{sample}}.nuclear.bam",
}
PRESEQ_TARGETS = (20000000, 30000000, 40000000, 100000000, 150000000, 250000000)


rule bam_qc:
    input:
        lambda w: BAMS[w.bam].format(sample=w.sample),
    output:
        stats=f"{OUTDIR}/{{sample}}/qc/{{sample}}.{{bam}}.stats.txt",
        flagstat=f"{OUTDIR}/{{sample}}/qc/{{sample}}.{{bam}}.flagstat.json",
    log:
        f"{LOGDIR}/bam_qc/{{sample}}.{{bam}}.txt",
    benchmark:
        f"{BENCHDIR}/bam_qc/{{sample}}.{{bam}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools stats -@ {threads} {input:q} > {output.stats:q}
        samtools flagstat -@ {threads} -O json {input:q} > {output.flagstat:q}

        sync {output:q}
        """


rule bamcounts:
    input:
        rules.merge.output,
    output:
        f"{OUTDIR}/{{sample}}/qc/{{sample}}.bamcounts.txt",
    log:
        f"{LOGDIR}/bamcounts/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bamcounts/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/bamcounts.py {input:q} {output:q}

        sync {output:q}
        """


rule dup_hist:
    input:
        rules.filter.output.bam,
    output:
        f"{OUTDIR}/{{sample}}/qc/{{sample}}.dup_hist.txt",
    log:
        f"{LOGDIR}/dup_hist/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/dup_hist/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/dup_hist.py {input:q} {output:q}

        sync {output:q}
        """


# The targets are stampipes' preseq_targets.sh: the first extrapolated depth with more
# distinct reads than each target, then the most distinct reads. Defects mode can write
# nan past some depth, which stampipes' awk compares as a string above every target.
rule preseq:
    input:
        rules.dup_hist.output,
    output:
        curve=f"{OUTDIR}/{{sample}}/qc/{{sample}}.preseq.txt",
        targets=f"{OUTDIR}/{{sample}}/qc/{{sample}}.preseq_targets.txt",
    params:
        defects=lambda w: ["-defects"] * SAMPLES[w.sample]["defects"],
        targets=PRESEQ_TARGETS,
    log:
        f"{LOGDIR}/preseq/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/preseq/{{sample}}.tsv"
    software:
        PRESEQ_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        preseq lc_extrap -hist {input:q} -extrap 1.001e9 -step 1e6 {params.defects:q} \
            -o {output.curve:q}
        for target in {params.targets}; do
            awk -v target="$target" \
                'NR > 1 && $2 != "nan" && $2 > target {{ printf "preseq-est-for-%d\t%d\n", target, $1; exit }}' \
                {output.curve:q}
        done > {output.targets:q}
        awk 'NR > 1 && $2 != "nan" {{ max = $2 }} END {{ printf "preseq-est-max\t%d\n", max }}' {output.curve:q} \
            >> {output.targets:q}

        sync {output:q}
        """


# Picard counted each non-duplicate pair once, by its |TLEN|. samtools stats' IS
# section stops near 99% of pairs, so the histogram comes from the first mates.
rule insert_sizes:
    input:
        rules.filter.output.bam,
    output:
        hist=f"{OUTDIR}/{{sample}}/qc/{{sample}}.insert_sizes.hist.txt",
        metrics=f"{OUTDIR}/{{sample}}/qc/{{sample}}.insert_sizes.tsv",
        plot=f"{OUTDIR}/{{sample}}/qc/{{sample}}.insert_sizes.svg",
        plot_png=f"{OUTDIR}/{{sample}}/qc/{{sample}}.insert_sizes.png",
    log:
        f"{LOGDIR}/insert_sizes/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/insert_sizes/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools view -@ {threads} -F 0x400 -f 0x41 {input:q} \
            | awk '{{ t = $9 < 0 ? -$9 : $9; if (t) n[t]++ }} END {{ for (t in n) print t "\t" n[t] }}' \
            | sort -k1,1n > {output.hist:q}
        python {workflow.basedir:q}/scripts/insert_sizes.py {output.hist:q} \
            --output {output.metrics:q} \
            --plot {output.plot:q}

        sync {output:q}
        """
