# Nuclear BAM to hotspot3 hotspots, peaks and density.


def hotspot3_calls(sample):
    """Return the hotspot and peak files hotspot3 writes for a sample.

    Parameters
    ----------
    sample : str
        Sample ID, or a `{sample}` wildcard pattern.

    Returns
    -------
    list of str
        bed.gz, its .tbi and the bigBed of hotspots and peaks at each FDR of
        `HOTSPOT_FDRS`.
    """
    return [
        f"{OUTDIR}/{sample}/hotspot3/fdr{fdr}/{sample}.{kind}.fdr{fdr}.{ext}"
        for fdr in HOTSPOT_FDRS
        for kind in ("hotspots", "peaks")
        for ext in ("bed.gz", "bed.gz.tbi", "bb")
    ]


# hotspot3 1.0.1 runs its cutcount script through /bin/sh with `set -o pipefail`, which
# dash rejects, so the rule runs the script with bash and passes the cutcounts in.
# Given cutcounts, hotspot3 also expects their total beside them. It only warns when
# bedToBigBed fails, so the bigBeds are declared. An FDR with no peaks writes no file
# and fails the job. debug/ and the p-values take about 6.5 GB on GRCh38 and nothing
# reads them.
rule hotspot3:
    input:
        bam=rules.filter.output.bam,
        bai=rules.filter.output.bai,
        sizes=rules.reference.output.sizes,
        mappable=rules.mappable_bases.output.bed,
        mappable_tbi=rules.mappable_bases.output.tbi,
    output:
        cutcounts=f"{OUTDIR}/{{sample}}/hotspot3/{{sample}}.cutcounts.bed.gz",
        cutcounts_tbi=f"{OUTDIR}/{{sample}}/hotspot3/{{sample}}.cutcounts.bed.gz.tbi",
        total=f"{OUTDIR}/{{sample}}/hotspot3/{{sample}}.total_cutcounts",
        fit_stats=f"{OUTDIR}/{{sample}}/hotspot3/{{sample}}.fit_stats.bed.gz",
        density=f"{OUTDIR}/{{sample}}/hotspot3/{{sample}}.normalized_density.bw",
        background=f"{OUTDIR}/{{sample}}/hotspot3/{{sample}}.background.bw",
        calls=hotspot3_calls("{sample}"),
    params:
        outdir=lambda w: f"{OUTDIR}/{w.sample}/hotspot3",
        fdrs=HOTSPOT_FDRS,
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/hotspot3/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/hotspot3/{{sample}}.tsv"
    software:
        HOTSPOT3_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        script=$(python -c 'import importlib.resources as r; print(r.files("hotspot3.scripts") / "extract_cutcounts.sh")')
        bash "$script" {input.bam:q} $(cut -f 1 {input.sizes:q}) | bgzip > {output.cutcounts:q}
        tabix -p bed {output.cutcounts:q}
        zcat {output.cutcounts:q} | awk '!/^#/ {{ s += $4 }} END {{ printf "%d\n", s }}' > {output.total:q}
        hotspot3 {wildcards.sample:q} \
            --cutcounts {output.cutcounts:q} \
            --chrom_sizes {input.sizes:q} \
            --mappable_bases {input.mappable:q} \
            --fdrs {params.fdrs} \
            --save_density \
            --cpus {threads} \
            --outdir {params.outdir:q} \
            --tempdir "$tmp"
        rm -rf {params.outdir:q}/debug {params.outdir:q}/{wildcards.sample:q}.pvals.parquet

        sync {output:q}
        """
