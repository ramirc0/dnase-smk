# Nuclear BAM and 5% hotspots to footprints, as ENCODE pipeline ENCPL336PGV's.


def footprint_files(sample):
    """Return the footprint files of a sample.

    Parameters
    ----------
    sample : str
        Sample ID, or a `{sample}` wildcard pattern.

    Returns
    -------
    list of str
        bed.gz and bigBed at each FDR of `FOOTPRINT_FDRS`.
    """
    return [
        f"{OUTDIR}/{sample}/footprints/{sample}.footprints.fdr{fdr}.{ext}"
        for fdr in FOOTPRINT_FDRS
        for ext in ("bed.gz", "bb")
    ]


rule dispersion_model:
    input:
        hotspots=f"{OUTDIR}/{{sample}}/hotspot3/fdr0.05/{{sample}}.hotspots.fdr0.05.bed.gz",
        bam=rules.filter.output.bam,
        bai=rules.filter.output.bai,
        fasta=rules.reference.output.fasta,
        fai=rules.reference.output.fai,
        bias=REFERENCES["bias_model"],
    output:
        f"{OUTDIR}/{{sample}}/footprints/{{sample}}.dm.json",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/dispersion_model/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/dispersion_model/{{sample}}.tsv"
    software:
        FOOTPRINTS_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        zcat {input.hotspots:q} | grep -v '^#' | cut -f 1-3 > "$tmp/hotspots.bed"
        ftd learn_dm "$tmp/hotspots.bed" {input.bam:q} {input.fasta:q} \
            --bias_model_file {input.bias:q} \
            --half_win_width 5 \
            --seed 12345 \
            --n_threads {threads} \
            --outfile {output:q}

        sync {output:q}
        """


# ENCODE's v3 settings. detect's shuffle count is passed because its default is now 100.
# Thresholding skips the bedgraph's # header, which ENCODE's awk let through, and clamps
# the padded start at 0. bedToBigBed reads its input twice, so it gets a file.
rule footprints:
    input:
        hotspots=rules.dispersion_model.input.hotspots,
        bam=rules.filter.output.bam,
        bai=rules.filter.output.bai,
        fasta=rules.reference.output.fasta,
        fai=rules.reference.output.fai,
        sizes=rules.reference.output.sizes,
        bias=REFERENCES["bias_model"],
        dm=rules.dispersion_model.output,
    output:
        calls=footprint_files("{sample}"),
        bedgraph=temp(f"{OUTDIR}/{{sample}}/footprints/{{sample}}.interval.all.bedgraph"),
    params:
        prefix=lambda w: f"{OUTDIR}/{w.sample}/footprints/{w.sample}.footprints",
        fdrs=FOOTPRINT_FDRS,
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/footprints/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/footprints/{{sample}}.tsv"
    software:
        FOOTPRINTS_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        zcat {input.hotspots:q} | grep -v '^#' | cut -f 1-3 > "$tmp/hotspots.bed"
        ftd detect "$tmp/hotspots.bed" {input.bam:q} {input.fasta:q} \
            --bias_model_file {input.bias:q} \
            --dispersion_model_file {input.dm:q} \
            --half_win_width 5 \
            --smooth_half_win_width 50 \
            --smooth_clip 0.01 \
            --fdr_shuffle_n 50 \
            --seed 12345 \
            --n_threads {threads} \
            --outprefix "$tmp/interval.all"
        mv "$tmp/interval.all.bedgraph" {output.bedgraph:q}
        for t in {params.fdrs}; do
            awk -F '\t' -v OFS='\t' -v t="$t" \
                '!/^#/ && $8 <= t {{ s = $2 - 3; print $1, (s < 0 ? 0 : s), $3 + 3 }}' \
                {output.bedgraph:q} \
                | sort-bed - \
                | bedops --merge - > "$tmp/fp.bed"
            bgzip -c "$tmp/fp.bed" > {params.prefix:q}.fdr"$t".bed.gz
            bedToBigBed -type=bed3 "$tmp/fp.bed" {input.sizes:q} {params.prefix:q}.fdr"$t".bb
        done

        sync {output:q}
        """
