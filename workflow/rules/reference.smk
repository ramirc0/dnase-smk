# Run-independent references, built once under results/refs/ and reused across runs.


# ftd and hotspot3 read the FASTA through its .fai, so both live in refs/. The chrom
# sizes keep ENCODE's hotspot3 order, sorted by name.
rule reference:
    input:
        REFERENCES["fasta"],
    output:
        fasta=f"{REFDIR}/reference.fa",
        fai=f"{REFDIR}/reference.fa.fai",
        sizes=f"{REFDIR}/nuclear.chrom.sizes",
        bed=f"{REFDIR}/nuclear.bed",
    params:
        chroms=lambda w: REFERENCES["nuclear_chroms"],
    log:
        f"{LOGDIR}/reference.txt",
    benchmark:
        f"{BENCHDIR}/reference.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cp {input:q} {output.fasta:q}
        samtools faidx {output.fasta:q}
        printf '%s\n' {params.chroms:q} \
            | awk 'BEGIN {{ FS = OFS = "\t" }} NR == FNR {{ keep[$1]; next }} $1 in keep {{ print $1, $2 }}' \
                - {output.fai:q} \
            | sort -k1,1 > {output.sizes:q}
        if [[ ! -s {output.sizes:q} ]]; then
            echo "references.nuclear_chroms names no contig of the FASTA" >&2
            exit 1
        fi
        awk 'BEGIN {{ OFS = "\t" }} {{ print $1, 0, $2 }}' {output.sizes:q} > {output.bed:q}

        sync {output:q}
        """


# hotspot3 reads the mappable bases with tabix, which needs BGZF.
rule mappable_bases:
    input:
        REFERENCES["mappable_bases"],
    output:
        bed=f"{REFDIR}/mappable_bases.bed.gz",
        tbi=f"{REFDIR}/mappable_bases.bed.gz.tbi",
    log:
        f"{LOGDIR}/mappable_bases.txt",
    benchmark:
        f"{BENCHDIR}/mappable_bases.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        rapidgzip -d -c {input:q} | bgzip -@ {threads} > {output.bed:q}
        tabix -p bed {output.bed:q}

        sync {output:q}
        """


rule bwa_mem3_index:
    input:
        rules.reference.output.fasta,
    output:
        multiext(BWA_INDEX, ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    params:
        prefix=BWA_INDEX,
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/bwa_mem3_index.txt",
    benchmark:
        f"{BENCHDIR}/bwa_mem3_index.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        bwa-mem3 index -p {params.prefix:q} -t {threads} --tmp-dir "$tmp" {input:q}

        sync {output:q}
        """
