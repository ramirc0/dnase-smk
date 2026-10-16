# Unfiltered BAM to the nuclear BAM: flag, mark duplicates, keep passing nuclear reads.


# flag_reads.py reads bwa's TLEN, so it runs before fixmate, whose -p keeps the 0x2 it
# set. markdup --include-fails puts QC-fail reads in duplicate sets, as Picard did.
# -d 100 is Picard's optical distance. Nothing removes duplicates: hotspot3 and ftd
# count them, as hotspot2 did.
rule filter:
    input:
        bam=rules.merge.output,
        bed=rules.reference.output.bed,
    output:
        bam=f"{OUTDIR}/{{sample}}/filter/{{sample}}.nuclear.bam",
        bai=f"{OUTDIR}/{{sample}}/filter/{{sample}}.nuclear.bam.bai",
        markdup=f"{OUTDIR}/{{sample}}/filter/{{sample}}.markdup.json",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/filter/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/filter/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        python {workflow.basedir:q}/scripts/flag_reads.py {input.bam:q} - \
            | samtools fixmate -@ {threads} -m -p -u - - \
            | samtools sort -@ {threads} -m 1536M -u -T "$tmp/sort" - \
            | samtools markdup -@ {threads} --include-fails -d 100 -T "$tmp/markdup" \
                -f {output.markdup:q} --json -u - - \
            | samtools view -@ {threads} -b -F 512 -L {input.bed:q} -o {output.bam:q} -
        samtools index -@ {threads} {output.bam:q}

        sync {output:q}
        """
