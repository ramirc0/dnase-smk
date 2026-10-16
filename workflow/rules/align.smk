# FASTQs to the unfiltered BAM. Each layout of a sample aligns on its own, as ENCODE's.


def cutadapt_flags(sample, layout):
    """Return the cutadapt trimming flags of a sample's layout.

    Parameters
    ----------
    sample : str
        Sample ID from the sheet.
    layout : {"pe", "se"}
        Layout.

    Returns
    -------
    list of str
        ENCODE's paired-end adapter trimming, then the crop to `read_length`.
        Interpolate them with `:q`.
    """
    s = SAMPLES[sample]
    flags = []
    if layout == "pe":
        r1, r2 = s["adapters"]
        flags += ["-a", r1, "-A", r2, "--pair-adapters", "-e", "0.1", "-m", "1", "--interleaved"]
    if s["read_length"]:
        flags += ["--length", str(s["read_length"])]
    return flags


# cutadapt crops after it trims adapters and before -m filters, as ENCODE's awk crop
# follows cutadapt. --discard-casava does what bwa aln -Y did. -K fixes bwa-mem3's
# batch size, so its output does not depend on the thread count. cutadapt cannot
# sniff the compression of a pipe, so it reads decompressed lanes.
rule align:
    input:
        r1=lambda w: LANES[w.sample, w.layout]["r1"],
        r2=lambda w: LANES[w.sample, w.layout]["r2"],
        index=rules.bwa_mem3_index.output,
    output:
        bam=temp(f"{OUTDIR}/{{sample}}/align/{{sample}}.{{layout}}.bam"),
        json=f"{OUTDIR}/{{sample}}/align/{{sample}}.{{layout}}.cutadapt.json",
    params:
        trim=lambda w: cutadapt_flags(w.sample, w.layout),
        pairs=lambda w: ["-p"] * (w.layout == "pe"),
        index=BWA_INDEX,
        rg=lambda w: rf"@RG\tID:{w.sample}\tSM:{w.sample}",
    log:
        f"{LOGDIR}/align/{{sample}}.{{layout}}.txt",
    benchmark:
        f"{BENCHDIR}/align/{{sample}}.{{layout}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        r1=({input.r1:q})
        r2=({input.r2:q})
        trim() {{ cutadapt -j {threads} --discard-casava {params.trim:q} --json {output.json:q} -o - "$@"; }}
        # A process substitution closes when its command ends, so each sits on trim's line.
        if (( ${{#r2[@]}} )); then
            trim <(cat "${{r1[@]}}" | rapidgzip -d -c) <(cat "${{r2[@]}}" | rapidgzip -d -c)
        else
            trim <(cat "${{r1[@]}}" | rapidgzip -d -c)
        fi | bwa-mem3 mem {params.pairs:q} -K 100000000 -R {params.rg:q} -t {threads} --bam=6 \
            -o {output.bam:q} {params.index:q} -

        sync {output:q}
        """


# samtools cat keeps aligner order, so a pair's records stay adjacent.
rule merge:
    input:
        lambda w: expand(
            f"{OUTDIR}/{{sample}}/align/{{sample}}.{{layout}}.bam",
            sample=w.sample,
            layout=LAYOUTS[w.sample],
        ),
    output:
        f"{OUTDIR}/{{sample}}/align/{{sample}}.unfiltered.bam",
    log:
        f"{LOGDIR}/merge/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/merge/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools cat -o {output:q} {input:q}

        sync {output:q}
        """
