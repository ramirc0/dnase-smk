#!/usr/bin/env bash
# Builds resources/fixture: the chr22 reads, and every tenth chrM read, of three ENCODE
# DNase-seq replicates, plus the chr22 + chrM genome and the mappable bases cut to chr22.
# Reads come from the raw FASTQs, picked by the names ENCODE's v3 unfiltered BAMs place
# there. A pair goes by its first mate.
#   UHT: ENCSR693UHT rep1, two PE 76 lanes. BAM ENCFF561MEY.
#   IEA: ENCSR426IEA rep1, four PE 36 lanes and one SE 36 lane. BAM ENCFF919CXW.
#   RWU: ENCSR420RWU rep1, two SE 36 lanes. BAM ENCFF021ZDN.
# Needs .tests/download.sh and the Weng Lab ENCODE mirror. Run from the repository root:
#   pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/integration/make_fixture.sh
set -euo pipefail
out=resources/fixture
refs=/data/zusers/ramirezc/dnase-smk-refs
encode=/data/projects/encode/data
fasta=/data/zusers/ramirezc/atac-smk/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
mkdir -p "$out"
tmp=$(mktemp -d -p "$out")
trap 'rm -rf "$tmp"' EXIT

samtools faidx "$fasta" chr22 chrM > "$out/genome.fa"
rapidgzip -d -c "$refs/GRCh38_no_alts.K36.n150.center_sites_and_extended_blacklist.bed.gz" \
    | awk '$1 == "chr22"' | bgzip > "$out/mappable_bases.bed.gz"
cp "$refs/vierstra_et_al.6mer-model.txt" "$out/"

names() {
    samtools view -@ 4 -F 0x980 "$encode/$1.bam" \
        | awk '$3 == "chr22" || ($3 == "chrM" && mito++ % 10 == 0) { print $1 }' > "$tmp/$2"
}
lane() {
    seqkit grep -j 4 -f "$tmp/$1" "$refs/fastq/$2.fastq.gz" -o "$out/$2.fastq.gz"
}
export -f names lane
export tmp out refs encode

parallel --colsep ' ' --halt now,fail=1 names <<'EOF'
ENCSR693UHT/ENCFF561MEY UHT
ENCSR426IEA/ENCFF919CXW IEA
ENCSR420RWU/ENCFF021ZDN RWU
EOF
parallel --colsep ' ' --halt now,fail=1 lane <<'EOF'
UHT ENCFF201UOC
UHT ENCFF149BGS
UHT ENCFF752UBC
UHT ENCFF559DWT
IEA ENCFF178WLB
IEA ENCFF645INN
IEA ENCFF377AOZ
IEA ENCFF536KHM
IEA ENCFF665FKP
IEA ENCFF580MND
IEA ENCFF901FSW
IEA ENCFF719IMH
IEA ENCFF302JOP
RWU ENCFF002DYE
RWU ENCFF002DYD
EOF
