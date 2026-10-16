#!/usr/bin/env bash
# Downloads what the integration and validation tests read into $refs: the raw FASTQs
# (about 28 GB), md5-checked against the ENCODE portal, the Altius mappable bases that
# hotspot3's README names, the footprint-tools bias model, and ENCODE's v3.0.0 image.
# Files already present are kept. Run from the repository root:
#   bash .tests/download.sh
set -euo pipefail
refs=/data/zusers/ramirezc/dnase-smk-refs
mkdir -p "$refs/fastq" "$refs/containers"

fastq() {
    local acc=$1 md5=$2 dest=$refs/fastq/$1.fastq.gz
    [[ -s "$dest" ]] && return
    curl -sSfL --retry 5 -o "$dest.part" "https://www.encodeproject.org/files/$acc/@@download/$acc.fastq.gz"
    if [[ "$(md5sum < "$dest.part" | cut -d' ' -f1)" != "$md5" ]]; then
        echo "$acc: md5 mismatch" >&2
        return 1
    fi
    mv "$dest.part" "$dest"
}
export -f fastq
export refs

# Accession and portal md5. ENCSR693UHT rep1, ENCSR426IEA reps 1 and 2, ENCSR420RWU rep1.
parallel --colsep ' ' -j 6 --halt now,fail=1 fastq <<'EOF'
ENCFF201UOC fe150ba56063307ec2a65074741f39ee
ENCFF149BGS 84673654e67e5c7d9b65c214790fb63b
ENCFF752UBC bcccccd17f984a33d885df5e758f4c54
ENCFF559DWT 8e3d36c99d97464cdd39ce38372d150c
ENCFF178WLB 6f33e443737d272c164a59f852e5e32e
ENCFF645INN a80d511100c9f6bbb1d15e5099be129f
ENCFF377AOZ f4678718b3a03b3badf067cf21b7905f
ENCFF536KHM 8b9fefa180981b3f1d6b1ec19eee26e1
ENCFF665FKP 06557db02d42769371a3a284813502f1
ENCFF580MND 85b275ce275671632b7ed1939f8fccb8
ENCFF901FSW 3437acbe4f9a66cdc7fc5a17f71db50e
ENCFF719IMH d33b88414a8efa3e1accf8b5a28311ec
ENCFF302JOP 3c0a04e82e73ba4297392d7de63014a6
ENCFF266WWA 1ad9fbaaab47eb84d6e1fd19ab2afa05
ENCFF314WHE 009e4ed53c6800900edfbc09d4251ab6
ENCFF355JSW 67bb958248f81bc8c2aea7e55b225fcc
ENCFF002DYE be48b51b79ef3a8a0d5b86884a0ef33c
ENCFF002DYD 400bfc5a6052ac3efcaa547ff082f009
EOF

mappable=$refs/GRCh38_no_alts.K36.n150.center_sites_and_extended_blacklist.bed.gz
[[ -s "$mappable" ]] || curl -sSfL --retry 5 -o "$mappable" \
    "https://resources.altius.org/~sabramov/files/hotspot3/${mappable##*/}"

bias=$refs/vierstra_et_al.6mer-model.txt
[[ -s "$bias" ]] || curl -sSfL --retry 5 -o "$bias" \
    "https://raw.githubusercontent.com/vierstralab/footprint-tools/e839b4e6d73dbd2cd0b2791b5f19fd2743ab147c/data/${bias##*/}"

sif=$refs/containers/dnase-seq-pipeline_v3.0.0.sif
[[ -s "$sif" ]] || apptainer pull "$sif" docker://encodedcc/dnase-seq-pipeline:v3.0.0
