#!/usr/bin/env bash
# Build TerL marker maximum-likelihood phylogenies for the Desulfovibrio-associated
# viral catalog, against an INPHARED reference backbone.
#
#   viral catalog contigs -> pyrodigal (meta) -> hmmsearch vs Pfam terminase
#   large-subunit profiles -> one marker protein per contig and per reference
#   genus -> MAFFT -> trimAl -> IQ-TREE (ModelFinder, UFBoot, SH-aLRT)
#
# Source archive only: this upstream reconstruction is not a default release step.
# All input resources must be supplied explicitly. Output must be a new directory.
set -euo pipefail

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${PC057_TERL_OUTPUT:?Set PC057_TERL_OUTPUT to a new output directory}"
CATALOG="${PC057_VIRAL_FASTA:?Set PC057_VIRAL_FASTA to the archived viral catalogue FASTA}"
DB="${PC057_DATABASE:?Set PC057_DATABASE to the annotation SQLite database}"
PFAM="${PFAM_HMM:?Set PFAM_HMM to the Pfam-A HMM file}"
INPHARED_DIR="${INPHARED_DIR:?Set INPHARED_DIR to the original reference data directory}"
INPHARED_FAA_GZ="${INPHARED_DIR}/7Apr2026_vConTACT2_proteins.faa.gz"
INPHARED_TSV="${INPHARED_DIR}/7Apr2026_data.tsv.gz"
[[ ! -e "${OUT}" ]] || { echo "Refuse to reuse an existing output directory" >&2; exit 1; }
THREADS="${THREADS:-10}"
SEED=20260917
PROFILES=(Terminase_6N GpA_ATPase)
# Pfam accessions of the terminase large-subunit families, searched together so
# that recovery can be reported per family before the tree profiles are chosen.
PFAM_ACC=(PF03354.22 PF20441.4 PF04466.20 PF03237.22 PF17289.9 \
          PF05876.19 PF20454.5 PF06056.19 PF11053.13 PF17288.9)

mkdir -p "${OUT}"
cd "${OUT}"

echo "[1/6] terminase HMM profiles"
if [[ ! -s terl_profiles.hmm ]]; then
    for acc in "${PFAM_ACC[@]}"; do hmmfetch "${PFAM}" "${acc}"; done > terl_profiles.hmm
fi

echo "[2/6] gene calling on the viral catalog"
if [[ ! -s viral_catalog.faa ]]; then
    pyrodigal -p meta -i "${CATALOG}" -a viral_catalog.faa -f gff -o viral_catalog.gff -j "${THREADS}"
fi

echo "[3/6] hmmsearch, queries and INPHARED references"
if [[ ! -s query_terl.domtbl ]]; then
    hmmsearch --cpu "${THREADS}" -E 1e-5 --domtblout query_terl.domtbl -o /dev/null \
        terl_profiles.hmm viral_catalog.faa
fi
if [[ ! -s ref_terl.domtbl ]]; then
    zcat "${INPHARED_FAA_GZ}" > inphared.faa
    hmmsearch --cpu "${THREADS}" -E 1e-5 --domtblout ref_terl.domtbl -o /dev/null \
        terl_profiles.hmm inphared.faa
fi
[[ -s inphared.faa ]] || zcat "${INPHARED_FAA_GZ}" > inphared.faa

echo "[4/6] marker selection"
if [[ ! -s terl_taxa_metadata.tsv ]]; then
    python3 "${SOURCE_DIR}/terl_select_seqs.py" \
        --query-faa viral_catalog.faa --query-domtbl query_terl.domtbl \
        --ref-faa inphared.faa --ref-domtbl ref_terl.domtbl \
        --inphared-tsv "${INPHARED_TSV}" --db "${DB}" --outdir "${OUT}" \
        --profiles "${PROFILES[@]}"
fi

echo "[5/6] alignment and trimming"
for P in "${PROFILES[@]}"; do
    [[ -s "${P}.aln.faa" ]]  || mafft --auto --thread "${THREADS}" --quiet "${P}.dom.faa" > "${P}.aln.faa"
    [[ -s "${P}.trim.faa" ]] || trimal -in "${P}.aln.faa" -out "${P}.trim.faa" -automated1
done

echo "[6/6] maximum-likelihood inference"
for P in "${PROFILES[@]}"; do
    # ModelFinder is restricted to LG and Q.pfam under gamma or free-rate
    # heterogeneity: Q.pfam is estimated on Pfam domain alignments, which is what
    # these are, and the full matrix sweep costs hours at this taxon count.
    # --nmax caps the tree search: at this taxon count small likelihood gains keep
    # resetting the default stopping rule long after UFBoot has converged.
    [[ -s "${P}_terl.contree" ]] || iqtree -s "${P}.trim.faa" -m MFP \
        -mset LG,Q.pfam -mrate G4,R4 --nmax 400 \
        -B 1000 -alrt 1000 -T "${THREADS}" -seed "${SEED}" -pre "${P}_terl" --quiet -redo
done

rm -f inphared.faa
echo "done -> ${OUT}"
