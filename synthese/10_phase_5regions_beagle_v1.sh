#!/usr/bin/env bash
# Script 10 — Phasage BEAGLE des 5 régions A sans VCF phasé existant
# Entrée  : analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/filtered_vcf/chrN.PASS_biallelic_snps.vcf.gz
# Sortie  : analyses/phasing_beagle/Phasing_Beagle_v1_5regions/
# Méthode : extraction ±1Mb autour de chaque région → phasage BEAGLE → index
#
# IMPORTANT : chemin absolu en dur ligne 8, à adapter si le dépôt est cloné ailleurs.
# Nécessite Beagle 5.5 (Browning et al. 2018, https://faculty.washington.edu/browning/beagle/beagle.html)
# installé sous BEAGLE_JAR ci-dessous — pas fourni dans ce dépôt.
set -euo pipefail
cd /home/tanguyruel/Bureau/genome_complet_Awassi  # chemin en dur à adapter

BEAGLE_JAR="tools/envs/beagle-awassi/share/beagle-5.5_27Feb25.75f-0/beagle.jar"  # chemin vers le .jar Beagle
VCF_DIR="analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/filtered_vcf"  # VCF déjà filtrés PASS+biallélique
OUTDIR="analyses/phasing_beagle/Phasing_Beagle_v1_5regions"
LOGDIR="logs"
mkdir -p "$OUTDIR/input" "$OUTDIR/phased"  # crée les sous-dossiers de sortie si absents

TIMESTAMP=$(date +%d_%m_%H%M)  # horodatage pour nommer les logs de façon unique

phase_region() {
    local LABEL="$1"   # ex: chr14_35.41_35.57Mb
    local CHR="$2"
    local START="$3"   # coordonnée start de la région A
    local END="$4"     # coordonnée end de la région A
    local FLANK=1000000  # marge de 1 Mb ajoutée de chaque côté (contexte utile au phasage statistique)

    local WIN_START=$(( START - FLANK ))
    local WIN_END=$(( END + FLANK ))
    # clamp à 1
    [[ $WIN_START -lt 1 ]] && WIN_START=1  # évite une coordonnée négative si la région est près du début du chromosome

    local SRC="${VCF_DIR}/chr${CHR}.PASS_biallelic_snps.vcf.gz"
    local UNPHASED="${OUTDIR}/input/${LABEL}.unphased.vcf.gz"
    local PHASED_PREFIX="${OUTDIR}/phased/${LABEL}.beagle_phased"

    echo "=== ${LABEL} (chr${CHR}:${WIN_START}-${WIN_END}) ==="

    echo "  1) Extraction bcftools..."
    bcftools view \
        -r "${CHR}:${WIN_START}-${WIN_END}" \
        "$SRC" \
        -Oz -o "$UNPHASED"  # extrait uniquement la fenêtre élargie, compresse en sortie (Oz)
    bcftools index -t "$UNPHASED"  # index tabix, nécessaire pour les accès par région ensuite
    echo "     SNPs extraits : $(bcftools view -H "$UNPHASED" | wc -l)"  # compte les lignes de variants (sans l'en-tête)

    echo "  2) Phasage BEAGLE..."
    java -Xmx4g -jar "$BEAGLE_JAR" \
        gt="$UNPHASED" \
        out="$PHASED_PREFIX" \
        impute=false \
        nthreads=8 \
        seed=12345 \
        > "${LOGDIR}/beagle_phase_${LABEL}_${TIMESTAMP}.log" 2>&1
        # gt=       : VCF non phasé en entrée
        # impute=false : phasage seul, pas d'imputation de génotypes manquants
        # nthreads=8   : parallélisation Beagle
        # seed=12345   : graine fixe, pour un résultat reproductible d'un run à l'autre

    bcftools index -t "${PHASED_PREFIX}.vcf.gz"  # index du VCF phasé produit par Beagle

    echo "  3) Vérification..."
    echo "     SNPs phasés   : $(bcftools view -H "${PHASED_PREFIX}.vcf.gz" | wc -l)"
    echo "     Génotypes |   : $(bcftools query -f '[%GT\n]' "${PHASED_PREFIX}.vcf.gz" | grep -c "|" || true)"  # "|" = allèles phasés
    echo "     Génotypes /   : $(bcftools query -f '[%GT\n]' "${PHASED_PREFIX}.vcf.gz" | grep -c "/" || true)"  # "/" = restés non phasés
    echo ""
}

# 5 régions A sans VCF phasé existant — un appel par région (label, chromosome, début, fin)
phase_region "chr14_35.41_35.57Mb"  14  35415001  35570000
phase_region "chr5_48.91_49.01Mb"    5  48910001  49015000
phase_region "chr17_34.16_34.22Mb"  17  34165001  34220000
phase_region "chr3_186.07_186.12Mb"  3 186075001 186125000
phase_region "chr9_77.59_77.69Mb"    9  77590001  77695000

echo "=== PHASAGE TERMINÉ ==="
echo "VCF phasés dans : ${OUTDIR}/phased/"
