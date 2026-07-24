#!/bin/bash
#
# PCA de référence : 25 chromosomes, tous les individus (387), filtrage + LD pruning complets
# (pas de sous-échantillonnage géographique/spatial - voir PCA A/B/C pour les variantes rapides).
# Entrée : $VCF_LIST, $METADATA_VCF_ORDER (config_project.sh).
# Sortie : analyses/pca/PCA_<date>/ (VCF filtré, PLINK, plots via 02_plot_pca_25chr.R).
# Usage  : ./01_run_pca_25chr_filtered.sh

set -euo pipefail  # arrête le script en cas d'erreur, de variable non définie ou d'échec dans un pipe

cd ~/Bureau/genome_complet_Awassi  # se place à la racine du projet
source config_project.sh  # charge les variables de configuration (VCF_LIST, PROJECT_DIR, etc.)

RUN_DATE=$(date +%d_%m_%Hh%M)  # horodatage utilisé pour nommer le dossier de ce run
D_RUN="$PROJECT_DIR/analyses/pca/PCA_${RUN_DATE}"  # dossier racine de sortie pour ce run

D_FILTERED="$D_RUN/filtered_vcf"  # sous-dossier des VCF filtrés
D_PLINK="$D_RUN/plink"  # sous-dossier des fichiers PLINK
D_PLOTS="$D_RUN/plots"  # sous-dossier des figures
D_CHECKS="$D_RUN/checks"  # sous-dossier des fichiers de vérification
D_RUN_LOGS="$D_RUN/logs"  # sous-dossier des logs

mkdir -p "$D_FILTERED" "$D_PLINK" "$D_PLOTS" "$D_CHECKS" "$D_RUN_LOGS"  # crée tous les sous-dossiers de sortie

cp "$PROJECT_DIR/config_project.sh" "$D_RUN/config_project_used.sh"  # copie la config utilisée pour traçabilité

echo "=== PCA 25 chromosomes filtrés avec PLINK 1.9 ==="
echo "Run : $D_RUN"
echo "VCF list : $VCF_LIST"
echo "Metadata : $METADATA_VCF_ORDER"
echo

command -v bcftools >/dev/null || { echo "ERREUR : bcftools introuvable"; exit 1; }  # vérifie que bcftools est installé
command -v plink >/dev/null || { echo "ERREUR : plink introuvable"; exit 1; }  # vérifie que plink est installé
command -v Rscript >/dev/null || { echo "ERREUR : Rscript introuvable"; exit 1; }  # vérifie que Rscript est installé

plink --version  # affiche la version de plink utilisée (traçabilité)

N_VCF=$(wc -l < "$VCF_LIST")  # compte le nombre de VCF listés
if [ "$N_VCF" -ne 25 ]; then  # vérifie qu'il y a bien 25 chromosomes dans la liste
  echo "ERREUR : la liste VCF ne contient pas 25 chromosomes : $N_VCF"
  exit 1
fi

echo
echo "1) Filtrage chromosome par chromosome"
echo "Filtres : SNP bialléliques, PASS, MAF >= 0.05, NS >= 349, donc missing variants <= 0.10"
echo

FILTERED_LIST="$D_FILTERED/list_filtered_vcf.txt"  # fichier qui listera les VCF filtrés
> "$FILTERED_LIST"  # crée/vide ce fichier avant la boucle

while read -r VCF; do  # boucle sur chaque VCF (un par chromosome) listé dans VCF_LIST
  base=$(basename "$VCF")  # nom du fichier sans le chemin
  chr=$(echo "$base" | sed -E 's/.*chr([0-9]+)\.vcf\.gz/\1/')  # extrait le numéro de chromosome du nom de fichier
  OUT="$D_FILTERED/awassi_basedata_chr${chr}.filtered.vcf.gz"  # chemin de sortie du VCF filtré pour ce chromosome

  echo "Filtrage chr${chr} : $base"

  # filtre SNP bialléliques PASS (-m2 -M2 -v snps -f PASS), calcule MAF et NS (+fill-tags),
  # puis ne garde que les variants avec MAF >= 0.05 et NS (nb échantillons génotypés) >= 349
  bcftools view -m2 -M2 -v snps -f PASS "$VCF" -Ou \
    | bcftools +fill-tags -Ou -- -t MAF,NS \
    | bcftools view -i 'MAF>=0.05 && NS>=349' -Oz -o "$OUT"

  bcftools index -t -f "$OUT"  # indexe le VCF filtré (index tabix)

  echo "$OUT" >> "$FILTERED_LIST"  # ajoute le chemin du VCF filtré à la liste
done < "$VCF_LIST"  # fin de boucle ; lit VCF_LIST en entrée

echo
echo "2) Concaténation des 25 chromosomes filtrés"

VCF_FILTERED_RAW="$D_FILTERED/awassi_basedata_25chr.filtered.raw_ids.vcf.gz"  # VCF concaténé, avant ajout des ID uniques
VCF_FILTERED="$D_FILTERED/awassi_basedata_25chr.filtered.PASS_biallelicSNP_maf005_miss010.vcf.gz"  # VCF filtré final (avec ID uniques)

bcftools concat -f "$FILTERED_LIST" -Oz -o "$VCF_FILTERED_RAW"  # concatène les VCF filtrés listés dans FILTERED_LIST
bcftools index -t -f "$VCF_FILTERED_RAW"  # indexe le VCF concaténé

echo
echo "3) Ajout d'identifiants uniques aux variants"

# ajoute un ID unique à chaque variant, au format CHROM:POS:REF:ALT
bcftools annotate \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  -Oz -o "$VCF_FILTERED" \
  "$VCF_FILTERED_RAW"

bcftools index -t -f "$VCF_FILTERED"  # indexe le VCF final

echo "VCF filtré final : $VCF_FILTERED"
echo

echo "4) Conversion PLINK 1.9"

# convertit le VCF filtré en PLINK binaire (bed/bim/fam) ;
# --double-id duplique l'ID échantillon en FID/IID, --allow-extra-chr autorise des noms de chromosomes non standards
plink \
  --vcf "$VCF_FILTERED" \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out "$D_PLINK/awassi_25chr_filtered"

echo
echo "5) Filtrage PLINK complémentaire"

# filtre PLINK supplémentaire : MAF >= 0.05, données manquantes <= 10% par SNP (geno) et par individu (mind)
plink \
  --bfile "$D_PLINK/awassi_25chr_filtered" \
  --allow-extra-chr \
  --maf 0.05 \
  --geno 0.10 \
  --mind 0.10 \
  --make-bed \
  --out "$D_PLINK/awassi_25chr_filtered_qc"

echo
echo "6) LD pruning avant PCA"

# LD pruning : fenêtre de 50 SNP, pas de 10, seuil r² 0.2 (retire les SNP trop corrélés entre eux)
plink \
  --bfile "$D_PLINK/awassi_25chr_filtered_qc" \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out "$D_PLINK/pruning_25chr"

echo
echo "7) PCA"

# calcule la PCA (20 composantes) sur les SNP indépendants issus du pruning
plink \
  --bfile "$D_PLINK/awassi_25chr_filtered_qc" \
  --allow-extra-chr \
  --extract "$D_PLINK/pruning_25chr.prune.in" \
  --pca 20 \
  --out "$D_PLINK/pca_25chr"

echo
echo "8) Plots PCA"

Rscript "$PROJECT_DIR/scripts/01_pca/02_plot_pca_25chr.R" "$D_RUN" "$METADATA_VCF_ORDER"  # génère les plots (script R)

echo
echo "=== PCA terminée ==="
echo "Résultats : $D_RUN"
echo "Plot principal : $D_PLOTS/PCA_PC1_PC2_geo_Awassi_highlight.pdf"
