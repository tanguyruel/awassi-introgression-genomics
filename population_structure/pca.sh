#!/bin/bash
#
# PCA de référence : 25 chromosomes, tous les individus (387), filtrage + LD pruning complets
# (pas de sous-échantillonnage géographique/spatial - voir PCA A/B/C pour les variantes rapides).
# Entrée : $VCF_LIST, $METADATA_VCF_ORDER (config_project.sh).
# Sortie : analyses/pca/PCA_<date>/ (VCF filtré, PLINK, plots via 02_plot_pca_25chr.R).
# Usage  : ./01_run_pca_25chr_filtered.sh
#
# IMPORTANT : chemin absolu en dur ligne ci-dessous, à adapter si le dépôt est cloné
# ailleurs. Dépend aussi de config_project.sh et scripts/01_pca/02_plot_pca_25chr.R,
# tous deux du dépôt de travail interne et NON inclus ici (script fourni à titre de
# référence de méthode, pas exécutable tel quel sans ces deux fichiers).

set -euo pipefail

cd ~/Bureau/genome_complet_Awassi  # chemin en dur à adapter
source config_project.sh  # charge les variables de configuration (VCF_LIST, PROJECT_DIR, etc.) — fichier non fourni, voir note ci-dessus

RUN_DATE=$(date +%d_%m_%Hh%M)
D_RUN="$PROJECT_DIR/analyses/pca/PCA_${RUN_DATE}"

D_FILTERED="$D_RUN/filtered_vcf"
D_PLINK="$D_RUN/plink"
D_PLOTS="$D_RUN/plots"
D_CHECKS="$D_RUN/checks"
D_RUN_LOGS="$D_RUN/logs"

mkdir -p "$D_FILTERED" "$D_PLINK" "$D_PLOTS" "$D_CHECKS" "$D_RUN_LOGS"

cp "$PROJECT_DIR/config_project.sh" "$D_RUN/config_project_used.sh"  # copie la config utilisée pour traçabilité

echo "=== PCA 25 chromosomes filtrés avec PLINK 1.9 ==="
echo "Run : $D_RUN"
echo "VCF list : $VCF_LIST"
echo "Metadata : $METADATA_VCF_ORDER"
echo

command -v bcftools >/dev/null || { echo "ERREUR : bcftools introuvable"; exit 1; }
command -v plink >/dev/null || { echo "ERREUR : plink introuvable"; exit 1; }
command -v Rscript >/dev/null || { echo "ERREUR : Rscript introuvable"; exit 1; }

plink --version  # traçabilité

N_VCF=$(wc -l < "$VCF_LIST")
if [ "$N_VCF" -ne 25 ]; then
  echo "ERREUR : la liste VCF ne contient pas 25 chromosomes : $N_VCF"
  exit 1
fi

echo
echo "1) Filtrage chromosome par chromosome"
echo "Filtres : SNP bialléliques, PASS, MAF >= 0.05, NS >= 349, donc missing variants <= 0.10"
echo

FILTERED_LIST="$D_FILTERED/list_filtered_vcf.txt"
> "$FILTERED_LIST"

while read -r VCF; do
  base=$(basename "$VCF")
  chr=$(echo "$base" | sed -E 's/.*chr([0-9]+)\.vcf\.gz/\1/')
  OUT="$D_FILTERED/awassi_basedata_chr${chr}.filtered.vcf.gz"

  echo "Filtrage chr${chr} : $base"

  # filtre SNP bialléliques PASS (-m2 -M2 -v snps -f PASS), calcule MAF et NS (+fill-tags),
  # puis ne garde que les variants avec MAF >= 0.05 et NS (nb échantillons génotypés) >= 349
  bcftools view -m2 -M2 -v snps -f PASS "$VCF" -Ou \
    | bcftools +fill-tags -Ou -- -t MAF,NS \
    | bcftools view -i 'MAF>=0.05 && NS>=349' -Oz -o "$OUT"

  bcftools index -t -f "$OUT"

  echo "$OUT" >> "$FILTERED_LIST"
done < "$VCF_LIST"

echo
echo "2) Concaténation des 25 chromosomes filtrés"

VCF_FILTERED_RAW="$D_FILTERED/awassi_basedata_25chr.filtered.raw_ids.vcf.gz"  # VCF concaténé, avant ajout des ID uniques
VCF_FILTERED="$D_FILTERED/awassi_basedata_25chr.filtered.PASS_biallelicSNP_maf005_miss010.vcf.gz"  # VCF filtré final (avec ID uniques)

bcftools concat -f "$FILTERED_LIST" -Oz -o "$VCF_FILTERED_RAW"
bcftools index -t -f "$VCF_FILTERED_RAW"

echo
echo "3) Ajout d'identifiants uniques aux variants"

bcftools annotate \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  -Oz -o "$VCF_FILTERED" \
  "$VCF_FILTERED_RAW"

bcftools index -t -f "$VCF_FILTERED"

echo "VCF filtré final : $VCF_FILTERED"
echo

echo "4) Conversion PLINK 1.9"

# --double-id duplique l'ID échantillon en FID/IID, --allow-extra-chr autorise des noms de chromosomes non standards
plink \
  --vcf "$VCF_FILTERED" \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out "$D_PLINK/awassi_25chr_filtered"

echo
echo "5) Filtrage PLINK complémentaire"

# MAF >= 0.05, données manquantes <= 10% par SNP (geno) et par individu (mind)
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

# fenêtre de 50 SNP, pas de 10, seuil r² 0.2 (retire les SNP trop corrélés entre eux)
plink \
  --bfile "$D_PLINK/awassi_25chr_filtered_qc" \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out "$D_PLINK/pruning_25chr"

echo
echo "7) PCA"

# sur les SNP indépendants issus du pruning
plink \
  --bfile "$D_PLINK/awassi_25chr_filtered_qc" \
  --allow-extra-chr \
  --extract "$D_PLINK/pruning_25chr.prune.in" \
  --pca 20 \
  --out "$D_PLINK/pca_25chr"

echo
echo "8) Plots PCA"

Rscript "$PROJECT_DIR/scripts/01_pca/02_plot_pca_25chr.R" "$D_RUN" "$METADATA_VCF_ORDER"

echo
echo "=== PCA terminée ==="
echo "Résultats : $D_RUN"
echo "Plot principal : $D_PLOTS/PCA_PC1_PC2_geo_Awassi_highlight.pdf"
