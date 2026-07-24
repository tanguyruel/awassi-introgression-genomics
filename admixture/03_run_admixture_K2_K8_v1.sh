#!/usr/bin/env bash
# Lance ADMIXTURE en cross-validation pour K=2 a K=8 sur le jeu LD-prune
# genome-wide (359 individus, tous groupes geo). Choix du meilleur K via
# l'erreur de validation croisee (CV error), la plus basse = meilleur K.
#
# IMPORTANT : nécessite le binaire ADMIXTURE (Alexander et al. 2009), pas fourni
# dans ce dépôt. À télécharger depuis https://dalexander.github.io/admixture/
# puis adapter la variable ADMIXTURE ci-dessous vers l'exécutable téléchargé.
set -euo pipefail  # arrête le script en cas d'erreur ou de variable non définie
cd ~/Bureau/genome_complet_Awassi/analyses/admixture_v2  # se place dans le dossier de l'analyse ADMIXTURE

ADMIXTURE=/home/tanguyruel/Bureau/3718_SNP_Awassi/test_3718SNPs/software/admixture_linux-1.4.0/admixture  # chemin local à adapter (voir note ci-dessus)
BED="genomewide_pruned.bed"  # fichier PLINK bed en entrée (SNP LD-prunés genome-wide)
THREADS=6  # nombre de threads utilisés par ADMIXTURE

mkdir -p admixture_runs  # crée le dossier de sortie des runs
cd admixture_runs  # se place dans ce dossier
ln -sf ../genomewide_pruned.bed .  # lien symbolique vers le bed (ADMIXTURE le lit en local)
ln -sf ../genomewide_pruned.bim .  # idem pour le bim
ln -sf ../genomewide_pruned.fam .  # idem pour le fam

> cv_errors.txt  # crée/vide le fichier récapitulatif des erreurs de validation croisée
for K in 2 3 4 5 6 7 8; do  # boucle sur chaque valeur de K testée
  echo "=== K=${K} ==="
  # lance ADMIXTURE pour ce K, validation croisée 10-fold (--cv=10), log redirigé dans un fichier
  "$ADMIXTURE" --cv=10 -j${THREADS} genomewide_pruned.bed ${K} \
    > log_K${K}.out 2>&1
  cv=$(grep "CV error" log_K${K}.out | awk '{print $NF}')  # extrait la valeur du CV error dans le log
  echo -e "${K}\t${cv}" >> cv_errors.txt  # ajoute K et son CV error au fichier récapitulatif
  echo "  K=${K} CV error = ${cv}"
done

echo "TERMINE"
