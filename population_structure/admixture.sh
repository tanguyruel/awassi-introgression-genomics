#!/usr/bin/env bash
# Lance ADMIXTURE en cross-validation pour K=2 a K=8 sur le jeu LD-prune
# genome-wide (359 individus, tous groupes geo). Choix du meilleur K via
# l'erreur de validation croisee (CV error), la plus basse = meilleur K.
#
# IMPORTANT : nécessite le binaire ADMIXTURE (Alexander et al. 2009), pas fourni
# dans ce dépôt. À télécharger depuis https://dalexander.github.io/admixture/,
# puis soit le mettre dans le PATH, soit régler la variable ADMIXTURE_BIN.
# Usage : AWASSI_PROJECT_DIR=/chemin/vers/le/projet ADMIXTURE_BIN=/chemin/vers/admixture ./run_admixture_K2_K8_v1.sh
set -euo pipefail
cd "${AWASSI_PROJECT_DIR:-$PWD}/analyses/admixture_v2"

ADMIXTURE="${ADMIXTURE_BIN:-admixture}"  # binaire ADMIXTURE : $ADMIXTURE_BIN, sinon cherché dans le PATH (voir note ci-dessus)
BED="genomewide_pruned.bed"
THREADS=6

mkdir -p admixture_runs
cd admixture_runs
ln -sf ../genomewide_pruned.bed .  # ADMIXTURE lit ses fichiers en local, d'où les liens symboliques
ln -sf ../genomewide_pruned.bim .
ln -sf ../genomewide_pruned.fam .

> cv_errors.txt
for K in 2 3 4 5 6 7 8; do
  echo "=== K=${K} ==="
  # --cv=10 : validation croisée 10-fold
  "$ADMIXTURE" --cv=10 -j${THREADS} genomewide_pruned.bed ${K} \
    > log_K${K}.out 2>&1
  cv=$(grep "CV error" log_K${K}.out | awk '{print $NF}')
  echo -e "${K}\t${cv}" >> cv_errors.txt
  echo "  K=${K} CV error = ${cv}"
done

echo "TERMINE"
