#!/usr/bin/env bash
#
# Rôle : version la plus récente du FST local genome-wide Awassi vs groupes,
#        fenêtres 20kb / step 5kb (plus fin que le script 20), Europe/America/
#        Australia séparés. Recrée les popfiles si absents (depuis les groupes
#        utilisés en analyse fd) et écrit une table brute + une table filtrée
#        par nombre minimal de SNP (MIN_SNPS).
# Entrée : VCF par chromosome, GROUP_SOURCE (analyses/fd/.../chr10_groups_used.tsv)
# Sortie : analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/merged/
#          FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_{RAW,MIN<n>SNP}.tsv
# Usage  : MIN_SNPS=10 MAX_JOBS=6 ./21_run_FST_local_20kb_step5kb_separate_EU_AM_AUS.sh
#
set -euo pipefail  # arrête le script en cas d'erreur, de variable non définie ou d'échec dans un pipe

cd ~/Bureau/genome_complet_Awassi  # se place à la racine du projet

OUT="analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1"  # dossier racine de sortie
SAMPLES="analyses/fst/popmaps_separees_v1"  # dossier des fichiers samples par groupe
GROUP_SOURCE="analyses/fd/fd_genomewide_clean7_26_06_11h21/results/chr10_groups_used.tsv"  # source des groupes (produite par l'analyse fd)

mkdir -p "$OUT"/{filtered_vcf,results,merged,logs}  # crée les sous-dossiers de sortie
mkdir -p "$SAMPLES"  # crée le dossier des fichiers samples

WINDOW=20000  # taille de fenêtre (pb)
STEP=5000  # pas entre fenêtres (pb)
MIN_SNPS="${MIN_SNPS:-10}"  # seuil minimal de SNP par fenêtre pour la table filtrée (paramétrable via l'environnement)
MAX_JOBS="${MAX_JOBS:-6}"  # nombre de chromosomes traités en parallèle (paramétrable via l'environnement)
GROUP_LIST="Africa Asia Europe America Australia MiddleEastNonAwassi"  # groupes comparés à Awassi

echo "============================================================"
echo "FST local 20kb step5kb séparé Europe / America / Australia"
echo "OUT      : $OUT"
echo "WINDOW   : $WINDOW"
echo "STEP     : $STEP"
echo "MIN_SNPS : $MIN_SNPS pour fichiers filtrés/ranking"
echo "MAX_JOBS : $MAX_JOBS"
echo "GROUPS   : $GROUP_LIST"
echo "============================================================"

command -v bcftools >/dev/null || { echo "ERREUR : bcftools absent"; exit 1; }  # vérifie que bcftools est installé
command -v vcftools >/dev/null || { echo "ERREUR : vcftools absent"; exit 1; }  # vérifie que vcftools est installé

# ------------------------------------------------------------
# 0. Créer les fichiers samples séparés si absents
# ------------------------------------------------------------
if [[ ! -s "$SAMPLES/Awassi.txt" || ! -s "$SAMPLES/Europe.txt" || ! -s "$SAMPLES/America.txt" || ! -s "$SAMPLES/Australia.txt" ]]; then  # si un des fichiers samples clés manque ou est vide
  echo "Création des fichiers samples depuis : $GROUP_SOURCE"

  if [[ ! -s "$GROUP_SOURCE" ]]; then  # vérifie que le fichier source des groupes existe
    echo "ERREUR : source groupes absente : $GROUP_SOURCE"
    exit 1
  fi

  python3 - <<'PY'
# imports : chemins (Path) et manipulation de tables (pandas)
from pathlib import Path
import pandas as pd

src = Path("analyses/fd/fd_genomewide_clean7_26_06_11h21/results/chr10_groups_used.tsv")  # table source des groupes/échantillons
outdir = Path("analyses/fst/popmaps_separees_v1")  # dossier de sortie des fichiers samples par groupe
outdir.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si besoin

# liste des groupes à conserver (Awassi + groupes de comparaison)
keep = [
    "Awassi",
    "MiddleEastNonAwassi",
    "Africa",
    "Asia",
    "Europe",
    "America",
    "Australia",
]

df = pd.read_csv(src, sep="\t")  # charge la table des groupes/échantillons
rows = []  # accumulateur pour construire la popmap globale

for _, r in df.iterrows():  # parcourt chaque groupe de la table
    g = str(r["group"])  # nom du groupe
    if g not in keep:  # ignore les groupes non retenus
        continue

    samples = str(r["samples"]).split(",")  # liste des échantillons du groupe (séparés par virgule)
    with open(outdir / f"{g}.txt", "w") as f:  # écrit un fichier .txt (un échantillon par ligne) pour ce groupe
        for s in samples:  # parcourt chaque échantillon
            s = s.strip()  # retire les espaces superflus
            if s:  # ignore les entrées vides
                f.write(s + "\n")  # écrit l'échantillon dans le fichier du groupe
                rows.append({"sample": s, "group": g})  # mémorise l'association échantillon/groupe

popmap = pd.DataFrame(rows)  # construit la popmap globale (échantillon + groupe)
# écrit la popmap globale en TSV, sans en-tête
popmap.to_csv(
    outdir / "popmap_Awassi_ME_Africa_Asia_Europe_America_Australia.tsv",
    sep="\t",
    index=False,
    header=False
)

print(popmap["group"].value_counts().reindex(keep))  # affiche l'effectif par groupe (ordre = keep)
PY
fi

for g in Awassi MiddleEastNonAwassi Africa Asia Europe America Australia; do  # vérifie que chaque fichier samples existe et n'est pas vide
  if [[ ! -s "$SAMPLES/${g}.txt" ]]; then
    echo "ERREUR : fichier samples absent ou vide : $SAMPLES/${g}.txt"
    exit 1
  fi
done

echo
echo "===== Effectifs samples ====="
for g in Awassi MiddleEastNonAwassi Africa Asia Europe America Australia; do  # affiche le nombre d'échantillons par groupe
  echo -e "$g\t$(wc -l < "$SAMPLES/${g}.txt")"
done

# ------------------------------------------------------------
# 1. Trouver VCF
# ------------------------------------------------------------
find_vcf() {  # fonction : cherche le VCF d'un chromosome parmi plusieurs emplacements possibles
  local chr="$1"  # numéro de chromosome passé en argument

  # liste des chemins possibles où chercher le VCF de ce chromosome
  local candidates=(
    "data/raw data_08_06/awassi_and_basedata_chr${chr}.vcf.gz"
    "data/raw_data_08_06/awassi_and_basedata_chr${chr}.vcf.gz"
    "/home/deschaoc/AWASSI/awassi_and_basedata_chr${chr}.vcf.gz"
    "data/awassi_and_basedata_chr${chr}.vcf.gz"
  )

  for f in "${candidates[@]}"; do  # teste chaque chemin candidat
    if [[ -f "$f" ]]; then  # si le fichier existe
      echo "$f"  # renvoie son chemin
      return 0  # trouvé : sort avec succès
    fi
  done

  return 1  # aucun VCF trouvé : sort en échec
}

# ------------------------------------------------------------
# 2. Run FST par chromosome
# ------------------------------------------------------------
run_one_chr() {  # fonction : traite un chromosome (filtrage puis FST vs chaque groupe)
  local chr="$1"  # numéro de chromosome passé en argument

  local raw  # contiendra le chemin du VCF brut trouvé
  # cherche le VCF du chromosome, sort en erreur si introuvable
  raw=$(find_vcf "$chr") || {
    echo "ERREUR chr${chr} : VCF introuvable" >&2
    return 1
  }

  local filt="$OUT/filtered_vcf/chr${chr}.PASS_biallelic_snps.vcf.gz"  # chemin du VCF filtré pour ce chromosome

  if [[ ! -s "$filt" ]]; then  # si le VCF filtré n'existe pas déjà
    echo "[chr${chr}] Filtrage PASS SNP bialléliques"
    bcftools view -f PASS -m2 -M2 -v snps "$raw" -Oz -o "$filt"  # ne garde que les SNP bialléliques PASS
    bcftools index -t "$filt"  # indexe le VCF filtré
  else
    echo "[chr${chr}] VCF filtré déjà présent"
  fi

  for g in $GROUP_LIST; do  # boucle sur chaque groupe comparé à Awassi
    local pop1="$SAMPLES/Awassi.txt"  # fichier samples Awassi (population 1)
    local pop2="$SAMPLES/${g}.txt"  # fichier samples du groupe comparé (population 2)
    local prefix="$OUT/results/chr${chr}_Awassi_vs_${g}_20kb_step5kb"  # préfixe des fichiers de sortie vcftools

    if [[ -s "${prefix}.windowed.weir.fst" ]]; then  # si le résultat existe déjà, on saute ce groupe
      echo "[chr${chr}] Awassi vs $g déjà fait"
      continue
    fi

    echo "[chr${chr}] FST Awassi vs $g"

    # calcule le FST (Weir & Cockerham) par fenêtre glissante entre Awassi et le groupe g ;
    # --weir-fst-pop (x2) définit les deux populations comparées, --fst-window-size/--fst-window-step définissent la fenêtre glissante
    vcftools \
      --gzvcf "$filt" \
      --weir-fst-pop "$pop1" \
      --weir-fst-pop "$pop2" \
      --fst-window-size "$WINDOW" \
      --fst-window-step "$STEP" \
      --out "$prefix" \
      > "$OUT/logs/chr${chr}_Awassi_vs_${g}.vcftools.log" 2>&1
  done
}

export -f find_vcf  # exporte la fonction pour les sous-shells lancés par xargs
export -f run_one_chr  # idem pour run_one_chr
export OUT SAMPLES WINDOW STEP GROUP_LIST  # exporte les variables utilisées par les fonctions dans les sous-shells

seq 1 25 | xargs -I{} -P "$MAX_JOBS" bash -c 'run_one_chr "$@"' _ {}  # lance run_one_chr en parallèle (MAX_JOBS jobs) pour les 25 chromosomes

# ------------------------------------------------------------
# 3. Fusion + version filtrée MIN_SNPS
# ------------------------------------------------------------
echo
echo "===== Fusion des résultats ====="

python3 - <<PY
# imports : chemins, tables, expressions régulières, calcul numérique
# (heredoc sans quotes : ${MIN_SNPS} ci-dessous est substitué par bash avant exécution)
from pathlib import Path
import pandas as pd
import re
import numpy as np

OUT = Path("analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1")  # dossier racine des résultats FST
MIN_SNPS = int("${MIN_SNPS}")  # seuil minimal de SNP par fenêtre, valeur substituée depuis bash

resdir = OUT / "results"  # dossier des résultats par chromosome/groupe
merged = OUT / "merged"  # dossier de sortie fusionné
merged.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si besoin

rows = []  # accumulateur des lignes FST fusionnées

for f in sorted(resdir.glob("chr*_Awassi_vs_*_20kb_step5kb.windowed.weir.fst")):  # parcourt tous les fichiers de résultats vcftools
    m = re.search(r"chr(\d+)_Awassi_vs_(.+)_20kb_step5kb\.windowed\.weir\.fst$", f.name)  # extrait chromosome et groupe depuis le nom du fichier
    if not m:  # ignore les fichiers qui ne correspondent pas au motif attendu
        continue

    chrom = int(m.group(1))  # numéro de chromosome
    group2 = m.group(2)  # nom du groupe comparé à Awassi

    try:
        df = pd.read_csv(f, sep="\t")  # charge le fichier de résultats FST de vcftools
    except Exception as e:
        print("Lecture impossible :", f, e)
        continue

    if df.empty:  # ignore les fichiers vides
        continue

    for _, r in df.iterrows():  # parcourt chaque fenêtre du fichier
        wf = pd.to_numeric(r.get("WEIGHTED_FST", np.nan), errors="coerce")  # FST pondéré de la fenêtre (NaN si absent/invalide)
        mf = pd.to_numeric(r.get("MEAN_FST", np.nan), errors="coerce")  # FST moyen de la fenêtre
        nv = pd.to_numeric(r.get("N_VARIANTS", np.nan), errors="coerce")  # nombre de variants dans la fenêtre

        # enregistre une ligne standardisée (position, groupes, FST, nb SNP) pour cette fenêtre
        rows.append({
            "chr": chrom,
            "start": int(r["BIN_START"]),
            "end": int(r["BIN_END"]),
            "group1": "Awassi",
            "group2": group2,
            "n_variants": int(nv) if pd.notna(nv) else 0,
            "weighted_fst": wf,
            "mean_fst": mf,
        })

out = pd.DataFrame(rows)  # assemble toutes les lignes en une table

if out.empty:  # arrête si aucun résultat n'a été fusionné
    raise SystemExit("Aucun résultat FST fusionné.")

out = out.sort_values(["chr", "start", "end", "group2"])  # trie par position génomique puis par groupe

raw_file = merged / "FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_RAW.tsv"  # chemin de la table brute (non filtrée)
filtered_file = merged / f"FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_MIN{MIN_SNPS}SNP.tsv"  # chemin de la table filtrée par MIN_SNPS

out.to_csv(raw_file, sep="\t", index=False)  # écrit la table brute complète

# ne garde que les fenêtres avec assez de SNP et un FST pondéré valide
out_f = out[
    (out["n_variants"] >= MIN_SNPS) &
    (out["weighted_fst"].notna())
].copy()

out_f.to_csv(filtered_file, sep="\t", index=False)  # écrit la table filtrée

print("RAW écrit      :", raw_file)
print("RAW lignes     :", len(out))
print()
print("Filtré écrit   :", filtered_file)
print("MIN_SNPS       :", MIN_SNPS)
print("Filtré lignes  :", len(out_f))
print()
print("Lignes filtrées par groupe :")
print(out_f.groupby("group2").size().sort_index())
print()
print("Résumé n_variants par groupe :")
print(out.groupby("group2")["n_variants"].describe()[["count","mean","min","25%","50%","75%","max"]])
PY

echo
echo "============================================================"
echo "TERMINÉ"
echo "Fichiers finaux :"
echo "$OUT/merged/FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_RAW.tsv"
echo "$OUT/merged/FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_MIN${MIN_SNPS}SNP.tsv"
echo "============================================================"
