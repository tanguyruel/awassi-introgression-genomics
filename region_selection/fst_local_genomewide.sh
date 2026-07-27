#!/usr/bin/env bash
#
# Rôle : version la plus récente du FST local genome-wide Awassi vs groupes,
#        fenêtres 20kb / step 5kb (plus fin qu'une version antérieure à 20kb simple), Europe/America/
#        Australia séparés. Recrée les popfiles si absents (depuis les groupes
#        utilisés en analyse fd) et écrit une table brute + une table filtrée
#        par nombre minimal de SNP (MIN_SNPS).
# Entrée : VCF par chromosome, GROUP_SOURCE (analyses/fd/.../chr10_groups_used.tsv)
# Sortie : analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/merged/
#          FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_{RAW,MIN<n>SNP}.tsv
# Usage  : MIN_SNPS=10 MAX_JOBS=6 ./21_run_FST_local_20kb_step5kb_separate_EU_AM_AUS.sh
#
set -euo pipefail

cd ~/Bureau/genome_complet_Awassi

OUT="analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1"
SAMPLES="analyses/fst/popmaps_separees_v1"
GROUP_SOURCE="analyses/fd/fd_genomewide_clean7_26_06_11h21/results/chr10_groups_used.tsv"  # source des groupes (produite par l'analyse fd)

mkdir -p "$OUT"/{filtered_vcf,results,merged,logs}
mkdir -p "$SAMPLES"

WINDOW=20000
STEP=5000
MIN_SNPS="${MIN_SNPS:-10}"  # paramétrable via l'environnement (voir Usage en tête de fichier)
MAX_JOBS="${MAX_JOBS:-6}"  # paramétrable via l'environnement (voir Usage en tête de fichier)
GROUP_LIST="Africa Asia Europe America Australia MiddleEastNonAwassi"

echo "============================================================"
echo "FST local 20kb step5kb séparé Europe / America / Australia"
echo "OUT      : $OUT"
echo "WINDOW   : $WINDOW"
echo "STEP     : $STEP"
echo "MIN_SNPS : $MIN_SNPS pour fichiers filtrés/ranking"
echo "MAX_JOBS : $MAX_JOBS"
echo "GROUPS   : $GROUP_LIST"
echo "============================================================"

command -v bcftools >/dev/null || { echo "ERREUR : bcftools absent"; exit 1; }
command -v vcftools >/dev/null || { echo "ERREUR : vcftools absent"; exit 1; }

# ------------------------------------------------------------
# 0. Créer les fichiers samples séparés si absents
# ------------------------------------------------------------
if [[ ! -s "$SAMPLES/Awassi.txt" || ! -s "$SAMPLES/Europe.txt" || ! -s "$SAMPLES/America.txt" || ! -s "$SAMPLES/Australia.txt" ]]; then
  echo "Création des fichiers samples depuis : $GROUP_SOURCE"

  if [[ ! -s "$GROUP_SOURCE" ]]; then
    echo "ERREUR : source groupes absente : $GROUP_SOURCE"
    exit 1
  fi

  python3 - <<'PY'
from pathlib import Path
import pandas as pd

src = Path("analyses/fd/fd_genomewide_clean7_26_06_11h21/results/chr10_groups_used.tsv")
outdir = Path("analyses/fst/popmaps_separees_v1")
outdir.mkdir(parents=True, exist_ok=True)

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

df = pd.read_csv(src, sep="\t")
rows = []

for _, r in df.iterrows():
    g = str(r["group"])
    if g not in keep:
        continue

    samples = str(r["samples"]).split(",")
    with open(outdir / f"{g}.txt", "w") as f:
        for s in samples:
            s = s.strip()
            if s:
                f.write(s + "\n")
                rows.append({"sample": s, "group": g})

popmap = pd.DataFrame(rows)
popmap.to_csv(
    outdir / "popmap_Awassi_ME_Africa_Asia_Europe_America_Australia.tsv",
    sep="\t",
    index=False,
    header=False
)

print(popmap["group"].value_counts().reindex(keep))  # ordre d'affichage = keep
PY
fi

for g in Awassi MiddleEastNonAwassi Africa Asia Europe America Australia; do
  if [[ ! -s "$SAMPLES/${g}.txt" ]]; then
    echo "ERREUR : fichier samples absent ou vide : $SAMPLES/${g}.txt"
    exit 1
  fi
done

echo
echo "===== Effectifs samples ====="
for g in Awassi MiddleEastNonAwassi Africa Asia Europe America Australia; do
  echo -e "$g\t$(wc -l < "$SAMPLES/${g}.txt")"
done

# ------------------------------------------------------------
# 1. Trouver VCF
# ------------------------------------------------------------
find_vcf() {
  local chr="$1"

  local candidates=(
    "data/raw data_08_06/awassi_and_basedata_chr${chr}.vcf.gz"
    "data/raw_data_08_06/awassi_and_basedata_chr${chr}.vcf.gz"
    "/home/deschaoc/AWASSI/awassi_and_basedata_chr${chr}.vcf.gz"
    "data/awassi_and_basedata_chr${chr}.vcf.gz"
  )

  for f in "${candidates[@]}"; do
    if [[ -f "$f" ]]; then
      echo "$f"
      return 0
    fi
  done

  return 1
}

# ------------------------------------------------------------
# 2. Run FST par chromosome
# ------------------------------------------------------------
run_one_chr() {  # traite un chromosome : filtrage VCF puis FST vs chaque groupe de GROUP_LIST
  local chr="$1"

  local raw
  raw=$(find_vcf "$chr") || {
    echo "ERREUR chr${chr} : VCF introuvable" >&2
    return 1
  }

  local filt="$OUT/filtered_vcf/chr${chr}.PASS_biallelic_snps.vcf.gz"

  if [[ ! -s "$filt" ]]; then
    echo "[chr${chr}] Filtrage PASS SNP bialléliques"
    bcftools view -f PASS -m2 -M2 -v snps "$raw" -Oz -o "$filt"
    bcftools index -t "$filt"
  else
    echo "[chr${chr}] VCF filtré déjà présent"
  fi

  for g in $GROUP_LIST; do
    local pop1="$SAMPLES/Awassi.txt"
    local pop2="$SAMPLES/${g}.txt"
    local prefix="$OUT/results/chr${chr}_Awassi_vs_${g}_20kb_step5kb"

    if [[ -s "${prefix}.windowed.weir.fst" ]]; then
      echo "[chr${chr}] Awassi vs $g déjà fait"
      continue
    fi

    echo "[chr${chr}] FST Awassi vs $g"

    # FST Weir & Cockerham par fenêtre glissante ; --weir-fst-pop (x2) définit les deux populations comparées,
    # --fst-window-size/--fst-window-step définissent la fenêtre glissante
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

export -f find_vcf  # exporte les fonctions et variables : nécessaire pour qu'elles soient visibles dans les sous-shells lancés par xargs
export -f run_one_chr
export OUT SAMPLES WINDOW STEP GROUP_LIST

seq 1 25 | xargs -I{} -P "$MAX_JOBS" bash -c 'run_one_chr "$@"' _ {}  # lance run_one_chr en parallèle (MAX_JOBS jobs) pour les 25 chromosomes

# ------------------------------------------------------------
# 3. Fusion + version filtrée MIN_SNPS
# ------------------------------------------------------------
echo
echo "===== Fusion des résultats ====="

python3 - <<PY
# (heredoc sans quotes : ${MIN_SNPS} ci-dessous est substitué par bash avant exécution)
from pathlib import Path
import pandas as pd
import re
import numpy as np

OUT = Path("analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1")
MIN_SNPS = int("${MIN_SNPS}")  # valeur substituée depuis bash

resdir = OUT / "results"
merged = OUT / "merged"
merged.mkdir(parents=True, exist_ok=True)

rows = []

for f in sorted(resdir.glob("chr*_Awassi_vs_*_20kb_step5kb.windowed.weir.fst")):
    m = re.search(r"chr(\d+)_Awassi_vs_(.+)_20kb_step5kb\.windowed\.weir\.fst$", f.name)
    if not m:
        continue

    chrom = int(m.group(1))
    group2 = m.group(2)

    try:
        df = pd.read_csv(f, sep="\t")
    except Exception as e:
        print("Lecture impossible :", f, e)
        continue

    if df.empty:
        continue

    for _, r in df.iterrows():
        wf = pd.to_numeric(r.get("WEIGHTED_FST", np.nan), errors="coerce")
        mf = pd.to_numeric(r.get("MEAN_FST", np.nan), errors="coerce")
        nv = pd.to_numeric(r.get("N_VARIANTS", np.nan), errors="coerce")

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

out = pd.DataFrame(rows)

if out.empty:
    raise SystemExit("Aucun résultat FST fusionné.")

out = out.sort_values(["chr", "start", "end", "group2"])

raw_file = merged / "FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_RAW.tsv"
filtered_file = merged / f"FST_genomewide_20kb_step5kb_long_sep_EU_AM_AUS_MIN{MIN_SNPS}SNP.tsv"

out.to_csv(raw_file, sep="\t", index=False)

# ne garde que les fenêtres avec assez de SNP et un FST pondéré valide
out_f = out[
    (out["n_variants"] >= MIN_SNPS) &
    (out["weighted_fst"].notna())
].copy()

out_f.to_csv(filtered_file, sep="\t", index=False)

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
