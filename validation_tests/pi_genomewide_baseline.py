#!/usr/bin/env python3
"""
Pi genome-wide "de fond" par groupe géo (7 groupes), pour servir de référence à
laquelle comparer le pi des 9 régions candidates (voir pi_regions_qc.py).

Méthode : plutôt que traiter les ~100M SNP du génome entier (VCF filtrés
~100 Go au total), échantillon aléatoire de 600 fenêtres de 20kb réparties
proportionnellement à la longueur de chaque chromosome (~12 Mb échantillonnés,
~0.5% du génome) — estimation non biaisée de la moyenne genome-wide (loi des
grands nombres), sans scanner tout le génome. Les 9 régions candidates (+
marge 200kb) sont exclues du tirage pour ne pas biaiser la référence de fond.

Même formule pi_site que pi_regions_qc.py : (2n/(2n-1))*2p(1-p), sommée puis divisée
par le nombre total de bp échantillonnés (convention vcftools --window-pi).

Source VCF : analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/filtered_vcf/
  (déjà filtrés PASS + biallélique + SNP only, réutilisés tels quels — pas de
  données brutes touchées)

Sortie : analyses/synthese_resultats/pi_genomewide_baseline/
  - windows_sampled.bed (fenêtres tirées, pour traçabilité/reproductibilité)
  - Pi_genomewide_baseline_by_group.tsv
  - pi_baseline_vs_regions.png/.pdf (comparaison avec les résultats de pi_regions_qc.py)

Script de référence pour "π — fond génomique" (cf. Annexe A du rapport de stage),
complémentaire du pi_regions_qc.py (π dans les 9 régions elles-mêmes).
Usage : python3 50_pi_genomewide_baseline_v1.py
"""

import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # racine du dépôt, pour importer _shared
from _shared import gt_to_alt_count, pi_site

random.seed(42)  # graine fixe pour un tirage reproductible des fenêtres

VCF_DIR = Path("analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/filtered_vcf")  # VCF déjà filtrés, par chromosome
POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # popmaps par groupe géo
OUTDIR = Path("analyses/synthese_resultats/pi_genomewide_baseline")
OUTDIR.mkdir(parents=True, exist_ok=True)

GROUPS = ["Awassi", "MiddleEastNonAwassi", "Africa", "Asia", "Europe", "America", "Australia"]
WINDOW = 20_000  # taille d'une fenêtre échantillonnée (bp)
N_WINDOWS_TOTAL = 600  # nombre total de fenêtres tirées sur tout le génome
EXCLUDE_MARGIN = 200_000  # marge autour de chaque région candidate à exclure du tirage

# Longueur (bp) de chaque chromosome, utilisée pour répartir les fenêtres proportionnellement
CHROM_LENGTHS = {
    1: 275406953, 2: 248966461, 3: 223996068, 4: 119216639, 5: 107836144,
    6: 116888256, 7: 100009711, 8: 90615088, 9: 94583238, 10: 86377204,
    11: 62170480, 12: 79028859, 13: 82951069, 14: 62568341, 15: 80783214,
    16: 71693149, 17: 72251135, 18: 68494538, 19: 60445663, 20: 51049468,
    21: 49987992, 22: 50780147, 23: 62282865, 24: 41976827, 25: 45223504,
}

# 9 régions candidates (à exclure + marge, pour une référence de fond non biaisée)
CANDIDATE_REGIONS = [
    ("2", 112785001, 112865000), ("3", 129220001, 129260000), ("3", 137390001, 137420000),
    ("3", 186075001, 186125000), ("5", 58070001, 58110000), ("6", 70240001, 70295000),
    ("10", 49360001, 49390000), ("17", 34160001, 34220000), ("20", 765001, 845000),
]
EXCLUDE_BY_CHROM = defaultdict(list)
for c, s, e in CANDIDATE_REGIONS:
    EXCLUDE_BY_CHROM[c].append((s - EXCLUDE_MARGIN, e + EXCLUDE_MARGIN))


def overlaps_excluded(chrom, start, end):
    """True si [start, end] chevauche une zone exclue (région candidate + marge) du chromosome `chrom`."""
    for es, ee in EXCLUDE_BY_CHROM.get(chrom, []):
        if start <= ee and es <= end:
            return True
    return False


# ── 1. Échantillonnage des fenêtres, proportionnel à la longueur des chromosomes ──
total_len = sum(CHROM_LENGTHS.values())
windows_by_chrom = {}
n_allocated = 0  # nombre de fenêtres déjà allouées (pour ajuster le dernier chromosome)
chrom_list = sorted(CHROM_LENGTHS)
for i, c in enumerate(chrom_list):
    if i == len(chrom_list) - 1:
        n = N_WINDOWS_TOTAL - n_allocated  # dernier chromosome : absorbe l'arrondi restant
    else:
        n = round(N_WINDOWS_TOTAL * CHROM_LENGTHS[c] / total_len)
    n_allocated += n
    windows_by_chrom[str(c)] = n

bed_rows = []  # liste des fenêtres tirées (chrom, start, end)
for c, n in windows_by_chrom.items():
    clen = CHROM_LENGTHS[int(c)]
    chosen = []
    attempts = 0  # compteur de tentatives (évite boucle infinie)
    while len(chosen) < n and attempts < n * 50:
        attempts += 1
        start = random.randint(1, clen - WINDOW)
        end = start + WINDOW - 1
        if overlaps_excluded(c, start, end):
            continue
        if any(start <= e2 and s2 <= end for s2, e2 in chosen):
            continue
        chosen.append((start, end))
    for s, e in sorted(chosen):
        bed_rows.append((c, s, e))

bed_path = OUTDIR / "windows_sampled.bed"  # fichier BED de traçabilité des fenêtres tirées
with open(bed_path, "w") as f:
    for c, s, e in bed_rows:
        f.write(f"{c}\t{s}\t{e}\n")

n_total_sampled = len(bed_rows)
bp_sampled = n_total_sampled * WINDOW
print(f"Fenêtres échantillonnées : {n_total_sampled} ({bp_sampled/1e6:.1f} Mb, {100*bp_sampled/total_len:.3f}% du génome)")
print(f"→ {bed_path}")

# ── 2. Pi pooled genome-wide, un chromosome à la fois (bcftools par chromosome) ──
groups_samples = {g: [x.strip() for x in open(POP_DIR / f"{g}.txt") if x.strip()] for g in GROUPS}
all_samples = sorted(set().union(*groups_samples.values()))


sum_pi = defaultdict(float)  # somme des pi de site, par groupe, sur toutes les fenêtres échantillonnées
n_sites_used = defaultdict(int)  # nombre de sites utilisés (pi défini), par groupe
n_sites_read = 0

chroms_with_windows = sorted({c for c, _, _ in bed_rows}, key=int)
for c in chroms_with_windows:  # traite un chromosome à la fois pour limiter la mémoire
    vcf = VCF_DIR / f"chr{c}.PASS_biallelic_snps.vcf.gz"
    sub_bed = OUTDIR / f"_tmp_chr{c}.bed"
    with open(sub_bed, "w") as f:
        for cc, s, e in bed_rows:
            if cc == c:
                f.write(f"{cc}\t{s}\t{e}\n")

    order = subprocess.check_output(
        ["bash", "-lc", f"bcftools view -R {sub_bed} -f PASS -m2 -M2 -v snps {vcf} | bcftools query -l"],
        text=True,
    ).splitlines()  # ordre des échantillons dans le VCF restreint aux fenêtres (région, PASS, biallélique, SNP)
    idx_map = {s: i for i, s in enumerate(order)}
    group_idx = {g: np.array([idx_map[s] for s in groups_samples[g] if s in idx_map], dtype=int) for g in GROUPS}

    cmd = f"bcftools view -R {sub_bed} -f PASS -m2 -M2 -v snps {vcf} | bcftools query -f '%CHROM[\\t%GT]\\n'"  # extrait CHROM + GT de chaque échantillon
    txt = subprocess.check_output(["bash", "-lc", cmd], text=True)
    sub_bed.unlink()

    n_lines_this_chr = 0
    for line in txt.splitlines():
        parts = line.rstrip("\n").split("\t")
        gts = parts[1:]
        alt_counts = np.array([gt_to_alt_count(gt) for gt in gts], dtype=float)
        for g in GROUPS:
            pi = pi_site(alt_counts, group_idx[g])
            if np.isfinite(pi):
                sum_pi[g] += pi
                n_sites_used[g] += 1
        n_lines_this_chr += 1
    n_sites_read += n_lines_this_chr
    print(f"  chr{c}: {windows_by_chrom[c]} fenêtres, {n_lines_this_chr} SNP lus")

print(f"\nTotal SNP lus (toutes fenêtres échantillonnées) : {n_sites_read}")

rows = []
for g in GROUPS:
    pi_per_bp = sum_pi[g] / bp_sampled  # pi moyen par paire de bases (convention vcftools --window-pi)
    rows.append({
        "group": g, "n_individuals": len(groups_samples[g]),
        "n_sites_used": n_sites_used[g], "bp_sampled": bp_sampled,
        "sum_pi_site": sum_pi[g], "pi_per_bp": pi_per_bp, "pi_per_kb": pi_per_bp * 1000,
    })
out = pd.DataFrame(rows).sort_values("pi_per_kb", ascending=False)
out_path = OUTDIR / "Pi_genomewide_baseline_by_group.tsv"
out.to_csv(out_path, sep="\t", index=False)
print(f"\n→ {out_path}")
print(out.to_string(index=False))
