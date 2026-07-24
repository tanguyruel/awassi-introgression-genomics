#!/usr/bin/env python3
"""
Hétérozygotie observée (Ho) et coefficient de
consanguinité (F = 1 - Ho/He, He basé sur LES PROPRES fréquences alléliques
du sous-groupe testé, jamais les fréquences globales — évite l'effet Wahlund) :
  - par zone géo (7) et par race (35), sur l'échantillon genome-wide de fond
    (mêmes 600 fenêtres que pi_genomewide_baseline.py/53)
  - par zone géo et par race, sur chacune des 9 régions candidates

v1 (vcftools --het --keep, 420 appels) beaucoup trop lent (>10 min, chaque
appel reparcourt tout le fichier) — remplacé par un calcul direct en une
seule lecture par fichier VCF (comme scripts 50/53/54 pour pi), agrégé
ensuite pour les 42 groupes à la fois.

F par individu = 1 - O(HOM)/E(HOM), E(HOM) = somme sur les sites appelés de
(p²+(1-p)²) où p = fréquence allélique du GROUPE à ce site. F du groupe =
moyenne des F individuels (même convention que vcftools --het).

Sources VCF (aucune donnée brute retouchée) :
  - analyses/synthese_resultats/pi_genomewide_baseline/sampled_windows_merged.vcf.gz
  - analyses/synthese_resultats/annotation_9regions/vcf/{region_id}.vcf.gz

Sorties : analyses/synthese_resultats/het_inbreeding/
  - Het_F_genomewide_by_zone.tsv, Het_F_genomewide_by_race.tsv
  - Het_F_9regions_by_zone.tsv, Het_F_9regions_by_race.tsv

Script de référence pour "F et Ho" (cf. Annexe A du rapport de stage).
Usage : python3 56_het_inbreeding_v1.py
(nécessite les sorties de 50_pi_genomewide_baseline_v1.py et 40_annotate_variants_9regions_v1.py)
"""

# subprocess : appelle bcftools en ligne de commande ; pathlib : gestion des chemins
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("analyses/synthese_resultats")
GENOMEWIDE_VCF = BASE / "pi_genomewide_baseline" / "sampled_windows_merged.vcf.gz"  # échantillon de fond (pi_genomewide_baseline.py)
REGION_VCF_DIR = BASE / "annotation_9regions" / "vcf"  # VCF par région candidate (annotation_variants_gff.py)
ZONE_POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # listes d'échantillons par zone géo
RACE_DIR = Path("analyses/LD/popmaps_races_v1")  # listes d'échantillons par race
OUTDIR = BASE / "het_inbreeding"
OUTDIR.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si absent

ZONES = ["Awassi", "MiddleEastNonAwassi", "Africa", "Asia", "Europe", "America", "Australia"]
race_map = pd.read_csv(RACE_DIR / "race_to_geo_group_mapping.tsv", sep="\t")  # table race -> zone géo
slug_to_geo = dict(zip(race_map["race_slug"], race_map["geo_group"]))  # dictionnaire race -> zone
RACES = sorted(p.stem for p in RACE_DIR.glob("*.txt"))  # une race par fichier liste trouvé

# Les 9 régions candidates avec leur P3 (groupe le plus proche localement) déjà identifié
REGIONS = {
    "chr2_112.8Mb_NIPA2_CYFIP1": "Australia", "chr3_129.2Mb_desert": "Europe",
    "chr3_137.4Mb_LOC101123547": "America", "chr3_186.1Mb_CCDC91": "Asia",
    "chr5_58.1Mb_ABLIM3_AFAP1L1": "Europe", "chr6_70.2Mb_KIT": "Africa",
    "chr10_49.4Mb_KLF12": "Asia", "chr17_34.2Mb_SPATA5": "Asia", "chr20_0.8Mb_KHDRBS2": "Europe",
}


def load_list(path):
    # lit un fichier texte, un identifiant d'échantillon par ligne
    return [x.strip() for x in open(path) if x.strip()]


# dictionnaire {nom_du_groupe: liste d'échantillons}, pour les 7 zones et pour chaque race
zone_samples = {z: load_list(ZONE_POP_DIR / f"{z}.txt") for z in ZONES}
race_samples = {r: load_list(RACE_DIR / f"{r}.txt") for r in RACES}
# GROUPS réunit zones et races dans un seul dict, clé = (type, nom) -> tous les groupes traités pareil
GROUPS = {**{("zone", z): zone_samples[z] for z in ZONES}, **{("race", r): race_samples[r] for r in RACES}}


def gt_alt_count(gt):
    # convertit un génotype VCF (ex "0/1", "1|1") en nombre d'allèles alternatifs (0, 1 ou 2)
    gt = gt.split(":")[0]  # ne garde que le champ GT (avant le premier ":")
    if gt in {"./.", ".|.", "."}:
        return np.nan  # génotype manquant
    sep = "|" if "|" in gt else "/"  # phasé ("|") ou non phasé ("/")
    parts = gt.split(sep)
    if len(parts) != 2 or "." in parts:
        return np.nan  # génotype incomplet/invalide
    try:
        return int(parts[0]) + int(parts[1])  # somme des 2 allèles = dosage 0/1/2
    except ValueError:
        return np.nan


def compute_het_f(vcf_path):
    """Une lecture du VCF -> matrice (n_sites, n_samples) alt_count, puis F/Ho
    pour chaque groupe défini dans GROUPS, en un seul passage vectorisé."""
    samples = subprocess.check_output(["bcftools", "query", "-l", str(vcf_path)], text=True).splitlines()  # liste des échantillons du VCF
    sample_idx = {s: i for i, s in enumerate(samples)}  # position de chaque échantillon dans les colonnes

    # une seule lecture du VCF : un génotype par échantillon et par site, sur toute la longueur du fichier
    txt = subprocess.check_output(
        ["bash", "-lc", f"bcftools query -f '[\\t%GT]\\n' {vcf_path}"], text=True,
    )
    rows = []
    for line in txt.splitlines():
        gts = line.split("\t")[1:]  # génotypes de tous les échantillons pour ce site
        rows.append([gt_alt_count(gt) for gt in gts])  # converti en dosage 0/1/2/NaN
    mat = np.array(rows, dtype=float)  # (n_sites, n_samples), valeurs 0/1/2/NaN

    is_het = mat == 1  # masque : génotype hétérozygote (dosage=1)
    is_called = ~np.isnan(mat)  # masque : génotype non manquant

    results = {}
    for (kind, name), sample_list in GROUPS.items():
        idx = np.array([sample_idx[s] for s in sample_list if s in sample_idx], dtype=int)  # colonnes du groupe
        if len(idx) == 0:
            results[(kind, name)] = (np.nan, np.nan, 0)  # groupe vide (aucun échantillon dans ce VCF)
            continue
        sub_called = is_called[:, idx]  # sous-matrice "appelé" restreinte au groupe
        sub_alt = mat[:, idx]  # sous-matrice dosages restreinte au groupe
        n_called_per_site = sub_called.sum(axis=1)  # nb d'individus appelés à chaque site, dans ce groupe

        # Fréquence allélique p du groupe à chaque site (>=2 appelés)
        alt_sum = np.where(sub_called, sub_alt, 0).sum(axis=1)
        two_n = 2 * n_called_per_site
        with np.errstate(invalid="ignore", divide="ignore"):
            p = np.where(n_called_per_site >= 2, alt_sum / np.where(two_n == 0, 1, two_n), np.nan)

        # Site retenu seulement s'il est POLYMORPHE au sein du groupe (0<p<1) —
        # convention vcftools --het : un site monomorphe dans le sous-groupe
        # n'apporte aucune information sur F/Ho et est exclu de N_SITES (pas
        # compté comme "homozygote", sinon ça dilue Ho vers le bas à tort).
        valid_sites = (n_called_per_site >= 2) & (p > 0) & (p < 1)
        # E(HOM) avec estimateur non biaisé (2n/(2n-1))*2pq pour E(HET), même
        # correction petit-échantillon que pi (scripts 28/50/53/54) — validé
        # empiriquement contre vcftools --het (E(HOM) 198.91 vs 198.9 attendu).
        with np.errstate(invalid="ignore", divide="ignore"):
            exp_hom_site = 1 - (two_n / (two_n - 1)) * 2 * p * (1 - p)

        # Par individu du groupe : O(HOM) et E(HOM) sommés sur les sites où il est appelé et p défini
        n_ind = len(idx)
        o_hom = np.zeros(n_ind)
        e_hom = np.zeros(n_ind)
        n_sites_ind = np.zeros(n_ind)
        for j in range(n_ind):
            called_j = sub_called[:, j] & valid_sites
            vals_j = sub_alt[called_j, j]
            o_hom[j] = np.sum((vals_j == 0) | (vals_j == 2))
            e_hom[j] = np.nansum(exp_hom_site[called_j])
            n_sites_ind[j] = called_j.sum()

        # F = 1 - O(HET)/E(HET), PAS 1-O(HOM)/E(HOM) malgré les noms de colonnes
        # vcftools --het — validé à la main sur 2 individus (CNOA-AM1 : O(HOM)=316,
        # E(HOM)=198.9, N=316 -> O(HET)=0, E(HET)=117.1 -> F=1-0/117.1=1.0, exact ;
        # CNOA-AM2 : O(HOM)=179, E(HOM)=189.2, N=301 -> F=-0.0912, exact).
        o_het = n_sites_ind - o_hom
        e_het = n_sites_ind - e_hom
        with np.errstate(invalid="ignore", divide="ignore"):
            ho_ind = 1 - o_hom / n_sites_ind
            f_ind = np.where(e_het > 0, 1 - o_het / e_het, np.nan)

        valid_mask = n_sites_ind > 0
        results[(kind, name)] = (
            np.nanmean(f_ind[valid_mask]) if valid_mask.any() else np.nan,
            np.nanmean(ho_ind[valid_mask]) if valid_mask.any() else np.nan,
            int(valid_mask.sum()),
        )
    return results


def rows_from_results(results, extra=None):
    # transforme le dict {(kind,name): (F,Ho,n)} en lignes de table prêtes pour un DataFrame
    rows = []
    for (kind, name), (F, Ho, n) in results.items():
        row = {"F": F, "Ho": Ho, "n_individuals_used": n}
        if kind == "zone":
            row["group"] = name
        else:
            row["race"] = name
            row["geo_group"] = slug_to_geo.get(name, "?")  # rattache la race à sa zone géo
        if extra:
            row.update(extra)  # ajoute des colonnes fixes (ex: region_id, P3_best)
        rows.append((kind, row))
    return rows


# ── 1. Genome-wide ────────────────────────────────────────────────────────────
print("Genome-wide : lecture + calcul (une seule passe)...")
res_gw = compute_het_f(GENOMEWIDE_VCF)  # F/Ho sur l'échantillon de fond, pour les 7 zones + toutes les races
rows = rows_from_results(res_gw)
zone_rows = [r for k, r in rows if k == "zone"]  # sépare les lignes "zone" des lignes "race"
race_rows = [r for k, r in rows if k == "race"]
pd.DataFrame(zone_rows).to_csv(OUTDIR / "Het_F_genomewide_by_zone.tsv", sep="\t", index=False, float_format="%.4f")
pd.DataFrame(race_rows).sort_values(["geo_group", "F"]).to_csv(  # trié par zone puis par F croissant
    OUTDIR / "Het_F_genomewide_by_race.tsv", sep="\t", index=False, float_format="%.4f")
print(pd.DataFrame(zone_rows).to_string(index=False))  # affiche le résumé par zone dans le terminal
print(f"→ {OUTDIR}/Het_F_genomewide_by_zone.tsv")
print(f"→ {OUTDIR}/Het_F_genomewide_by_race.tsv")

# ── 2. 9 régions ──────────────────────────────────────────────────────────────
print("\n9 régions candidates...")
all_zone_rows, all_race_rows = [], []
for region_id, p3 in REGIONS.items():  # une passe complète par région (VCF différent à chaque fois)
    vcf = REGION_VCF_DIR / f"{region_id}.vcf.gz"
    res = compute_het_f(vcf)
    rows = rows_from_results(res, extra={"region_id": region_id, "P3_best": p3})
    for k, r in rows:
        if k == "zone":
            all_zone_rows.append(r)
        else:
            r["is_P3_zone"] = r["geo_group"] == p3  # marque si cette race appartient au groupe P3 de la région
            all_race_rows.append(r)
    print(f"  {region_id} : terminé")

pd.DataFrame(all_zone_rows).to_csv(OUTDIR / "Het_F_9regions_by_zone.tsv", sep="\t", index=False, float_format="%.4f")
pd.DataFrame(all_race_rows).to_csv(OUTDIR / "Het_F_9regions_by_race.tsv", sep="\t", index=False, float_format="%.4f")
print(f"→ {OUTDIR}/Het_F_9regions_by_zone.tsv")
print(f"→ {OUTDIR}/Het_F_9regions_by_race.tsv")
