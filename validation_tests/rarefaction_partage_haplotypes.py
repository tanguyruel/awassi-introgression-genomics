#!/usr/bin/env python3
"""
Contrôle par raréfaction du partage d'haplotypes Awassi/P3.

Problème signalé : le partage d'haplotypes brut (calculé en amont, dans le dépôt de
travail) dépend du nombre
d'haplotypes disponibles dans le groupe P3 — un grand groupe (Africa n=334
haplotypes) a mécaniquement plus de chances de contenir par hasard un
haplotype proche qu'un petit groupe (America n=16). Pour comparer les 9
régions à armes égales, on sous-échantillonne chaque groupe P3 à un nombre
d'haplotypes commun (= la plus petite taille de P3 parmi les 9 régions,
America n=16) et on répète le tirage N_ITER fois.

Réutilise les allele_string déjà calculées par le heatmap_haplotypes_arbre_ld.py (*_row_order.tsv),
pas de nouvelle extraction VCF.

Sortie : analyses/synthese_resultats/rarefaction_haplotype_sharing_9regions/
  Rarefaction_haplotype_sharing_9regions.tsv

Script de référence pour "Raréfaction" (cf. Annexe A du rapport de stage).
Usage : python3 34_rarefaction_haplotype_sharing_9regions_v1.py
(nécessite les sorties de 26_heatmap_LD_arbre_9regions_v2.py)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # racine du dépôt, pour importer _shared
from _shared import default_project

PROJECT = default_project()  # racine des données : --project, sinon $AWASSI_PROJECT_DIR, sinon cwd
HEATMAP_DIR = PROJECT / "analyses/haplotype_heatmap/Awassi_haplo/results/figures_finales/regions_A_9regions_v2"  # dossier des sorties du heatmap_haplotypes_arbre_ld.py (row_order.tsv)
OUTDIR = PROJECT / "analyses/synthese_resultats/rarefaction_haplotype_sharing_9regions"
OUTDIR.mkdir(parents=True, exist_ok=True)

MISMATCH_THRESHOLD = 0.01  # seuil de fraction de mésappariement pour un match "quasi identique" (<1%)
N_ITER = 500
RNG_SEED = 12345  # graine du générateur aléatoire, pour reproductibilité

REGIONS = [  # (préfixe de fichiers de la région, groupe P3 associé) pour les 9 régions
    ("chr2_112.785_112.865Mb_Australia", "Australia"),
    ("chr3_129.220_129.260Mb_Europe", "Europe"),
    ("chr3_137.390_137.420Mb_America", "America"),
    ("chr3_186.075_186.125Mb_Asia", "Asia"),
    ("chr5_58.070_58.110Mb_Europe", "Europe"),
    ("chr6_70.240_70.295Mb_Africa", "Africa"),
    ("chr10_49.360_49.390Mb_Asia", "Asia"),
    ("chr17_34.160_34.220Mb_Asia", "Asia"),
    ("chr20_0.765_0.845Mb_Europe", "Europe"),
]


def string_to_array(s):
    """Convertit une allele_string ("0"/"1"/"N") en tableau numpy (0/1/-1 pour manquant)."""
    return np.array([{"0": 0, "1": 1, "N": -1}.get(c, -1) for c in s], dtype=np.int8)


def min_mismatch_fracs(A, B):
    """Pour chaque haplotype de A, fraction de mésappariement avec son meilleur match dans B."""
    # comparaison croisée A x B x SNP : True si les 2 allèles sont valides (non manquants) à ce SNP
    valid = (A[:, None, :] != -1) & (B[None, :, :] != -1)
    # True si les allèles diffèrent ET sont tous deux valides (mésappariement comptabilisable)
    mismatch = (A[:, None, :] != B[None, :, :]) & valid
    n_valid = valid.sum(axis=2)  # nb de SNP comparables pour chaque paire (haplotype A, haplotype B)
    n_mismatch = mismatch.sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        # fraction de mésappariement par paire = n_mismatch / n_valid ; NaN si aucun SNP comparable
        frac = np.where(n_valid > 0, n_mismatch / np.maximum(n_valid, 1), np.nan)
    return np.nanmin(frac, axis=1)


def load_matrices(prefix, p3):
    """Charge les matrices haplotypes x SNP d'Awassi et du groupe P3 depuis le row_order.tsv déjà calculé.

    Parameters
    ----------
    prefix : str
        Préfixe de fichiers de la région (voir REGIONS).
    p3 : str
        Nom du groupe P3 associé à cette région.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        (A, B) : matrices d'haplotypes (0/1/-1) d'Awassi et de P3.
    """
    df = pd.read_csv(HEATMAP_DIR / f"{prefix}_row_order.tsv", sep="\t")
    awassi = df[df["group"] == "Awassi"]
    p3_df = df[df["group"] == p3]
    A = np.stack([string_to_array(s) for s in awassi["allele_string"]])  # matrice haplotypes Awassi x SNP
    B = np.stack([string_to_array(s) for s in p3_df["allele_string"]])  # matrice haplotypes P3 x SNP
    return A, B


if __name__ == "__main__":
    # 1) taille P3 commune = la plus petite parmi les 9 régions
    sizes = {}
    for prefix, p3 in REGIONS:
        _, B = load_matrices(prefix, p3)  # seul B nous intéresse ici, pour connaître sa taille
        sizes[prefix] = B.shape[0]
    common_n = min(sizes.values())  # plus petit effectif P3 parmi toutes les régions -> taille de raréfaction
    print(f"Tailles P3 (haplotypes) par région : {sizes}")
    print(f"-> taille commune de raréfaction : {common_n} haplotypes\n")

    rng = np.random.default_rng(RNG_SEED)
    rows = []

    for prefix, p3 in REGIONS:
        A, B = load_matrices(prefix, p3)
        n_awassi, n_p3 = A.shape[0], B.shape[0]

        pct_exact_list, pct_lt1_list, mean_min_frac_list = [], [], []
        for _ in range(N_ITER):
            # tire common_n haplotypes P3 sans remise (ou tous les B si déjà <= common_n)
            idx = rng.choice(n_p3, size=common_n, replace=False) if n_p3 > common_n else np.arange(n_p3)
            Bsub = B[idx]
            min_frac = min_mismatch_fracs(A, Bsub)
            pct_exact_list.append(100 * np.mean(min_frac == 0))
            pct_lt1_list.append(100 * np.mean(min_frac < MISMATCH_THRESHOLD))
            mean_min_frac_list.append(np.nanmean(min_frac))

        row = {
            "region": prefix, "P3": p3, "n_awassi_haplotypes": n_awassi,
            "n_P3_haplotypes_total": n_p3, "n_P3_haplotypes_rarefied": common_n,
            "n_iter": N_ITER,
            "pct_exact_mean": round(float(np.mean(pct_exact_list)), 1),
            "pct_exact_sd": round(float(np.std(pct_exact_list)), 1),
            "pct_lt1pct_mean": round(float(np.mean(pct_lt1_list)), 1),
            "pct_lt1pct_sd": round(float(np.std(pct_lt1_list)), 1),
            "mean_min_mismatch_frac": round(float(np.mean(mean_min_frac_list)), 4),
        }
        rows.append(row)
        print(f"{prefix:35s} P3={p3:10s} (n_P3_total={n_p3:4d} -> raréfié à {common_n}) : "
              f"exact={row['pct_exact_mean']:.1f}%±{row['pct_exact_sd']:.1f}  "
              f"<1%={row['pct_lt1pct_mean']:.1f}%±{row['pct_lt1pct_sd']:.1f}")

    out = pd.DataFrame(rows)
    out_tsv = OUTDIR / "Rarefaction_haplotype_sharing_9regions.tsv"
    out.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nTable écrite : {out_tsv}")
