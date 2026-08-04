#!/usr/bin/env python3
"""
LD (r²) recalculé séparément pour Awassi seul, P3 seul et ME seul,
en plus du LD combiné Awassi+ME+P3 déjà utilisé dans les heatmaps (heatmap_haplotypes_arbre_ld.py).

Objectif (retour utilisateur) : le LD combiné peut être gonflé par la
structure ENTRE groupes (si Awassi/ME/P3 ont des fréquences alléliques très
différentes à plusieurs SNPs, le LD calculé sur le mélange peut sembler élevé
même sans vrai bloc haplotypique à l'intérieur de chaque groupe — effet
Wahlund-like). Recalculer le LD DANS chaque groupe séparément permet de
vérifier si le bloc haplotypique existe réellement à l'intérieur d'Awassi
et/ou de P3, ou si c'est surtout un artefact de mélange de groupes.

Réutilise les fonctions de heatmap_haplotypes_arbre_ld.py (import direct, même VCF phasés).

Sortie : analyses/synthese_resultats/ld_per_subgroup_9regions/
  LD_per_subgroup_9regions.tsv

Script de référence pour "LD par sous-groupe" (appelé "ld_par_sousgroupe.py" dans l'Annexe A du
rapport de stage). Usage : python3 ld_par_sousgroupe.py  (nécessite heatmap_haplotypes_arbre_ld.py
dans le même dossier, importé directement)
"""

import importlib.util
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # racine du dépôt, pour importer _shared
from _shared import default_project

# ── Import dynamique de heatmap_haplotypes_arbre_ld.py, pour réutiliser
# directement ses fonctions (load_haplotypes, restrict_positions,
# filter_snps, compute_r2) sans dupliquer le code ──
SCRIPT_HEATMAP = Path(__file__).with_name("heatmap_haplotypes_arbre_ld.py")  # même dossier
spec = importlib.util.spec_from_file_location("heatmap_ld", SCRIPT_HEATMAP)
heatmap_ld = importlib.util.module_from_spec(spec)  # instancie le module (vide pour l'instant)
spec.loader.exec_module(heatmap_ld)  # exécute heatmap_haplotypes_arbre_ld.py -> ses fonctions deviennent accessibles via heatmap_ld.xxx

PROJECT = default_project()  # racine des données : --project, sinon $AWASSI_PROJECT_DIR, sinon cwd
POP_DIR = PROJECT / "analyses/fst/popmaps_separees_v1"  # dossier des popmaps (listes d'individus par groupe)
OUTDIR = PROJECT / "analyses/synthese_resultats/ld_per_subgroup_9regions"
OUTDIR.mkdir(parents=True, exist_ok=True)

MAX_MISSING = 0.05
MIN_MAF = 0.05

# (chr, start, end, P3, vcf) — mêmes fenêtres et VCF phasés que heatmap_haplotypes_arbre_ld.py
REGIONS = [
    {"region_id": "chr3_129.2Mb_desert", "chr": 3, "start": 129220001, "end": 129260000, "P3": "Europe",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr3_129.22_129.26Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr17_SPATA5", "chr": 17, "start": 34160001, "end": 34220000, "P3": "Asia",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v1_5regions/phased/chr17_34.16_34.22Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr3_CCDC91", "chr": 3, "start": 186075001, "end": 186125000, "P3": "Asia",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v1_5regions/phased/chr3_186.07_186.12Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr20_KHDRBS2", "chr": 20, "start": 765001, "end": 845000, "P3": "Europe",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr20_0.765_0.845Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr6_KIT", "chr": 6, "start": 70240001, "end": 70295000, "P3": "Africa",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_12_06_14h28/phased/chr6_68.86_70.95Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr2_NIPA2_CYFIP1", "chr": 2, "start": 112785001, "end": 112865000, "P3": "Australia",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr2_112.79_112.87Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr5_ABLIM3_AFAP1L1", "chr": 5, "start": 58070001, "end": 58110000, "P3": "Europe",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr5_58.07_58.11Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr10_KLF12", "chr": 10, "start": 49360001, "end": 49390000, "P3": "Asia",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v3_2regions_new/phased/chr10_49.36_49.39Mb.beagle_phased.vcf.gz"},
    {"region_id": "chr3_LOC101123547", "chr": 3, "start": 137390001, "end": 137420000, "P3": "America",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v3_2regions_new/phased/chr3_137.39_137.42Mb.beagle_phased.vcf.gz"},
]


def ld_summary_for_subset(H, positions, mask):
    """Calcule le LD (r²) et ses statistiques résumées pour un sous-ensemble d'haplotypes.

    Parameters
    ----------
    H : numpy.ndarray
        Matrice haplotypes x SNP de la région.
    positions : numpy.ndarray
        Positions (bp) des SNP, alignées sur les colonnes de H.
    mask : numpy.ndarray
        Masque booléen des haplotypes du sous-groupe à analyser.

    Returns
    -------
    dict | None
        n_haplotypes, n_snps_used, mean_r2, median_r2, prop_r2_ge_0.5,
        prop_r2_ge_0.8 ; None si le sous-groupe est trop petit ou n'a pas
        assez de SNP après filtre qualité.
    """
    Hs = H[mask, :]
    if Hs.shape[0] < 4:
        return None  # trop peu d'haplotypes pour un LD fiable
    try:
        Hf, posf = heatmap_ld.filter_snps(Hs, positions.copy(), MAX_MISSING, MIN_MAF)  # filtre qualité SNP (fonction de heatmap_haplotypes_arbre_ld.py)
    except RuntimeError:
        return None  # moins de 2 SNP après filtre -> sous-groupe ignoré
    r2, ld_pos = heatmap_ld.compute_r2(Hf, posf)  # LD (r²) calculé dans ce sous-groupe uniquement (fonction de heatmap_haplotypes_arbre_ld.py)
    tri = np.tril_indices(len(ld_pos), k=-1)  # indices du triangle inférieur (paires uniques, hors diagonale)
    vals = r2[tri]
    off = vals[np.isfinite(vals)]
    if len(off) == 0:
        return None
    return {
        "n_haplotypes": Hs.shape[0], "n_snps_used": len(ld_pos),
        "mean_r2": float(np.nanmean(off)), "median_r2": float(np.nanmedian(off)),
        "prop_r2_ge_0.5": float(np.mean(off >= 0.5)), "prop_r2_ge_0.8": float(np.mean(off >= 0.8)),
    }


def compute_region(reg):
    """Calcule le LD par sous-groupe (Awassi seul, P3 seul, ME seul, combiné) pour une région candidate.

    Parameters
    ----------
    reg : dict
        Une entrée de REGIONS (region_id, chr, start, end, P3, vcf).

    Returns
    -------
    list[dict]
        Une ligne par sous-groupe évalué (voir ld_summary_for_subset).
    """
    region_id, chrom, start, end, p3, vcf = (
        reg["region_id"], reg["chr"], reg["start"], reg["end"], reg["P3"], reg["vcf"]
    )
    print(f"\n{'='*70}\n{region_id} | P3={p3}")

    group_order = ["Awassi", "MiddleEastNonAwassi", p3]
    H, positions, haplotypes, individuals, hap_index, groups, selected, s2g = heatmap_ld.load_haplotypes(
        str(PROJECT / vcf), str(POP_DIR), group_order
    )  # charge la matrice haplotypes x SNP depuis le VCF phasé (fonction de heatmap_haplotypes_arbre_ld.py)
    H, positions = heatmap_ld.restrict_positions(H, positions, start, end)  # restreint à la fenêtre de la région (fonction de heatmap_haplotypes_arbre_ld.py)

    rows = []
    subsets = {  # masques booléens définissant chaque sous-groupe à comparer
        "Awassi_seul": groups == "Awassi",
        "P3_seul": groups == p3,
        "ME_seul": groups == "MiddleEastNonAwassi",
        "Awassi+ME+P3_combine": np.isin(groups, ["Awassi", "MiddleEastNonAwassi", p3]),
    }
    for label, mask in subsets.items():
        res = ld_summary_for_subset(H, positions, mask)
        if res is None:
            print(f"  {label:22s} : insuffisant (trop peu d'haplotypes/SNPs après filtre)")
            continue
        print(f"  {label:22s} : n_hap={res['n_haplotypes']:4d} n_snps={res['n_snps_used']:4d} "
              f"mean_r2={res['mean_r2']:.3f} prop>=0.8={res['prop_r2_ge_0.8']:.3f}")
        rows.append({"region_id": region_id, "P3": p3, "subset": label, **res})

    return rows


if __name__ == "__main__":
    all_rows = []
    for reg in REGIONS:
        all_rows += compute_region(reg)

    df = pd.DataFrame(all_rows)
    out_tsv = OUTDIR / "LD_per_subgroup_9regions.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)  # sauvegarde la table finale (une ligne par région x sous-groupe)
    print(f"\nTable écrite : {out_tsv}")
