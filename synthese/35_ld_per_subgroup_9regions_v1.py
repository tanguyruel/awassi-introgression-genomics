#!/usr/bin/env python3
"""
Script 35 — LD (r²) recalculé séparément pour Awassi seul, P3 seul et ME seul,
en plus du LD combiné Awassi+ME+P3 déjà utilisé dans les heatmaps (script 26).

Objectif (retour utilisateur) : le LD combiné peut être gonflé par la
structure ENTRE groupes (si Awassi/ME/P3 ont des fréquences alléliques très
différentes à plusieurs SNPs, le LD calculé sur le mélange peut sembler élevé
même sans vrai bloc haplotypique à l'intérieur de chaque groupe — effet
Wahlund-like). Recalculer le LD DANS chaque groupe séparément permet de
vérifier si le bloc haplotypique existe réellement à l'intérieur d'Awassi
et/ou de P3, ou si c'est surtout un artefact de mélange de groupes.

Réutilise les fonctions du script 26 (import direct, même VCF phasés).

Sortie : analyses/synthese_resultats/ld_per_subgroup_9regions/
  LD_per_subgroup_9regions.tsv

Script de référence pour "LD par sous-groupe" (cf. Annexe A du rapport de stage).
Usage : python3 35_ld_per_subgroup_9regions_v1.py  (nécessite 26_heatmap_LD_arbre_9regions_v2.py
dans le même dossier, importé directement)
"""

# ── Imports : import dynamique de module, chemins, calcul, tables ──
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

# ── Import dynamique du script 26 (26_heatmap_LD_arbre_9regions_v2.py), pour
# réutiliser directement ses fonctions (load_haplotypes, restrict_positions,
# filter_snps, compute_r2) sans dupliquer le code ──
SCRIPT26 = Path(__file__).with_name("26_heatmap_LD_arbre_9regions_v2.py")  # chemin du script 26 (même dossier)
spec = importlib.util.spec_from_file_location("heatmap26", SCRIPT26)  # spec du module à charger dynamiquement
heatmap26 = importlib.util.module_from_spec(spec)  # instancie le module (vide pour l'instant)
spec.loader.exec_module(heatmap26)  # exécute le script 26 -> ses fonctions deviennent accessibles via heatmap26.xxx

POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # dossier des popmaps (listes d'individus par groupe)
OUTDIR = Path("analyses/synthese_resultats/ld_per_subgroup_9regions")  # dossier de sortie de la table
OUTDIR.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie s'il n'existe pas

MAX_MISSING = 0.05  # seuil max de génotypes manquants toléré par SNP
MIN_MAF = 0.05  # fréquence allélique mineure minimale pour garder un SNP

# (chr, start, end, P3, vcf) — mêmes fenêtres et VCF phasés que script 26
REGIONS = [  # une entrée par région, mêmes 9 régions/VCF que le script 26
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
    Hs = H[mask, :]  # sous-matrice des haplotypes appartenant au sous-groupe (mask booléen)
    if Hs.shape[0] < 4:
        return None  # trop peu d'haplotypes pour un LD fiable
    try:
        Hf, posf = heatmap26.filter_snps(Hs, positions.copy(), MAX_MISSING, MIN_MAF)  # filtre qualité SNP (fonction du script 26)
    except RuntimeError:
        return None  # moins de 2 SNP après filtre -> sous-groupe ignoré
    r2, ld_pos = heatmap26.compute_r2(Hf, posf)  # LD (r²) calculé dans ce sous-groupe uniquement (fonction du script 26)
    tri = np.tril_indices(len(ld_pos), k=-1)  # indices du triangle inférieur (paires uniques, hors diagonale)
    vals = r2[tri]  # valeurs r² de toutes les paires de SNP
    off = vals[np.isfinite(vals)]  # exclut les valeurs non finies
    if len(off) == 0:
        return None
    return {
        "n_haplotypes": Hs.shape[0], "n_snps_used": len(ld_pos),
        "mean_r2": float(np.nanmean(off)), "median_r2": float(np.nanmedian(off)),
        "prop_r2_ge_0.5": float(np.mean(off >= 0.5)), "prop_r2_ge_0.8": float(np.mean(off >= 0.8)),  # proportions de paires en fort LD
    }


def compute_region(reg):
    region_id, chrom, start, end, p3, vcf = (
        reg["region_id"], reg["chr"], reg["start"], reg["end"], reg["P3"], reg["vcf"]
    )  # déballe les paramètres de la région
    print(f"\n{'='*70}\n{region_id} | P3={p3}")

    group_order = ["Awassi", "MiddleEastNonAwassi", p3]
    H, positions, haplotypes, individuals, hap_index, groups, selected, s2g = heatmap26.load_haplotypes(
        vcf, str(POP_DIR), group_order
    )  # charge la matrice haplotypes x SNP depuis le VCF phasé (fonction du script 26)
    H, positions = heatmap26.restrict_positions(H, positions, start, end)  # restreint à la fenêtre de la région (fonction du script 26)

    rows = []
    subsets = {  # masques booléens définissant chaque sous-groupe à comparer
        "Awassi_seul": groups == "Awassi",
        "P3_seul": groups == p3,
        "ME_seul": groups == "MiddleEastNonAwassi",
        "Awassi+ME+P3_combine": np.isin(groups, ["Awassi", "MiddleEastNonAwassi", p3]),
    }
    for label, mask in subsets.items():
        res = ld_summary_for_subset(H, positions, mask)  # LD résumé pour ce sous-groupe
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
        all_rows += compute_region(reg)  # calcule le LD par sous-groupe pour chaque région

    df = pd.DataFrame(all_rows)
    out_tsv = OUTDIR / "LD_per_subgroup_9regions.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)  # sauvegarde la table finale (une ligne par région x sous-groupe)
    print(f"\nTable écrite : {out_tsv}")
