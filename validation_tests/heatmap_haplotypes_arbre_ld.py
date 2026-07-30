#!/usr/bin/env python3
"""
Heatmap génotypes + arbre + LD pour les 9 régions candidates
(7 régions A v10b + 2 régions nouvelles KLF12/chr10 et LOC101123547/chr3).

Reprend une version antérieure du dépôt de travail (7 régions) à l'identique, avec
2 ajouts :
  1) 2 régions supplémentaires (VCF phasés séparément, voir phasing_beagle.sh pour la méthode).
  2) Annotation du/des gène(s) d'intérêt dans chaque fenêtre : bande verte
     semi-transparente + traits pointillés à la position exacte (GFF Oar_v4.0,
     coords vérifiées directement dans le GFF, cf. AWASSI_AGENT_LOG.md) sur la
     heatmap ET sur le triangle LD (même axe X), nom du gène affiché au-dessus.

Sorties : analyses/haplotype_heatmap/Awassi_haplo/results/figures_finales/regions_A_9regions_v2/

Script de référence pour "Heatmap + arbre + LD" (cf. Annexe A du rapport de stage,
et skill .claude/skills/heatmap-ld-region/). Remplace 15_heatmap_LD_arbre_7regions_v1.py
et les anciennes tentatives dans scripts/04_heatmap_ggtree/ et 04_heatmap_haplotype/
(exploratoires, non retenues). Importé directement par le ld_par_sousgroupe.py (LD par sous-groupe).
Usage : python3 26_heatmap_LD_arbre_9regions_v2.py
"""

import os
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram

from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from matplotlib.collections import PolyCollection
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # racine du dépôt, pour importer _shared
from _shared import run as run_cmd, read_list

POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # dossier des popmaps (listes d'individus par groupe)
OUTDIR  = Path("analyses/haplotype_heatmap/Awassi_haplo/results/figures_finales/regions_A_9regions_v2")
OUTDIR.mkdir(parents=True, exist_ok=True)

MAX_MISSING = 0.05
MIN_MAF     = 0.05
LABEL_FONTSIZE = 2.6
TREE_WIDTH     = 5.2

GENE_COLOR = "#2ca02c"  # couleur verte des annotations de gène

# ── Régions (7 régions A v10b + 2 nouvelles) — VCF phasés déjà calculés ─────
# "genes" : gène(s) d'intérêt de la région, coords exactes GFF Oar_v4.0
# (vérifiées par requête directe dans le GFF, cf. AWASSI_AGENT_LOG.md 07/07).
REGIONS = [
    {"chr": 3,  "start": 129220001, "end": 129260000, "P3": "Europe",  # région 1 : chr3, P3=Europe
     "label": "chr3 129.220–129.260 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr3_129.22_129.26Mb.beagle_phased.vcf.gz",
     "genes": []},  # désert génique confirmé (aucun gène nommé dans la fenêtre)
    {"chr": 17, "start": 34160001,  "end": 34220000,  "P3": "Asia",  # région 2 : chr17, gène SPATA5
     "label": "chr17 34.160–34.220 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v1_5regions/phased/chr17_34.16_34.22Mb.beagle_phased.vcf.gz",
     "genes": [{"name": "SPATA5", "start": 34076424, "end": 34424110}]},
    {"chr": 3,  "start": 186075001, "end": 186125000, "P3": "Asia",  # région 3 : chr3, gène CCDC91
     "label": "chr3 186.075–186.125 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v1_5regions/phased/chr3_186.07_186.12Mb.beagle_phased.vcf.gz",
     "genes": [{"name": "CCDC91", "start": 185706850, "end": 186137244}]},
    {"chr": 20, "start": 765001,    "end": 845000,    "P3": "Europe",  # région 4 : chr20, gène KHDRBS2
     "label": "chr20 0.765–0.845 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr20_0.765_0.845Mb.beagle_phased.vcf.gz",
     "genes": [{"name": "KHDRBS2", "start": 198543, "end": 1026276}]},
    {"chr": 6,  "start": 70240001,  "end": 70295000,  "P3": "Africa",  # région 5 : chr6, désert génique près de KIT
     "label": "chr6 70.240–70.295 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_12_06_14h28/phased/chr6_68.86_70.95Mb.beagle_phased.vcf.gz",
     "genes": []},  # désert génique ; KIT est ~98kb en amont, hors fenêtre
    {"chr": 2,  "start": 112785001, "end": 112865000, "P3": "Australia",  # région 6 : chr2, gènes NIPA2/CYFIP1
     "label": "chr2 112.785–112.865 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr2_112.79_112.87Mb.beagle_phased.vcf.gz",
     "genes": [{"name": "NIPA2", "start": 112788268, "end": 112810606},
               {"name": "CYFIP1", "start": 112811686, "end": 112885989}]},
    {"chr": 5,  "start": 58070001,  "end": 58110000,  "P3": "Europe",  # région 7 : chr5, gènes ABLIM3/AFAP1L1
     "label": "chr5 58.070–58.110 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v2_5regions_v10b/phased/chr5_58.07_58.11Mb.beagle_phased.vcf.gz",
     "genes": [{"name": "ABLIM3", "start": 57956105, "end": 58090022},
               {"name": "AFAP1L1", "start": 58109106, "end": 58177188}]},
    # ── 2 régions nouvelles (méthode FST/fd séparés, cf. AWASSI_AGENT_LOG.md) ──
    {"chr": 10, "start": 49360001,  "end": 49390000,  "P3": "Asia",  # région 8 : chr10, gène KLF12
     "label": "chr10 49.360–49.390 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v3_2regions_new/phased/chr10_49.36_49.39Mb.beagle_phased.vcf.gz",
     "genes": [{"name": "KLF12", "start": 48919222, "end": 49664271}]},
    {"chr": 3,  "start": 137390001, "end": 137420000, "P3": "America",  # région 9 : chr3, gène LOC101123547
     "label": "chr3 137.390–137.420 Mb",
     "vcf": "analyses/phasing_beagle/Phasing_Beagle_v3_2regions_new/phased/chr3_137.39_137.42Mb.beagle_phased.vcf.gz",
     "genes": [{"name": "LOC101123547", "start": 137397357, "end": 137398295}]},
]

# ── Couleurs : mêmes que les scripts 13b_plot_fst_fd_ld_aligned_v*.py ───────
BASE_COLORS = {
    "Awassi":               "#E41A1C",
    "MiddleEastNonAwassi":  "#1f77b4",
    "Africa":               "#e07b39",
    "Asia":                 "#6a5acd",
    "Australia":            "#2ca02c",
    "America":              "#d62728",
    "Europe":               "#8c564b",
}

DISPLAY_LABELS = {
    "Awassi": "Awassi",
    "MiddleEastNonAwassi": "ME (Moyen-Orient)",
    "Africa": "Africa",
    "Asia": "Asia",
    "Australia": "Australia",
    "America": "America",
    "Europe": "Europe",
}


def display_label(g):
    """Libellé d'affichage d'un groupe (DISPLAY_LABELS, ou le nom brut avec "_" remplacés par des espaces)."""
    return DISPLAY_LABELS.get(g, g.replace("_", " "))


def group_color(g):
    """Couleur associée à un groupe (BASE_COLORS, gris par défaut si inconnu)."""
    return BASE_COLORS.get(g, "#999999")


def safe_name(x):
    """Nettoie une chaîne en un nom sûr pour fichiers (alphanumériques/._- uniquement, "_" dédupliqués)."""
    x = re.sub(r"[^A-Za-z0-9._-]+", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def parse_gt(gt_field):
    """Décompose un champ GT VCF en ses deux allèles (float, NaN si manquant/non phasé sur 2 allèles)."""
    gt = gt_field.split(":")[0]
    if gt in ["./.", ".|."]:
        return np.nan, np.nan
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)
    if len(parts) != 2:
        return np.nan, np.nan

    def conv(x):
        """Convertit un allèle texte ("0"/"1"/".") en float (NaN si manquant)."""
        if x == ".":
            return np.nan
        return float(int(x))

    return conv(parts[0]), conv(parts[1])


def load_haplotypes(vcf, pop_files_dir, groups):
    """Charge la matrice haplotypes x SNP d'un VCF phasé, restreinte aux individus des groupes demandés.

    Chaque individu diploïde donne 2 haplotypes (h1/h2), déduits des GT phasés.

    Parameters
    ----------
    vcf : str
        Chemin du VCF phasé.
    pop_files_dir : str
        Dossier des popmaps (un fichier liste par groupe, "{group}.txt").
    groups : list[str]
        Groupes à inclure ; en cas d'individu présent dans plusieurs popmaps,
        le premier groupe de la liste rencontré est prioritaire.

    Returns
    -------
    tuple
        (H, positions, hap_names, hap_individuals, hap_index, hap_groups,
        selected_samples, sample_to_group). H est la matrice (n_haplotypes,
        n_snps) de valeurs 0/1/NaN.
    """
    vcf_samples = run_cmd(["bcftools", "query", "-l", vcf]).splitlines()

    sample_to_group = {}
    for g in groups:
        p = os.path.join(pop_files_dir, f"{g}.txt")
        if not os.path.exists(p):
            print(f"Attention : pop file absent, ignoré : {p}")
            continue
        for s in read_list(p):
            if s not in sample_to_group:
                # le premier groupe rencontré est prioritaire (pas d'écrasement)
                sample_to_group[s] = g

    selected_samples = [s for s in vcf_samples if s in sample_to_group]
    if len(selected_samples) == 0:
        raise RuntimeError("Aucun individu sélectionné dans le VCF.")

    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        for s in selected_samples:
            tmp.write(s + "\n")
        sample_file = tmp.name

    try:
        q_vcf = shlex.quote(vcf)
        q_sample = shlex.quote(sample_file)
        cmd = [
            "bash", "-lc",
            f"bcftools view -S {q_sample} -m2 -M2 -v snps -Ou {q_vcf} | "
            "bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'"
        ]
        txt = run_cmd(cmd)
    finally:
        os.remove(sample_file)

    hap_names, hap_individuals, hap_index, hap_groups = [], [], [], []
    for s in selected_samples:
        # chaque individu diploïde donne 2 haplotypes (h1 et h2)
        hap_names += [s + "_h1", s + "_h2"]
        hap_individuals += [s, s]
        hap_index += ["h1", "h2"]
        hap_groups += [sample_to_group[s], sample_to_group[s]]

    positions, rows_by_snp = [], []
    for line in txt.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        positions.append(int(parts[1]))
        one_snp = []
        for gt in parts[2:]:
            a1, a2 = parse_gt(gt)
            one_snp.append(a1)
            one_snp.append(a2)
        rows_by_snp.append(one_snp)  # une ligne = un SNP, colonnes = haplotypes

    H = np.array(rows_by_snp, dtype=float).T  # transposée : lignes=haplotypes, colonnes=SNP
    positions = np.array(positions)

    return (H, positions, np.array(hap_names), np.array(hap_individuals),
            np.array(hap_index), np.array(hap_groups), selected_samples, sample_to_group)


def restrict_positions(H, positions, pos_min, pos_max):
    """Restreint H et positions aux SNP dans [pos_min, pos_max] (bornes incluses, None = pas de borne)."""
    keep = np.ones(len(positions), dtype=bool)
    if pos_min is not None:
        keep &= positions >= pos_min
    if pos_max is not None:
        keep &= positions <= pos_max
    return H[:, keep], positions[keep]


def filter_snps(H, positions, max_missing, min_maf):
    """Filtre les SNP sur missingness et fréquence de l'allèle mineur.

    Parameters
    ----------
    H : numpy.ndarray
        Matrice haplotypes x SNP (0/1/NaN).
    positions : numpy.ndarray
        Positions (bp) alignées sur les colonnes de H.
    max_missing : float
        Fraction maximale de valeurs manquantes tolérée par SNP.
    min_maf : float
        Fréquence minimale de l'allèle mineur tolérée par SNP.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        (H, positions) restreints aux SNP retenus.

    Raises
    ------
    RuntimeError
        S'il reste moins de 2 SNP après filtre.
    """
    keep = []
    for j in range(H.shape[1]):
        col = H[:, j]
        miss = np.mean(np.isnan(col))
        if miss > max_missing:
            continue
        vals = col[~np.isnan(col)]
        if len(vals) == 0:
            continue
        if len(np.unique(vals)) < 2:
            continue
        n0 = np.sum(vals == 0)
        n1 = np.sum(vals == 1)
        total = n0 + n1
        if total == 0:
            continue
        maf = min(n0, n1) / total  # fréquence de l'allèle mineur
        if maf < min_maf:
            continue
        keep.append(j)

    keep = np.array(keep, dtype=int)
    if len(keep) < 2:
        raise RuntimeError("Moins de 2 SNPs après filtre.")
    return H[:, keep], positions[keep]


def allele_string(row):
    """Convertit une ligne d'allèles en chaîne "0"/"1"/"N" (N = manquant), sert d'identité exacte."""
    return "".join("N" if np.isnan(x) else str(int(x)) for x in row)


def assign_exact_ids(H):
    """Attribue un identifiant (H001, H002, ...) à chaque séquence d'allèles distincte de H.

    Deux haplotypes ayant exactement la même allele_string reçoivent le même
    identifiant (identité exacte SNP par SNP, pas de tolérance).

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        (exact_ids, allele_strings) pour chaque haplotype (ligne de H).
    """
    strings = [allele_string(H[i, :]) for i in range(H.shape[0])]
    str_to_id = {}  # associe chaque chaîne unique à un identifiant H001, H002, ...
    exact_ids = []
    for s in strings:
        if s not in str_to_id:
            str_to_id[s] = f"H{len(str_to_id) + 1:03d}"
        exact_ids.append(str_to_id[s])
    return np.array(exact_ids), np.array(strings)


def make_tree_order(H):
    """Clustering hiérarchique (UPGMA, distance de Hamming) des haplotypes de H.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        (Z, order) : matrice de linkage scipy et ordre des feuilles du
        dendrogramme (indices dans H).
    """
    X = H.copy()
    X[np.isnan(X)] = 2  # code les valeurs manquantes en 2 (catégorie à part pour la distance)
    D = pdist(X, metric="hamming")  # distance de Hamming entre chaque paire d'haplotypes
    Z = linkage(D, method="average")  # clustering hiérarchique (average linkage / UPGMA)
    d = dendrogram(Z, no_plot=True)
    return Z, np.array(d["leaves"], dtype=int)


def matrix_for_plot(H):
    """Copie H avec les valeurs manquantes recodées à 2, pour un affichage à 3 couleurs (0/1/manquant)."""
    M = H.copy()
    M[np.isnan(M)] = 2  # code le manquant en 2 pour un affichage à 3 couleurs (0/1/manquant)
    return M


def compute_r2(H, positions):
    """Calcule le LD (r²) par paire de SNP, à partir de la matrice haplotypes x SNP.

    Le manquant est imputé par la fréquence allélique moyenne du SNP avant
    corrélation. Les SNP restés sans variance après imputation sont exclus.

    Parameters
    ----------
    H : numpy.ndarray
        Matrice haplotypes (lignes) x SNP (colonnes), valeurs 0/1/NaN.
    positions : numpy.ndarray
        Positions (bp) alignées sur les colonnes de H.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        (r2, positions) : matrice carrée de r² entre SNP retenus, et leurs positions.
    """
    X = H.copy()
    means = np.nanmean(X, axis=0)  # fréquence allélique moyenne de chaque SNP (NaN ignorés)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(means, inds[1])  # impute le manquant par la moyenne du SNP correspondant
    X = X - X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    valid = sd > 0  # exclut les SNP sans variance (devenus constants après imputation)
    X = X[:, valid]
    positions = positions[valid]
    sd = X.std(axis=0, ddof=1)  # recalcule l'écart-type sur les SNP restants uniquement
    X = X / sd
    n = X.shape[0]
    corr = (X.T @ X) / (n - 1)  # matrice de corrélation de Pearson entre SNP (produit croisé normalisé)
    corr = np.clip(corr, -1, 1)  # sécurité numérique : corrélation strictement bornée à [-1, 1]
    r2 = corr ** 2  # r² = carré du coefficient de corrélation de Pearson = mesure du LD
    np.fill_diagonal(r2, 1.0)  # r² d'un SNP avec lui-même fixé à 1 par convention
    return r2, positions


def gene_index_span(positions, gene_start, gene_end):
    """Convertit les bornes bp d'un gène en indices SNP (espace de la heatmap),
    en clippant aux bornes de la fenêtre affichée. None si hors fenêtre."""
    lo = max(gene_start, int(positions[0]))
    hi = min(gene_end, int(positions[-1]))
    if lo > hi:
        return None
    i0 = int(np.searchsorted(positions, lo, side="left"))
    i1 = int(np.searchsorted(positions, hi, side="right"))
    if i1 <= i0:
        i1 = i0 + 1  # garantit une largeur minimale d'1 SNP pour que ce soit visible
    return i0, i1


def annotate_genes(ax_mat, ax_ld, n_rows, n_snps, genes, positions):
    """Marque la position du/des gène(s) d'intérêt SANS teinter la heatmap ni le
    triangle LD : une ligne verte fine pile sur le bord supérieur de la heatmap
    (+ nom du gène au-dessus) et une ligne verte sur l'axe X du panneau LD.
    Épaisseurs de trait et décalages du texte en points fixes (linewidth,
    "offset points") : taille visuelle strictement constante quel que soit
    n_rows, jamais de chevauchement avec le titre."""
    # transformations mixtes : X en coordonnées données (SNP), Y en coordonnées axes (0-1)
    trans_mat = blended_transform_factory(ax_mat.transData, ax_mat.transAxes)
    trans_ld = blended_transform_factory(ax_ld.transData, ax_ld.transAxes)

    for gene in genes:
        span = gene_index_span(positions, gene["start"], gene["end"])
        if span is None:
            continue
        i0, i1 = span

        # Ligne verte fine pile sur le bord supérieur de la heatmap
        ax_mat.plot([i0, i1], [1.0, 1.0], transform=trans_mat, color=GENE_COLOR,
                    linewidth=4, solid_capstyle="butt", clip_on=False, zorder=6)
        # Nom du gène, décalé d'un nombre fixe de points au-dessus de la ligne
        ax_mat.annotate(gene["name"], xy=((i0 + i1) / 2, 1.0), xycoords=trans_mat,
                         xytext=(0, 8), textcoords="offset points",
                         ha="center", va="bottom", fontsize=17, fontweight="bold",
                         color=GENE_COLOR, clip_on=False, zorder=7)

        # Ligne verte sur l'axe X du panneau LD (bord bas)
        ax_ld.plot([i0, i1], [0.0, 0.0], transform=trans_ld, color=GENE_COLOR,
                   linewidth=4, solid_capstyle="butt", clip_on=False, zorder=6)


def draw_ld_down_triangle(ax, r2, positions):
    """Dessine le triangle LD (r², pointe vers le bas) aligné sur l'axe X en indices SNP de la heatmap.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axe cible.
    r2 : numpy.ndarray
        Matrice carrée de r² (voir compute_r2).
    positions : numpy.ndarray
        Positions (bp) des SNP, pour les ticks de l'axe X.

    Returns
    -------
    matplotlib.collections.PolyCollection
        La collection de polygones tracée (pour construire la colorbar).
    """
    n = len(positions)
    # palette blanc -> rouge foncé pour représenter r² croissant (0 à 1)
    cmap = LinearSegmentedColormap.from_list(
        "ld_red", ["#FFFFFF", "#FFF5EB", "#FDD0A2", "#FC8D59", "#D7301F", "#7F0000"]
    )
    verts, values = [], []  # verts : polygones (losanges) ; values : r² associé à chaque losange
    for i in range(1, n):
        for j in range(i):
            # construit un losange centré sur la position de la paire de SNP (i,j)
            x_old = (i + j + 1) / 2.0
            scale_x = n / (n - 1) if n > 1 else 1.0  # facteur d'échelle pour aligner sur l'axe de la heatmap
            x = (x_old - 0.5) * scale_x
            dx = 0.5 * scale_x
            y = (i - j) / 2.0  # hauteur du losange = distance (en indices) entre les 2 SNP
            verts.append([(x - dx, y), (x, y + 0.5), (x + dx, y), (x, y - 0.5)])
            values.append(r2[i, j])

    values = np.array(values)
    # collection de polygones colorés selon r² (triangle LD "pointe vers le bas")
    pc = PolyCollection(verts, array=values, cmap=cmap, edgecolors="none", linewidths=0.0)
    pc.set_clim(0, 1)  # échelle de couleur fixée entre 0 et 1 (bornes possibles de r²)
    ax.add_collection(pc)
    ax.set_xlim(0, n)
    ax.set_ylim(n / 2.0 + 0.5, -0.1)  # axe Y inversé : sommet du triangle en haut
    ax.set_aspect("auto")
    # contour du triangle (bord supérieur + les deux côtés obliques)
    ax.plot([0, n, n / 2.0, 0], [0, 0, n / 2.0, 0], color="black", linewidth=0.35, zorder=5)

    tick_idx = np.linspace(0, n - 1, min(9, n)).astype(int)
    tick_lab = [f"{positions[i] / 1e6:.3f}" for i in tick_idx]
    ax.set_xticks(tick_idx + 0.5)
    ax.set_xticklabels(tick_lab, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("Position génomique (Mb)", fontsize=12)
    ax.set_yticks([])
    for spine in ["left", "right", "top"]:
        ax.spines[spine].set_visible(False)  # masque les bordures non utiles
    ax.spines["bottom"].set_linewidth(0.8)
    return pc  # renvoyé pour construire la colorbar dans plot_combined


def write_exact_sharing_table(H, haplotypes, individuals, hap_index, groups, out_tsv):
    """Regroupe les haplotypes par identité exacte et écrit un résumé de partage inter-groupes en TSV.

    Parameters
    ----------
    H : numpy.ndarray
        Matrice haplotypes x SNP.
    haplotypes, individuals, hap_index, groups : numpy.ndarray
        Métadonnées alignées sur les lignes de H (voir load_haplotypes).
    out_tsv : pathlib.Path
        Chemin de sortie du tableau (une ligne par identifiant d'haplotype exact).
    """
    exact_ids, allele_strings = assign_exact_ids(H)
    df = pd.DataFrame({
        "exact_haplotype_id": exact_ids, "haplotype": haplotypes,
        "individual": individuals, "hap_index": hap_index,
        "group": groups, "allele_string": allele_strings,
    })
    rows = []
    for hid, sub in df.groupby("exact_haplotype_id", sort=False):
        # sub = tous les haplotypes partageant exactement le même identifiant (même séquence d'allèles)
        group_counts = sub["group"].value_counts().to_dict()
        rows.append({
            "exact_haplotype_id": hid, "n_haplotypes": len(sub),
            "n_individuals": sub["individual"].nunique(), "n_groups": sub["group"].nunique(),
            "groups": ",".join(sorted(sub["group"].unique())),
            "individuals": ",".join(sorted(sub["individual"].unique())),
            "haplotypes": ",".join(sub["haplotype"].tolist()),
            **{f"n_{g}": n for g, n in group_counts.items()}
        })
    out = pd.DataFrame(rows).fillna(0)  # NaN -> 0 pour les groupes absents d'une ligne
    out = out.sort_values(["n_groups", "n_haplotypes"], ascending=[False, False])
    out.to_csv(out_tsv, sep="\t", index=False)


def plot_combined(region_name, H, positions, haplotypes, individuals, hap_index, groups,
                   group_order, genes, out_png, out_pdf, out_rows, out_sharing, out_snps,
                   out_ld_matrix, out_ld_summary):
    """Construit et sauvegarde la figure combinée (arbre + heatmap génotypes + triangle LD) d'une région.

    Écrit aussi les tables associées : ordre des lignes, partage
    d'haplotypes exact, SNP utilisés, matrice LD et résumé LD.

    Parameters
    ----------
    region_name : str
        Libellé de la région (titre de la figure).
    H : numpy.ndarray
        Matrice haplotypes x SNP (0/1/NaN), déjà restreinte/filtrée.
    positions : numpy.ndarray
        Positions (bp) des SNP, alignées sur les colonnes de H.
    haplotypes, individuals, hap_index, groups : numpy.ndarray
        Métadonnées alignées sur les lignes de H (voir load_haplotypes).
    group_order : list[str]
        Ordre d'affichage préféré des groupes dans la légende/bande couleur.
    genes : list[dict]
        Gènes d'intérêt à annoter (voir REGIONS / annotate_genes).
    out_png, out_pdf : pathlib.Path
        Chemins de sortie de la figure.
    out_rows, out_sharing, out_snps, out_ld_matrix, out_ld_summary : pathlib.Path
        Chemins de sortie des tables associées.
    """
    Z, order = make_tree_order(H)
    H_ord = H[order, :]
    M_ord = matrix_for_plot(H_ord)
    haplotypes_ord = haplotypes[order]
    individuals_ord = individuals[order]
    hap_index_ord = hap_index[order]
    groups_ord = groups[order]

    exact_ids, allele_strings = assign_exact_ids(H)
    exact_ids_ord = exact_ids[order]
    allele_strings_ord = allele_strings[order]

    r2, ld_positions = compute_r2(H, positions)  # LD recalculé sur les haplotypes affichés

    n_rows, n_snps = M_ord.shape
    n_ld = len(ld_positions)
    if n_ld != n_snps:
        # compute_r2 a pu retirer des SNP (variance nulle) : on resynchronise heatmap et LD
        pos_keep = set(ld_positions)
        idx_keep = [i for i, p in enumerate(positions) if p in pos_keep]
        H_ord = H_ord[:, idx_keep]
        M_ord = matrix_for_plot(H_ord)
        positions = positions[idx_keep]
        n_rows, n_snps = M_ord.shape

    fig_w = 34
    heat_h = max(9.0, min(31.0, n_rows * 0.070))  # hauteur heatmap proportionnelle au nb d'haplotypes, bornée
    ld_h = max(4.2, min(8.5, n_snps * 0.030))  # hauteur triangle LD proportionnelle au nb de SNP, bornée
    fig_h = heat_h + ld_h + 3.0

    fig = plt.figure(figsize=(fig_w, fig_h))
    # grille 3 lignes x 5 colonnes : heatmap/LD/colorbar en lignes, arbre/groupe/matrice/labels/légende en colonnes
    gs = GridSpec(3, 5, height_ratios=[heat_h, ld_h, 0.65],
                  width_ratios=[TREE_WIDTH, 0.45, 18.0, 5.0, 7.0],
                  hspace=0.00, wspace=0.045)

    ax_tree = fig.add_subplot(gs[0, 0])  # panneau arbre (dendrogramme)
    ax_group = fig.add_subplot(gs[0, 1])  # bande couleur = groupe d'origine
    ax_mat = fig.add_subplot(gs[0, 2])  # heatmap génotypes
    ax_lab = fig.add_subplot(gs[0, 3])  # étiquettes d'identifiants d'haplotypes
    ax_leg = fig.add_subplot(gs[0, 4])  # légende
    ax_ld = fig.add_subplot(gs[1, 2], sharex=ax_mat)  # triangle LD (même axe X que la heatmap)
    ax_cbar = fig.add_subplot(gs[2, 2])  # colorbar du LD

    dendrogram(Z, orientation="left", no_labels=True, color_threshold=0,
              above_threshold_color="black", ax=ax_tree)
    ax_tree.set_title("Arbre", fontsize=20, fontweight="bold")
    ax_tree.set_xticks([]); ax_tree.set_yticks([])
    ax_tree.set_xlabel("Distance 0/1", fontsize=10)

    present_groups = [g for g in group_order if g in groups_ord]
    for g in groups_ord:
        if g not in present_groups:
            present_groups.append(g)  # ajoute les groupes imprévus en fin de liste

    g_to_i = {g: i for i, g in enumerate(present_groups)}
    g_vals = np.array([g_to_i[g] for g in groups_ord]).reshape(-1, 1)
    g_cols = [group_color(g) for g in present_groups]

    # bande verticale colorée = groupe d'origine de chaque haplotype (même ordre que la heatmap)
    ax_group.imshow(g_vals, aspect="auto", interpolation="nearest",
                    cmap=ListedColormap(g_cols), origin="lower",
                    extent=[0, 1, 0, n_rows * 10])
    ax_group.set_xticks([]); ax_group.set_yticks([])
    ax_group.set_title("Origine", fontsize=20, fontweight="bold")

    cmap_h = ListedColormap(["#F7F7F7", "#D95F02", "#111111"])  # gris clair=REF, orange=ALT, noir=manquant
    norm_h = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap_h.N)  # associe les valeurs 0/1/2 aux 3 couleurs
    # heatmap principale : lignes = haplotypes (ordre arbre), colonnes = SNP
    ax_mat.imshow(M_ord, aspect="auto", interpolation="nearest", cmap=cmap_h,
                 norm=norm_h, origin="lower", extent=[0, n_snps, 0, n_rows * 10])
    title_pad = 42 if genes else 24  # plus de marge sous le titre si des gènes sont annotés au-dessus
    ax_mat.set_title(f"{region_name} | {n_rows} haplotypes ; {n_snps} SNPs",
                     fontsize=20, fontweight="bold", pad=title_pad)
    ax_mat.set_yticks([]); ax_mat.set_xlim(0, n_snps); ax_mat.set_xticks([])

    ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, n_rows * 10); ax_lab.axis("off")
    ax_lab.set_title("ID haplotype", fontsize=14, fontweight="bold")
    for i in range(n_rows):
        y = i * 10 + 5  # centre vertical de la ligne i
        lab = f"{exact_ids_ord[i]}  {individuals_ord[i]}_{hap_index_ord[i]}"  # ex: "H001  Ind1_h1"
        ax_lab.text(0, y, lab, va="center", ha="left", fontsize=LABEL_FONTSIZE,
                   color=group_color(groups_ord[i]))

    for ax in [ax_tree, ax_group, ax_mat, ax_lab]:
        ax.set_ylim(0, n_rows * 10)  # aligne verticalement les 4 panneaux

    pc = draw_ld_down_triangle(ax_ld, r2, positions)
    ax_ld.set_xlim(0, n_snps)

    annotate_genes(ax_mat, ax_ld, n_rows, n_snps, genes, positions)

    cbar = fig.colorbar(pc, cax=ax_cbar, orientation="horizontal")
    cbar.set_label("LD pairwise recalculé uniquement sur les haplotypes affichés — r² entre SNPs", fontsize=14)
    cbar.ax.tick_params(labelsize=12)

    ax_leg.axis("off")
    y = 0.96  # position verticale de départ (coordonnées axes 0-1)
    ax_leg.text(0, y, "Allèles SNPs", fontweight="bold", fontsize=18, transform=ax_leg.transAxes)
    y -= 0.08
    for color, label in [("#F7F7F7", "0 = REF"), ("#D95F02", "1 = ALT"), ("#111111", "missing")]:
        # carré de couleur + texte pour chaque catégorie d'allèle
        ax_leg.add_patch(plt.Rectangle((0.02, y - 0.022), 0.15, 0.050, facecolor=color,
                                       edgecolor="black", linewidth=1.0,
                                       transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.23, y, label, fontsize=14.0, va="center", transform=ax_leg.transAxes)
        y -= 0.075

    y -= 0.035
    ax_leg.text(0, y, "Groupes", fontweight="bold", fontsize=19, transform=ax_leg.transAxes)
    y -= 0.075
    for g in present_groups:
        n_hap = int(np.sum(groups_ord == g))
        n_ind = len(set(individuals_ord[groups_ord == g]))
        ax_leg.add_patch(plt.Rectangle((0.02, y - 0.022), 0.15, 0.050, facecolor=group_color(g),
                                       edgecolor="black", linewidth=1.0,
                                       transform=ax_leg.transAxes, clip_on=False))
        ax_leg.text(0.23, y, f"{display_label(g)} (n={n_ind} ind.; {n_hap} hap.)",
                   fontsize=14.0, va="center", transform=ax_leg.transAxes)
        y -= 0.064

    if genes:
        y -= 0.035
        ax_leg.text(0, y, "Gène(s) d'intérêt", fontweight="bold", fontsize=19, transform=ax_leg.transAxes)
        y -= 0.075
        for gene in genes:
            ax_leg.add_patch(plt.Rectangle((0.02, y - 0.022), 0.15, 0.050, facecolor=GENE_COLOR,
                                           edgecolor="black", linewidth=1.0, alpha=0.5,
                                           transform=ax_leg.transAxes, clip_on=False))
            ax_leg.text(0.23, y, f"{gene['name']} ({gene['start']/1e6:.3f}-{gene['end']/1e6:.3f}Mb)",
                       fontsize=13.0, va="center", transform=ax_leg.transAxes)
            y -= 0.064

    fig.subplots_adjust(left=0.018, right=0.992, top=0.940, bottom=0.055)
    plt.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02)
    plt.close()

    rows = pd.DataFrame({
        "plot_row_bottom_to_top_1based": np.arange(1, n_rows + 1),
        "exact_haplotype_id": exact_ids_ord, "haplotype": haplotypes_ord,
        "individual": individuals_ord, "hap_index": hap_index_ord,
        "group": groups_ord, "allele_string": allele_strings_ord,
    })
    rows.to_csv(out_rows, sep="\t", index=False)
    write_exact_sharing_table(H, haplotypes, individuals, hap_index, groups, out_sharing)
    pd.DataFrame({"snp_index_1based": np.arange(1, len(positions) + 1),
                  "position": positions, "position_Mb": positions / 1e6}).to_csv(out_snps, sep="\t", index=False)
    pd.DataFrame(r2, index=[str(x) for x in positions],
                columns=[str(x) for x in positions]).to_csv(out_ld_matrix, sep="\t")

    tri = np.tril_indices(len(positions), k=-1)  # indices du triangle inférieur (paires uniques, hors diagonale)
    vals = r2[tri]
    off = vals[np.isfinite(vals)]
    pd.DataFrame({
        "region": [region_name], "n_haplotypes_total": [H.shape[0]], "n_snps_used": [len(positions)],
        "mean_r2_offdiag": [float(np.nanmean(off))], "median_r2_offdiag": [float(np.nanmedian(off))],
        "prop_pairs_r2_ge_0.5": [float(np.mean(off >= 0.5))], "prop_pairs_r2_ge_0.8": [float(np.mean(off >= 0.8))],
        "note": ["LD recalculé sur les haplotypes affichés (Awassi+ME+P3_best) ; triangle inversé, même axe X que la matrice"]
    }).to_csv(out_ld_summary, sep="\t", index=False)


def run_region(reg):
    """Charge, filtre et trace la heatmap + arbre + LD pour une région candidate (voir REGIONS).

    Parameters
    ----------
    reg : dict
        Une entrée de REGIONS (chr, start, end, P3, label, vcf, genes).

    Returns
    -------
    list[pathlib.Path]
        Fichiers générés (figure PNG/PDF + tables associées), liste vide si
        le VCF phasé de la région est introuvable.
    """
    chrom, start, end, p3, label, vcf = reg["chr"], reg["start"], reg["end"], reg["P3"], reg["label"], reg["vcf"]
    genes = reg.get("genes", [])
    print(f"\n{'='*70}\n{label}  |  P3 = {p3}")

    if not Path(vcf).exists():
        print(f"  [!] VCF phasé introuvable : {vcf} — région ignorée")
        return []

    group_order = ["Awassi", "MiddleEastNonAwassi", p3]

    H, positions, haplotypes, individuals, hap_index, groups, selected, s2g = load_haplotypes(
        vcf, str(POP_DIR), group_order
    )
    H, positions = restrict_positions(H, positions, start, end)
    H, positions = filter_snps(H, positions, max_missing=MAX_MISSING, min_maf=MIN_MAF)

    start_mb, end_mb = start / 1e6, end / 1e6
    prefix = f"chr{chrom}_{start_mb:.3f}_{end_mb:.3f}Mb_{p3}"

    out_png = OUTDIR / f"{prefix}.png"
    out_pdf = OUTDIR / f"{prefix}.pdf"
    out_rows = OUTDIR / f"{prefix}_row_order.tsv"
    out_sharing = OUTDIR / f"{prefix}_haplotype_sharing.tsv"
    out_snps = OUTDIR / f"{prefix}_snps_used.tsv"
    out_ld_matrix = OUTDIR / f"{prefix}_LD_matrix_r2.tsv"
    out_ld_summary = OUTDIR / f"{prefix}_LD_summary.tsv"

    plot_combined(label, H, positions, haplotypes, individuals, hap_index, groups, group_order, genes,
                 out_png, out_pdf, out_rows, out_sharing, out_snps, out_ld_matrix, out_ld_summary)

    print(f"  {H.shape[0]} haplotypes ; {len(positions)} SNPs")
    print(f"  → {out_png.name}")

    return [out_png, out_pdf, out_rows, out_sharing, out_snps, out_ld_matrix, out_ld_summary]


if __name__ == "__main__":
    generated = []
    for reg in REGIONS:
        try:
            generated += run_region(reg)
        except Exception as e:
            print(f"  [!] Erreur région {reg['label']} : {e}")  # continue même en cas d'échec d'une région

    print(f"\nTerminé — fichiers dans {OUTDIR}")
    for f in generated:
        print(f"  {f}")
