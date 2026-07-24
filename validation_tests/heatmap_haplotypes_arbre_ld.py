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

# ── Imports : stdlib (fichiers, regex, appels externes), calcul (numpy/pandas),
# tracé (matplotlib), clustering hiérarchique (scipy) ──
import os
import re
import shlex
import subprocess
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

POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # dossier des popmaps (listes d'individus par groupe)
OUTDIR  = Path("analyses/haplotype_heatmap/Awassi_haplo/results/figures_finales/regions_A_9regions_v2")  # dossier de sortie figures/tables
OUTDIR.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie s'il n'existe pas

MAX_MISSING = 0.05  # taux max de génotypes manquants toléré par SNP
MIN_MAF     = 0.05  # fréquence allélique mineure minimale pour garder un SNP
LABEL_FONTSIZE = 2.6  # taille de police des étiquettes d'haplotypes
TREE_WIDTH     = 5.2  # largeur relative du panneau arbre dans la figure

GENE_COLOR = "#2ca02c"  # couleur verte des annotations de gène

# ── Régions (7 régions A v10b + 2 nouvelles) — VCF phasés déjà calculés ─────
# "genes" : gène(s) d'intérêt de la région, coords exactes GFF Oar_v4.0
# (vérifiées par requête directe dans le GFF, cf. AWASSI_AGENT_LOG.md 07/07).
REGIONS = [  # une entrée (dict) par région candidate à tracer
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
BASE_COLORS = {  # couleur associée à chaque groupe pour les graphiques
    "Awassi":               "#E41A1C",
    "MiddleEastNonAwassi":  "#1f77b4",
    "Africa":               "#e07b39",
    "Asia":                 "#6a5acd",
    "Australia":            "#2ca02c",
    "America":              "#d62728",
    "Europe":               "#8c564b",
}

DISPLAY_LABELS = {  # libellé affiché pour chaque groupe (légende de la figure)
    "Awassi": "Awassi",
    "MiddleEastNonAwassi": "ME (Moyen-Orient)",
    "Africa": "Africa",
    "Asia": "Asia",
    "Australia": "Australia",
    "America": "America",
    "Europe": "Europe",
}


def display_label(g):
    # renvoie le libellé d'affichage du groupe, ou nom brut (sans "_") si absent du dict
    return DISPLAY_LABELS.get(g, g.replace("_", " "))


def group_color(g):
    # renvoie la couleur du groupe, gris par défaut si groupe inconnu
    return BASE_COLORS.get(g, "#999999")


def run_cmd(cmd):
    # exécute une commande externe (liste ou chaîne shell) et renvoie sa sortie texte
    return subprocess.check_output(cmd, text=True)


def read_list(path):
    # lit un fichier texte, renvoie la liste des lignes non vides (espaces retirés)
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def safe_name(x):
    # remplace tout caractère non alphanumérique (hors ._-) par "_"
    x = re.sub(r"[^A-Za-z0-9._-]+", "_", x)
    # regroupe les underscores consécutifs en un seul
    x = re.sub(r"_+", "_", x)
    # retire les underscores en début/fin de chaîne
    return x.strip("_")


def parse_gt(gt_field):
    # extrait le champ GT (avant les ":") du champ VCF complet
    gt = gt_field.split(":")[0]
    if gt in ["./.", ".|."]:
        # génotype manquant pour les deux allèles
        return np.nan, np.nan
    # séparateur "|" si phasé, sinon "/" (non phasé)
    sep = "|" if "|" in gt else "/"
    # sépare les deux allèles du génotype
    parts = gt.split(sep)
    if len(parts) != 2:
        # format inattendu (pas exactement 2 allèles) -> considéré manquant
        return np.nan, np.nan

    def conv(x):
        # convertit un allèle "." en NaN, sinon en float (0 ou 1)
        if x == ".":
            return np.nan
        return float(int(x))

    # renvoie (allèle de l'haplotype 1, allèle de l'haplotype 2)
    return conv(parts[0]), conv(parts[1])


def load_haplotypes(vcf, pop_files_dir, groups):
    # liste des échantillons présents dans le VCF (bcftools query -l)
    vcf_samples = run_cmd(["bcftools", "query", "-l", vcf]).splitlines()

    sample_to_group = {}  # associe chaque échantillon à son groupe d'origine
    for g in groups:
        p = os.path.join(pop_files_dir, f"{g}.txt")  # fichier popmap du groupe g
        if not os.path.exists(p):
            print(f"Attention : pop file absent, ignoré : {p}")
            continue
        for s in read_list(p):
            if s not in sample_to_group:
                # le premier groupe rencontré est prioritaire (pas d'écrasement)
                sample_to_group[s] = g

    # échantillons du VCF qui appartiennent à un des groupes demandés
    selected_samples = [s for s in vcf_samples if s in sample_to_group]
    if len(selected_samples) == 0:
        raise RuntimeError("Aucun individu sélectionné dans le VCF.")

    # écrit la liste des échantillons sélectionnés dans un fichier temporaire
    # (nécessaire pour l'option -S de bcftools view)
    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        for s in selected_samples:
            tmp.write(s + "\n")
        sample_file = tmp.name

    try:
        q_vcf = shlex.quote(vcf)  # échappement shell du chemin du VCF
        q_sample = shlex.quote(sample_file)  # échappement shell du fichier d'échantillons
        cmd = [
            "bash", "-lc",
            # filtre : échantillons sélectionnés, SNPs bialléliques uniquement,
            # puis extrait CHROM/POS/GT de chaque individu
            f"bcftools view -S {q_sample} -m2 -M2 -v snps -Ou {q_vcf} | "
            "bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'"
        ]
        txt = run_cmd(cmd)  # sortie texte : une ligne par SNP
    finally:
        os.remove(sample_file)  # nettoyage du fichier temporaire

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
        positions.append(int(parts[1]))  # position du SNP (colonne POS)
        one_snp = []
        for gt in parts[2:]:
            # parcourt les GT de chaque individu et les éclate en 2 haplotypes
            a1, a2 = parse_gt(gt)
            one_snp.append(a1)
            one_snp.append(a2)
        rows_by_snp.append(one_snp)  # une ligne = un SNP, colonnes = haplotypes

    H = np.array(rows_by_snp, dtype=float).T  # transposée : lignes=haplotypes, colonnes=SNP
    positions = np.array(positions)

    # renvoie la matrice haplotypes x SNP + toutes les métadonnées associées
    return (H, positions, np.array(hap_names), np.array(hap_individuals),
            np.array(hap_index), np.array(hap_groups), selected_samples, sample_to_group)


def restrict_positions(H, positions, pos_min, pos_max):
    keep = np.ones(len(positions), dtype=bool)  # masque booléen, tout gardé au départ
    if pos_min is not None:
        keep &= positions >= pos_min  # exclut les SNP avant la borne basse
    if pos_max is not None:
        keep &= positions <= pos_max  # exclut les SNP après la borne haute
    return H[:, keep], positions[keep]  # sous-matrice et positions restreintes à la fenêtre


def filter_snps(H, positions, max_missing, min_maf):
    keep = []  # indices des SNP conservés
    for j in range(H.shape[1]):
        col = H[:, j]  # colonne = génotypes de tous les haplotypes pour ce SNP
        miss = np.mean(np.isnan(col))  # proportion de valeurs manquantes
        if miss > max_missing:
            continue  # trop de données manquantes -> SNP rejeté
        vals = col[~np.isnan(col)]  # valeurs non manquantes
        if len(vals) == 0:
            continue
        if len(np.unique(vals)) < 2:
            continue  # SNP non variable (monomorphe) -> rejeté
        n0 = np.sum(vals == 0)  # nb d'allèles REF
        n1 = np.sum(vals == 1)  # nb d'allèles ALT
        total = n0 + n1
        if total == 0:
            continue
        maf = min(n0, n1) / total  # fréquence de l'allèle mineur
        if maf < min_maf:
            continue  # MAF trop faible -> SNP rejeté
        keep.append(j)

    keep = np.array(keep, dtype=int)
    if len(keep) < 2:
        raise RuntimeError("Moins de 2 SNPs après filtre.")
    return H[:, keep], positions[keep]  # matrice et positions filtrées


def allele_string(row):
    # convertit une ligne d'allèles en chaîne "0/1/N" (N = manquant), sert d'identité exacte
    return "".join("N" if np.isnan(x) else str(int(x)) for x in row)


def assign_exact_ids(H):
    # chaîne d'allèles de chaque haplotype (identité exacte SNP par SNP)
    strings = [allele_string(H[i, :]) for i in range(H.shape[0])]
    str_to_id = {}  # associe chaque chaîne unique à un identifiant H001, H002, ...
    exact_ids = []
    for s in strings:
        if s not in str_to_id:
            str_to_id[s] = f"H{len(str_to_id) + 1:03d}"  # nouvel identifiant séquentiel
        exact_ids.append(str_to_id[s])
    return np.array(exact_ids), np.array(strings)


def make_tree_order(H):
    X = H.copy()
    X[np.isnan(X)] = 2  # code les valeurs manquantes en 2 (catégorie à part pour la distance)
    D = pdist(X, metric="hamming")  # distance de Hamming entre chaque paire d'haplotypes
    Z = linkage(D, method="average")  # clustering hiérarchique (average linkage / UPGMA)
    d = dendrogram(Z, no_plot=True)  # calcule l'ordre des feuilles sans tracer de figure
    return Z, np.array(d["leaves"], dtype=int)  # matrice de linkage + ordre des feuilles


def matrix_for_plot(H):
    M = H.copy()
    M[np.isnan(M)] = 2  # code le manquant en 2 pour un affichage à 3 couleurs (0/1/manquant)
    return M


def compute_r2(H, positions):
    # H : matrice haplotypes (lignes) x SNP (colonnes), valeurs 0/1/NaN
    X = H.copy()
    means = np.nanmean(X, axis=0)  # fréquence allélique moyenne de chaque SNP (NaN ignorés)
    inds = np.where(np.isnan(X))  # coordonnées des valeurs manquantes
    X[inds] = np.take(means, inds[1])  # impute le manquant par la moyenne du SNP correspondant
    X = X - X.mean(axis=0)  # centre chaque colonne (SNP) sur sa moyenne
    sd = X.std(axis=0, ddof=1)  # écart-type de chaque SNP (variance d'échantillon, ddof=1)
    valid = sd > 0  # exclut les SNP sans variance (devenus constants après imputation)
    X = X[:, valid]
    positions = positions[valid]
    sd = X.std(axis=0, ddof=1)  # recalcule l'écart-type sur les SNP restants uniquement
    X = X / sd  # réduit chaque colonne (SNP centré-réduit, moyenne 0 écart-type 1)
    n = X.shape[0]  # nombre d'haplotypes (observations)
    corr = (X.T @ X) / (n - 1)  # matrice de corrélation de Pearson entre SNP (produit croisé normalisé)
    corr = np.clip(corr, -1, 1)  # sécurité numérique : corrélation strictement bornée à [-1, 1]
    r2 = corr ** 2  # r² = carré du coefficient de corrélation de Pearson = mesure du LD
    np.fill_diagonal(r2, 1.0)  # r² d'un SNP avec lui-même fixé à 1 par convention
    return r2, positions  # matrice r² (SNP x SNP) et positions correspondantes (après filtrage variance)


def gene_index_span(positions, gene_start, gene_end):
    """Convertit les bornes bp d'un gène en indices SNP (espace de la heatmap),
    en clippant aux bornes de la fenêtre affichée. None si hors fenêtre."""
    lo = max(gene_start, int(positions[0]))  # borne basse, clippée au premier SNP affiché
    hi = min(gene_end, int(positions[-1]))  # borne haute, clippée au dernier SNP affiché
    if lo > hi:
        return None  # le gène est entièrement hors de la fenêtre affichée
    i0 = int(np.searchsorted(positions, lo, side="left"))  # indice SNP correspondant au début du gène
    i1 = int(np.searchsorted(positions, hi, side="right"))  # indice SNP correspondant à la fin du gène
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
        span = gene_index_span(positions, gene["start"], gene["end"])  # indices SNP occupés par le gène
        if span is None:
            continue  # gène hors fenêtre, rien à tracer
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
    n = len(positions)  # nombre de SNP affichés
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
            values.append(r2[i, j])  # valeur r² de la paire de SNP (i,j)

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

    tick_idx = np.linspace(0, n - 1, min(9, n)).astype(int)  # au plus 9 graduations, réparties régulièrement
    tick_lab = [f"{positions[i] / 1e6:.3f}" for i in tick_idx]  # labels en Mb
    ax.set_xticks(tick_idx + 0.5)
    ax.set_xticklabels(tick_lab, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("Position génomique (Mb)", fontsize=12)
    ax.set_yticks([])
    for spine in ["left", "right", "top"]:
        ax.spines[spine].set_visible(False)  # masque les bordures non utiles
    ax.spines["bottom"].set_linewidth(0.8)
    return pc  # renvoyé pour construire la colorbar dans plot_combined


def write_exact_sharing_table(H, haplotypes, individuals, hap_index, groups, out_tsv):
    exact_ids, allele_strings = assign_exact_ids(H)  # identifiant d'haplotype exact + séquence d'allèles
    df = pd.DataFrame({
        "exact_haplotype_id": exact_ids, "haplotype": haplotypes,
        "individual": individuals, "hap_index": hap_index,
        "group": groups, "allele_string": allele_strings,
    })
    rows = []
    for hid, sub in df.groupby("exact_haplotype_id", sort=False):
        # sub = tous les haplotypes partageant exactement le même identifiant (même séquence d'allèles)
        group_counts = sub["group"].value_counts().to_dict()  # nb d'occurrences par groupe
        rows.append({
            "exact_haplotype_id": hid, "n_haplotypes": len(sub),
            "n_individuals": sub["individual"].nunique(), "n_groups": sub["group"].nunique(),
            "groups": ",".join(sorted(sub["group"].unique())),
            "individuals": ",".join(sorted(sub["individual"].unique())),
            "haplotypes": ",".join(sub["haplotype"].tolist()),
            **{f"n_{g}": n for g, n in group_counts.items()}  # une colonne de comptage par groupe
        })
    out = pd.DataFrame(rows).fillna(0)  # NaN -> 0 pour les groupes absents d'une ligne
    out = out.sort_values(["n_groups", "n_haplotypes"], ascending=[False, False])  # les plus partagés en premier
    out.to_csv(out_tsv, sep="\t", index=False)


def plot_combined(region_name, H, positions, haplotypes, individuals, hap_index, groups,
                   group_order, genes, out_png, out_pdf, out_rows, out_sharing, out_snps,
                   out_ld_matrix, out_ld_summary):
    Z, order = make_tree_order(H)  # clustering hiérarchique + ordre des lignes selon l'arbre
    H_ord = H[order, :]  # matrice réordonnée selon l'arbre
    M_ord = matrix_for_plot(H_ord)  # version affichable (manquant codé en 2)
    haplotypes_ord = haplotypes[order]
    individuals_ord = individuals[order]
    hap_index_ord = hap_index[order]
    groups_ord = groups[order]

    exact_ids, allele_strings = assign_exact_ids(H)  # identifiants d'haplotypes exacts (ordre original)
    exact_ids_ord = exact_ids[order]  # réordonnés selon l'arbre
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

    fig_w = 34  # largeur de figure fixe (pouces)
    heat_h = max(9.0, min(31.0, n_rows * 0.070))  # hauteur heatmap proportionnelle au nb d'haplotypes, bornée
    ld_h = max(4.2, min(8.5, n_snps * 0.030))  # hauteur triangle LD proportionnelle au nb de SNP, bornée
    fig_h = heat_h + ld_h + 3.0  # hauteur totale de la figure

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

    # dessine le dendrogramme, feuilles orientées vers la gauche, sans labels de feuilles
    dendrogram(Z, orientation="left", no_labels=True, color_threshold=0,
              above_threshold_color="black", ax=ax_tree)
    ax_tree.set_title("Arbre", fontsize=20, fontweight="bold")
    ax_tree.set_xticks([]); ax_tree.set_yticks([])  # pas de graduations
    ax_tree.set_xlabel("Distance 0/1", fontsize=10)

    present_groups = [g for g in group_order if g in groups_ord]  # groupes présents, dans l'ordre attendu
    for g in groups_ord:
        if g not in present_groups:
            present_groups.append(g)  # ajoute les groupes imprévus en fin de liste

    g_to_i = {g: i for i, g in enumerate(present_groups)}  # index numérique par groupe
    g_vals = np.array([g_to_i[g] for g in groups_ord]).reshape(-1, 1)  # colonne d'index pour chaque haplotype
    g_cols = [group_color(g) for g in present_groups]  # couleurs correspondant à chaque groupe

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
                   color=group_color(groups_ord[i]))  # couleur du texte = groupe de l'haplotype

    for ax in [ax_tree, ax_group, ax_mat, ax_lab]:
        ax.set_ylim(0, n_rows * 10)  # aligne verticalement les 4 panneaux

    pc = draw_ld_down_triangle(ax_ld, r2, positions)  # trace le triangle LD
    ax_ld.set_xlim(0, n_snps)

    annotate_genes(ax_mat, ax_ld, n_rows, n_snps, genes, positions)  # ajoute les repères de gène(s)

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
        y -= 0.075  # descend pour l'entrée suivante

    y -= 0.035
    ax_leg.text(0, y, "Groupes", fontweight="bold", fontsize=19, transform=ax_leg.transAxes)
    y -= 0.075
    for g in present_groups:
        n_hap = int(np.sum(groups_ord == g))  # nb d'haplotypes de ce groupe affichés
        n_ind = len(set(individuals_ord[groups_ord == g]))  # nb d'individus distincts de ce groupe
        # carré de couleur + texte récapitulatif pour chaque groupe
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
            # carré de couleur + texte pour chaque gène annoté
            ax_leg.add_patch(plt.Rectangle((0.02, y - 0.022), 0.15, 0.050, facecolor=GENE_COLOR,
                                           edgecolor="black", linewidth=1.0, alpha=0.5,
                                           transform=ax_leg.transAxes, clip_on=False))
            ax_leg.text(0.23, y, f"{gene['name']} ({gene['start']/1e6:.3f}-{gene['end']/1e6:.3f}Mb)",
                       fontsize=13.0, va="center", transform=ax_leg.transAxes)
            y -= 0.064

    fig.subplots_adjust(left=0.018, right=0.992, top=0.940, bottom=0.055)  # ajuste les marges globales
    plt.savefig(out_png, dpi=300, bbox_inches="tight", pad_inches=0.02)  # export PNG haute résolution
    plt.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02)  # export PDF vectoriel
    plt.close()  # libère la mémoire de la figure

    rows = pd.DataFrame({
        "plot_row_bottom_to_top_1based": np.arange(1, n_rows + 1),
        "exact_haplotype_id": exact_ids_ord, "haplotype": haplotypes_ord,
        "individual": individuals_ord, "hap_index": hap_index_ord,
        "group": groups_ord, "allele_string": allele_strings_ord,
    })
    rows.to_csv(out_rows, sep="\t", index=False)  # table de correspondance ligne <-> haplotype
    write_exact_sharing_table(H, haplotypes, individuals, hap_index, groups, out_sharing)  # table de partage d'haplotypes exacts
    pd.DataFrame({"snp_index_1based": np.arange(1, len(positions) + 1),
                  "position": positions, "position_Mb": positions / 1e6}).to_csv(out_snps, sep="\t", index=False)  # liste des SNP utilisés
    pd.DataFrame(r2, index=[str(x) for x in positions],
                columns=[str(x) for x in positions]).to_csv(out_ld_matrix, sep="\t")  # matrice r² complète

    tri = np.tril_indices(len(positions), k=-1)  # indices du triangle inférieur (paires uniques, hors diagonale)
    vals = r2[tri]  # valeurs r² de toutes les paires de SNP
    off = vals[np.isfinite(vals)]  # exclut d'éventuelles valeurs non finies
    pd.DataFrame({
        "region": [region_name], "n_haplotypes_total": [H.shape[0]], "n_snps_used": [len(positions)],
        "mean_r2_offdiag": [float(np.nanmean(off))], "median_r2_offdiag": [float(np.nanmedian(off))],
        "prop_pairs_r2_ge_0.5": [float(np.mean(off >= 0.5))], "prop_pairs_r2_ge_0.8": [float(np.mean(off >= 0.8))],
        "note": ["LD recalculé sur les haplotypes affichés (Awassi+ME+P3_best) ; triangle inversé, même axe X que la matrice"]
    }).to_csv(out_ld_summary, sep="\t", index=False)  # résumé statistique du LD de la région


def run_region(reg):
    chrom, start, end, p3, label, vcf = reg["chr"], reg["start"], reg["end"], reg["P3"], reg["label"], reg["vcf"]  # déballe les paramètres de la région
    genes = reg.get("genes", [])  # liste des gènes à annoter (vide si absent)
    print(f"\n{'='*70}\n{label}  |  P3 = {p3}")

    if not Path(vcf).exists():
        print(f"  [!] VCF phasé introuvable : {vcf} — région ignorée")
        return []  # rien à générer si le VCF n'existe pas

    group_order = ["Awassi", "MiddleEastNonAwassi", p3]  # ordre d'affichage des groupes

    H, positions, haplotypes, individuals, hap_index, groups, selected, s2g = load_haplotypes(
        vcf, str(POP_DIR), group_order
    )  # charge la matrice haplotypes x SNP depuis le VCF phasé
    H, positions = restrict_positions(H, positions, start, end)  # restreint aux SNP de la fenêtre de la région
    H, positions = filter_snps(H, positions, max_missing=MAX_MISSING, min_maf=MIN_MAF)  # filtre qualité des SNP

    start_mb, end_mb = start / 1e6, end / 1e6
    prefix = f"chr{chrom}_{start_mb:.3f}_{end_mb:.3f}Mb_{p3}"  # préfixe commun des fichiers de sortie

    out_png = OUTDIR / f"{prefix}.png"
    out_pdf = OUTDIR / f"{prefix}.pdf"
    out_rows = OUTDIR / f"{prefix}_row_order.tsv"
    out_sharing = OUTDIR / f"{prefix}_haplotype_sharing.tsv"
    out_snps = OUTDIR / f"{prefix}_snps_used.tsv"
    out_ld_matrix = OUTDIR / f"{prefix}_LD_matrix_r2.tsv"
    out_ld_summary = OUTDIR / f"{prefix}_LD_summary.tsv"

    plot_combined(label, H, positions, haplotypes, individuals, hap_index, groups, group_order, genes,
                 out_png, out_pdf, out_rows, out_sharing, out_snps, out_ld_matrix, out_ld_summary)  # génère la figure et les tables

    print(f"  {H.shape[0]} haplotypes ; {len(positions)} SNPs")
    print(f"  → {out_png.name}")

    return [out_png, out_pdf, out_rows, out_sharing, out_snps, out_ld_matrix, out_ld_summary]  # fichiers générés pour cette région


if __name__ == "__main__":
    generated = []  # liste cumulée de tous les fichiers générés
    for reg in REGIONS:
        try:
            generated += run_region(reg)  # traite chaque région, une par une
        except Exception as e:
            print(f"  [!] Erreur région {reg['label']} : {e}")  # continue même en cas d'échec d'une région

    print(f"\nTerminé — fichiers dans {OUTDIR}")
    for f in generated:
        print(f"  {f}")  # liste finale des fichiers produits
