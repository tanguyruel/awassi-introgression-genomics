#!/usr/bin/env python3
"""
FST SNP par SNP (per-SNP, vcftools --weir-fst-pop SANS fenêtre) + fd + LD, sur une région donnée.

Différence majeure vs v8 :
- Le panel FST n'utilise plus les fenêtres 20kb/step5kb chevauchantes.
  Il charge/calcule le FST Weir & Cockerham SNP par SNP (.weir.fst, pas
  .windowed.weir.fst) pour chaque comparaison Awassi vs {Africa, Asia,
  Europe, America, Australia, MiddleEastNonAwassi}, restreint aux SNPs
  de la région affichée.
- Pas d'interpolation, pas de moyenne de fenêtre : scatter par SNP +
  médiane glissante indicative (optionnelle, en tirets fins, non
  ajoutée à la légende) sur ROLLING_WINDOW_SNPS SNPs.
- Calcul : extraction régionale une fois par région (bcftools view -r,
  rapide car indexé) puis vcftools --weir-fst-pop sur ce petit VCF
  régional (pas de scan du chromosome entier). Résultats mis en cache
  dans analyses/fst/per_snp_regions_v9/<region>/ et réutilisés si déjà
  présents.
- fd et LD : inchangés par rapport à v8 (mêmes fonctions, même
  mécanisme d'alignement de largeur avec le triangle LD).
- xlim des 3 panels = bornes exactes de la région (start/end), plus le
  léger débordement du triangle LD nécessaire à l'alignement pixel.
Sorties  : analyses/LD/figures_aligned_v9/ — suffixe _v9_perSNPfst
À la fin : tentative d'ouverture des PNG, puis
           Enter = conserver | a+Enter = supprimer les nouveaux fichiers.

Script de référence pour "FST par SNP" (cf. Annexe A du rapport de stage) — la version
la plus élevée (v9) de la famille 13b_plot_fst_fd_ld_aligned_v1...v15 est la version
citée ; les autres numéros de version sont des étapes intermédiaires conservées telles
quelles (historique de mise au point), pas des alternatives à utiliser.
"""

# Imports : appels système (bcftools/vcftools), chemins, calcul, dataframes, graphiques
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # backend non interactif (génère les PNG/PDF sans écran)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import PolyCollection

SCRIPT_NAME = Path(__file__).name  # nom du script, affiché en filigrane sur les figures

# ── Régions A v10b ────────────────────────────────────────────────────────────
# Liste des régions candidates à tracer : coordonnées, groupe P3 au signal fd
# le plus fort, label d'affichage, chemins vers les matrices LD déjà calculées.
REGIONS = [
    {"chr": 3,  "start": 129220001, "end": 129260000, "P3": "Europe",
     "label": "chr3 129.220–129.260 Mb",
     "ld_matrix": "analyses/LD/LD_zones_geo_v1/chr3_129.220_129.260Mb_combined_LD_matrix_r2.tsv",
     "ld_snps":   "analyses/LD/LD_zones_geo_v1/chr3_129.220_129.260Mb_combined_LD_snps_used.tsv"},
    {"chr": 17, "start": 34160001,  "end": 34220000,  "P3": "Asia",
     "label": "chr17 34.160–34.220 Mb",
     "ld_matrix": "analyses/LD/LD_zones_geo_v1/chr17_34.165_34.220Mb/combined/chr17_34.165_34.220Mb_combined_LD_matrix_r2.tsv",
     "ld_snps":   "analyses/LD/LD_zones_geo_v1/chr17_34.165_34.220Mb/combined/chr17_34.165_34.220Mb_combined_LD_snps_used.tsv"},
    {"chr": 3,  "start": 186075001, "end": 186125000, "P3": "Asia",
     "label": "chr3 186.075–186.125 Mb",
     "ld_matrix": "analyses/LD/LD_zones_geo_v1/chr3_186.075_186.125Mb/combined/chr3_186.075_186.125Mb_combined_LD_matrix_r2.tsv",
     "ld_snps":   "analyses/LD/LD_zones_geo_v1/chr3_186.075_186.125Mb/combined/chr3_186.075_186.125Mb_combined_LD_snps_used.tsv"},
    {"chr": 20, "start": 765001,    "end": 845000,    "P3": "Europe",
     "label": "chr20 0.765–0.845 Mb",
     "ld_matrix": "analyses/LD/LD_zones_geo_v1/chr20_0.765_0.845Mb_combined_LD_matrix_r2.tsv",
     "ld_snps":   "analyses/LD/LD_zones_geo_v1/chr20_0.765_0.845Mb_combined_LD_snps_used.tsv"},
    {"chr": 6,  "start": 70240001,  "end": 70295000,  "P3": "Africa",
     "label": "chr6 70.240–70.295 Mb",
     "ld_matrix": "analyses/LD/LD_zones_geo_v1/chr6_70.240_70.295Mb_combined_LD_matrix_r2.tsv",
     "ld_snps":   "analyses/LD/LD_zones_geo_v1/chr6_70.240_70.295Mb_combined_LD_snps_used.tsv"},
    {"chr": 2,  "start": 112785001, "end": 112865000, "P3": "Australia",
     "label": "chr2 112.785–112.865 Mb",
     "ld_matrix": "analyses/LD/LD_zones_geo_v1/chr2_112.785_112.865Mb_combined_LD_matrix_r2.tsv",
     "ld_snps":   "analyses/LD/LD_zones_geo_v1/chr2_112.785_112.865Mb_combined_LD_snps_used.tsv"},
    {"chr": 5,  "start": 58070001,  "end": 58110000,  "P3": "Europe",
     "label": "chr5 58.070–58.110 Mb",
     "ld_matrix": "analyses/LD/LD_zones_geo_v1/chr5_58.070_58.110Mb_combined_LD_matrix_r2.tsv",
     "ld_snps":   "analyses/LD/LD_zones_geo_v1/chr5_58.070_58.110Mb_combined_LD_snps_used.tsv"},
]

VCF_DIR        = Path("analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/filtered_vcf")  # VCF filtrés par chromosome
POP_DIR        = Path("analyses/fst/popmaps_separees_v1")  # listes d'individus par groupe
FD_DIR         = Path("analyses/fd/fd_genomewide_20kb_step5kb_sep_groups_v1_01_07_15h10/results")  # fd genome-wide déjà calculé
FST_PERSNP_DIR = Path("analyses/fst/per_snp_regions_v9")  # cache des FST par SNP calculés à la volée
OUTDIR         = Path("analyses/LD/figures_aligned_v9")  # dossier de sortie des figures
OUTDIR.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si absent
FST_PERSNP_DIR.mkdir(parents=True, exist_ok=True)  # crée le dossier cache si absent

# Ordre d'affichage / légende (couleurs identiques aux scripts précédents)
FST_GROUPS = ["Africa", "Asia", "Europe", "America", "Australia", "MiddleEastNonAwassi"]  # groupes comparés à Awassi pour le FST
ALL_P3     = ["Africa", "Asia", "Australia", "America", "Europe"]  # groupes candidats P3 pour le fd (hors ME)
N_TICKS    = 8  # nombre de graduations sur l'axe des positions
LD_MAX_SNPS = 250  # nombre max de SNP affichés dans le triangle LD
ROLLING_WINDOW_SNPS = 25  # taille (en nb de SNP) de la médiane glissante indicative

GEO_COLORS = {  # couleur attribuée à chaque groupe géographique (légende commune aux 3 panels)
    "Africa":               "#e07b39",
    "Asia":                 "#6a5acd",
    "Australia":            "#2ca02c",
    "America":              "#d62728",
    "Europe":               "#8c564b",
    "MiddleEastNonAwassi":  "#1f77b4",
}
GEO_LABELS = {  # libellé court pour affichage
    "MiddleEastNonAwassi": "ME",
}

# palette blanc -> rouge foncé pour la couleur du r² dans le triangle LD
CMAP_LD = LinearSegmentedColormap.from_list(
    "ld_red", ["#ffffff", "#ffe0e0", "#ff9090", "#e00000", "#7a0000"]
)


def region_safe_name(chrom, start, end):
    return f"chr{chrom}_{start}_{end}"  # nom sûr (sans espaces/points) pour fichiers et dossiers


def extract_region_vcf(chrom, start, end, outdir):
    """bcftools view -r (indexé, rapide) : extrait une seule fois le VCF
    régional, réutilisé pour les 6 comparaisons de groupes."""
    region_vcf = outdir / "region.PASS_biallelic_snps.vcf.gz"  # chemin du VCF régional en cache
    if region_vcf.exists():
        return region_vcf, None  # déjà extrait -> réutilise le cache, pas de recalcul

    src = VCF_DIR / f"chr{chrom}.PASS_biallelic_snps.vcf.gz"  # VCF complet du chromosome
    if not src.exists():
        return None, f"VCF chromosome introuvable : {src}"

    # extrait uniquement la région d'intérêt (rapide car VCF source indexé)
    cmd_view = ["bcftools", "view", "-r", f"{chrom}:{start}-{end}",
                str(src), "-Oz", "-o", str(region_vcf)]
    r = subprocess.run(cmd_view, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0 or not region_vcf.exists():
        return None, f"bcftools view a échoué : {r.stderr[-300:]}"

    cmd_idx = ["bcftools", "index", "-t", str(region_vcf)]  # indexe le VCF régional (nécessaire pour vcftools/tabix)
    r = subprocess.run(cmd_idx, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        return None, f"bcftools index a échoué : {r.stderr[-300:]}"

    return region_vcf, None  # succès : chemin du VCF régional prêt à l'emploi


def run_persnp_fst(region_vcf, group, outdir):
    """vcftools --weir-fst-pop SANS --fst-window-size/--fst-window-step
    -> fichier .weir.fst SNP par SNP (pas de fenêtre)."""
    prefix = outdir / f"Awassi_vs_{group}"  # préfixe des fichiers produits par vcftools
    fst_file = Path(str(prefix) + ".weir.fst")  # fichier FST par SNP attendu
    if fst_file.exists():
        return fst_file, None  # déjà calculé -> réutilise le cache

    pop_awassi = POP_DIR / "Awassi.txt"  # liste des individus Awassi
    pop_group  = POP_DIR / f"{group}.txt"  # liste des individus du groupe comparé
    if not pop_awassi.exists() or not pop_group.exists():
        return None, f"popmap manquant ({pop_awassi} / {pop_group})"

    # commande vcftools --weir-fst-pop SANS fenêtre -> un FST par SNP (pas de moyenne)
    cmd = [
        "vcftools", "--gzvcf", str(region_vcf),
        "--weir-fst-pop", str(pop_awassi),
        "--weir-fst-pop", str(pop_group),
        "--out", str(prefix),
    ]
    log_file = Path(str(prefix) + ".vcftools.log")  # log complet de vcftools
    with open(log_file, "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)  # exécute vcftools, sortie redirigée dans le log

    if r.returncode != 0 or not fst_file.exists():
        return None, f"vcftools a échoué (voir {log_file})"

    return fst_file, None  # succès : chemin du fichier FST par SNP


def load_persnp_fst(fst_file):
    """Charge CHROM/POS/WEIR_AND_COCKERHAM_FST ; '-nan' -> NaN."""
    df = pd.read_csv(fst_file, sep="\t", na_values=["-nan", "nan", "NaN"])  # charge le fichier .weir.fst
    df.columns = [c.strip() for c in df.columns]  # nettoie les espaces dans les noms de colonnes
    fst_col = next(c for c in df.columns if "FST" in c.upper())  # trouve la colonne FST (nom variable selon vcftools)
    df = df.rename(columns={fst_col: "FST", "POS": "pos"})  # renomme en noms standard
    df["FST"] = pd.to_numeric(df["FST"], errors="coerce")  # force en numérique, valeurs invalides -> NaN
    return df[["pos", "FST"]]  # ne garde que position et FST


def load_fd_per_snp(chrom, snp_positions_bp):
    """Inchangé vs v8 : fd assigné aux positions SNP LD depuis la fenêtre 20kb."""
    fd_file = FD_DIR / f"fd_chr{chrom}_all_awassi_pairs_20kb_step5kb.tsv.gz"  # fichier fd du chromosome
    if not fd_file.exists():
        return pd.DataFrame(), 0  # pas de fd dispo pour ce chromosome
    fd = pd.read_csv(fd_file, sep="\t",
                     usecols=["chrom", "start", "end", "P1", "P2", "P3", "fd", "n_snps"])  # charge les colonnes utiles
    # filtre : P1=ME fixé, P2=Awassi, P3 parmi les groupes candidats, fd valide dans (0,1], >=10 SNP/fenêtre
    fd = fd[(fd["P1"] == "MiddleEastNonAwassi") &
            (fd["P2"] == "Awassi") &
            (fd["P3"].isin(ALL_P3)) &
            (fd["fd"] > 0) & (fd["fd"] <= 1) &
            (fd["n_snps"] >= 10)]
    if fd.empty:
        return pd.DataFrame(), 0

    rows = []  # accumulateur : une ligne (P3, position SNP, fd) par SNP couvert
    for p3, grp in fd.groupby("P3"):  # traite chaque groupe P3 séparément
        starts = grp["start"].values  # débuts des fenêtres fd
        ends   = grp["end"].values  # fins des fenêtres fd
        fds    = grp["fd"].values  # valeurs fd des fenêtres
        for pos in snp_positions_bp:  # pour chaque SNP du panel LD
            mask = (starts <= pos) & (ends >= pos)  # fenêtres fd qui couvrent ce SNP
            if mask.any():
                rows.append({"P3": p3, "pos_bp": pos, "fd": fds[mask].mean()})  # fd moyen des fenêtres couvrantes

    df = pd.DataFrame(rows)  # tableau fd par SNP et par P3
    n_snps_fd = df["pos_bp"].nunique() if not df.empty else 0  # nombre de SNP distincts avec un fd assigné
    return df, n_snps_fd


def load_ld_matrix(mat_path, snps_path):
    """Inchangé vs v8."""
    if not Path(mat_path).exists():
        return None, None  # matrice LD absente pour cette région
    mat = pd.read_csv(mat_path, sep="\t", index_col=0)  # charge la matrice r² (SNP x SNP)
    r2  = mat.values.astype(float)  # matrice numpy des r²
    np.fill_diagonal(r2, np.nan)  # diagonale (SNP vs lui-même) mise à NaN
    snps = pd.read_csv(snps_path, sep="\t")  # positions des SNP utilisés dans la matrice
    positions = snps["position"].values.astype(int)  # positions en bp
    if len(positions) > LD_MAX_SNPS:
        # sous-échantillonne à LD_MAX_SNPS SNP répartis régulièrement (lisibilité du triangle)
        idx = np.linspace(0, len(positions) - 1, LD_MAX_SNPS).round().astype(int)
        positions = positions[idx]
        r2 = r2[np.ix_(idx, idx)]
    return r2, positions


def mb_to_idx(mb_values, pos_mb):
    """Interpole Mb -> index SNP LD, avec extrapolation linéaire aux bords
    (pente locale) pour pouvoir caler l'axe X exactement sur les bornes de
    région (qui débordent légèrement des positions SNP LD extrêmes)."""
    mb_values = np.atleast_1d(np.asarray(mb_values, dtype=float))  # force en tableau numpy 1D
    idx = np.arange(len(pos_mb), dtype=float)  # indices 0..n-1 des positions SNP LD
    out = np.interp(mb_values, pos_mb, idx)  # interpolation linéaire Mb -> index SNP

    left_slope  = (idx[1] - idx[0]) / (pos_mb[1] - pos_mb[0])  # pente locale au premier SNP (extrapolation gauche)
    right_slope = (idx[-1] - idx[-2]) / (pos_mb[-1] - pos_mb[-2])  # pente locale au dernier SNP (extrapolation droite)

    left_mask  = mb_values < pos_mb[0]  # valeurs avant le premier SNP LD
    right_mask = mb_values > pos_mb[-1]  # valeurs après le dernier SNP LD
    out[left_mask]  = idx[0]  + (mb_values[left_mask]  - pos_mb[0])  * left_slope  # extrapole à gauche
    out[right_mask] = idx[-1] + (mb_values[right_mask] - pos_mb[-1]) * right_slope  # extrapole à droite
    return out


def make_figure(reg):
    # extraction des paramètres de la région à traiter
    chrom  = reg["chr"]
    start  = reg["start"]
    end    = reg["end"]
    p3     = reg["P3"]
    label  = reg["label"]

    print(f"\n{'='*70}\n{label}")

    region_safe = region_safe_name(chrom, start, end)  # nom sûr pour les chemins de cache
    persnp_dir  = FST_PERSNP_DIR / region_safe  # dossier cache spécifique à cette région
    persnp_dir.mkdir(parents=True, exist_ok=True)  # crée le dossier cache si absent

    vcf_chrom_path = VCF_DIR / f"chr{chrom}.PASS_biallelic_snps.vcf.gz"  # VCF du chromosome complet (pour info)
    print(f"VCF régional : {vcf_chrom_path}  (restreint à {chrom}:{start}-{end})")
    print("FST source   : per-SNP .weir.fst (vcftools --weir-fst-pop, sans --fst-window-size/--fst-window-step)")

    region_vcf, err = extract_region_vcf(chrom, start, end, persnp_dir)  # extrait (ou réutilise) le VCF régional
    if err:
        print(f"  [!] extraction VCF régional impossible : {err} — région ignorée")
        return None

    # ── LD (inchangé) ─────────────────────────────────────────────────────────
    r2, positions = load_ld_matrix(reg["ld_matrix"], reg["ld_snps"])  # matrice LD déjà calculée pour la région
    if r2 is None:
        print("  LD absent — skip")
        return None
    n_ld = len(positions)  # nombre de SNP LD
    pos_mb = positions / 1e6  # positions LD converties en Mb

    # ── FST per-SNP pour les 6 groupes ──────────────────────────────────────
    fst_series = {}  # FST par SNP, un DataFrame par groupe comparé à Awassi
    n_fst_unique = {}  # nombre de SNP avec FST valide, par groupe
    n_snps_region_total = None  # nombre total de SNP dans la région (fixé au 1er groupe traité)
    for grp in FST_GROUPS:  # calcule/charge le FST par SNP pour chaque groupe comparé
        fst_file, err = run_persnp_fst(region_vcf, grp, persnp_dir)  # lance (ou réutilise) vcftools
        if err:
            print(f"  [!] Awassi_vs_{grp} : {err}")
            fst_series[grp] = pd.DataFrame(columns=["pos", "FST"])  # groupe ignoré : série vide
            n_fst_unique[grp] = 0
            continue
        df = load_persnp_fst(fst_file)  # charge les FST par SNP calculés
        if n_snps_region_total is None:
            n_snps_region_total = len(df)  # nombre de SNP de référence (identique pour tous les groupes)
        fst_series[grp] = df
        n_valid = int(df["FST"].notna().sum())  # compte les FST non manquants
        n_fst_unique[grp] = n_valid

    print(f"  SNPs disponibles dans la région (vcftools) : {n_snps_region_total}")
    for grp in FST_GROUPS:  # affiche le nombre de SNP FST valides par groupe
        lbl = "ME" if grp == "MiddleEastNonAwassi" else grp
        print(f"  Awassi_vs_{lbl:10s} n={n_fst_unique.get(grp, 0)}")

    # ── fd (inchangé, positions = SNP LD) ────────────────────────────────────
    fd_data, n_snps_fd = load_fd_per_snp(chrom, positions)  # fd assigné aux positions des SNP LD
    print(f"  fd n={n_snps_fd}")
    print(f"  LD n={n_ld}")

    # ── Axe X : bornes exactes de la région pour les 3 panels ───────────────
    start_mb, end_mb = start / 1e6, end / 1e6  # bornes de la région en Mb
    xlim_mb  = (start_mb, end_mb)  # limites d'affichage des panels FST/fd
    tick_mb  = np.linspace(start_mb, end_mb, N_TICKS)  # graduations régulières entre les bornes
    tick_idx = mb_to_idx(tick_mb, pos_mb)  # graduations converties en index SNP LD (panel triangle)
    xlim_idx = (float(mb_to_idx([start_mb], pos_mb)[0]),
                float(mb_to_idx([end_mb],   pos_mb)[0]))  # bornes région converties en index SNP LD

    # ── Figure ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(15, 16))  # figure avec 3 panels empilés
    gs  = fig.add_gridspec(3, 1, height_ratios=[2.2, 2.2, 4.0], hspace=0.20)  # grille : FST, fd, triangle LD (plus haut)
    ax1 = fig.add_subplot(gs[0])  # panel FST
    ax2 = fig.add_subplot(gs[1])  # panel fd
    ax3 = fig.add_subplot(gs[2])  # panel triangle LD

    fig.suptitle(f"{label}  |  P3_best = {p3}", fontsize=13, fontweight="bold", y=0.997)  # titre général
    subtitle = (f"FST per-SNP : n uniques par groupe  |  "
                f"fd : n={n_snps_fd} affichés  |  LD : n={n_ld} utilisés")
    fig.text(0.5, 0.978, subtitle, ha="center", fontsize=9, color="#333333")  # sous-titre récapitulatif
    breakdown = ", ".join(f"{GEO_LABELS.get(g, g)}={n_fst_unique.get(g, 0)}" for g in FST_GROUPS)  # détail n par groupe
    fig.text(0.5, 0.969, f"FST per-SNP n unique : {breakdown}",
              ha="center", fontsize=7.5, color="#555555")

    fig.text(0.006, 0.5, SCRIPT_NAME, rotation=90, va="center", ha="left",
              fontsize=7, color="#888888", alpha=0.85)  # nom du script en filigrane, bord gauche

    # ── Panel FST : scatter SNP par SNP (+ médiane glissante indicative) ────
    for grp in FST_GROUPS:  # un nuage de points par groupe comparé à Awassi
        df = fst_series[grp]
        if df.empty:
            continue  # rien à tracer pour ce groupe
        sub = df.dropna(subset=["FST"]).copy()
        if sub.empty:
            continue
        sub["FST"] = sub["FST"].clip(lower=0)  # FST négatif ramené à 0 (valeurs négatives = bruit d'estimation)
        sub = sub.sort_values("pos")  # trie par position pour un tracé propre

        is_best = (grp == p3)  # ce groupe est-il le meilleur P3 de la région ?
        is_me   = (grp == "MiddleEastNonAwassi")  # ce groupe est-il le groupe de référence ME ?
        lbl     = GEO_LABELS.get(grp, grp) + (" ★" if is_best else "")  # étoile si groupe P3 retenu

        # nuage de points FST par SNP ; points plus gros/opaques et devant pour P3 et ME
        ax1.scatter(sub["pos"] / 1e6, sub["FST"],
                    color=GEO_COLORS[grp],
                    s=16 if (is_best or is_me) else 9,
                    alpha=0.55 if (is_best or is_me) else 0.35,
                    linewidths=0,
                    zorder=4 if (is_best or is_me) else 2,
                    label=lbl)

        if len(sub) >= ROLLING_WINDOW_SNPS:  # assez de SNP pour une médiane glissante
            roll = sub["FST"].rolling(ROLLING_WINDOW_SNPS, center=True,
                                      min_periods=max(5, ROLLING_WINDOW_SNPS // 3)).median()  # médiane glissante indicative
            ax1.plot(sub["pos"] / 1e6, roll,
                     color=GEO_COLORS[grp],
                     lw=1.6 if (is_best or is_me) else 0.9,
                     alpha=0.9 if (is_best or is_me) else 0.45,
                     ls="--", zorder=5, label="_nolegend_")  # tracé en tirets, hors légende

    ax1.axhline(0.15, color="black", lw=0.8, ls="--", alpha=0.45, label="seuil 0.15")  # seuil FST de référence
    for x in (start_mb, end_mb):
        ax1.axvline(x, color="dimgray", lw=0.9, ls=":", alpha=0.7, zorder=6)  # lignes verticales aux bornes de région
    ax1.set_xlim(xlim_mb)  # cale l'axe X sur les bornes de la région
    ax1.set_xticks(tick_mb)
    ax1.set_xticklabels([f"{x:.3f}" for x in tick_mb],
                        rotation=35, ha="right", fontsize=9)  # étiquettes des graduations en Mb
    ax1.set_ylim(bottom=0)  # FST toujours affiché à partir de 0
    ax1.set_ylabel("FST Awassi vs groupe\n(Weir & Cockerham, per-SNP)", fontsize=10)
    ax1.set_title("FST SNP par SNP — Awassi vs groupes géographiques  (pas de moyennage fenêtre)",
                  fontsize=9, loc="left", pad=3)
    ax1.text(0.004, 0.965,
             "tirets = médiane glissante indicative (25 SNPs) — pas un lissage de fenêtre génomique",
             transform=ax1.transAxes, fontsize=6.5, color="#777777", style="italic", va="top")  # note explicative
    ax1.legend(fontsize=8, loc="upper right", framealpha=0.85, ncol=3)
    for sp in ["top", "right"]: ax1.spines[sp].set_visible(False)  # retire le cadre haut/droite

    # ── Panel fd : tous P3 par SNP (inchangé vs v8) ──────────────────────────
    if fd_data.empty:
        ax2.text(0.5, 0.5, "fd non disponible", ha="center", va="center",
                 transform=ax2.transAxes, color="gray", fontsize=10)  # message si fd indisponible
    else:
        for grp_p3, df_p3 in fd_data.groupby("P3"):  # un nuage de points par groupe P3 candidat
            is_best = (grp_p3 == p3)  # ce P3 est-il le meilleur de la région ?
            df_p3   = df_p3.sort_values("pos_bp")  # trie par position
            # nuage de points fd par SNP ; points plus gros/opaques pour le meilleur P3
            ax2.scatter(df_p3["pos_bp"] / 1e6, df_p3["fd"],
                        color=GEO_COLORS.get(grp_p3, "#aaa"),
                        s=22 if is_best else 9,
                        alpha=1.0 if is_best else 0.50,
                        zorder=4 if is_best else 2,
                        linewidths=0,
                        label=f"P3={grp_p3}" + (" ★" if is_best else ""))
    ax2.axhline(0, color="black", lw=0.5, ls=":", alpha=0.4)  # ligne de référence fd=0
    for x in (start_mb, end_mb):
        ax2.axvline(x, color="dimgray", lw=0.9, ls=":", alpha=0.7, zorder=6)  # lignes verticales aux bornes de région
    ax2.set_xlim(xlim_mb)  # cale l'axe X sur les bornes de la région
    ax2.set_xticks(tick_mb)
    ax2.set_xticklabels([f"{x:.3f}" for x in tick_mb],
                        rotation=35, ha="right", fontsize=9)
    ax2.set_ylim(-0.05, 1.05)  # fd borné entre 0 et 1, légère marge d'affichage
    ax2.set_ylabel("fd par SNP (Martin 2015)\nP1=ME fixé, tous P3", fontsize=10)
    ax2.set_title(f"fd — n SNPs avec valeur = {n_snps_fd}  |  fd assigné depuis fenêtre 20kb",
                  fontsize=9, loc="left", pad=3)
    ax2.legend(fontsize=8, loc="upper right", framealpha=0.85, ncol=2)
    for sp in ["top", "right"]: ax2.spines[sp].set_visible(False)  # retire le cadre haut/droite

    # ── Panel LD triangle inversé — INCHANGÉ (aspect="equal" conservé) ───────
    verts, colors = [], []  # sommets des losanges et couleurs (r²) du triangle LD
    for i in range(n_ld):  # boucle sur chaque paire de SNP (i,j), i<j
        for j in range(i + 1, n_ld):
            cx = (i + j) / 2.0  # centre X du losange (milieu des deux SNP)
            cy = -(j - i) / 2.0  # centre Y du losange (distance entre les deux SNP, triangle inversé)
            verts.append([(cx - 0.5, cy), (cx, cy - 0.5),
                          (cx + 0.5, cy), (cx, cy + 0.5)])  # 4 sommets du losange représentant la paire
            v = r2[i, j]  # valeur r² de la paire
            colors.append(0.0 if np.isnan(v) else v)  # NaN affiché comme 0 (blanc)

    # collection de losanges colorés selon r² (plus rapide qu'un scatter individuel)
    pc = PolyCollection(verts, array=np.array(colors), cmap=CMAP_LD,
                        clim=(0, 1), linewidths=0)
    ax3.add_collection(pc)  # ajoute la collection au panel
    ax3.set_xlim(xlim_idx)  # cale l'axe X sur les bornes de région (en index SNP)
    ax3.set_ylim(-(n_ld / 2.0 + 0.45), 0.5)  # limite verticale du triangle
    ax3.set_aspect("equal")  # losanges non déformés (nécessaire pour l'alignement pixel)

    ax3.set_xticks(tick_idx)
    ax3.set_xticklabels([f"{x:.3f}" for x in tick_mb],
                        rotation=35, ha="right", fontsize=9)
    ax3.set_yticks([])  # pas de graduation verticale (axe sans signification directe)
    ax3.set_xlabel("Position génomique (Mb)", fontsize=11)
    ax3.set_title(f"LD combiné (Awassi + ME + {p3}) — n SNPs = {n_ld}",
                  fontsize=9, loc="left", pad=3)
    for sp in ["left", "right", "top"]: ax3.spines[sp].set_visible(False)  # retire ces cadres
    ax3.spines["bottom"].set_linewidth(0.8)  # cadre bas plus épais (axe X)
    for xi in (xlim_idx[0], xlim_idx[1]):
        pass  # bornes région = bornes xlim déjà -> pas de ligne verticale redondante
    idx_s = float(mb_to_idx([start_mb], pos_mb)[0])  # borne début région en index SNP
    idx_e = float(mb_to_idx([end_mb], pos_mb)[0])  # borne fin région en index SNP
    for xi in (idx_s, idx_e):
        ax3.axvline(xi, color="dimgray", lw=0.9, ls=":", alpha=0.7, zorder=5)  # lignes verticales aux bornes de région

    cbar = fig.colorbar(pc, ax=ax3, fraction=0.018, pad=0.01)  # barre de couleur du r²
    cbar.set_label("r²", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    # ── Alignement largeur : FST/fd calés sur la largeur RÉELLE de LD ───────
    box1_orig = ax1.get_position()  # position actuelle du panel FST (avant réalignement)
    box2_orig = ax2.get_position()  # position actuelle du panel fd
    fig.canvas.draw()  # nécessaire pour que apply_aspect ait fixé la boîte de ax3
    box3 = ax3.get_position()  # position réelle du panel LD après rendu (aspect equal appliqué)
    ax1.set_position([box3.x0, box1_orig.y0, box3.width, box1_orig.height])  # recale FST sur la largeur de LD
    ax2.set_position([box3.x0, box2_orig.y0, box3.width, box2_orig.height])  # recale fd sur la largeur de LD

    safe    = f"chr{chrom}_{start_mb:.3f}_{end_mb:.3f}Mb_{p3}"  # nom de fichier basé sur région et P3
    out_png = OUTDIR / f"{safe}_v9_perSNPfst.png"
    out_pdf = OUTDIR / f"{safe}_v9_perSNPfst.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")  # sauvegarde en PNG haute résolution
    plt.savefig(out_pdf,           bbox_inches="tight")  # sauvegarde en PDF (vectoriel)
    plt.close()  # ferme la figure (libère la mémoire)
    print(f"  → {out_png.name}")
    return out_png, out_pdf


if __name__ == "__main__":
    if subprocess.run(["bash", "-lc", "command -v vcftools"],
                       stdout=subprocess.DEVNULL).returncode != 0:
        raise SystemExit("ERREUR : vcftools introuvable")  # vérifie que vcftools est installé
    if subprocess.run(["bash", "-lc", "command -v bcftools"],
                       stdout=subprocess.DEVNULL).returncode != 0:
        raise SystemExit("ERREUR : bcftools introuvable")  # vérifie que bcftools est installé

    generated = []  # liste des fichiers produits (png + pdf)
    for reg in REGIONS:  # traite chaque région de la liste
        result = make_figure(reg)  # génère la figure pour cette région
        if result:
            generated.extend(result)  # ajoute les fichiers produits à la liste

    print(f"\nTerminé — {len(generated)//2} figures dans {OUTDIR}")
    print("Fichiers générés :")
    for f in generated:
        print(f"  {f}")

    for f in generated:  # ouvre chaque PNG généré pour vérification visuelle
        if f.suffix == ".png":
            subprocess.run(["bash", "-lc", f"xdg-open '{f}' >/dev/null 2>&1 || true"])  # tente l'ouverture système

    print("\nEntrée = conserver les fichiers  |  a + Entrée = supprimer les nouveaux fichiers")
    try:
        choice = input("> ").strip().lower()  # lit le choix utilisateur
    except EOFError:
        choice = ""  # pas de terminal interactif -> conserve par défaut

    if choice == "a":  # l'utilisateur a demandé la suppression
        for f in generated:
            Path(f).unlink(missing_ok=True)  # supprime les fichiers générés
        print("Fichiers supprimés.")
    else:
        print("Fichiers conservés.")
