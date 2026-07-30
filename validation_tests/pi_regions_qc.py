#!/usr/bin/env python3
"""
Pi sur les 9 régions candidates, fenêtres 20kb/pas 5kb, AVEC filtres QC explicites
et 1 plot par région (corrige une version antérieure du dépôt de travail : légende absente
sur la figure combinée, pas de filtre missingness, pas de flag qualité par palier de SNPs).

Filtres SNP appliqués (demande utilisateur 07/07) :
  - PASS uniquement (bcftools -f PASS)
  - Biallélique uniquement (-m2 -M2 -v snps)
  - PAS de filtre MAF strict (pas de seuil 0.05 — volontairement absent)
  - Missing par SNP raisonnable : max MAX_MISSING_SITE (15%, entre les 10-20%
    demandés) calculé sur l'ensemble des individus utilisés (les 7 groupes
    réunis) — site rejeté entièrement au-delà, pour tous les groupes.
  - Par groupe, au moins MIN_CALLED_FRAC_GROUP (75%, entre les 70-80% demandés)
    des individus du groupe doivent être génotypés à ce site, sinon ce site
    est traité comme manquant pour CE groupe seulement (site gardé pour les
    autres groupes si leur propre taux de génotypage est suffisant).

Seuils de fiabilité par fenêtre (par groupe, sur n_snps_used après QC) :
  - < 10 SNPs   : pi = NA, quality_flag = "NA_too_few_snps"
  - 10-19 SNPs  : pi gardé, quality_flag = "low_snp_count"
  - >= 20 SNPs  : quality_flag = "reliable"

Même astuce P3 que scripts 28/29 : pi calculé pour les 7 groupes à chaque
fenêtre, jamais de branche par région ; `is_P3_best`/`P3_best` = flag d'affichage.

Sortie :
  analyses/synthese_resultats/pi_9regions_20kb_step5kb_QC/Pi_9regions_20kb_step5kb_QC_by_group.tsv
  analyses/synthese_resultats/pi_9regions_20kb_step5kb_QC/plots/{region_id}_pi_20kb_step5kb.png/.pdf  (1 par région)

Script de référence pour "π avec contrôle qualité" (cf. Annexe A du rapport de stage) —
remplace les scripts 28 (valeur unique par région) et 29 (fenêtré sans QC explicite).
Usage : python3 30_pi_9regions_20kb_step5kb_QC_v1.py
"""

from pathlib import Path
from collections import defaultdict
import subprocess
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # racine du dépôt, pour importer _shared
from _shared import read_list, gt_to_alt_count, pi_from_values as pi_site_value

POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # dossier des popmaps par groupe
VCF_DIR = Path("data/raw data_08_06")  # dossier des VCF bruts par chromosome
OUTDIR = Path("analyses/synthese_resultats/pi_9regions_20kb_step5kb_QC")
PLOTDIR = OUTDIR / "plots"
OUTDIR.mkdir(parents=True, exist_ok=True)
PLOTDIR.mkdir(parents=True, exist_ok=True)

GROUPS = ["Awassi", "MiddleEastNonAwassi", "Africa", "Asia", "Europe", "America", "Australia"]

WINDOW = 20_000  # taille de fenêtre glissante (bp)
STEP = 5_000     # pas d'avancement entre deux fenêtres (bp)

MAX_MISSING_SITE = 0.15        # site rejeté si > 15% de missing tous groupes confondus
MIN_CALLED_FRAC_GROUP = 0.75   # site traité comme manquant POUR CE GROUPE si < 75% génotypés

N_SNPS_NA = 10        # < 10 : pi = NA, pas interprétable
N_SNPS_RELIABLE = 20  # >= 20 : fenêtre fiable ; entre les deux : "low_snp_count"

BASE_COLORS = {
    "Awassi":               "#E41A1C",
    "MiddleEastNonAwassi":  "#1f77b4",
    "Africa":               "#e07b39",
    "Asia":                 "#6a5acd",
    "Australia":            "#2ca02c",
    "America":              "#d62728",
    "Europe":               "#8c564b",
}

# Liste des 9 régions candidates : coordonnées et groupe P3 identifié comme meilleur donneur
REGIONS = [
    {"region_id": "chr2_112.8Mb_NIPA2_CYFIP1",   "chr": "2",  "start": 112785001, "end": 112865000, "P3_best": "Australia"},
    {"region_id": "chr3_129.2Mb_desert",           "chr": "3",  "start": 129220001, "end": 129260000, "P3_best": "Europe"},
    {"region_id": "chr3_137.4Mb_LOC101123547",     "chr": "3",  "start": 137390001, "end": 137420000, "P3_best": "America"},
    {"region_id": "chr3_186.1Mb_CCDC91",           "chr": "3",  "start": 186075001, "end": 186125000, "P3_best": "Asia"},
    {"region_id": "chr5_58.1Mb_ABLIM3_AFAP1L1",    "chr": "5",  "start": 58070001,  "end": 58110000,  "P3_best": "Europe"},
    {"region_id": "chr6_70.2Mb_KIT",               "chr": "6",  "start": 70240001,  "end": 70295000,  "P3_best": "Africa"},
    {"region_id": "chr10_49.4Mb_KLF12",            "chr": "10", "start": 49360001,  "end": 49390000,  "P3_best": "Asia"},
    {"region_id": "chr17_34.2Mb_SPATA5",           "chr": "17", "start": 34160001,  "end": 34220000,  "P3_best": "Asia"},
    {"region_id": "chr20_0.8Mb_KHDRBS2",           "chr": "20", "start": 765001,    "end": 845000,    "P3_best": "Europe"},
]


def quality_flag(n_snps):
    """Attribue un niveau de fiabilité à une fenêtre selon le nombre de SNPs utilisés."""
    if n_snps < N_SNPS_NA:
        return "NA_too_few_snps"
    if n_snps < N_SNPS_RELIABLE:
        return "low_snp_count"  # nombre de SNPs limité, pi gardé mais à interpréter avec prudence
    return "reliable"


def compute_region(reg):
    """Calcule pi en fenêtres glissantes (20kb/pas 5kb), pour les 7 groupes, sur une région candidate.

    Applique les filtres QC du module (missingness globale par site,
    couverture minimale par groupe) puis agrège par fenêtre avec un
    quality_flag (voir quality_flag).

    Parameters
    ----------
    reg : dict
        Une entrée de REGIONS (region_id, chr, start, end, P3_best).

    Returns
    -------
    list[dict]
        Une ligne par fenêtre x groupe.
    """
    region_id, chrom, start, end = reg["region_id"], reg["chr"], reg["start"], reg["end"]
    vcf = str(VCF_DIR / f"awassi_and_basedata_chr{chrom}.vcf.gz")
    print(f"\n{'='*70}\n{region_id} | chr{chrom}:{start}-{end}")

    groups = {g: read_list(POP_DIR / f"{g}.txt") for g in GROUPS}
    samples = sorted(set().union(*groups.values()))
    sample_file = OUTDIR / f"samples_{region_id}.txt"  # fichier temporaire pour bcftools -S
    sample_file.write_text("\n".join(samples) + "\n")

    order_txt = subprocess.check_output(
        ["bash", "-lc", f"bcftools view -S {sample_file} -Ou '{vcf}' | bcftools query -l"],
        text=True,
    )  # récupère l'ordre des échantillons dans le VCF sous-échantillonné
    order = order_txt.splitlines()
    idx_map = {s: i for i, s in enumerate(order)}
    group_idx = {g: np.array([idx_map[s] for s in groups[g] if s in idx_map], dtype=int) for g in GROUPS}
    n_total_group = {g: len(group_idx[g]) for g in GROUPS}
    n_total_all = len(order)

    cmd = (
        f"bcftools view -r {chrom}:{start}-{end} -f PASS -m2 -M2 -v snps "
        f"-S {sample_file} -Ou '{vcf}' | bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'"
    )  # filtre : région, PASS, biallélique, SNP uniquement, échantillons choisis ; extrait POS + GT
    txt = subprocess.check_output(["bash", "-lc", cmd], text=True)
    sample_file.unlink(missing_ok=True)

    positions = []  # positions des SNPs retenus après filtre missing global
    pi_by_group_site = defaultdict(list)
    n_read = 0
    n_dropped_missing_site = 0

    for line in txt.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        n_read += 1
        pos = int(parts[1])
        gts = parts[2:]
        alt_counts = np.array([gt_to_alt_count(gt) for gt in gts], dtype=float)

        missing_frac_all = np.isnan(alt_counts).sum() / n_total_all
        if missing_frac_all > MAX_MISSING_SITE:
            n_dropped_missing_site += 1
            continue

        positions.append(pos)
        for g in GROUPS:
            vals = alt_counts[group_idx[g]]
            vals = vals[~np.isnan(vals)]
            called_frac_g = len(vals) / n_total_group[g] if n_total_group[g] > 0 else 0.0
            if called_frac_g < MIN_CALLED_FRAC_GROUP:
                pi_by_group_site[g].append(np.nan)  # site traité comme manquant pour ce groupe seulement
            else:
                pi_by_group_site[g].append(pi_site_value(vals))

    positions = np.array(positions)
    for g in GROUPS:
        pi_by_group_site[g] = np.array(pi_by_group_site[g])

    print(f"  SNPs lus : {n_read} ; rejetés (missing site > {MAX_MISSING_SITE:.0%}) : {n_dropped_missing_site} ; retenus : {len(positions)}")

    rows = []
    for wstart in range(start, end - WINDOW + 2, STEP):
        wend = wstart + WINDOW - 1
        wmid = (wstart + wend) // 2
        in_win = (positions >= wstart) & (positions <= wend) if len(positions) else np.array([], dtype=bool)

        for g in GROUPS:
            if len(positions):
                vals = pi_by_group_site[g][in_win]
                vals = vals[np.isfinite(vals)]
            else:
                vals = np.array([])
            n_snps_used = len(vals)
            qflag = quality_flag(n_snps_used)
            pi_per_bp = vals.sum() / WINDOW if n_snps_used > 0 else np.nan  # pi moyen par paire de bases (somme / taille fenêtre)
            pi_per_kb = pi_per_bp * 1000 if np.isfinite(pi_per_bp) else np.nan
            if qflag == "NA_too_few_snps":
                pi_per_kb = np.nan  # non interprétable, cf. consigne utilisateur

            rows.append({
                "region_id": region_id, "chrom": chrom,
                "window_start": wstart, "window_end": wend, "window_mid": wmid,
                "group": g, "n_snps_used": n_snps_used, "quality_flag": qflag,
                "pi_per_kb": pi_per_kb,
                "is_P3_best": (g == reg["P3_best"]), "P3_best": reg["P3_best"],
            })

    return rows


def make_region_plot(df_region, region_id, p3, out_png, out_pdf):
    """Trace la courbe pi en fonction de la position pour les 7 groupes d'une région, sauvegarde PNG+PDF."""
    fig, ax = plt.subplots(figsize=(10, 6.5))

    for g in GROUPS:
        gsub = df_region[df_region["group"] == g].sort_values("window_mid")
        if gsub.empty:
            continue
        lw = 3.0 if g in ("Awassi", p3) else 1.4  # épaisseur de ligne accentuée pour Awassi et P3
        alpha = 1.0 if g in ("Awassi", p3) else 0.75
        x = gsub["window_mid"].to_numpy() / 1e6
        y = gsub["pi_per_kb"].to_numpy()

        ax.plot(x, y, color=BASE_COLORS[g], linewidth=lw, alpha=alpha, label=g, zorder=3)

        reliable = gsub["quality_flag"].eq("reliable").to_numpy()
        low = gsub["quality_flag"].eq("low_snp_count").to_numpy()
        ax.scatter(x[reliable], y[reliable], color=BASE_COLORS[g], s=28, zorder=4, edgecolor="none")
        ax.scatter(x[low], y[low], facecolor="none", edgecolor=BASE_COLORS[g], s=40, zorder=4,
                   linewidth=1.3, marker="o")

    ax.set_title(f"{region_id}  (P3 = {p3})\npi — fenêtres 20kb / pas 5kb (cercle plein = fiable ≥20 SNPs, "
                f"cercle vide = 10-19 SNPs, trou = <10 SNPs/NA)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Position (Mb)", fontsize=11)
    ax.set_ylabel("pi (x10⁻³ / kb)", fontsize=11)
    ax.legend(loc="best", fontsize=10, ncol=2, framealpha=0.9, title="Groupe")
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    all_rows = []
    for reg in REGIONS:
        rows = compute_region(reg)
        all_rows += rows

        df_region = pd.DataFrame(rows)
        p3 = reg["P3_best"]
        out_png = PLOTDIR / f"{reg['region_id']}_pi_20kb_step5kb.png"
        out_pdf = PLOTDIR / f"{reg['region_id']}_pi_20kb_step5kb.pdf"
        make_region_plot(df_region, reg["region_id"], p3, out_png, out_pdf)
        print(f"  Figure région écrite : {out_png}")

    df = pd.DataFrame(all_rows)
    out_tsv = OUTDIR / "Pi_9regions_20kb_step5kb_QC_by_group.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)

    n_na = (df["quality_flag"] == "NA_too_few_snps").sum()
    n_low = (df["quality_flag"] == "low_snp_count").sum()
    n_ok = (df["quality_flag"] == "reliable").sum()
    print(f"\nTable écrite : {out_tsv}")
    print(f"Lignes : reliable={n_ok}  low_snp_count={n_low}  NA_too_few_snps={n_na}  (total={len(df)})")
    print(f"9 figures écrites dans {PLOTDIR}/")
