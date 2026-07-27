#!/usr/bin/env python3
"""
Sélection des 50 meilleures régions candidates + sélection stricte, à partir du
ranking FST et fd déjà calculés genome-wide (score géométrique combiné des deux).

Point de méthode important : filtre P1=MiddleEastNonAwassi lors du
chargement du fd. Tout le reste (FST source, critères stricts v9, fusion,
comptage SNPs) est identique.

Vérification en tête : fenêtres FST et fd toutes deux 20kb/5kb.

Sorties dans analyses/synthese_resultats/relecture_fst_fd/v10b_top50_P1_ME_strict/ :
  ALL_windows_FST_fd_ranked_v10b.tsv
  TOP50_candidates_introgression_v10b.tsv
  strict_regions_v10b.tsv
  strict_regions_v10b.html

Script de référence pour la sélection des régions candidates (méthode "score combiné",
cf. Annexe A du rapport de stage). Usage : python3 15_rebuild_v10b.py (aucun argument,
chemins en dur ci-dessous).
"""

import subprocess
from pathlib import Path
import pandas as pd
import numpy as np

# ── Chemins ───────────────────────────────────────────────────────────────────
BASE    = Path("analyses")
FST_IN  = BASE / "synthese_resultats/relecture_fst_fd/v7_FST20step5_sep_groups/FST20step5_sep_all_ranked_windows.tsv"  # FST par fenêtre déjà classé
FD_DIR  = BASE / "fd/fd_genomewide_20kb_step5kb_sep_groups_v1_01_07_15h10/results"  # fd genome-wide déjà calculé
VCF_DIR = BASE / "fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/filtered_vcf"  # VCF filtrés par chromosome (comptage SNPs)
OUT_DIR = BASE / "synthese_resultats/relecture_fst_fd/v10b_top50_P1_ME_strict"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Paramètres ────────────────────────────────────────────────────────────────
MIN_SNPS    = 20
TOP_WINDOWS = 1000  # nombre de meilleures fenêtres (score combiné) fusionnées en régions
TOP_REGIONS = 50
MERGE_GAP   = 250_000  # distance max (bp) pour fusionner deux régions strictes voisines
STRICT = dict(  # seuils du filtre de sélection stricte (cœur de la méthode)
    n_windows       = 5,     # au moins 5 fenêtres 20kb dans la région fusionnée
    FST_ME          = 0.15,  # FST Awassi vs Middle-East minimum
    FST_score       = 0.15,  # score FST (écart Awassi/ME vs groupe le plus proche) minimum
    FST_closest_max = 0.05,  # FST vs groupe géographique le plus proche : maximum toléré
    fd_max          = 0.85,  # fraction d'introgression fd minimum
    D_min           = 0.20,  # statistique D (signe de l'introgression) minimum
)

# ── 0. Vérification cohérence FST / fd (fenêtres 20kb, step 5kb) ─────────────
print("0. Vérification cohérence FST / fd...")

fst_check = pd.read_csv(FST_IN, sep="\t", usecols=["chr", "start", "end"], nrows=500)
win_fst  = int((fst_check["end"] - fst_check["start"]).median())
starts_f = fst_check[fst_check["chr"] == fst_check["chr"].iloc[0]]["start"].sort_values()
diffs    = starts_f.diff().dropna()
step_fst = int(diffs[diffs > 0].median())

fd_sample = pd.read_csv(
    sorted(FD_DIR.glob("fd_chr1_*.tsv.gz"))[0],
    sep="\t", usecols=["start", "end"], nrows=200
)
win_fd  = int((fd_sample["end"] - fd_sample["start"]).median())
uniq_starts = fd_sample["start"].drop_duplicates().sort_values()
step_fd = int(uniq_starts.diff().dropna().iloc[0]) if len(uniq_starts) > 1 else 0

print(f"   FST : fenêtre ~{win_fst+1} bp, step ~{step_fst} bp")
print(f"   fd  : fenêtre ~{win_fd+1} bp, step ~{step_fd} bp")

if abs(win_fst - win_fd) > 200 or abs(step_fst - step_fd) > 200:
    raise RuntimeError(
        f"INCOHÉRENCE : FST {win_fst+1}bp/{step_fst}bp ≠ fd {win_fd+1}bp/{step_fd}bp"
    )
print("   → Cohérent : FST et fd sur mêmes fenêtres 20kb/5kb ✅")

# ── 1. FST v7 ────────────────────────────────────────────────────────────────
print("\n1. Chargement FST v7...")
fst = pd.read_csv(FST_IN, sep="\t", usecols=[
    "chr", "start", "end",
    "FST_ME", "FST_score", "closest_group", "FST_closest_group", "rank_FST"
])
fst["chr"] = fst["chr"].astype(str)  # chr en texte pour une jointure fiable avec fd
print(f"   {len(fst):,} fenêtres FST")

# ── 2. fd : P1=ME fixé, P3 varie (Africa/Asia/Australia/America/Europe) ───────
print("\n2. Chargement fd 20kb/5kb (filtre P1=MiddleEastNonAwassi)...")
chunks = []  # accumulateur : une meilleure ligne fd par fenêtre, par chromosome
for gz in sorted(FD_DIR.glob("fd_chr*_all_awassi_pairs_20kb_step5kb.tsv.gz")):
    df = pd.read_csv(gz, sep="\t",
                     usecols=["chrom", "start", "end", "P1", "P2", "P3",
                               "D", "fd", "n_snps"])
    df = df[
        (df["P1"] == "MiddleEastNonAwassi") &   # ← correction v10b
        (df["P2"] == "Awassi") &
        (df["P3"] != "MiddleEastNonAwassi") &
        (df["D"]  >  0) &
        (df["fd"] >  0) &
        (df["fd"] <= 1) &
        (df["n_snps"] >= MIN_SNPS)
    ]  # filtre : trio P1=ME/P2=Awassi/P3≠ME, signe D positif, fd valide, assez de SNP
    if df.empty:
        continue
    best = df.loc[df.groupby(["chrom", "start", "end"])["fd"].idxmax()].copy()  # garde le P3 au fd max par fenêtre
    best = best.rename(columns={
        "chrom": "chr", "fd": "fd_max",
        "D": "D_at_max_fd", "P3": "P3_best", "n_snps": "n_snps_fd"
    })
    chunks.append(best[["chr", "start", "end", "fd_max", "D_at_max_fd", "P3_best", "n_snps_fd"]])

fd = pd.concat(chunks, ignore_index=True)
fd["chr"] = fd["chr"].astype(str)
print(f"   {len(fd):,} fenêtres fd valides (P1=ME, P3≠ME, D>0, fd∈(0,1], n_snps≥{MIN_SNPS})")

# ── 3. Jointure FST + fd ──────────────────────────────────────────────────────
print("\n3. Jointure (chr, start, end)...")
merged = fst.merge(fd, on=["chr", "start", "end"], how="inner")
N = len(merged)
print(f"   {N:,} fenêtres communes")

# ── 4. Scores ────────────────────────────────────────────────────────────────
print("\n4. Calcul scores...")
merged = merged.sort_values("fd_max", ascending=False).reset_index(drop=True)
merged["rank_fd"]  = merged.index + 1  # rang fd (1 = meilleur fd_max)
merged["rank_sum"] = merged["rank_FST"] + merged["rank_fd"]  # somme des rangs (indicatif)

# score combiné = moyenne géométrique des percentiles FST et fd (cœur de la méthode)
merged["pct_FST"] = (N + 1 - merged["rank_FST"]) / N  # percentile FST (1 = meilleur)
merged["pct_fd"]  = (N + 1 - merged["rank_fd"])  / N  # percentile fd (1 = meilleur)
merged["geom_pct"] = np.sqrt(
    merged["pct_FST"].clip(lower=0) * merged["pct_fd"].clip(lower=0)
)  # moyenne géométrique des 2 percentiles : score final de classement des fenêtres

fst_mean, fst_std = merged["FST_score"].mean(), merged["FST_score"].std()
fd_mean,  fd_std  = merged["fd_max"].mean(),    merged["fd_max"].std()
merged["z_FST"] = (merged["FST_score"] - fst_mean) / fst_std  # score Z du FST (indicatif)
merged["z_fd"]  = (merged["fd_max"]    - fd_mean)  / fd_std  # score Z du fd (indicatif)
merged["z_sum"] = merged["z_FST"] + merged["z_fd"]  # somme des Z-scores (indicatif, utilisé pour départager les fusions)

all_out = OUT_DIR / "ALL_windows_FST_fd_ranked_v10b.tsv"
merged.sort_values("geom_pct", ascending=False).rename(
    columns={"geom_pct": "score_geom_pct"}
).to_csv(all_out, sep="\t", index=False)
print(f"   → {all_out.name}")

# ── 5. Fusion TOP 1000 → régions ─────────────────────────────────────────────
print(f"\n5. Fusion des {TOP_WINDOWS} meilleures fenêtres...")
top_win = merged.sort_values("geom_pct", ascending=False).head(TOP_WINDOWS).copy()
top_win["chr_num"] = top_win["chr"].apply(lambda x: int(x) if x.isdigit() else 99)  # chr numérique pour trier (X/Y -> 99)
top_win = top_win.sort_values(["chr_num", "start"]).reset_index(drop=True)

regions = []
cur = None
for _, row in top_win.iterrows():  # parcourt les fenêtres triées par position, fusionne les fenêtres chevauchantes
    if cur is None:
        cur = row.to_dict(); cur["n_windows"] = 1  # première fenêtre : démarre une région
    elif row["chr"] == cur["chr"] and row["start"] <= cur["end"]:
        # fenêtre suivante chevauche la région en cours -> l'étend
        cur["end"] = max(cur["end"], row["end"])
        cur["n_windows"] += 1
        if row["geom_pct"] > cur["geom_pct"]:
            for k in ("geom_pct", "pct_FST", "pct_fd", "z_FST", "z_fd", "z_sum"):
                cur[k] = row[k]  # garde les valeurs du pic de score combiné de la région
        if row["FST_score"] > cur["FST_score"]:
            for k in ("FST_score", "FST_ME", "FST_closest_group",
                      "closest_group", "rank_FST"):
                cur[k] = row[k]  # garde les valeurs du pic FST de la région
        if row["fd_max"] > cur["fd_max"]:
            for k in ("fd_max", "D_at_max_fd", "P3_best", "rank_fd"):
                cur[k] = row[k]  # garde les valeurs du pic fd de la région
        cur["rank_sum"] = cur["rank_FST"] + cur["rank_fd"]  # recalcule la somme des rangs après mise à jour
    else:
        regions.append(cur); cur = row.to_dict(); cur["n_windows"] = 1  # pas de chevauchement -> nouvelle région
if cur:
    regions.append(cur)

reg = pd.DataFrame(regions)
reg = reg.sort_values("geom_pct", ascending=False).reset_index(drop=True)
reg["region_rank"]      = reg.index + 1  # rang final de la région (1 = meilleure)
reg["region_length_kb"] = ((reg["end"] - reg["start"]) / 1000).round(1)
reg = reg.rename(columns={"geom_pct": "score_geom_pct"})

cols = ["region_rank", "chr", "start", "end", "region_length_kb", "n_windows",
        "FST_ME", "FST_score", "closest_group", "FST_closest_group",
        "fd_max", "D_at_max_fd", "P3_best",
        "pct_FST", "pct_fd", "score_geom_pct",
        "rank_FST", "rank_fd", "rank_sum", "z_FST", "z_fd", "z_sum"]
reg = reg[[c for c in cols if c in reg.columns]]

top50 = reg.head(TOP_REGIONS)
out50 = OUT_DIR / "TOP50_candidates_introgression_v10b.tsv"
top50.to_csv(out50, sep="\t", index=False)
print(f"   {len(reg)} régions fusionnées → TOP50 : {out50.name}")

print("\nAperçu TOP 10 :")
view = ["region_rank", "chr", "start", "end", "n_windows",
        "FST_score", "fd_max", "P3_best", "score_geom_pct"]
print(top50[view].head(10).to_string(index=False))

# ── 6. Sélection stricte (critères v9) ───────────────────────────────────────
print("\n6. Sélection stricte (critères v9)...")
df = top50.copy()

# masque booléen : une région est "stricte" si elle passe les 6 seuils simultanément
strict_mask = (
    (df["n_windows"]         >= STRICT["n_windows"])       &
    (df["FST_ME"]            >= STRICT["FST_ME"])          &
    (df["FST_score"]         >= STRICT["FST_score"])       &
    (df["FST_closest_group"] <= STRICT["FST_closest_max"]) &
    (df["fd_max"]            >= STRICT["fd_max"])          &
    (df["D_at_max_fd"]       >  STRICT["D_min"])
)
strict     = df[strict_mask].copy()  # régions classe A (signal strict)
non_strict = df[~strict_mask].copy()  # régions classe B (au moins un critère raté)
print(f"   Strictes : {len(strict)}  |  Non strictes (B) : {len(non_strict)}")

# ── Fusion régions strictes proches (gap < 250 kb) ────────────────────────────
strict = strict.sort_values(["chr", "start"]).reset_index(drop=True)
fused = []
cur = None
for _, row in strict.iterrows():  # fusionne les régions strictes séparées de moins de MERGE_GAP bp
    if cur is None:
        cur = row.to_dict()
        cur["ranks_included"]  = [int(row["region_rank"])]  # rangs TOP50 d'origine inclus dans la fusion
        cur["n_windows_total"] = int(row["n_windows"])
    elif row["chr"] == cur["chr"] and (row["start"] - cur["end"]) < MERGE_GAP:
        # région suivante proche (gap < 250kb) -> fusionne avec la région en cours
        cur["end"] = max(cur["end"], row["end"])
        cur["ranks_included"].append(int(row["region_rank"]))
        cur["n_windows_total"] += int(row["n_windows"])
        cur["z_sum"] = max(cur["z_sum"], row["z_sum"])
        if row["FST_ME"]    > cur["FST_ME"]:    cur["FST_ME"]    = row["FST_ME"]
        if row["FST_score"] > cur["FST_score"]: cur["FST_score"] = row["FST_score"]
        if row["FST_closest_group"] < cur["FST_closest_group"]:
            cur["FST_closest_group"] = row["FST_closest_group"]  # garde le FST_closest_group min (plus strict)
            cur["closest_group"]     = row["closest_group"]
        if row["fd_max"] > cur["fd_max"]:
            cur["fd_max"]      = row["fd_max"]  # garde le fd_max le plus élevé
            cur["D_at_max_fd"] = row["D_at_max_fd"]
            cur["P3_best"]     = row["P3_best"]
    else:
        fused.append(cur)  # pas assez proche -> clôture la région en cours
        cur = row.to_dict()
        cur["ranks_included"]  = [int(row["region_rank"])]
        cur["n_windows_total"] = int(row["n_windows"])
if cur:
    fused.append(cur)

ms = pd.DataFrame(fused)
ms["ranks_included"]   = ms["ranks_included"].apply(lambda x: ",".join(map(str, x)))  # liste -> chaîne pour export
ms["region_length_kb"] = ((ms["end"] - ms["start"]) / 1000).round(1)
print(f"   Régions strictes après fusion : {len(ms)}")

# ── Comptage SNPs tabix ───────────────────────────────────────────────────────
def count_snps(chrom, start, end):
    vcf = VCF_DIR / f"chr{chrom}.PASS_biallelic_snps.vcf.gz"
    if not vcf.exists():
        return -1  # VCF absent -> comptage impossible
    try:
        r = subprocess.run(
            ["tabix", str(vcf), f"{chrom}:{start}-{end}"],
            capture_output=True, text=True, timeout=60
        )
        return r.stdout.count("\n")  # 1 ligne de sortie = 1 variant dans la région
    except Exception:
        return -1

print("   Comptage SNPs (tabix)...")
ms["n_SNPs_region"] = ms.apply(
    lambda r: count_snps(r["chr"], int(r["start"]), int(r["end"])), axis=1
)
ms["priorite"] = ms["n_SNPs_region"].apply(lambda n: "A" if n >= 50 else "B")  # A = signal strict + assez de SNP

def comment_strict(row):
    parts = []
    if row["n_SNPs_region"] < 50: parts.append("peu de SNPs")  # signale une région pauvre en SNP malgré le filtre strict
    if row["FST_ME"] >= 0.40:     parts.append("FST_ME très élevé")
    if row["fd_max"] >= 0.95:     parts.append("fd quasi-maximal")
    nb = len(row["ranks_included"].split(","))
    if nb > 1: parts.append(f"fusion de {nb} régions")
    return "; ".join(parts) if parts else "signal strict propre"

ms["commentaire"] = ms.apply(comment_strict, axis=1)

# Non-strictes : classe B
non_strict["n_SNPs_region"] = non_strict.apply(
    lambda r: count_snps(r["chr"], int(r["start"]), int(r["end"])), axis=1
)
non_strict["priorite"]        = "B"
non_strict["ranks_included"]  = non_strict["region_rank"].astype(str)  # pas de fusion, un seul rang d'origine
non_strict["n_windows_total"] = non_strict["n_windows"]
# commentaire : liste les critères stricts non atteints, pour chaque région
non_strict["commentaire"] = non_strict.apply(lambda r: "; ".join(filter(None, [
    f"n_windows<{STRICT['n_windows']}"         if r["n_windows"]         < STRICT["n_windows"]       else "",
    f"FST_ME<{STRICT['FST_ME']}"               if r["FST_ME"]            < STRICT["FST_ME"]          else "",
    f"FST_score<{STRICT['FST_score']}"         if r["FST_score"]         < STRICT["FST_score"]       else "",
    f"FST_closest>{STRICT['FST_closest_max']}" if r["FST_closest_group"] > STRICT["FST_closest_max"] else "",
    f"fd_max<{STRICT['fd_max']}"               if r["fd_max"]            < STRICT["fd_max"]          else "",
    f"D<{STRICT['D_min']}"                     if r["D_at_max_fd"]       <= STRICT["D_min"]          else "",
])), axis=1)

# ── 7. Assembler et sauvegarder ───────────────────────────────────────────────
# colonnes finales exportées, classes A (strictes fusionnées) et B (non strictes)
COLS = [
    "priorite", "chr", "start", "end", "region_length_kb",
    "ranks_included", "n_windows_total",
    "FST_ME", "FST_closest_group", "FST_score", "closest_group",
    "fd_max", "D_at_max_fd", "P3_best",
    "n_SNPs_region", "score_geom_pct", "z_sum", "commentaire"
]
final = pd.concat(
    [ms[[c for c in COLS if c in ms.columns]],
     non_strict[[c for c in COLS if c in non_strict.columns]]],
    ignore_index=True
).sort_values(["priorite", "score_geom_pct"], ascending=[True, False]).reset_index(drop=True)  # A avant B, puis par score

tsv_out = OUT_DIR / "strict_regions_v10b.tsv"
final.to_csv(tsv_out, sep="\t", index=False, float_format="%.4f")
print(f"\n→ {tsv_out.name}")

# ── 8. Export HTML ────────────────────────────────────────────────────────────
def render_html(df):
    rows_html = []
    for _, row in df.iterrows():
        cells = []
        for col in df.columns:
            val = row[col]
            if col == "commentaire" and row["priorite"] == "B":
                cells.append(f'<td style="color:red">{val}</td>')
            elif col == "priorite" and val == "B":
                cells.append(f'<td style="color:red;font-weight:bold">{val}</td>')
            elif isinstance(val, float):
                cells.append(f"<td>{val:.4f}</td>")
            else:
                cells.append(f"<td>{val}</td>")
        rows_html.append("<tr>" + "".join(cells) + "</tr>")
    headers = "".join(f"<th>{c}</th>" for c in df.columns)
    return f"""<!DOCTYPE html>
<html lang="fr"><head><meta charset="utf-8">
<title>Sélection stricte v10b</title>
<style>
  body {{ font-family: monospace; font-size: 12px; padding: 20px; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th {{ background: #3a3a3a; color: white; padding: 6px 8px; text-align: left; }}
  td {{ padding: 4px 8px; border-bottom: 1px solid #ddd; }}
  tr:nth-child(even) {{ background: #f5f5f5; }}
  tr:hover {{ background: #e8f4fd; }}
</style></head><body>
<h2>Sélection stricte v10b — correction P1=MiddleEastNonAwassi</h2>
<p>A = signal strict + ≥50 SNPs &nbsp;|&nbsp;
<span style="color:red">B = critère non atteint (commentaire en rouge)</span></p>
<table><thead><tr>{headers}</tr></thead>
<tbody>{"".join(rows_html)}</tbody></table>
</body></html>"""

html_out = OUT_DIR / "strict_regions_v10b.html"
html_out.write_text(render_html(final), encoding="utf-8")
print(f"→ {html_out.name}")

print(f"\nClasse A : {(final['priorite']=='A').sum()}")
print(f"Classe B : {(final['priorite']=='B').sum()}")
print(f"\nTerminé → {OUT_DIR}")
