#!/usr/bin/env python3
"""
Script 23 — Manhattan plots FST et fd, genome entier (25 chromosomes),
panels empilés et alignés sur le même axe X (position cumulée), pour
regarder séparément où chaque statistique ressort fort — au lieu du score
combiné (moyenne géométrique) du script 15_rebuild_v10b.py.

Demande de François (via l'utilisateur) : identifier les régions qui
ressortent indépendamment dans FST ET dans fd. Méthode proche de Hu et al.
2019 (Manhattan plot FST fig. 6, Manhattan plot fd fig. 5A, seuils en
pointillés).

v2 — chr25 était coupé (xlim trop court) -> corrigé ; étiquettes de région
passées à l'horizontale (avant : rotation 90°, illisible), remplacées par le
nom du gène réel annoté depuis le GFF Oar_v4.0 (au lieu du simple nom de
région) — mêmes gènes que script 21, étendu aux 25 chromosomes.

v3 — retour demandé à 1 point rouge par FENÊTRE brute de 20kb (pas 1 par
région fusionnée) : le chevauchement visuel (fenêtres qui se recouvrent tous
les 5kb) n'est pas un problème, c'est la réalité du pas de la fenêtre
glissante — gardé tel quel, juste indiqué explicitement sur le plot
(fenêtres 20kb / pas 5kb). Les régions fusionnées ne servent plus qu'à
positionner UNE étiquette de gène par région (pas à agréger les points).

v4 — correction du filtre "strict" : ne vérifiait que 2 critères (FST_ME et
fd_max), alors que le filtre STRICT de référence (script 15) en a 5 (hors
n_windows, volontairement omis ici) : FST_ME≥0.15, FST_score≥0.15,
FST_closest_group≤0.05, fd_max≥0.85, D_at_max_fd>0.20. Avec seulement 2
critères, chr13:37Mb (CSRP2BP) apparaissait à tort comme strict alors qu'elle
échoue sur FST_score (0.1499) et D_at_max_fd (0.127). Corrigé : les 5
critères complets sont maintenant appliqués → 9 régions strictes (pas 10).

Entrées (déjà construites, script 22, rien recalculé depuis les VCF) :
  analyses/synthese_resultats/relecture_fst_fd/genomewide_separate_rank/
    FST_all_windows_ranked_genomewide.tsv  (267 067 fenêtres)
    fd_all_windows_ranked_genomewide.tsv   (317 729 fenêtres valides)
  data/reference/Oar_v4.0/GCF_000298735.2_Oar_v4.0_genomic.gff.gz (annotation gènes)

Sortie : analyses/synthese_resultats/manhattan_fst_fd_v1/
  manhattan_fst_fd_genomewide.png / .pdf
  regions_top_intersection_K{K}.tsv (régions fusionnées de l'intersection,
  avec gène annoté, pour retrouver le détail derrière chaque marqueur rouge)
"""

# Imports : lecture GFF compressé, regex, chemins, dataframes, calcul, graphiques
import gzip
import re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # backend non interactif (génère les PNG/PDF sans écran)
import matplotlib.pyplot as plt

BASE   = Path("analyses/synthese_resultats/relecture_fst_fd/genomewide_separate_rank")  # tables FST/fd genome-wide (script 22)
OUTDIR = Path("analyses/synthese_resultats/manhattan_fst_fd_v1")  # dossier de sortie
OUTDIR.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si absent
GFF = Path("data/reference/Oar_v4.0/GCF_000298735.2_Oar_v4.0_genomic.gff.gz")  # annotation des gènes

K = 2000       # top-K par statistique, séparément (stable à partir de K~1500-2000, voir échange précédent)
FLANK = 50_000  # marge pour associer un gène flanquant si rien dans la région

CHR_COLORS = ["#3b5b92", "#8faadc"]  # alternance 2 couleurs par chromosome


def load_data():
    fst = pd.read_csv(BASE / "FST_all_windows_ranked_genomewide.tsv", sep="\t",
                       usecols=["chr", "start", "end", "FST_ME", "FST_score", "FST_closest_group", "rank_FST"])  # FST classé genome-wide
    fd = pd.read_csv(BASE / "fd_all_windows_ranked_genomewide.tsv", sep="\t",
                      usecols=["chr", "start", "end", "fd_max", "D_at_max_fd", "rank_fd", "P3_best"])  # fd classé genome-wide
    fst["chr"] = fst["chr"].astype(int)  # chr en entier pour trier/comparer
    fd["chr"] = fd["chr"].astype(int)
    return fst, fd


def build_chrom_offsets(fst, fd):
    """Offset cumulé par chromosome, basé sur la position max observée dans
    les deux jeux de données (pour que les 2 panels restent alignés)."""
    max_end = pd.concat([
        fst.groupby("chr")["end"].max(),
        fd.groupby("chr")["end"].max(),
    ], axis=1).max(axis=1)  # position max par chromosome, sur FST et fd réunis
    chroms = sorted(max_end.index)  # liste des chromosomes triée
    offsets = {}  # décalage cumulé (bp) à ajouter à chaque chromosome pour l'axe X commun
    cum = 0
    for c in chroms:  # construit l'offset cumulé chromosome par chromosome
        offsets[c] = cum
        cum += int(max_end[c]) + 1_000_000  # marge entre chromosomes
    return offsets, chroms, max_end


def add_cum_pos(df, offsets):
    df = df.copy()
    df["cum_pos"] = df["chr"].map(offsets) + (df["start"] + df["end"]) / 2  # position cumulée = offset + milieu de fenêtre
    return df


def load_genes_by_chr():
    """Genes du GFF Oar_v4.0 pour les 25 chromosomes (mêmes accessions NC_
    que script 21, étendues : NC_019458.2=chr1 ... NC_019482.2=chr25)."""
    chr_map = {f"NC_0194{58 + i}.2": str(i + 1) for i in range(25)}  # accession NCBI -> numéro de chromosome (1-25)
    genes_by_chr = {c: [] for c in chr_map.values()}  # liste de gènes (start,end,name) par chromosome
    with gzip.open(GFF, "rt") as f:  # lit le GFF compressé ligne par ligne (fichier volumineux)
        for line in f:
            if line.startswith("#"):
                continue  # ignore les lignes d'en-tête/commentaire GFF
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue  # ne garde que les entités de type "gene"
            rseq = parts[0]
            if rseq not in chr_map:
                continue  # ignore les scaffolds hors chromosomes 1-25
            chrom = chr_map[rseq]
            start, end = int(parts[3]), int(parts[4])
            m_name = re.search(r"(?:^|;)Name=([^;]+)", parts[8])  # extrait l'attribut Name= de la colonne 9
            name = m_name.group(1) if m_name else "?"
            genes_by_chr[chrom].append((start, end, name))  # ajoute le gène à la liste du chromosome
    return genes_by_chr


def best_gene(genes_by_chr, chrom, start, end):
    chrom = str(chrom)
    # gènes situés dans la région ou dans la marge FLANK autour (candidats)
    hits = [g for g in genes_by_chr.get(chrom, []) if g[1] >= start - FLANK and g[0] <= end + FLANK]
    if not hits:
        return "(désert génique)"  # aucun gène, même avec la marge
    in_region = [g for g in hits if g[0] <= end and g[1] >= start]  # gènes réellement dans la région (sans marge)
    named = [g for g in in_region if not g[2].startswith("LOC")]  # gènes avec un nom réel (pas un LOC générique)
    if named:
        return "+".join(sorted(set(g[2] for g in named)))  # combine les noms si plusieurs gènes nommés
    if in_region:
        return in_region[0][2]  # sinon renvoie le premier gène (même LOC) trouvé dans la région
    hits.sort(key=lambda g: min(abs(g[0] - end), abs(g[1] - start)))  # trie par distance au bord de région
    return hits[0][2] + " (flanc)"  # aucun gène dans la région : renvoie le plus proche, marqué "flanc"


def fuse_intersection(fst, fd, inter_keys):
    """Fusionne les fenêtres en intersection en régions ; garde le pic
    FST_ME/FST_score (max) et le P3_best associé au pic fd_max (pas juste la
    première fenêtre rencontrée, contrairement au 1er jet)."""
    sub = fst.merge(fd, on=["chr", "start", "end"], how="inner")  # fenêtres présentes dans FST et fd
    sub = sub[sub.apply(lambda r: (r["chr"], r["start"]) in inter_keys, axis=1)]  # ne garde que l'intersection top-K/top-K
    sub = sub.sort_values(["chr", "start"]).reset_index(drop=True)  # ordre génomique

    regions = []  # régions fusionnées de l'intersection
    cur = None  # région en cours de construction
    for _, row in sub.iterrows():  # fusionne les fenêtres chevauchantes en régions
        if cur is None or not (row["chr"] == cur["chr"] and row["start"] <= cur["end"]):
            if cur is not None:
                regions.append(cur)  # clôture la région précédente
            cur = {"chr": row["chr"], "start": row["start"], "end": row["end"], "n_windows": 1,
                   "FST_ME": row["FST_ME"], "FST_score": row["FST_score"],
                   "FST_closest_group": row["FST_closest_group"],
                   "fd_max": row["fd_max"], "D_at_max_fd": row["D_at_max_fd"], "P3_best": row["P3_best"]}  # nouvelle région
        else:
            cur["end"] = max(cur["end"], row["end"])  # étend la région en cours
            cur["n_windows"] += 1
            if row["FST_score"] > cur["FST_score"]:
                cur["FST_score"] = row["FST_score"]  # garde le pic FST_score de la région
                cur["FST_ME"] = row["FST_ME"]
                cur["FST_closest_group"] = row["FST_closest_group"]
            if row["fd_max"] > cur["fd_max"]:
                cur["fd_max"] = row["fd_max"]  # garde le pic fd_max de la région, et son P3 associé
                cur["D_at_max_fd"] = row["D_at_max_fd"]
                cur["P3_best"] = row["P3_best"]
    if cur is not None:
        regions.append(cur)  # ajoute la dernière région en cours
    return pd.DataFrame(regions)


def main():
    print("Chargement des tables genome-wide (déjà construites, script 22)...")
    fst, fd = load_data()
    print(f"  FST : {len(fst):,} fenêtres")
    print(f"  fd  : {len(fd):,} fenêtres")

    offsets, chroms, max_end = build_chrom_offsets(fst, fd)  # décalage cumulé par chromosome (axe X commun)
    fst = add_cum_pos(fst, offsets)  # ajoute la position cumulée à chaque fenêtre FST
    fd  = add_cum_pos(fd,  offsets)  # ajoute la position cumulée à chaque fenêtre fd

    thresh_fst = fst["FST_score"].quantile(0.99)  # seuil FST du top 1% (repère visuel)
    thresh_fd  = fd["fd_max"].quantile(0.99)  # seuil fd du top 1% (repère visuel)
    print(f"  Seuil top 1% FST_score : {thresh_fst:.4f}")
    print(f"  Seuil top 1% fd_max    : {thresh_fd:.4f}")

    top_fst = fst.nsmallest(K, "rank_FST")  # les K meilleures fenêtres FST (rang le plus petit = meilleur)
    top_fd  = fd.nsmallest(K, "rank_fd")  # les K meilleures fenêtres fd
    inter_keys = set(zip(top_fst["chr"], top_fst["start"])) & set(zip(top_fd["chr"], top_fd["start"]))  # fenêtres fortes dans les 2 statistiques
    print(f"  Top {K} FST seul ∩ top {K} fd seul : {len(inter_keys)} fenêtres (avant fusion)")

    fst["is_top_both"] = list(zip(fst["chr"], fst["start"]))
    fst["is_top_both"] = fst["is_top_both"].isin(inter_keys)  # marque les fenêtres FST dans l'intersection
    fd["is_top_both"] = list(zip(fd["chr"], fd["start"]))
    fd["is_top_both"] = fd["is_top_both"].isin(inter_keys)  # marque les fenêtres fd dans l'intersection

    reg = fuse_intersection(fst, fd, inter_keys)  # fusionne les fenêtres de l'intersection en régions
    reg["cum_pos"] = reg["chr"].map(offsets) + (reg["start"] + reg["end"]) / 2  # position cumulée pour le placement des étiquettes
    # 5 critères complets du filtre STRICT (script 15), hors n_windows (voir
    # échange précédent : n_windows est un critère de largeur, pas de force
    # de signal, volontairement omis dans cette exploration "pics")
    reg["strict"] = (  # une région est stricte si elle passe les 5 seuils simultanément
        (reg["FST_ME"] >= 0.15) &
        (reg["FST_score"] >= 0.15) &
        (reg["FST_closest_group"] <= 0.05) &
        (reg["fd_max"] >= 0.85) &
        (reg["D_at_max_fd"] > 0.20)
    )
    print(f"  → {len(reg)} régions après fusion (pour les étiquettes gènes), dont {reg['strict'].sum()} strictes "
          f"(FST_ME≥0.15, FST_score≥0.15, FST_closest≤0.05, fd_max≥0.85, D>0.20 — n_windows exclu)")

    print("Annotation gènes (GFF Oar_v4.0, ±50kb)...")
    genes_by_chr = load_genes_by_chr()  # charge tous les gènes du GFF, groupés par chromosome
    reg["gene"] = reg.apply(lambda r: best_gene(genes_by_chr, r["chr"], r["start"], r["end"]), axis=1)  # annote chaque région

    # chr6 région stricte = désert génique mais connue comme étant en aval de
    # KIT (chr6:70,068,578-70,152,325) -> label plus informatif que "désert génique"
    is_kit_region = (reg["chr"] == 6) & (reg["start"] >= 70_200_000) & (reg["end"] <= 70_300_000)  # repère la région près de KIT
    reg.loc[is_kit_region, "gene"] = "aval de KIT"  # remplace le label par défaut

    out_tsv = OUTDIR / f"regions_top_intersection_K{K}.tsv"
    reg.sort_values(["chr", "start"]).to_csv(out_tsv, sep="\t", index=False)  # exporte les régions annotées
    print(f"  → {out_tsv}")

    # ── Figure : 2 panels empilés, même axe X (position cumulée) ────────────
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(32, 12), sharex=True,
                                    gridspec_kw={"hspace": 0.10})  # panel FST (ax1) au-dessus du panel fd (ax2)

    for i, c in enumerate(chroms):  # nuage de fond, tous les points, 2 couleurs en alternance par chromosome
        color = CHR_COLORS[i % 2]
        subf = fst[fst["chr"] == c]
        ax1.scatter(subf["cum_pos"], subf["FST_score"], s=3, color=color, alpha=0.5, linewidths=0)  # tous les FST_score
        subd = fd[fd["chr"] == c]
        ax2.scatter(subd["cum_pos"], subd["fd_max"], s=3, color=color, alpha=0.5, linewidths=0)  # tous les fd_max

    # 1 point par FENÊTRE brute de 20kb (pas par région) — le chevauchement
    # visuel entre fenêtres voisines (step 5kb) est normal, pas un défaut.
    reg_keys_strict = set(zip(reg[reg["strict"]]["chr"], reg[reg["strict"]]["start"], reg[reg["strict"]]["end"]))  # régions strictes

    def in_strict_region(row):
        return any(row["chr"] == c and s <= row["start"] and row["end"] <= e for c, s, e in reg_keys_strict)  # fenêtre incluse dans une région stricte ?

    hi_fst = fst[fst["is_top_both"]].copy()  # fenêtres FST dans l'intersection top-K/top-K
    hi_fd  = fd[fd["is_top_both"]].copy()  # fenêtres fd dans l'intersection top-K/top-K
    hi_fst["strict"] = hi_fst.apply(in_strict_region, axis=1)  # marque celles situées dans une région stricte
    hi_fd["strict"]  = hi_fd.apply(in_strict_region, axis=1)

    # points rouges de l'intersection : plus petits/pâles hors région stricte, plus gros/opaques dans une région stricte
    for df_, col, size, alpha in [(hi_fst[~hi_fst["strict"]], "FST_score", 12, 0.55),
                                   (hi_fst[hi_fst["strict"]],  "FST_score", 22, 0.95)]:
        ax1.scatter(df_["cum_pos"], df_[col], s=size, color="#d62728", alpha=alpha, zorder=5, linewidths=0)
    for df_, col, size, alpha in [(hi_fd[~hi_fd["strict"]], "fd_max", 12, 0.55),
                                   (hi_fd[hi_fd["strict"]],  "fd_max", 22, 0.95)]:
        ax2.scatter(df_["cum_pos"], df_[col], s=size, color="#d62728", alpha=alpha, zorder=5, linewidths=0)

    ax1.scatter([], [], s=22, color="#d62728", label=f"fenêtre 20kb, région stricte (n régions={reg['strict'].sum()})")  # entrée légende (points vides)
    ax1.scatter([], [], s=12, color="#d62728", alpha=0.55, label="fenêtre 20kb, région non stricte")  # entrée légende (points vides)

    ax1.axhline(thresh_fst, color="black", lw=0.9, ls="--", alpha=0.6, label="seuil top 1%")  # ligne = seuil top 1% FST
    ax2.axhline(thresh_fd,  color="black", lw=0.9, ls="--", alpha=0.6, label="seuil top 1%")  # ligne = seuil top 1% fd
    ax1.axhline(0.15, color="#b3121f", lw=1.1, ls="-.", alpha=0.7, label="seuil strict (FST_ME et FST_score ≥ 0.15)")  # seuil strict FST
    ax2.axhline(0.85, color="#b3121f", lw=1.1, ls="-.", alpha=0.7, label="seuil strict (fd_max ≥ 0.85)")  # seuil strict fd

    # étiquettes horizontales, JUSTE AU-DESSUS du pic — uniquement pour les
    # régions strictes (les 41 autres ne sont plus annotées, pour la lisibilité)
    ax1.set_ylim(0, 0.50)  # marge en haut pour l'empilement des étiquettes

    strict_only = reg[reg["strict"]].sort_values("cum_pos").reset_index(drop=True)  # régions strictes, ordre génomique
    genome_len = offsets[chroms[-1]] + int(max_end[chroms[-1]])  # longueur totale du génome (position cumulée)
    min_gap = genome_len / 90  # distance horizontale minimale entre deux étiquettes d'un même niveau
    n_tiers = 5  # nombre de niveaux verticaux disponibles pour empiler les étiquettes
    dy_step = 0.028  # décalage vertical entre niveaux
    tier_last_x = [-np.inf] * n_tiers  # dernière position utilisée sur chaque niveau

    for _, row in strict_only.iterrows():  # place chaque étiquette de gène sans chevauchement horizontal
        tier = n_tiers - 1  # niveau par défaut si tous les niveaux sont trop proches
        for t in range(n_tiers):
            if row["cum_pos"] - tier_last_x[t] >= min_gap:  # niveau libre (assez loin de la dernière étiquette)
                tier = t
                break
        tier_last_x[tier] = row["cum_pos"]
        y = row["FST_score"] + 0.012 + tier * dy_step  # hauteur de l'étiquette = juste au-dessus du pic + niveau
        ax1.text(row["cum_pos"], y, row["gene"], rotation=0, ha="center", va="bottom",
                  fontsize=8.5, color="#b3121f", fontweight="bold")

    ax1.set_ylabel("FST_score\n(Awassi vs ME − vs groupe le plus proche)", fontsize=13)
    ax1.set_title("Manhattan FST_score — genome entier, 25 chromosomes (fenêtres 20kb/step5kb)",
                   fontsize=13, loc="left")
    ax1.legend(fontsize=13, markerscale=2.2, loc="upper left", bbox_to_anchor=(1.005, 1.0), borderaxespad=0)
    ax1.tick_params(axis="both", labelsize=11)
    for sp in ["top", "right"]: ax1.spines[sp].set_visible(False)  # retire le cadre haut/droite

    ax2.set_ylabel("fd_max\n(P1=ME, P2=Awassi, meilleur P3)", fontsize=13)
    ax2.set_title("Manhattan fd_max — genome entier, 25 chromosomes (fenêtres 20kb/step5kb)",
                   fontsize=13, loc="left")
    ax2.set_ylim(0, 1.05)
    ax2.legend(fontsize=13, markerscale=2.2, loc="upper left", bbox_to_anchor=(1.005, 1.0), borderaxespad=0)
    ax2.tick_params(axis="both", labelsize=11)
    for sp in ["top", "right"]: ax2.spines[sp].set_visible(False)  # retire le cadre haut/droite

    # ticks = milieu de chaque chromosome ; xlim = genome ENTIER (fix bug chr25 tronqué)
    tick_pos = [offsets[c] + int(max_end[c]) / 2 for c in chroms]  # position cumulée du milieu de chaque chromosome
    ax2.set_xticks(tick_pos)
    ax2.set_xticklabels([str(c) for c in chroms], fontsize=11)  # étiquettes = numéro de chromosome
    ax2.set_xlabel("Chromosome", fontsize=13)
    last_c = chroms[-1]
    ax2.set_xlim(0, offsets[last_c] + int(max_end[last_c]) + 500_000)  # borne droite = fin du dernier chromosome

    fig.text(-0.01, 0.5,
              f"Intersection top{K} FST seul ∩ top{K} fd seul : {len(inter_keys)} fenêtres 20kb "
              f"→ {len(reg)} régions fusionnées, dont {reg['strict'].sum()} strictes (FST_ME≥0.15 & fd_max≥0.85)",
              rotation=90, va="center", ha="left", fontsize=13, color="#444444")  # texte vertical, bord gauche

    fig.suptitle("FST et fd regardés séparément (pas de score combiné) — recoupement en rouge, nom du gène annoté",
                 fontsize=13, fontweight="bold", y=1.0)

    out_png = OUTDIR / "manhattan_fst_fd_genomewide.png"
    out_pdf = OUTDIR / "manhattan_fst_fd_genomewide.pdf"
    plt.savefig(out_png, dpi=200, bbox_inches="tight")  # sauvegarde en PNG
    plt.savefig(out_pdf, bbox_inches="tight")  # sauvegarde en PDF (vectoriel)
    plt.close()  # ferme la figure (libère la mémoire)
    print(f"\n→ {out_png}")
    print(f"→ {out_pdf}")


if __name__ == "__main__":
    main()
