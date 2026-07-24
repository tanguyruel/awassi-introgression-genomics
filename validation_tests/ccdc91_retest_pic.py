#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CCDC91 : reprise du test de partage haplotypique sur le cœur du pic.

Contexte
--------
La région chr3:186,075,001–186,125,000 est une fusion de fenêtres glissantes.
Le détail fenêtre par fenêtre montre que le signal n'y est pas homogène :

    186,075–186,095   fd 0,852   Asia
    186,080–186,100   fd 0,912   Asia
    186,085–186,105   fd 0,933   Asia     <- pic
    186,090–186,110   fd 0,929   Asia
    186,095–186,115   fd 0,908   Asia
    186,100–186,120   fd 0,898   Europe   <- la queue bascule
    186,105–186,125   fd 0,884   Africa

Les tests publiés ont porté sur la région entière, qui contient le pic (vérifié
séparément par une analyse fenêtre par fenêtre du fd). Mais elle contient aussi
une queue de 20 kb où le groupe candidat n'est plus l'Asie. Si cette queue dilue
le signal, le verdict négatif obtenu sur la région entière pourrait être un
artefact de fenêtrage.

Ce script rejoue donc le test sur trois fenêtres emboîtées, avec exactement la
même méthode que partage_haplotypes_tous_groupes.py (même raréfaction à 16
haplotypes, même bootstrap sur les 22 individus Awassi) :

    A. région entière      186,075–186,125   (ce qui a été publié)
    B. cœur Asia           186,075–186,115   (fenêtres dont le P3 est l'Asie)
    C. pic strict          186,085–186,105   (la seule fenêtre du maximum de fd)

Si le verdict est identique sur les trois, il ne dépend pas du fenêtrage et la
conclusion sur CCDC91 est solide. S'il change, la région doit être rouverte.

Sortie : analyses/synthese_resultats/CCDC91_peak_core_retest.tsv

Script de référence pour "CCDC91 — reprise sur le cœur du pic" (§3.6, appelé
"ccdc91_retest_pic.py" dans l'Annexe A du rapport de stage) — cas particulier, pas une
méthode générique réutilisable telle quelle pour d'autres régions.
IMPORTANT : importe partage_haplotypes_tous_groupes.py (import dynamique par chemin de
fichier, voir plus bas) — ce fichier doit rester présent dans le même dossier.
Chemin PROJ ci-dessous en dur, à adapter si relancé ailleurs (sert uniquement à localiser
le VCF déjà phasé, hors de ce dépôt).
Usage : python3 ccdc91_retest_pic.py
"""

# ── Imports : import dynamique de module, combinatoire, chemins, calcul ──
import importlib.util
from math import comb
from pathlib import Path

import numpy as np

PROJ = Path("/home/tanguyruel/Bureau/genome_complet_Awassi")  # chemin en dur à adapter (dossier de travail, hors dépôt)

# ── Import dynamique de partage_haplotypes_tous_groupes.py (même dossier que ce
# script), pour réutiliser ses fonctions (load_pop, read_haplotypes, filter_snps,
# mismatch_matrix) sans dupliquer le code ──
SCRIPT_PARTAGE = Path(__file__).with_name("partage_haplotypes_tous_groupes.py")
spec = importlib.util.spec_from_file_location("partage_haplotypes", SCRIPT_PARTAGE)
m64 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m64)  # exécute le script -> ses fonctions deviennent accessibles via m64.xxx

VCF = (PROJ / "analyses/phasing_beagle/Phasing_Beagle_v1_5regions/phased/"
       "chr3_186.07_186.12Mb.beagle_phased.vcf.gz")  # VCF phasé de la région CCDC91
OUT = PROJ / "analyses/synthese_resultats/CCDC91_peak_core_retest.tsv"  # fichier de sortie

K, N_BOOT, THRESH = 16, 2000, 0.01  # taille de tirage (raréfaction), nb de bootstraps, seuil "jumeau"
CANDIDATES = ["Africa", "Asia", "Europe", "America", "Australia"]  # groupes P3 candidats testés
P3_REGION = "Asia"  # groupe P3 retenu pour cette région (à tester contre les autres)
RNG = np.random.default_rng(20260723)  # générateur aléatoire reproductible (graine fixe)

WINDOWS = [  # les 3 fenêtres emboîtées à comparer (tag, start, end, description)
    ("A_region_entiere", 186075001, 186125000, "région publiée"),
    ("B_coeur_Asia",     186075001, 186115000, "fenêtres dont le P3 est l'Asie"),
    ("C_pic_strict",     186085001, 186105000, "fenêtre du maximum de fd"),
]


def rarefied_score(D_by_group, group, k=K):
    """Probabilité moyenne qu'un tirage de k haplotypes du groupe contienne
    un jumeau (< THRESH de sites divergents) d'un haplotype Awassi donné.
    Calcul hypergéométrique exact, identique au rarefaction_partage_haplotypes.py."""
    D = D_by_group[group]  # matrice de mésappariement Awassi x groupe, pour ce groupe
    n = D.shape[1]  # effectif du groupe (nb d'haplotypes)
    if n < k:
        return float("nan")  # groupe trop petit pour tirer k haplotypes
    probs = []
    for i in range(D.shape[0]):
        d = D[i]  # mésappariements de l'haplotype Awassi i avec tous les haplotypes du groupe
        m = int(np.sum(d[~np.isnan(d)] < THRESH))     # nb de jumeaux disponibles
        if m == 0:
            probs.append(0.0)  # aucun jumeau -> probabilité nulle pour cet haplotype
        else:
            # P(au moins un jumeau dans un tirage de k parmi n)
            probs.append(1.0 - comb(n - m, k) / comb(n, k) if n - m >= k else 1.0)
    return float(np.mean(probs)) * 100  # moyenne sur tous les haplotypes Awassi, en %


def main():
    pops = m64.load_pop()  # groupes -> ensembles d'échantillons (fonction de partage_haplotypes_tous_groupes.py)
    rows = []  # résultats accumulés (une ligne par fenêtre)

    for tag, start, end, note in WINDOWS:
        H, samples = m64.read_haplotypes(VCF, "3", start, end)  # charge la matrice haplotypes x SNP de la fenêtre (fonction de partage_haplotypes_tous_groupes.py)
        H, n_snps = m64.filter_snps(H)  # filtre qualité des SNP (fonction de partage_haplotypes_tous_groupes.py)

        # index des haplotypes par groupe (2 haplotypes par individu)
        idx = {}
        for g, members in pops.items():
            cols = [2 * i + h for i, s in enumerate(samples)
                    if s in members for h in (0, 1)]  # indices des 2 haplotypes de chaque individu du groupe g
            if cols:
                idx[g] = np.array(cols)
        A = H[idx["Awassi"]]  # sous-matrice des haplotypes Awassi
        D_by_group = {g: m64.mismatch_matrix(A, H[idx[g]])
                      for g in CANDIDATES if g in idx}  # mésappariement Awassi x chaque groupe candidat (fonction de partage_haplotypes_tous_groupes.py)

        scores = {g: rarefied_score(D_by_group, g) for g in D_by_group}  # score de raréfaction exact par groupe
        best_other = max((g for g in scores if g != P3_REGION),
                         key=lambda g: scores[g])  # meilleur groupe concurrent du P3 retenu (Asia)
        marge = scores[P3_REGION] - scores[best_other]  # marge observée du P3 sur son meilleur concurrent

        # bootstrap sur les 22 individus Awassi (les 2 haplotypes restent solidaires)
        n_ind = A.shape[0] // 2  # nombre d'individus Awassi
        boot = []
        for _ in range(N_BOOT):
            pick = RNG.integers(0, n_ind, n_ind)  # tirage avec remise de n_ind individus
            hap = np.concatenate([[2 * p, 2 * p + 1] for p in pick])  # les 2 haplotypes de chaque individu tiré
            Db = {g: D_by_group[g][hap] for g in D_by_group}  # mésappariements restreints à cet échantillon bootstrap
            sc = {g: rarefied_score(Db, g) for g in Db}  # score de chaque groupe sur cet échantillon
            bo = max((g for g in sc if g != P3_REGION), key=lambda g: sc[g])  # meilleur concurrent sur ce tirage
            boot.append(sc[P3_REGION] - sc[bo])  # marge sur ce tirage bootstrap
        lo, hi = np.percentile(boot, [2.5, 97.5])  # intervalle de confiance à 95% de la marge (percentiles bootstrap)

        verdict = "oui" if lo > 0 else ("non" if hi < 0 else "indécis")  # marge signif. positive / négative / incertaine

        print(f"\n=== {tag}  chr3:{start/1e6:.3f}-{end/1e6:.3f} Mb  ({note}) ===")
        print(f"  SNPs retenus : {n_snps}")
        for g in sorted(scores, key=lambda g: -scores[g]):
            star = "  <- P3 retenu" if g == P3_REGION else ""
            print(f"    {g:<10s} {scores[g]:6.2f}{star}")
        print(f"  marge Asia - {best_other} = {marge:+.2f}  "
              f"IC95 [{lo:+.2f} ; {hi:+.2f}]  -> {verdict}")

        rows.append([tag, f"chr3:{start/1e6:.3f}-{end/1e6:.3f}", str(n_snps),
                     f"{scores[P3_REGION]:.2f}", best_other,
                     f"{scores[best_other]:.2f}", f"{marge:+.2f}",
                     f"[{lo:+.2f} ; {hi:+.2f}]", verdict])  # ligne de résultat pour cette fenêtre

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("fenetre\tcoordonnees\tn_snps\tscore_Asia\tmeilleur_autre\t"
                 "score_meilleur_autre\tmarge\tIC95\tverdict\n")  # en-tête du TSV
        for r in rows:
            fh.write("\t".join(r) + "\n")  # une ligne par fenêtre testée
    print(f"\nÉcrit : {OUT}")


if __name__ == "__main__":
    main()
