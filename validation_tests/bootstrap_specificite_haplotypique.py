#!/usr/bin/env python3
"""
Spécificité haplotypique du P3, testée proprement.

Deux corrections par rapport au partage_haplotypes_tous_groupes.py :

1. L'espérance de la raréfaction est calculée EXACTEMENT, sans Monte-Carlo.
   Pour l'haplotype Awassi i et le groupe g comptant n_g haplotypes dont m_ig
   sont des « jumeaux » (< 1 % de sites différents), la probabilité qu'un tirage
   de k=16 haplotypes contienne au moins un jumeau vaut :

       p_ig = 1 − C(n_g − m_ig, k) / C(n_g, k)

   Le score du groupe est la moyenne des p_ig sur les 44 haplotypes Awassi.
   L'écart-type entre tirages du partage_haplotypes_tous_groupes.py n'est que du bruit de simulation :
   il ne mesure PAS l'incertitude du résultat.

2. L'incertitude réelle vient du petit nombre d'individus Awassi (22). On la
   mesure par bootstrap sur ces individus (2000 rééchantillonnages avec remise),
   ce qui donne : la probabilité que le P3 soit le meilleur groupe, et
   l'intervalle de confiance de sa marge sur le meilleur des autres.

Sortie : analyses/synthese_resultats/haplotype_sharing_all_groups/
         Haplotype_specificity_bootstrap.tsv

Script de référence pour "Bootstrap de spécificité" (appelé "bootstrap_specificite_haplotypique.py" dans l'Annexe A du
rapport de stage) — remplace le script "64" (Monte-Carlo, ici partage_haplotypes_tous_groupes.py)
pour le résultat rapporté, mais l'importe encore pour les données de partage : IMPORTANT,
partage_haplotypes_tous_groupes.py doit rester présent dans le même dossier (import dynamique
par chemin de fichier, voir plus bas).
Usage : python3 bootstrap_specificite_haplotypique.py [--project <dossier>]

Le dossier racine des données (sert uniquement à localiser les résultats déjà calculés,
hors de ce dépôt) se règle via --project, sinon la variable d'environnement
AWASSI_PROJECT_DIR, sinon le répertoire courant.
"""
import argparse
import csv
import importlib.util
import os
from math import comb
from pathlib import Path
import numpy as np

SCRIPT_PARTAGE = Path(__file__).with_name("partage_haplotypes_tous_groupes.py")


def _load_partage_module():
    """Importe dynamiquement partage_haplotypes_tous_groupes.py (même dossier que ce script).

    Réutilise ainsi ses fonctions et données (load_pop, read_haplotypes,
    filter_snps, mismatch_matrix, REGIONS, GROUPS, PHASE) sans dupliquer le
    code. Doit être appelée après avoir positionné la variable d'environnement
    AWASSI_PROJECT_DIR, dont ce module dépend pour localiser ses propres données.
    """
    spec = importlib.util.spec_from_file_location("partage_haplotypes", SCRIPT_PARTAGE)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


K, N_BOOT, THRESH = 16, 2000, 0.01  # taille de tirage (raréfaction), nb de bootstraps, seuil "jumeau" (<1% de sites différents)
CANDIDATES = ["Africa", "Asia", "Europe", "America", "Australia"]  # groupes P3 candidats testés (ME exclu)


def p_at_least_one(n, m, k):
    """Proba qu'un tirage sans remise de k haplotypes parmi n en contienne >= 1 des m jumeaux."""
    if m == 0:
        return 0.0
    if n - m < k:
        return 1.0  # moins de k non-jumeaux disponibles -> un jumeau est certain dans le tirage
    # probabilité complémentaire : 1 - P(aucun jumeau tiré)
    # P(aucun jumeau) = C(n-m, k) / C(n, k)  (tirages de k parmi les n-m non-jumeaux, sur tous les tirages de k parmi n)
    return 1.0 - comb(n - m, k) / comb(n, k)


def main():
    """Point d'entrée CLI : bootstrap de spécificité haplotypique du P3 pour chaque région de REGIONS."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--project",
        default=os.environ.get("AWASSI_PROJECT_DIR", str(Path.cwd())),
        help="Dossier racine des données (contient analyses/). "
             "Par défaut : variable d'environnement AWASSI_PROJECT_DIR, sinon le répertoire courant.",
    )
    args = ap.parse_args()

    os.environ["AWASSI_PROJECT_DIR"] = args.project  # lu par partage_haplotypes_tous_groupes.py au chargement
    m64 = _load_partage_module()
    out_dir = Path(args.project) / "analyses/synthese_resultats/haplotype_sharing_all_groups"

    pop = m64.load_pop()  # groupes -> ensembles d'échantillons (fonction de partage_haplotypes_tous_groupes.py)
    rng = np.random.default_rng(7)  # générateur aléatoire reproductible (graine fixe) pour le bootstrap
    rows = []

    for rid, rel, chrom, start, end, p3best in m64.REGIONS:  # REGIONS défini dans partage_haplotypes_tous_groupes.py
        vcf = m64.PHASE / rel  # chemin complet du VCF phasé (PHASE défini dans partage_haplotypes_tous_groupes.py)
        H, samples = m64.read_haplotypes(vcf, chrom, start, end)  # charge la matrice haplotypes x SNP (fonction de partage_haplotypes_tous_groupes.py)
        idx = {g: np.array([j for i, s in enumerate(samples) if s in pop[g] for j in (2 * i, 2 * i + 1)])
               for g in m64.GROUPS}  # indices d'haplotypes par groupe (GROUPS défini dans partage_haplotypes_tous_groupes.py)
        used = np.concatenate([idx[g] for g in m64.GROUPS])  # tous les haplotypes utilisés, dans l'ordre des groupes
        H = H[used]
        off, pos = 0, {}
        for g in m64.GROUPS:
            pos[g] = np.arange(off, off + len(idx[g])); off += len(idx[g])  # plage d'indices de chaque groupe
        H, n_snps = m64.filter_snps(H)  # filtre qualité des SNP (fonction de partage_haplotypes_tous_groupes.py)
        A = H[pos["Awassi"]]                       # 44 haplotypes = 22 individus × 2
        n_ind = A.shape[0] // 2

        # m[i, g] = nb de jumeaux de l'haplotype Awassi i dans le groupe g
        M, NG = {}, {}
        for g in CANDIDATES:
            B = H[pos[g]]
            D = m64.mismatch_matrix(A, B)  # fraction de mésappariement Awassi x g (fonction de partage_haplotypes_tous_groupes.py)
            M[g] = (D < THRESH).sum(axis=1)  # nb de "jumeaux" (< 1% mésappariement) par haplotype Awassi, dans g
            NG[g] = B.shape[0]  # effectif total du groupe g

        # p[i, g] : espérance exacte, puis score du groupe
        P = {g: np.array([p_at_least_one(NG[g], int(mi), K) for mi in M[g]]) for g in CANDIDATES}  # proba exacte par haplotype Awassi et par groupe
        score = {g: 100 * P[g].mean() for g in CANDIDATES}  # score du groupe = moyenne des proba sur tous les haplotypes Awassi (en %)

        # bootstrap sur les 22 individus Awassi (les 2 haplotypes d'un individu restent ensemble)
        wins, margins = 0, []
        for _ in range(N_BOOT):
            ind = rng.integers(0, n_ind, n_ind)  # tirage avec remise de n_ind individus (indices d'individus)
            hap = np.concatenate([[2 * i, 2 * i + 1] for i in ind])
            sc = {g: 100 * P[g][hap].mean() for g in CANDIDATES}
            best_other = max((g for g in CANDIDATES if g != p3best), key=lambda g: sc[g])  # meilleur concurrent du P3 retenu
            margins.append(sc[p3best] - sc[best_other])  # marge du P3 retenu sur son meilleur concurrent
            if sc[p3best] >= max(sc.values()):
                wins += 1  # compte les tirages où le P3 retenu est (ex-aequo) le meilleur groupe

        margins = np.array(margins)
        lo, hi = np.percentile(margins, [2.5, 97.5])  # intervalle de confiance à 95% de la marge (percentiles bootstrap)
        best_obs = max(score, key=score.get)  # groupe avec le meilleur score observé (sans bootstrap)
        rows.append({
            "region_id": rid, "P3_best": p3best, "n_snps": n_snps,
            "score_P3": round(score[p3best], 1),
            "best_group": best_obs, "score_best": round(score[best_obs], 1),
            "margin_P3_vs_best_other": round(score[p3best] - max(score[g] for g in CANDIDATES if g != p3best), 1),  # marge observée (sans bootstrap)
            "boot_margin_lo": round(float(lo), 1), "boot_margin_hi": round(float(hi), 1),
            "P_P3_is_best": round(wins / N_BOOT, 3),
            **{f"score_{g}": round(score[g], 1) for g in CANDIDATES},
        })
        print(f"{rid:28s} P3={p3best:10s} score_P3={score[p3best]:5.1f}  meilleur={best_obs:10s} "
              f"marge={rows[-1]['margin_P3_vs_best_other']:+6.1f} [{lo:+.1f};{hi:+.1f}]  P(P3 best)={wins/N_BOOT:.2f}")

    f = out_dir / "Haplotype_specificity_bootstrap.tsv"
    with open(f, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("\nÉcrit :", f)


if __name__ == "__main__":
    main()
