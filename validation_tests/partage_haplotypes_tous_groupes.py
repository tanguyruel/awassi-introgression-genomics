#!/usr/bin/env python3
"""
Partage d'haplotypes Awassi ↔ CHAQUE groupe (pas seulement le P3 retenu).

Pour chaque région et chaque groupe G, on calcule le pourcentage d'haplotypes Awassi
qui possèdent dans G un haplotype à moins de 1 % de sites différents (et, séparément,
strictement identique).

Deux versions sont produites :
  - BRUTE      : tous les haplotypes de G. Comparable ENTRE régions pour un même G,
                 mais PAS entre groupes de tailles différentes (Africa=334 haplotypes
                 contre America=16 : le gros groupe gagne mécaniquement).
  - RARÉFIÉE   : G sous-échantillonné à 16 haplotypes (taille d'America), 500 tirages.
                 Comparable entre groupes ET entre régions.

Ascertainment des SNP : un seul jeu par région, filtré sur les 7 groupes réunis
(missing ≤ 5 %, MAF ≥ 5 %), donc indépendant du groupe testé. Ce choix diffère des
scripts 33/34 qui filtraient sur Awassi+P3_best+ME — les valeurs ne sont donc pas
strictement superposables aux leurs.

Sortie : analyses/synthese_resultats/haplotype_sharing_all_groups/

Pas cité directement dans l'Annexe A (le résultat rapporté vient du bootstrap_specificite_haplotypique.py), mais
requis pour l'exécuter : 65 (et 69) importent ce fichier par nom de module, il doit
rester dans le même dossier. Usage : python3 64_haplotype_sharing_all_groups_v1.py [--project <dossier>]

Le dossier racine des données (hors dépôt, contient analyses/) se règle via --project
en exécution directe ; en import dynamique (par bootstrap_specificite_haplotypique.py et
ccdc91_retest_pic.py), il est lu depuis la variable d'environnement AWASSI_PROJECT_DIR
(sinon le répertoire courant) au moment de l'import.
"""
import argparse
import subprocess, csv, sys
import os
from pathlib import Path
import numpy as np

def _default_project():
    """Dossier racine des données par défaut : $AWASSI_PROJECT_DIR, sinon le répertoire courant."""
    return Path(os.environ.get("AWASSI_PROJECT_DIR", str(Path.cwd())))

PROJ = _default_project()
POP = PROJ / "analyses/fst/popmaps_separees_v1"  # dossier des popmaps (listes d'individus par groupe)
PHASE = PROJ / "analyses/phasing_beagle"
OUT = PROJ / "analyses/synthese_resultats/haplotype_sharing_all_groups"

MAX_MISSING, MIN_MAF = 0.05, 0.05
MISMATCH_THRESHOLD = 0.01  # seuil de fraction de mésappariement pour un match "quasi identique"
N_ITER, RAREFY_TO = 500, 16  # nb de tirages de raréfaction, taille cible (= plus petit groupe, America)
GROUPS = ["Awassi", "MiddleEastNonAwassi", "Africa", "Asia", "Europe", "America", "Australia"]
CANDIDATES = ["Africa", "Asia", "Europe", "America", "Australia"]   # P3 possibles (ME exclu)

REGIONS = [  # (id région, chemin VCF relatif, chr, start, end, P3 retenu) — 10 régions
    ("chr2_112.8Mb_NIPA2_CYFIP1",  "Phasing_Beagle_v2_5regions_v10b/phased/chr2_112.79_112.87Mb.beagle_phased.vcf.gz",  "2",  112785001, 112865000, "Australia"),
    ("chr3_129.2Mb_desert",        "Phasing_Beagle_v2_5regions_v10b/phased/chr3_129.22_129.26Mb.beagle_phased.vcf.gz",  "3",  129220001, 129260000, "Europe"),
    ("chr3_137.4Mb_LOC101123547",  "Phasing_Beagle_v3_2regions_new/phased/chr3_137.39_137.42Mb.beagle_phased.vcf.gz",   "3",  137390001, 137420000, "America"),
    ("chr3_186.1Mb_CCDC91",        "Phasing_Beagle_v1_5regions/phased/chr3_186.07_186.12Mb.beagle_phased.vcf.gz",       "3",  186075001, 186125000, "Asia"),
    ("chr5_58.1Mb_ABLIM3_AFAP1L1", "Phasing_Beagle_v2_5regions_v10b/phased/chr5_58.07_58.11Mb.beagle_phased.vcf.gz",    "5",  58070001,  58110000,  "Europe"),
    ("chr6_70.2Mb_KIT",            "Phasing_Beagle_12_06_14h28/phased/chr6_68.86_70.95Mb.beagle_phased.vcf.gz",         "6",  70240001,  70295000,  "Africa"),
    ("chr10_49.4Mb_KLF12",         "Phasing_Beagle_v3_2regions_new/phased/chr10_49.36_49.39Mb.beagle_phased.vcf.gz",    "10", 49360001,  49390000,  "Asia"),
    ("chr16_31.1Mb_NNT",           "Phasing_Beagle_NNT/phased/chr16_31.10_31.22Mb.beagle_phased.vcf.gz",                "16", 31100000,  31220000,  "Asia"),
    ("chr17_34.2Mb_SPATA5",        "Phasing_Beagle_v1_5regions/phased/chr17_34.16_34.22Mb.beagle_phased.vcf.gz",        "17", 34160001,  34220000,  "Asia"),
    ("chr20_0.8Mb_KHDRBS2",        "Phasing_Beagle_v2_5regions_v10b/phased/chr20_0.765_0.845Mb.beagle_phased.vcf.gz",   "20", 765001,    845000,    "Europe"),
]


def load_pop():
    """Construit {groupe: set des noms d'échantillons} depuis les fichiers popmap sous POP."""
    return {g: {l.split()[0] for l in (POP / f"{g}.txt").read_text().split("\n") if l.strip()}
            for g in GROUPS}


def read_haplotypes(vcf, chrom, start, end):
    """Retourne (H, samples) : H = matrice (2*n_samples) × n_snps, valeurs 0/1/-1."""
    samples = subprocess.run(["bcftools", "query", "-l", str(vcf)],
                             capture_output=True, text=True, check=True).stdout.split()
    res = subprocess.run(["bcftools", "query", "-r", f"{chrom}:{start}-{end}",
                          "-f", "[%GT ]\n", str(vcf)],
                         capture_output=True, text=True, check=True).stdout.strip("\n")
    if not res:
        raise RuntimeError(f"aucun variant dans {chrom}:{start}-{end}")
    rows = []
    for line in res.split("\n"):
        gts = line.split()
        hap = np.empty(2 * len(gts), dtype=np.int8)  # 2 allèles (haplotypes) par individu
        for i, gt in enumerate(gts):
            sep = "|" if "|" in gt else "/"  # séparateur phasé "|" ou non phasé "/"
            a, _, b = gt.partition(sep)
            hap[2 * i]     = -1 if a in (".", "") else int(a)  # -1 = allèle manquant
            hap[2 * i + 1] = -1 if b in (".", "") else int(b)
        rows.append(hap)  # une ligne = un SNP (2*n_samples valeurs)
    return np.array(rows, dtype=np.int8).T, samples          # (n_hap, n_snps)


def filter_snps(H):
    """Filtre les SNP sur missingness (MAX_MISSING) et fréquence allélique mineure (MIN_MAF).

    Parameters
    ----------
    H : numpy.ndarray
        Matrice haplotypes x SNP, valeurs 0/1/-1 (-1 = manquant).

    Returns
    -------
    tuple[numpy.ndarray, int]
        Matrice filtrée et nombre de SNP retenus.
    """
    valid = H >= 0  # True si l'allèle n'est pas manquant (-1)
    miss = 1 - valid.mean(axis=0)  # proportion de valeurs manquantes par SNP
    keep = miss <= MAX_MISSING
    with np.errstate(invalid="ignore"):
        af = np.where(valid, H, 0).sum(axis=0) / np.maximum(valid.sum(axis=0), 1)  # fréquence de l'allèle 1 par SNP
    maf = np.minimum(af, 1 - af)
    keep &= maf >= MIN_MAF
    return H[:, keep], int(keep.sum())  # matrice filtrée + nombre de SNP retenus


def mismatch_matrix(A, B):
    """D[i, j] = fraction de sites différents entre l'haplotype Awassi i et l'haplotype j de B.

    Calculée une seule fois : la raréfaction se ramène ensuite à un tirage de colonnes,
    ce qui évite de recalculer les distances à chaque itération.
    """
    D = np.empty((A.shape[0], B.shape[0]))
    for i in range(A.shape[0]):
        valid = (A[i] >= 0) & (B >= 0)  # SNP valides (non manquants) pour la paire (haplotype i, chaque j de B)
        diff = (A[i] != B) & valid
        n_valid = valid.sum(axis=1)  # nb de SNP comparables pour chaque j
        D[i] = np.where(n_valid > 0, diff.sum(axis=1) / np.maximum(n_valid, 1), np.nan)
    return D


def main():
    """Point d'entrée CLI : calcule le partage d'haplotypes Awassi/groupe pour chaque région de REGIONS."""
    global PROJ, POP, PHASE, OUT

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--project",
        default=str(PROJ),
        help="Dossier racine des données (contient analyses/). "
             "Par défaut : variable d'environnement AWASSI_PROJECT_DIR, sinon le répertoire courant.",
    )
    args = ap.parse_args()

    PROJ = Path(args.project)
    POP = PROJ / "analyses/fst/popmaps_separees_v1"
    PHASE = PROJ / "analyses/phasing_beagle"
    OUT = PROJ / "analyses/synthese_resultats/haplotype_sharing_all_groups"

    OUT.mkdir(parents=True, exist_ok=True)
    pop = load_pop()
    rng = np.random.default_rng(12345)  # générateur aléatoire reproductible (graine fixe)
    rows = []

    for rid, rel, chrom, start, end, p3best in REGIONS:
        vcf = PHASE / rel
        if not vcf.exists():
            print(f"  !! VCF absent : {vcf}", file=sys.stderr); continue
        H, samples = read_haplotypes(vcf, chrom, start, end)
        # indices de colonnes (2 haplotypes par individu) appartenant à chaque groupe
        idx = {g: np.array([j for i, s in enumerate(samples) if s in pop[g] for j in (2 * i, 2 * i + 1)])
               for g in GROUPS}
        used = np.concatenate([idx[g] for g in GROUPS])  # tous les haplotypes utilisés, dans l'ordre des groupes
        H = H[used]
        off, pos = 0, {}
        for g in GROUPS:
            pos[g] = np.arange(off, off + len(idx[g])); off += len(idx[g])  # plage d'indices de chaque groupe dans H réordonné

        H, n_snps = filter_snps(H)  # filtre qualité des SNP (sur les 7 groupes réunis)
        A = H[pos["Awassi"]]
        print(f"{rid}: {n_snps} SNP retenus, {A.shape[0]} haplotypes Awassi")

        for g in GROUPS:
            if g == "Awassi":
                continue
            B = H[pos[g]]
            D = mismatch_matrix(A, B)
            mm = np.nanmin(D, axis=1)  # meilleur match dans G pour chaque haplotype Awassi (version BRUTE)
            raw_exact = 100 * np.mean(mm == 0)
            raw_lt1 = 100 * np.mean(mm < MISMATCH_THRESHOLD)

            if B.shape[0] > RAREFY_TO:
                # raréfaction : G a plus d'haplotypes que la taille cible -> sous-échantillonnage répété N_ITER fois
                ex, lt = [], []
                for _ in range(N_ITER):
                    sub = rng.choice(B.shape[0], RAREFY_TO, replace=False)
                    m = np.nanmin(D[:, sub], axis=1)  # meilleur match dans ce sous-échantillon, pour chaque haplotype Awassi
                    ex.append(100 * np.mean(m == 0)); lt.append(100 * np.mean(m < MISMATCH_THRESHOLD))
                r_ex, r_ex_sd, r_lt, r_lt_sd = np.mean(ex), np.std(ex), np.mean(lt), np.std(lt)
            else:
                # G déjà <= taille cible : raréfaction impossible, on reprend les valeurs brutes
                r_ex, r_ex_sd, r_lt, r_lt_sd = raw_exact, 0.0, raw_lt1, 0.0

            rows.append({"region_id": rid, "P3_best": p3best, "group": g,
                         "is_P3_best": g == p3best, "is_candidate": g in CANDIDATES,
                         "n_snps": n_snps, "n_hap_group": B.shape[0],
                         "raw_exact_pct": round(raw_exact, 1), "raw_lt1pct_pct": round(raw_lt1, 1),
                         "rarefied_exact_pct": round(float(r_ex), 1), "rarefied_exact_sd": round(float(r_ex_sd), 1),
                         "rarefied_lt1pct_pct": round(float(r_lt), 1), "rarefied_lt1pct_sd": round(float(r_lt_sd), 1)})

    f = OUT / "Haplotype_sharing_all_groups.tsv"
    with open(f, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("\nÉcrit :", f, f"({len(rows)} lignes)")


if __name__ == "__main__":
    main()
