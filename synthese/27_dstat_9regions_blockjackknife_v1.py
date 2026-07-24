#!/usr/bin/env python3
"""
Script 27 — D-stat + Z-score (block-jackknife) sur les 9 régions candidates,
fenêtre élargie à 150kb (centrée sur chaque région) pour un jackknife robuste.

Reprend le même code que le dossier dstat (scripts/05_dstat/03_dstat_small_windows_zscore.py
et 06_dstat_zoom_regions_zscore.py) : ABBA-BABA polarisé sur Ovis_canadensis,
D = (ABBA-BABA)/(ABBA+BABA), SE/Z par block-jackknife delete-one-block.

Objectif : les 9 régions ont été sélectionnées sur un pic de fd dans une fenêtre
étroite (30-90kb) — ici on vérifie que le signal D reste significatif (|Z|>=3
en général) sur une fenêtre plus large (150kb, block_size 15kb → 10 blocs),
ce qui écarte l'hypothèse d'un simple pic bruité sur quelques SNPs.

P1=MiddleEastNonAwassi, P2=Awassi, OUT=Ovis_canadensis, P3=meilleur groupe déjà
identifié par région (pas de boucle sur plusieurs P3 ici, un seul par région).

Sortie : analyses/synthese_resultats/dstat_9regions_150kb/Dstat_9regions_150kb_blockjackknife.tsv

Script de référence pour "D-stat 150kb + jackknife" (cf. Annexe A du rapport de stage).
Usage : python3 27_dstat_9regions_blockjackknife_v1.py
"""

# Bibliothèque standard : écriture TSV, calcul (racine carrée), appel de bcftools, chemins, dict à valeur par défaut
import csv
import math
import subprocess
from pathlib import Path
from collections import defaultdict

POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # dossier des popmaps par groupe (1 fichier liste d'échantillons par groupe)
POPMAP_MAIN = Path("data/popmap_main5.tsv")  # seule source pour Ovis_canadensis (absent de POP_DIR)
VCF_DIR = Path("data/raw data_08_06")  # dossier des VCF bruts par chromosome
OUTDIR = Path("analyses/synthese_resultats/dstat_9regions_150kb")  # dossier de sortie des résultats
OUTDIR.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si absent

P1 = "MiddleEastNonAwassi"  # population 1 (référence, non-Awassi Moyen-Orient)
P2 = "Awassi"               # population 2 (testée pour introgression)
OUT = "Ovis_canadensis"     # groupe extérieur (outgroup) servant à polariser les allèles

WINDOW = 150_000     # taille de la fenêtre autour de chaque région (bp)
BLOCK_SIZE = 15_000  # taille d'un bloc jackknife (bp) -> ~10 blocs par fenêtre
MIN_SNPS = 15        # nombre minimum de SNPs utilisables pour garder une région
MIN_BLOCKS = 5       # nombre minimum de blocs jackknife pour calculer un SE fiable

OUT_LOW = 0.10   # fréquence alt outgroup en dessous : considéré comme allèle ancestral
OUT_HIGH = 0.90  # fréquence alt outgroup au-dessus : considéré comme allèle dérivé (polarité inversée)

# ── 9 régions strictes, fenêtre stricte FST/fd d'origine + P3 déjà identifié ──
# (coords des fenêtres strictes, cf. AWASSI_AGENT_LOG.md / project_awassi_candidats)
# Liste des 9 régions candidates : coordonnées de la fenêtre stricte d'origine + groupe P3 associé
REGIONS = [
    {"region_id": "chr2_112.8Mb_NIPA2_CYFIP1",   "chr": "2",  "strict_start": 112785001, "strict_end": 112865000, "P3": "Australia"},
    {"region_id": "chr3_129.2Mb_desert",           "chr": "3",  "strict_start": 129220001, "strict_end": 129260000, "P3": "Europe"},
    {"region_id": "chr3_137.4Mb_LOC101123547",     "chr": "3",  "strict_start": 137390001, "strict_end": 137420000, "P3": "America"},
    {"region_id": "chr3_186.1Mb_CCDC91",           "chr": "3",  "strict_start": 186075001, "strict_end": 186125000, "P3": "Asia"},
    {"region_id": "chr5_58.1Mb_ABLIM3_AFAP1L1",    "chr": "5",  "strict_start": 58070001,  "strict_end": 58110000,  "P3": "Europe"},
    {"region_id": "chr6_70.2Mb_KIT",               "chr": "6",  "strict_start": 70240001,  "strict_end": 70295000,  "P3": "Africa"},
    {"region_id": "chr10_49.4Mb_KLF12",            "chr": "10", "strict_start": 49360001,  "strict_end": 49390000,  "P3": "Asia"},
    {"region_id": "chr17_34.2Mb_SPATA5",           "chr": "17", "strict_start": 34160001,  "strict_end": 34220000,  "P3": "Asia"},
    {"region_id": "chr20_0.8Mb_KHDRBS2",           "chr": "20", "strict_start": 765001,    "strict_end": 845000,    "P3": "Europe"},
]


# Lit un fichier texte et renvoie la liste des lignes non vides (sans espaces de bord)
def read_list(path):
    with open(path) as f:  # ouvre le fichier en lecture
        return [x.strip() for x in f if x.strip()]  # une entrée par ligne non vide


# Extrait la liste des échantillons appartenant à `group` depuis un popmap générique (colonnes échantillon/groupe)
def read_group_from_popmap(path, group):
    samples = []  # accumulateur des échantillons du groupe
    with open(path) as f:
        f.readline()  # saute la ligne d'en-tête
        for line in f:  # parcourt chaque échantillon
            if not line.strip():
                continue  # ignore les lignes vides
            s, g = line.rstrip("\n").split("\t")[:2]  # colonnes échantillon (s) et groupe (g)
            if g == group:
                samples.append(s)  # garde l'échantillon si son groupe correspond
    return samples


def read_group(group):
    """P1/P2/P3 viennent des popmaps séparés (POP_DIR) ; OUT (Ovis_canadensis)
    n'existe que dans POPMAP_MAIN."""
    p = POP_DIR / f"{group}.txt"  # chemin du popmap dédié à ce groupe
    if p.exists():
        return read_list(p)  # cas normal : fichier liste dédié au groupe
    return read_group_from_popmap(POPMAP_MAIN, group)  # sinon (ex: outgroup) : popmap général


# Calcule la fréquence de l'allèle alternatif pour un sous-ensemble d'échantillons à un site donné
def alt_freq(gts, idxs):
    alt = 0  # compteur d'allèles alternatifs (=1)
    n = 0    # compteur total d'allèles appelés
    for i in idxs:  # parcourt les index d'échantillons du groupe
        gt = gts[i].replace("|", "/")  # uniformise le séparateur phasé/non phasé
        if gt in {".", "./.", ".|."}:
            continue  # génotype manquant, ignoré
        for a in gt.split("/"):  # parcourt les deux allèles du génotype
            if a == ".":
                continue  # allèle manquant, ignoré
            if a == "0":
                n += 1  # allèle référence : compte dans le total
            elif a == "1":
                alt += 1  # allèle alternatif : compte dans alt
                n += 1    # et dans le total
    return None if n == 0 else alt / n  # fréquence alt, ou None si aucun allèle appelé


# Récupère l'ordre des échantillons tel qu'il apparaît dans le VCF restreint au fichier d'échantillons
def get_sample_order(vcf, sample_file):
    p1 = subprocess.Popen(["bcftools", "view", "-S", str(sample_file), "-Ou", vcf], stdout=subprocess.PIPE)  # sous-sélectionne les échantillons, sortie BCF non compressée
    p2 = subprocess.run(["bcftools", "query", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True, check=True)  # liste les noms d'échantillons dans l'ordre du VCF
    p1.stdout.close()  # ferme le pipe intermédiaire
    p1.wait()  # attend la fin du premier processus
    return p2.stdout.strip().splitlines()  # une ligne = un nom d'échantillon


# Calcule SE (erreur-type) et Z-score du D-stat par block-jackknife delete-one-block
def jackknife(num, den, blocks):
    if den == 0:
        return None, None, 0  # dénominateur nul : D indéfini
    D = num / den  # D-statistic global = (ABBA-BABA)/(ABBA+BABA)
    values = []  # pseudo-valeurs D obtenues en retirant un bloc à la fois
    for bnum, bden, bn in blocks.values():  # parcourt chaque bloc jackknife
        num_j = num - bnum  # numérateur total moins celui du bloc retiré
        den_j = den - bden  # dénominateur total moins celui du bloc retiré
        if den_j != 0:
            values.append(num_j / den_j)  # D recalculé sans ce bloc (delete-one-block)
    B = len(values)  # nombre de blocs effectivement utilisés
    if B < MIN_BLOCKS:
        return None, None, B  # pas assez de blocs pour un SE fiable
    mean_j = sum(values) / B  # moyenne des pseudo-valeurs D
    var = ((B - 1) / B) * sum((x - mean_j) ** 2 for x in values)  # variance jackknife (formule delete-one-block)
    SE = math.sqrt(var) if var > 0 else None  # erreur-type = racine carrée de la variance jackknife
    Z = D / SE if SE and SE > 0 else None  # score Z = D / SE, teste si D diffère significativement de 0
    return SE, Z, B


# Calcule le D-stat + jackknife pour une région donnée (fenêtre 150kb autour du centre de la région stricte)
def run_region(reg):
    region_id, chrom, strict_start, strict_end, P3 = (
        reg["region_id"], reg["chr"], reg["strict_start"], reg["strict_end"], reg["P3"]
    )  # déballe les champs du dict région
    mid = (strict_start + strict_end) // 2  # centre de la région stricte d'origine
    win_start = mid - WINDOW // 2  # début de la fenêtre élargie 150kb
    win_end = win_start + WINDOW - 1  # fin de la fenêtre élargie 150kb

    vcf = str(VCF_DIR / f"awassi_and_basedata_chr{chrom}.vcf.gz")  # VCF source du chromosome concerné
    print(f"\n{'='*70}\n{region_id} | chr{chrom}:{win_start}-{win_end} (150kb) | P3={P3}")  # log de progression
    print(f"  (fenêtre stricte d'origine : {strict_start}-{strict_end})")  # rappel de la fenêtre stricte

    groups = {g: read_group(g) for g in [P1, P2, P3, OUT]}  # charge la liste d'échantillons de chaque groupe
    samples = sorted(set(groups[P1] + groups[P2] + groups[P3] + groups[OUT]))  # union unique de tous les échantillons utiles
    sample_file = OUTDIR / f"samples_{region_id}.txt"  # fichier temporaire liste d'échantillons pour bcftools -S
    sample_file.write_text("\n".join(samples) + "\n")  # écrit la liste, un échantillon par ligne

    order = get_sample_order(vcf, sample_file)  # ordre des échantillons dans le VCF sous-échantillonné
    idx = {s: i for i, s in enumerate(order)}  # position de chaque échantillon dans les colonnes GT
    idx1 = [idx[s] for s in groups[P1] if s in idx]  # indices colonnes des échantillons P1
    idx2 = [idx[s] for s in groups[P2] if s in idx]  # indices colonnes des échantillons P2
    idx3 = [idx[s] for s in groups[P3] if s in idx]  # indices colonnes des échantillons P3
    idxO = [idx[s] for s in groups[OUT] if s in idx]  # indices colonnes des échantillons outgroup
    print(f"  P1={len(idx1)} P2={len(idx2)} P3={len(idx3)} OUT={len(idxO)}")  # log des effectifs par groupe

    ABBA = BABA = num_tot = den_tot = 0.0  # accumulateurs globaux ABBA/BABA et numérateur/dénominateur du D-stat
    n_snps = 0  # nombre de SNPs polarisés utilisés
    blocks = defaultdict(lambda: [0.0, 0.0, 0])  # accumulateurs par bloc jackknife : [num, den, n_snps]

    pview = subprocess.Popen(
        ["bcftools", "view", "-r", f"{chrom}:{win_start}-{win_end}", "-f", "PASS",
         "-m2", "-M2", "-v", "snps", "-S", str(sample_file), "-Ou", vcf],
        stdout=subprocess.PIPE,
    )  # filtre : région, variants PASS, biallélique, SNP uniquement, échantillons sélectionnés
    pquery = subprocess.Popen(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n"],
        stdin=pview.stdout, stdout=subprocess.PIPE, text=True,
    )  # extrait CHROM, POS et les génotypes de chaque échantillon
    pview.stdout.close()  # ferme le pipe intermédiaire côté parent

    n_read = 0  # nombre total de lignes (SNPs) lues
    for line in pquery.stdout:  # parcourt chaque SNP renvoyé par bcftools
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue  # ligne incomplète, ignorée
        pos = int(parts[1])  # position du SNP
        gts = parts[2:]  # génotypes de tous les échantillons à ce site
        n_read += 1

        o = alt_freq(gts, idxO)  # fréquence allélique alt dans l'outgroup
        p1a = alt_freq(gts, idx1)  # fréquence allélique alt dans P1
        p2a = alt_freq(gts, idx2)  # fréquence allélique alt dans P2
        p3a = alt_freq(gts, idx3)  # fréquence allélique alt dans P3
        if o is None or p1a is None or p2a is None or p3a is None:
            continue  # donnée manquante dans un des groupes, SNP écarté

        if o <= OUT_LOW:
            p1, p2, p3 = p1a, p2a, p3a  # outgroup proche de l'ancestral : pas d'inversion
        elif o >= OUT_HIGH:
            p1, p2, p3 = 1 - p1a, 1 - p2a, 1 - p3a  # outgroup proche du dérivé : polarité inversée
        else:
            continue  # outgroup ambigu (ni fixé ni quasi-fixé), SNP écarté

        abba = (1 - p1) * p2 * p3  # probabilité du pattern ABBA (P1 ancestral, P2 et P3 dérivés)
        baba = p1 * (1 - p2) * p3  # probabilité du pattern BABA (P2 ancestral, P1 et P3 dérivés)
        num = abba - baba  # contribution de ce SNP au numérateur du D-stat
        den = abba + baba  # contribution de ce SNP au dénominateur du D-stat

        ABBA += abba  # cumul ABBA sur la fenêtre
        BABA += baba  # cumul BABA sur la fenêtre
        num_tot += num  # cumul numérateur
        den_tot += den  # cumul dénominateur
        n_snps += 1  # un SNP polarisé de plus

        block_id = int((pos - win_start) // BLOCK_SIZE)  # numéro du bloc jackknife contenant ce SNP
        blocks[block_id][0] += num  # cumul numérateur du bloc
        blocks[block_id][1] += den  # cumul dénominateur du bloc
        blocks[block_id][2] += 1    # cumul nombre de SNPs du bloc

    pquery.wait()  # attend la fin du process bcftools query
    pview.wait()   # attend la fin du process bcftools view

    print(f"  SNPs lus : {n_read} ; SNPs polarisés/utilisés : {n_snps} ; blocs : {len(blocks)}")  # log récapitulatif

    if n_snps < MIN_SNPS or den_tot == 0:
        print(f"  [!] Trop peu de SNPs ({n_snps}) ou den=0 — région ignorée")
        return None  # région non exploitable, aucun résultat renvoyé

    D = num_tot / den_tot  # D-stat final sur toute la fenêtre = (ABBA-BABA)/(ABBA+BABA)
    SE, Z, n_blocks = jackknife(num_tot, den_tot, blocks)  # erreur-type et Z-score par block-jackknife
    print(f"  D={D:.4f}  SE={SE}  Z={Z}  n_blocks={n_blocks}")  # log du résultat

    # construit la ligne de résultat pour cette région
    return {
        "region_id": region_id, "chrom": chrom,
        "strict_start": strict_start, "strict_end": strict_end,
        "window_start": win_start, "window_end": win_end, "window_size": WINDOW,
        "block_size": BLOCK_SIZE,
        "P1": P1, "P2": P2, "P3": P3, "O": OUT,
        "D": D, "SE": "" if SE is None else SE, "Z": "" if Z is None else Z,
        "ABBA": ABBA, "BABA": BABA, "n_snps": n_snps, "n_blocks": n_blocks,
    }


if __name__ == "__main__":
    results = []  # résultats accumulés pour toutes les régions
    for reg in REGIONS:  # traite chaque région candidate
        r = run_region(reg)
        if r is not None:
            results.append(r)  # ignore les régions écartées (trop peu de SNPs)

    out = OUTDIR / "Dstat_9regions_150kb_blockjackknife.tsv"  # fichier de sortie final
    cols = ["region_id", "chrom", "strict_start", "strict_end", "window_start", "window_end",
            "window_size", "block_size", "P1", "P2", "P3", "O", "D", "SE", "Z",
            "ABBA", "BABA", "n_snps", "n_blocks"]  # colonnes du TSV de sortie

    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")  # écrivain TSV basé sur les colonnes définies
        writer.writeheader()  # écrit la ligne d'en-tête
        for r in results:
            writer.writerow(r)  # écrit une ligne par région

    print(f"\n{'='*70}\nTerminé — table écrite : {out}")
    print(f"Régions traitées : {len(results)}/{len(REGIONS)}")
