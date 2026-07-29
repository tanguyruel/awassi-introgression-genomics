#!/usr/bin/env python3
"""
D-stat + Z-score (block-jackknife) sur les 9 régions candidates,
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

import csv
import math
import subprocess
from pathlib import Path
from collections import defaultdict

POP_DIR = Path("analyses/fst/popmaps_separees_v1")  # dossier des popmaps par groupe (1 fichier liste d'échantillons par groupe)
POPMAP_MAIN = Path("data/popmap_main5.tsv")  # seule source pour Ovis_canadensis (absent de POP_DIR)
VCF_DIR = Path("data/raw data_08_06")
OUTDIR = Path("analyses/synthese_resultats/dstat_9regions_150kb")
OUTDIR.mkdir(parents=True, exist_ok=True)

P1 = "MiddleEastNonAwassi"  # population 1 (référence, non-Awassi Moyen-Orient)
P2 = "Awassi"               # population 2 (testée pour introgression)
OUT = "Ovis_canadensis"     # groupe extérieur (outgroup) servant à polariser les allèles

WINDOW = 150_000
BLOCK_SIZE = 15_000  # taille d'un bloc jackknife (bp) -> ~10 blocs par fenêtre
MIN_SNPS = 15
MIN_BLOCKS = 5       # nombre minimum de blocs jackknife pour calculer un SE fiable

OUT_LOW = 0.10   # fréquence alt outgroup en dessous : considéré comme allèle ancestral
OUT_HIGH = 0.90  # fréquence alt outgroup au-dessus : considéré comme allèle dérivé (polarité inversée)

# ── 9 régions strictes, fenêtre stricte FST/fd d'origine + P3 déjà identifié ──
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


def read_list(path):
    """Lit un fichier liste (un échantillon par ligne, lignes vides ignorées)."""
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def read_group_from_popmap(path, group):
    """Extrait la liste des échantillons appartenant à `group` depuis un popmap générique (colonnes échantillon/groupe)."""
    samples = []
    with open(path) as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            s, g = line.rstrip("\n").split("\t")[:2]  # colonnes échantillon (s) et groupe (g)
            if g == group:
                samples.append(s)
    return samples


def read_group(group):
    """P1/P2/P3 viennent des popmaps séparés (POP_DIR) ; OUT (Ovis_canadensis)
    n'existe que dans POPMAP_MAIN."""
    p = POP_DIR / f"{group}.txt"
    if p.exists():
        return read_list(p)
    return read_group_from_popmap(POPMAP_MAIN, group)


def alt_freq(gts, idxs):
    """Fréquence de l'allèle alternatif pour un sous-ensemble d'échantillons à un site donné (None si aucun génotype exploitable)."""
    alt = 0
    n = 0
    for i in idxs:
        gt = gts[i].replace("|", "/")
        if gt in {".", "./.", ".|."}:
            continue
        for a in gt.split("/"):
            if a == ".":
                continue
            if a == "0":
                n += 1
            elif a == "1":
                alt += 1
                n += 1
    return None if n == 0 else alt / n


def get_sample_order(vcf, sample_file):
    """Récupère l'ordre des échantillons tel qu'il apparaît dans le VCF restreint au fichier d'échantillons."""
    p1 = subprocess.Popen(["bcftools", "view", "-S", str(sample_file), "-Ou", vcf], stdout=subprocess.PIPE)
    p2 = subprocess.run(["bcftools", "query", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True, check=True)
    p1.stdout.close()
    p1.wait()
    return p2.stdout.strip().splitlines()


def jackknife(num, den, blocks):
    """Calcule SE (erreur-type) et Z-score du D-stat par block-jackknife delete-one-block.

    Parameters
    ----------
    num, den : float
        Numérateur et dénominateur du D-stat cumulés sur toute la fenêtre.
    blocks : dict[int, list]
        Par identifiant de bloc jackknife, [num, den, n_snps] accumulés sur ce bloc.

    Returns
    -------
    tuple[float | None, float | None, int]
        (SE, Z, nombre de blocs utilisés). SE et Z valent None si le
        dénominateur global est nul ou s'il n'y a pas assez de blocs
        (< MIN_BLOCKS) pour un jackknife fiable.
    """
    if den == 0:
        return None, None, 0
    D = num / den  # D-statistic global = (ABBA-BABA)/(ABBA+BABA)
    values = []  # pseudo-valeurs D obtenues en retirant un bloc à la fois
    for bnum, bden, bn in blocks.values():
        num_j = num - bnum  # numérateur total moins celui du bloc retiré
        den_j = den - bden  # dénominateur total moins celui du bloc retiré
        if den_j != 0:
            values.append(num_j / den_j)  # D recalculé sans ce bloc (delete-one-block)
    B = len(values)
    if B < MIN_BLOCKS:
        return None, None, B  # pas assez de blocs pour un SE fiable
    mean_j = sum(values) / B
    var = ((B - 1) / B) * sum((x - mean_j) ** 2 for x in values)  # variance jackknife (formule delete-one-block)
    SE = math.sqrt(var) if var > 0 else None
    Z = D / SE if SE and SE > 0 else None  # score Z = D / SE, teste si D diffère significativement de 0
    return SE, Z, B


def run_region(reg):
    """Calcule le D-stat + jackknife pour une région (fenêtre 150kb centrée sur le milieu de la région stricte).

    Parameters
    ----------
    reg : dict
        Une entrée de REGIONS (region_id, chr, strict_start, strict_end, P3).

    Returns
    -------
    dict | None
        Ligne de résultats (D, SE, Z, ABBA, BABA, n_snps, n_blocks, ...), ou
        None si trop peu de SNPs utilisables dans la fenêtre.
    """
    region_id, chrom, strict_start, strict_end, P3 = (
        reg["region_id"], reg["chr"], reg["strict_start"], reg["strict_end"], reg["P3"]
    )
    mid = (strict_start + strict_end) // 2
    win_start = mid - WINDOW // 2
    win_end = win_start + WINDOW - 1

    vcf = str(VCF_DIR / f"awassi_and_basedata_chr{chrom}.vcf.gz")
    print(f"\n{'='*70}\n{region_id} | chr{chrom}:{win_start}-{win_end} (150kb) | P3={P3}")
    print(f"  (fenêtre stricte d'origine : {strict_start}-{strict_end})")

    groups = {g: read_group(g) for g in [P1, P2, P3, OUT]}
    samples = sorted(set(groups[P1] + groups[P2] + groups[P3] + groups[OUT]))
    sample_file = OUTDIR / f"samples_{region_id}.txt"
    sample_file.write_text("\n".join(samples) + "\n")

    order = get_sample_order(vcf, sample_file)
    idx = {s: i for i, s in enumerate(order)}
    idx1 = [idx[s] for s in groups[P1] if s in idx]
    idx2 = [idx[s] for s in groups[P2] if s in idx]
    idx3 = [idx[s] for s in groups[P3] if s in idx]
    idxO = [idx[s] for s in groups[OUT] if s in idx]
    print(f"  P1={len(idx1)} P2={len(idx2)} P3={len(idx3)} OUT={len(idxO)}")

    ABBA = BABA = num_tot = den_tot = 0.0
    n_snps = 0
    blocks = defaultdict(lambda: [0.0, 0.0, 0])  # accumulateurs par bloc jackknife : [num, den, n_snps]

    pview = subprocess.Popen(
        ["bcftools", "view", "-r", f"{chrom}:{win_start}-{win_end}", "-f", "PASS",
         "-m2", "-M2", "-v", "snps", "-S", str(sample_file), "-Ou", vcf],
        stdout=subprocess.PIPE,
    )
    pquery = subprocess.Popen(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n"],
        stdin=pview.stdout, stdout=subprocess.PIPE, text=True,
    )
    pview.stdout.close()

    n_read = 0
    for line in pquery.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        pos = int(parts[1])
        gts = parts[2:]
        n_read += 1

        o = alt_freq(gts, idxO)
        p1a = alt_freq(gts, idx1)
        p2a = alt_freq(gts, idx2)
        p3a = alt_freq(gts, idx3)
        if o is None or p1a is None or p2a is None or p3a is None:
            continue

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

        ABBA += abba
        BABA += baba
        num_tot += num
        den_tot += den
        n_snps += 1

        block_id = int((pos - win_start) // BLOCK_SIZE)  # numéro du bloc jackknife contenant ce SNP
        blocks[block_id][0] += num
        blocks[block_id][1] += den
        blocks[block_id][2] += 1

    pquery.wait()
    pview.wait()

    print(f"  SNPs lus : {n_read} ; SNPs polarisés/utilisés : {n_snps} ; blocs : {len(blocks)}")

    if n_snps < MIN_SNPS or den_tot == 0:
        print(f"  [!] Trop peu de SNPs ({n_snps}) ou den=0 — région ignorée")
        return None

    D = num_tot / den_tot  # D-stat final sur toute la fenêtre = (ABBA-BABA)/(ABBA+BABA)
    SE, Z, n_blocks = jackknife(num_tot, den_tot, blocks)
    print(f"  D={D:.4f}  SE={SE}  Z={Z}  n_blocks={n_blocks}")

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
    results = []
    for reg in REGIONS:
        r = run_region(reg)
        if r is not None:
            results.append(r)

    out = OUTDIR / "Dstat_9regions_150kb_blockjackknife.tsv"
    cols = ["region_id", "chrom", "strict_start", "strict_end", "window_start", "window_end",
            "window_size", "block_size", "P1", "P2", "P3", "O", "D", "SE", "Z",
            "ABBA", "BABA", "n_snps", "n_blocks"]

    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"\n{'='*70}\nTerminé — table écrite : {out}")
    print(f"Régions traitées : {len(results)}/{len(REGIONS)}")
