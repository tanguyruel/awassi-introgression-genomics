#!/usr/bin/env python3
"""Génère un petit jeu de données synthétique avec un signal d'introgression connu.

But : permettre de tester les scripts region_selection/fd_genomewide.py et
validation_tests/dxy.py (et tout autre script suivant la même convention
--project) sans avoir besoin des données réelles du stage.

Le jeu simulé contient un seul "chromosome" (TEST1, 300 kb) et 5 groupes
domestiques (Awassi + 4 groupes partenaires) plus un outgroup sauvage
(Ovis_canadensis), suffisant pour reconstituer le quadruplet P1/P2=Awassi/P3/O
utilisé par les tests d'introgression du dépôt.

Modèle :
- Fréquence ancestrale par SNP tirée dans Beta(0.4, 0.4).
- Fréquence par groupe domestique dérivée par un modèle de Balding-Nichols
  (dérive génétique, Fst ~ 0.06 — différenciation modeste, réaliste entre
  races domestiques).
- Dans une fenêtre choisie (140 000-175 000 pb), la fréquence du groupe
  Awassi est mélangée avec celle du groupe "Asia" (fraction ~0.75) pour
  simuler un bloc introgressé : Awassi y devient anormalement proche
  d'Asia, comme attendu après un flux génique local.
- L'outgroup Ovis_canadensis est simulé fixé (tous les individus homozygotes
  pour un allèle tiré aléatoirement par site), ce qui garantit une
  polarisation ancestral/dérivé toujours possible (hypothèse simplificatrice
  : pas de tri de lignées incomplet modélisé).
- ~1% de génotypes manquants (./.) par site, pour rester réaliste sans nuire
  à la détection du signal.

Ce jeu ne modélise pas de déséquilibre de liaison (chaque SNP est simulé
indépendamment) : il convient aux méthodes basées sur des fréquences
alléliques par fenêtre (fd, D-stat, dXY, FST) mais pas aux méthodes qui
ont besoin d'haplotypes réalistes (LD, phasage, partage d'haplotypes).

Sortie (sous ce dossier) :
  data/raw_data_08_06/awassi_and_basedata_chrTEST1.vcf.gz(.tbi)
  analyses/haplotype_heatmap/Awassi_haplo/data/metadata/sample_metadata_387_FST_groups.tsv
  regions_test.tsv

Usage : python3 generate_test_dataset.py
"""
from pathlib import Path
import gzip
import subprocess

import numpy as np
import pandas as pd

SEED = 42
CHROM = "TEST1"
CHR_LEN = 300_000
N_SNPS = 1_800

INTROGRESSED_START = 140_000
INTROGRESSED_END = 175_000
INTROGRESSION_FRACTION_MEAN = 0.75
INTROGRESSION_FRACTION_SD = 0.05

FST_DOMESTIC = 0.06
MISSING_RATE = 0.01

GROUPS = {
    "Awassi": 10,
    "Asia": 10,               # vrai donneur de l'introgression simulée
    "Africa": 10,
    "Europe": 10,
    "MiddleEastNonAwassi": 10,
}
OUTGROUP_N = 4  # Ovis_canadensis

HERE = Path(__file__).resolve().parent
VCF_DIR = HERE / "data" / "raw_data_08_06"
METADATA_PATH = HERE / "analyses" / "haplotype_heatmap" / "Awassi_haplo" / "data" / "metadata" / "sample_metadata_387_FST_groups.tsv"
REGIONS_PATH = HERE / "regions_test.tsv"

NUCLEOTIDES = ["A", "C", "G", "T"]


def balding_nichols(rng, p_anc, fst):
    """Tire une fréquence de groupe à partir d'une fréquence ancestrale (dérive génétique, modèle de Balding-Nichols)."""
    a = p_anc * (1 - fst) / fst
    b = (1 - p_anc) * (1 - fst) / fst
    return rng.beta(a, b)


def build_samples():
    """Construit la liste ordonnée des échantillons et leur groupe (domestiques puis outgroup)."""
    samples = []
    for group, n in GROUPS.items():
        for i in range(1, n + 1):
            samples.append((f"{group}_{i:02d}", group))
    for i in range(1, OUTGROUP_N + 1):
        samples.append((f"OCAN_{i:02d}", "Ovis_canadensis"))
    return samples


def simulate_genotypes(rng, samples):
    """Simule les génotypes de tous les SNP pour tous les échantillons.

    Returns
    -------
    positions : np.ndarray (N_SNPS,)
    ref, alt : listes de N_SNPS caractères
    genotypes : np.ndarray (N_SNPS, n_samples) de chaînes "0/0"/"0/1"/"1/1"/"./."
    """
    positions = np.sort(rng.choice(np.arange(1, CHR_LEN + 1), size=N_SNPS, replace=False))

    domestic_groups = [g for g in GROUPS]
    sample_groups = [g for _, g in samples]
    n_samples = len(samples)

    genotypes = np.empty((N_SNPS, n_samples), dtype=object)
    ref_alleles = []
    alt_alleles = []

    for i, pos in enumerate(positions):
        ref, alt = rng.choice(NUCLEOTIDES, size=2, replace=False)
        ref_alleles.append(ref)
        alt_alleles.append(alt)

        p_anc = np.clip(rng.beta(0.4, 0.4), 0.02, 0.98)

        group_freq = {}
        for g in domestic_groups:
            group_freq[g] = np.clip(balding_nichols(rng, p_anc, FST_DOMESTIC), 0.01, 0.99)

        in_window = INTROGRESSED_START <= pos <= INTROGRESSED_END
        if in_window:
            m = np.clip(rng.normal(INTROGRESSION_FRACTION_MEAN, INTROGRESSION_FRACTION_SD), 0.5, 0.95)
            group_freq["Awassi"] = (1 - m) * group_freq["Awassi"] + m * group_freq["Asia"]

        outgroup_ancestral_allele = rng.integers(0, 2)  # 0 = REF fixé, 1 = ALT fixé chez l'outgroup

        for j, (sample, group) in enumerate(samples):
            if group == "Ovis_canadensis":
                a1 = a2 = outgroup_ancestral_allele
            else:
                p = group_freq[group]
                a1 = int(rng.random() < p)
                a2 = int(rng.random() < p)

            if rng.random() < MISSING_RATE:
                genotypes[i, j] = "./."
            else:
                genotypes[i, j] = f"{a1}/{a2}"

    return positions, ref_alleles, alt_alleles, genotypes


def write_vcf(samples, positions, ref_alleles, alt_alleles, genotypes):
    """Écrit le VCF (texte), le bgzippe et l'indexe avec tabix."""
    VCF_DIR.mkdir(parents=True, exist_ok=True)
    vcf_plain = VCF_DIR / f"awassi_and_basedata_chr{CHROM}.vcf"
    vcf_gz = VCF_DIR.with_suffix("")  # placeholder, non utilisé

    sample_names = [s for s, _ in samples]

    header = [
        "##fileformat=VCFv4.2",
        "##source=generate_test_dataset.py (jeu de test synthétique, non biologique)",
        f"##contig=<ID={CHROM},length={CHR_LEN}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names),
    ]

    with open(vcf_plain, "w") as f:
        f.write("\n".join(header) + "\n")
        for i, pos in enumerate(positions):
            row = [
                CHROM, str(int(pos)), ".", ref_alleles[i], alt_alleles[i],
                "99", "PASS", ".", "GT",
            ] + list(genotypes[i])
            f.write("\t".join(row) + "\n")

    out_gz = VCF_DIR / f"awassi_and_basedata_chr{CHROM}.vcf.gz"
    subprocess.run(["bgzip", "-f", str(vcf_plain)], check=True)
    (vcf_plain.parent / (vcf_plain.name + ".gz")).rename(out_gz)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(out_gz)], check=True)
    return out_gz


def write_metadata(samples):
    """Écrit la table de metadata (sample_id, fst_group, is_awassi)."""
    METADATA_PATH.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"sample_id": s, "fst_group": g, "is_awassi": (g == "Awassi")}
        for s, g in samples
    ]
    pd.DataFrame(rows).to_csv(METADATA_PATH, sep="\t", index=False)


def write_regions():
    """Écrit le fichier de régions pour dxy.py : la fenêtre introgressée + une région témoin
    (hors introgression, pour comparaison — même logique que la validation multi-critères
    du pipeline principal)."""
    pd.DataFrame([
        {
            "region_id": "test_region_awassi_asia",
            "chr": CHROM,
            "start": INTROGRESSED_START,
            "end": INTROGRESSED_END,
            "highlight": "yes",
            "title": "Signal test introgression Awassi x Asia (synthétique)",
        },
        {
            "region_id": "control_region",
            "chr": CHROM,
            "start": 20_000,
            "end": 55_000,
            "highlight": "no",
            "title": "Région témoin, hors introgression (synthétique)",
        },
    ]).to_csv(REGIONS_PATH, sep="\t", index=False)


def main():
    rng = np.random.default_rng(SEED)
    samples = build_samples()
    positions, ref_alleles, alt_alleles, genotypes = simulate_genotypes(rng, samples)
    vcf_path = write_vcf(samples, positions, ref_alleles, alt_alleles, genotypes)
    write_metadata(samples)
    write_regions()

    print(f"VCF     : {vcf_path}")
    print(f"Metadata: {METADATA_PATH}")
    print(f"Regions : {REGIONS_PATH}")
    print(f"{N_SNPS} SNPs sur {CHR_LEN} pb, fenêtre introgressée : "
          f"{INTROGRESSED_START}-{INTROGRESSED_END} (Awassi <- Asia)")


if __name__ == "__main__":
    main()
