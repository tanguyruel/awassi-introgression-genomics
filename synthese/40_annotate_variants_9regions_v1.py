#!/usr/bin/env python3
"""
Script 40 — Annotation fonctionnelle des SNP des 9 régions candidates (étape 1
de l'annotation, zoom géo/race/individu à venir dans les scripts suivants).

Généralise la méthode déjà établie pour KIT (scripts/kit/02_annotate_KIT_variants_gene_features.py) :
parsing du GFF Oar_v4.0, classification de chaque SNP par feature chevauchante
(CDS > exon > mRNA/transcript (intron) > gene sans transcript > intergénique),
gène le plus proche + distance si intergénique (recherche dans une fenêtre
élargie ±300kb autour de la région stricte).

Sources (aucune donnée brute modifiée) :
  - VCF par région : analyses/synthese_resultats/annotation_9regions/vcf/{region_id}.vcf.gz
    (extraits par bcftools view -r depuis les VCF filtrés PASS+biallélique déjà utilisés
    dans tout le pipeline : analyses/fst/local_20kb_step5kb_sep_EU_AM_AUS_v1/filtered_vcf/)
  - GFF : data/reference/Oar_v4.0/GCF_000298735.2_Oar_v4.0_genomic.gff.gz

Sorties : analyses/synthese_resultats/annotation_9regions/tables/
  - variants_all_annotated.tsv       (tous les SNP des 9 régions, 1 ligne par SNP)
  - variants_summary_by_region.tsv   (résumé feature_class x gène x région)
  - variants_CDS_exonic.tsv          (sous-ensemble CDS/exonique, si présent)

Script de référence pour "Annotation fonctionnelle (GFF)" (cf. Annexe A du rapport de
stage) — généralise scripts/kit/02_annotate_KIT_variants_gene_features.py aux 9 régions.
Usage : python3 40_annotate_variants_9regions_v1.py
"""

# gzip : lit le GFF compressé ; subprocess : appelle bcftools
import gzip
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("analyses/synthese_resultats/annotation_9regions")
VCF_DIR = BASE / "vcf"
OUT_DIR = BASE / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si absent

GFF = Path("data/reference/Oar_v4.0/GCF_000298735.2_Oar_v4.0_genomic.gff.gz")
FLANK = 300_000  # fenêtre élargie pour trouver le gène le plus proche si intergénique

# les 9 régions candidates avec leurs coordonnées exactes et leur P3 déjà identifié
REGIONS = [
    {"region_id": "chr2_112.8Mb_NIPA2_CYFIP1",  "chr": "2",  "start": 112785001, "end": 112865000, "P3_best": "Australia"},
    {"region_id": "chr3_129.2Mb_desert",         "chr": "3",  "start": 129220001, "end": 129260000, "P3_best": "Europe"},
    {"region_id": "chr3_137.4Mb_LOC101123547",   "chr": "3",  "start": 137390001, "end": 137420000, "P3_best": "America"},
    {"region_id": "chr3_186.1Mb_CCDC91",         "chr": "3",  "start": 186075001, "end": 186125000, "P3_best": "Asia"},
    {"region_id": "chr5_58.1Mb_ABLIM3_AFAP1L1",  "chr": "5",  "start": 58070001,  "end": 58110000,  "P3_best": "Europe"},
    {"region_id": "chr6_70.2Mb_KIT",             "chr": "6",  "start": 70240001,  "end": 70295000,  "P3_best": "Africa"},
    {"region_id": "chr10_49.4Mb_KLF12",          "chr": "10", "start": 49360001,  "end": 49390000,  "P3_best": "Asia"},
    {"region_id": "chr17_34.2Mb_SPATA5",         "chr": "17", "start": 34160001,  "end": 34220000,  "P3_best": "Asia"},
    {"region_id": "chr20_0.8Mb_KHDRBS2",         "chr": "20", "start": 765001,    "end": 845000,    "P3_best": "Europe"},
]

# Reprise telle quelle de scripts/02_fst/09_annotation/09_extract_genes_candidate_regions_OarV4.py
# table de correspondance numéro de chromosome -> identifiant RefSeq utilisé dans le GFF
CHR_TO_ACCESSION = {
    "1": "NC_019458.2", "2": "NC_019459.2", "3": "NC_019460.2", "4": "NC_019461.2",
    "5": "NC_019462.2", "6": "NC_019463.2", "7": "NC_019464.2", "8": "NC_019465.2",
    "9": "NC_019466.2", "10": "NC_019467.2", "11": "NC_019468.2", "12": "NC_019469.2",
    "13": "NC_019470.2", "14": "NC_019471.2", "15": "NC_019472.2", "16": "NC_019473.2",
    "17": "NC_019474.2", "18": "NC_019475.2", "19": "NC_019476.2", "20": "NC_019477.2",
    "21": "NC_019478.2", "22": "NC_019479.2", "23": "NC_019480.2", "24": "NC_019481.2",
    "25": "NC_019482.2", "26": "NC_019483.2", "X": "NC_019484.2",
}


def parse_attr(attr):
    # transforme la colonne 9 du GFF ("ID=x;gene=y;...") en dictionnaire {clé: valeur}
    d = {}
    for item in attr.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


# ── 1. Charger les features GFF utiles (uniquement les 7 chromosomes, fenêtre élargie) ──
print("Chargement GFF (chromosomes des 9 régions uniquement)...")
windows = []  # une fenêtre (accession, début-FLANK, fin+FLANK) par région, pour filtrer le GFF
for r in REGIONS:
    acc = CHR_TO_ACCESSION[r["chr"]]
    windows.append((acc, r["start"] - FLANK, r["end"] + FLANK))

features = []
with gzip.open(GFF, "rt") as f:  # lecture directe du GFF gzippé, ligne par ligne (fichier volumineux)
    for line in f:
        if line.startswith("#"):
            continue  # ligne d'en-tête/commentaire GFF
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue  # ligne mal formée, ignorée
        seqid, source, ftype, start, end, score, strand, phase, attr = p  # 9 colonnes standard du GFF
        if ftype == "region":
            continue  # feature structurelle couvrant tout le chromosome (GFF RefSeq), pas un gène
        start, end = int(start), int(end)
        for acc, wstart, wend in windows:
            if seqid == acc and end >= wstart and start <= wend:  # la feature chevauche une des 9 fenêtres élargies
                a = parse_attr(attr)
                features.append({
                    "seqid": seqid, "feature": ftype, "start": start, "end": end,
                    "strand": strand, "gene": a.get("gene", a.get("Name", "")),  # nom du gène si présent
                })
                break  # inutile de tester les autres fenêtres, la feature est déjà retenue

features_df = pd.DataFrame(features)
genes_df = features_df[features_df["feature"] == "gene"].copy()  # sous-table des seules features "gene"
print(f"  {len(features_df)} features chargées, {len(genes_df)} gènes.")


def nearest_gene(seqid, pos):
    # cherche le gène le plus proche d'une position, sur le même chromosome
    g = genes_df[genes_df["seqid"] == seqid]
    if g.empty:
        return "", np.nan, ""  # aucun gène connu sur ce chromosome dans la fenêtre chargée
    inside = g[(g["start"] <= pos) & (g["end"] >= pos)]
    if not inside.empty:
        row = inside.iloc[0]
        return row["gene"], 0, "inside"  # la position est à l'intérieur d'un gène, distance nulle
    dist_up = pos - g["end"]  # distance si le gène est en amont (avant) de la position
    dist_down = g["start"] - pos  # distance si le gène est en aval (après) de la position
    # ne garde que la distance positive pertinente (le gène est bien d'un seul côté)
    g = g.assign(dist=np.where(dist_up > 0, dist_up, np.where(dist_down > 0, dist_down, np.nan)))
    g = g.dropna(subset=["dist"])
    if g.empty:
        return "", np.nan, ""
    row = g.loc[g["dist"].idxmin()]  # le gène avec la plus petite distance
    side = "upstream" if row["start"] > pos else "downstream"  # position du gène par rapport au SNP
    return row["gene"], int(row["dist"]), side


def annotate_pos(seqid, pos):
    # classe une position du génome : quelle(s) feature(s) GFF la couvrent exactement ?
    ov = features_df[(features_df["seqid"] == seqid) & (features_df["start"] <= pos) & (features_df["end"] >= pos)]
    if ov.empty:
        # aucune feature ne couvre cette position -> intergénique, on cherche le gène le plus proche
        gene, dist, side = nearest_gene(seqid, pos)
        return pd.Series({
            "feature_class": "intergenic", "gene": gene,
            "nearest_gene_distance_bp": dist, "nearest_gene_side": side,
        })
    genes_here = sorted(set(g for g in ov["gene"] if g))  # gène(s) chevauchant cette position
    feats_here = set(ov["feature"])  # types de features chevauchant (CDS, exon, mRNA, gene...)
    # priorité décroissante : CDS > exon > intron/transcrit > gène sans transcrit > autre
    if "CDS" in feats_here:
        fclass = "CDS"
    elif "exon" in feats_here:
        fclass = "exon_non_CDS_or_UTR"
    elif {"mRNA", "transcript"} & feats_here:
        fclass = "intron_or_non_exonic_transcript"
    elif "gene" in feats_here:
        fclass = "gene_region_no_transcript_feature"
    else:
        fclass = ",".join(sorted(feats_here))  # cas rare, type de feature non prévu ci-dessus
    return pd.Series({
        "feature_class": fclass, "gene": ",".join(genes_here),
        "nearest_gene_distance_bp": 0, "nearest_gene_side": "inside",
    })


# ── 2. Annoter les SNP de chaque région ──────────────────────────────────────
all_rows = []
for r in REGIONS:
    acc = CHR_TO_ACCESSION[r["chr"]]
    vcf = VCF_DIR / f'{r["region_id"]}.vcf.gz'
    out = subprocess.run(  # liste chrom/position/ref/alt de tous les SNP du VCF de la région
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf)],
        capture_output=True, text=True, check=True,
    ).stdout
    print(f'  {r["region_id"]} : annotation en cours...')
    for line in out.splitlines():
        chrom, pos, ref, alt = line.split("\t")
        pos = int(pos)
        ann = annotate_pos(acc, pos)  # classification fonctionnelle de ce SNP
        all_rows.append({
            "region_id": r["region_id"], "chr": r["chr"], "pos": pos,
            "ref": ref, "alt": alt, "P3_best": r["P3_best"],
            **ann.to_dict(),
        })

df = pd.DataFrame(all_rows)  # une ligne par SNP annoté, toutes régions confondues
out_all = OUT_DIR / "variants_all_annotated.tsv"
df.to_csv(out_all, sep="\t", index=False)

# résumé : nombre de SNP par région x classe fonctionnelle x gène, triés du plus au moins nombreux
summary = (
    df.groupby(["region_id", "feature_class", "gene"])
    .size().reset_index(name="n_variants")
    .sort_values(["region_id", "n_variants"], ascending=[True, False])
)
out_sum = OUT_DIR / "variants_summary_by_region.tsv"
summary.to_csv(out_sum, sep="\t", index=False)

exonic = df[df["feature_class"].isin(["CDS", "exon_non_CDS_or_UTR"])].copy()  # sous-ensemble codant/exonique
out_cds = OUT_DIR / "variants_CDS_exonic.tsv"
exonic.to_csv(out_cds, sep="\t", index=False)

print(f"\nTotal SNP annotés : {len(df)}")
print(f"→ {out_all}")
print(f"→ {out_sum}")
print(f"→ {out_cds} ({len(exonic)} variants CDS/exoniques)")
print()
print("Résumé par région (top feature_class) :")
print(summary.groupby("region_id").head(3).to_string(index=False))  # aperçu : 3 classes principales par région
