#!/usr/bin/env python3
#
# Calcule dXY (Awassi vs chaque groupe cible, dont Ovis_canadensis) en
# fenêtres 20kb/step5kb, sur une liste de régions candidates fournie en entrée.
# Entrée : --regions <tsv avec colonnes region_id/chr/start/end/highlight/title>,
#          VCF correspondants, metadata
# Sortie : <outdir>/<region_id>_dxy_20kb_step5kb.tsv, _dxy_by_site.tsv.gz,
#          _dxy_ranked_by_lowest.tsv
# Usage  : python3 02_compute_dxy_local_candidates.py --regions <regions.tsv> --outdir <outdir>
# IMPORTANT : PROJECT ci-dessous est un chemin absolu en dur (machine d'origine) — à adapter
# si le dépôt est cloné ailleurs.
#
# imports : chemins, arguments CLI, appels externes (bcftools), fichiers temporaires,
# expressions régulières, calcul numérique et tables
from pathlib import Path
import argparse
import subprocess
import tempfile
import os
import re
import numpy as np
import pandas as pd

PROJECT = Path("/home/tanguyruel/Bureau/genome_complet_Awassi")  # chemin en dur à adapter

# liste des chemins de metadata possibles (le premier trouvé est utilisé)
METADATA_CANDIDATES = [
    PROJECT / "analyses/haplotype_heatmap/Awassi_haplo/data/metadata/sample_metadata_387_FST_groups.tsv",
    PROJECT / "analyses/synthese_resultats/nnt/06_NNT_indiv_pairwise_base_metadata/results/metadata_used.tsv",
]

WINDOW = 20_000  # taille de fenêtre (pb)
STEP = 5_000  # pas entre fenêtres (pb)
MIN_SNPS = 20  # nombre minimal de SNP requis dans une fenêtre pour la garder
MIN_CALLED_FRAC_GROUP = 0.50  # fraction minimale d'individus génotypés requise par groupe

# groupes comparés à Awassi pour le calcul de dXY
TARGET_GROUPS = [
    "Asia",
    "Africa",
    "MiddleEastNonAwassi",
    "Europe",
    "America",
    "Australia",
    "Ovis_canadensis",
]

# exécute une commande externe et renvoie sa sortie texte
def run(cmd):
    return subprocess.check_output(cmd, text=True)

# convertit une valeur texte en booléen (true/1/yes/oui -> Vrai)
def parse_bool(x):
    return str(x).strip().lower() in {"true", "1", "yes", "oui"}

def load_metadata():  # charge la table de metadata des échantillons (premier fichier trouvé)
    for p in METADATA_CANDIDATES:  # parcourt les chemins de metadata candidats
        if p.exists():  # si le fichier existe
            meta = pd.read_csv(p, sep="\t")  # charge la table
            print(f"Metadata utilisée : {p}")
            break  # arrête la recherche dès qu'un fichier est trouvé
    else:
        raise FileNotFoundError("Aucune metadata trouvée.")  # aucun fichier trouvé : erreur

    if "sample_id" in meta.columns:  # harmonise le nom de colonne "sample_id" en "sample"
        meta = meta.rename(columns={"sample_id": "sample"})
    if "fst_group" in meta.columns:  # harmonise le nom de colonne "fst_group" en "group"
        meta = meta.rename(columns={"fst_group": "group"})

    meta["sample"] = meta["sample"].astype(str)  # force le type texte pour les ID échantillons
    meta["group"] = meta["group"].astype(str)  # force le type texte pour les groupes

    if "is_awassi" in meta.columns:  # si la colonne existe déjà, la convertit en booléen
        meta["is_awassi"] = meta["is_awassi"].apply(parse_bool)
    else:
        meta["is_awassi"] = meta["group"].eq("Awassi")  # sinon la déduit du nom de groupe "Awassi"

    return meta.drop_duplicates("sample")  # renvoie la metadata sans doublons d'échantillon

def get_vcf(chrom):  # trouve le chemin du VCF pour un chromosome donné
    # chemins possibles du VCF pour ce chromosome
    candidates = [
        PROJECT / f"data/raw data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
        PROJECT / f"data/raw_data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
    ]
    for p in candidates:  # teste chaque chemin
        if p.exists():
            return p  # renvoie le premier trouvé
    raise FileNotFoundError(f"VCF introuvable pour chr{chrom}")  # aucun VCF trouvé

def get_groups(meta, vcf_samples):  # construit le dictionnaire groupe -> échantillons (Awassi + groupes cibles)
    vcf_set = set(vcf_samples)  # ensemble des échantillons présents dans le VCF (recherche rapide)
    groups = {}  # dictionnaire de sortie groupe -> liste d'échantillons

    # échantillons Awassi (metadata) présents dans le VCF
    groups["Awassi"] = [
        s for s in meta.loc[meta["is_awassi"], "sample"].tolist()
        if s in vcf_set
    ]

    for g in TARGET_GROUPS:  # parcourt chaque groupe cible
        if g == "Ovis_canadensis":  # cas particulier : outgroup identifié par motif dans le nom
            # échantillons dont le nom ressemble à l'outgroup
            groups[g] = [
                s for s in vcf_samples
                if re.search(r"OCAN|canad|bighorn|ovis_canadensis", s, flags=re.IGNORECASE)
            ]
        else:
            # échantillons du groupe g (metadata) présents dans le VCF
            groups[g] = [
                s for s in meta.loc[meta["group"].eq(g), "sample"].tolist()
                if s in vcf_set
            ]

    print("Groupes utilisés :")
    for g, s in groups.items():
        print(f"  {g}: {len(s)} individus")

    if len(groups["Awassi"]) == 0:  # vérifie qu'il y a au moins un Awassi
        raise RuntimeError("Awassi vide.")
    if len(groups["Ovis_canadensis"]) == 0:  # vérifie qu'il y a au moins un outgroup
        raise RuntimeError("Ovis_canadensis vide.")

    return {g: s for g, s in groups.items() if len(s) >= 2}  # ne garde que les groupes avec au moins 2 individus

def gt_to_alt_count(gt):  # convertit un génotype VCF en nombre d'allèles alternatifs (0,1,2)
    gt = str(gt).split(":")[0]  # garde seulement le champ GT

    if gt in {"./.", ".|.", "."}:  # génotype manquant
        return np.nan

    sep = "|" if "|" in gt else "/"  # détecte le séparateur (phasé "|" ou non-phasé "/")
    parts = gt.split(sep)  # sépare les deux allèles

    if len(parts) != 2 or "." in parts:  # génotype invalide ou partiellement manquant
        return np.nan

    try:
        return int(parts[0]) + int(parts[1])  # somme des deux allèles (0/1/2 copies de l'allèle alternatif)
    except Exception:
        return np.nan

def allele_freq_alt(counts, idx):  # calcule la fréquence de l'allèle alternatif pour un groupe d'individus
    vals = counts[idx]  # comptes d'allèles alternatifs des individus du groupe
    vals = vals[~np.isnan(vals)]  # retire les valeurs manquantes

    n_total = len(idx)  # nombre total d'individus dans le groupe
    n_called = len(vals)  # nombre d'individus effectivement génotypés

    if n_total == 0:  # groupe vide
        return np.nan, n_called

    if n_called / n_total < MIN_CALLED_FRAC_GROUP:  # pas assez d'individus génotypés dans ce groupe
        return np.nan, n_called

    return vals.sum() / (2 * n_called), n_called  # fréquence alt = somme des comptes / (2 x nb génotypés)

def compute_region(region, meta, outdir):  # calcule dXY (Awassi vs chaque groupe cible) pour une région candidate
    region_id = region["region_id"]  # identifiant de la région
    chrom = str(region["chr"])  # chromosome de la région
    start = int(region["start"])  # début de la région
    end = int(region["end"])  # fin de la région

    print()
    print("============================================================")
    print(f"Région : {region_id}")
    print(f"chr{chrom}:{start}-{end}")
    print("============================================================")

    vcf = get_vcf(chrom)  # localise le VCF du chromosome
    vcf_samples = run(["bcftools", "query", "-l", str(vcf)]).splitlines()  # liste des échantillons présents dans le VCF

    groups = get_groups(meta, vcf_samples)  # détermine les groupes (Awassi + cibles) disponibles pour ce VCF

    selected = set()  # ensemble de tous les échantillons utilisés
    for samples in groups.values():  # parcourt les échantillons de chaque groupe
        selected.update(samples)

    selected_order = [s for s in vcf_samples if s in selected]  # ordre des échantillons sélectionnés, aligné sur le VCF
    sample_to_idx = {s: i for i, s in enumerate(selected_order)}  # index de chaque échantillon dans cet ordre

    group_idx = {}  # indices (dans selected_order) des échantillons de chaque groupe
    for g, samples in groups.items():  # parcourt chaque groupe
        group_idx[g] = np.array([sample_to_idx[s] for s in samples if s in sample_to_idx], dtype=int)  # indices pour ce groupe

    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:  # fichier temporaire listant les échantillons à extraire
        for s in selected_order:
            tmp.write(s + "\n")
        sample_file = tmp.name  # chemin du fichier temporaire

    region_str = f"{chrom}:{start}-{end}"  # région au format bcftools (chr:start-end)

    # commande bcftools : filtre SNP bialléliques PASS, restreint aux échantillons sélectionnés
    # et à la région -r, extrait CHROM/POS/GT par site
    cmd = (
        f"bcftools view "
        f"-f PASS -m2 -M2 -v snps "
        f"-S {sample_file} "
        f"-r {region_str} "
        f"'{vcf}' -Ou | "
        f"bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'"
    )

    print("Extraction :")
    print(cmd)

    try:
        txt = run(["bash", "-lc", cmd])  # exécute bcftools et récupère toute la sortie (région limitée, pas besoin de flux)
    finally:
        os.remove(sample_file)  # supprime le fichier temporaire dans tous les cas

    site_rows = []  # accumulateur des dXY par site et par comparaison
    n_total = 0  # compteur de sites lus
    n_skip_awassi = 0  # compteur de sites où Awassi n'est pas assez informatif

    targets = [g for g in groups if g != "Awassi"]  # tous les groupes cibles (hors Awassi)

    for line in txt.splitlines():  # parcourt chaque ligne (site) de la sortie bcftools
        parts = line.rstrip("\n").split("\t")  # sépare CHROM, POS et les génotypes
        if len(parts) < 3:  # ligne incomplète, on ignore
            continue

        n_total += 1
        pos = int(parts[1])  # position du site
        gts = parts[2:]  # génotypes bruts pour ce site

        alt_counts = np.array([gt_to_alt_count(gt) for gt in gts], dtype=float)  # convertit chaque génotype en nb d'allèles alternatifs

        pA, nA = allele_freq_alt(alt_counts, group_idx["Awassi"])  # fréquence alt et effectif chez Awassi
        if not np.isfinite(pA):  # Awassi non informatif à ce site
            n_skip_awassi += 1
            continue

        for g in targets:  # pour chaque groupe cible, calcule dXY au site avec Awassi
            pB, nB = allele_freq_alt(alt_counts, group_idx[g])  # fréquence alt et effectif chez le groupe cible
            if not np.isfinite(pB):  # groupe cible non informatif à ce site
                continue

            dxy_site = pA * (1 - pB) + (1 - pA) * pB  # divergence nucléotidique au site entre Awassi et le groupe cible

            # enregistre la ligne de dXY au site pour cette comparaison
            site_rows.append({
                "chr": chrom,
                "pos": pos,
                "comparison": f"Awassi_vs_{g}",
                "group2": g,
                "p_awassi_alt": pA,
                "p_group_alt": pB,
                "dxy_site": dxy_site,
                "n_called_awassi": nA,
                "n_called_group": nB,
            })

    site_df = pd.DataFrame(site_rows)  # assemble tous les sites en une table

    print(f"SNPs lus : {n_total}")
    print(f"SNPs ignorés car Awassi insuffisant : {n_skip_awassi}")
    print(f"Lignes site dXY : {len(site_df)}")

    if site_df.empty:  # aucun SNP utilisable pour cette région
        raise RuntimeError(f"Aucun SNP utilisable pour {region_id}")

    win_rows = []  # accumulateur des résultats par fenêtre

    for wstart in range(start, end - WINDOW + 2, STEP):  # parcourt la région par fenêtres glissantes (pas = STEP)
        wend = wstart + WINDOW - 1  # fin de la fenêtre
        wmid = (wstart + wend) // 2  # milieu de la fenêtre

        subw = site_df[(site_df["pos"] >= wstart) & (site_df["pos"] <= wend)]  # sites dont la position tombe dans cette fenêtre

        if subw.empty:  # fenêtre sans site, on l'ignore
            continue

        for comp, subc in subw.groupby("comparison"):  # regroupe par comparaison (Awassi vs chaque groupe)
            if len(subc) < MIN_SNPS:  # pas assez de SNP dans la fenêtre pour cette comparaison
                continue

            g = subc["group2"].iloc[0]  # nom du groupe cible de cette comparaison

            # enregistre le dXY moyen de la fenêtre pour cette comparaison
            win_rows.append({
                "chr": chrom,
                "window_start": wstart,
                "window_end": wend,
                "window_mid": wmid,
                "comparison": comp,
                "group2": g,
                "n_snps": len(subc),
                "dxy": subc["dxy_site"].mean(),
            })

    win_df = pd.DataFrame(win_rows)  # assemble les résultats par fenêtre en une table

    if win_df.empty:  # aucune fenêtre utilisable pour cette région
        raise RuntimeError(f"Aucune fenêtre dXY utilisable pour {region_id}")

    out_table = outdir / f"{region_id}_dxy_20kb_step5kb.tsv"  # chemin de la table par fenêtre
    out_sites = outdir / f"{region_id}_dxy_by_site.tsv.gz"  # chemin de la table détaillée par site
    ranked_path = outdir / f"{region_id}_dxy_ranked_by_lowest.tsv"  # chemin de la table triée du dXY le plus bas au plus haut

    win_df.to_csv(out_table, sep="\t", index=False)  # écrit la table par fenêtre
    site_df.to_csv(out_sites, sep="\t", index=False, compression="gzip")  # écrit la table détaillée par site (compressée)
    win_df.sort_values(["dxy", "group2"]).to_csv(ranked_path, sep="\t", index=False)  # trie par dXY croissant et écrit le classement

    print(f"Table fenêtre : {out_table}")
    print(f"Table site    : {out_sites}")
    print(f"Ranking bas   : {ranked_path}")

def main():  # point d'entrée du script
    ap = argparse.ArgumentParser()  # parseur d'arguments en ligne de commande
    ap.add_argument("--regions", required=True)  # fichier TSV des régions candidates
    ap.add_argument("--outdir", required=True)  # dossier de sortie
    args = ap.parse_args()  # parse les arguments fournis

    regions = pd.read_csv(args.regions, sep="\t")  # charge la table des régions candidates
    outdir = Path(args.outdir)  # dossier de sortie en objet Path
    outdir.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si besoin

    required = {"region_id", "chr", "start", "end", "highlight", "title"}  # colonnes obligatoires attendues
    missing = required - set(regions.columns)  # colonnes manquantes éventuelles
    if missing:  # arrête si des colonnes obligatoires manquent
        raise RuntimeError(f"Colonnes manquantes dans régions : {missing}. Colonnes vues : {regions.columns.tolist()}")

    meta = load_metadata()  # charge la metadata des échantillons

    for _, row in regions.iterrows():  # traite chaque région candidate une par une
        compute_region(row, meta, outdir)

if __name__ == "__main__":  # point d'entrée du script en exécution directe
    main()
