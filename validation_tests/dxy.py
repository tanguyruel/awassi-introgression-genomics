#!/usr/bin/env python3
#
# Calcule dXY (Awassi vs chaque groupe cible, dont Ovis_canadensis) en
# fenêtres 20kb/step5kb, sur une liste de régions candidates fournie en entrée.
# Entrée : --regions <tsv avec colonnes region_id/chr/start/end/highlight/title>,
#          VCF correspondants, metadata
# Sortie : <outdir>/<region_id>_dxy_20kb_step5kb.tsv, _dxy_by_site.tsv.gz,
#          _dxy_ranked_by_lowest.tsv
# Usage  : python3 02_compute_dxy_local_candidates.py --regions <regions.tsv> --outdir <outdir> [--project <dossier>]
# Le dossier racine des données (contenant data/ et analyses/, hors dépôt) se règle
# via --project, sinon la variable d'environnement AWASSI_PROJECT_DIR, sinon le
# répertoire courant.
#
from pathlib import Path
import argparse
import subprocess
import tempfile
import os
import sys
import re
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # racine du dépôt, pour importer _shared
from _shared import run, parse_bool, load_metadata, get_vcf, gt_to_alt_count, allele_freq_alt

WINDOW = 20_000
STEP = 5_000
MIN_SNPS = 20

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

def get_groups(meta, vcf_samples):
    """Construit le dictionnaire groupe -> échantillons (Awassi + groupes cibles).

    Parameters
    ----------
    meta : pandas.DataFrame
        Metadata des échantillons (issue de load_metadata).
    vcf_samples : list[str]
        Échantillons présents dans le VCF.

    Returns
    -------
    dict[str, list[str]]
        Groupes retenus (au moins 2 individus), restreints aux échantillons du VCF.
    """
    vcf_set = set(vcf_samples)
    groups = {}

    groups["Awassi"] = [
        s for s in meta.loc[meta["is_awassi"], "sample"].tolist()
        if s in vcf_set
    ]

    for g in TARGET_GROUPS:
        if g == "Ovis_canadensis":  # cas particulier : outgroup identifié par motif dans le nom
            groups[g] = [
                s for s in vcf_samples
                if re.search(r"OCAN|canad|bighorn|ovis_canadensis", s, flags=re.IGNORECASE)
            ]
        else:
            groups[g] = [
                s for s in meta.loc[meta["group"].eq(g), "sample"].tolist()
                if s in vcf_set
            ]

    print("Groupes utilisés :")
    for g, s in groups.items():
        print(f"  {g}: {len(s)} individus")

    if len(groups["Awassi"]) == 0:
        raise RuntimeError("Awassi vide.")
    if len(groups["Ovis_canadensis"]) == 0:
        raise RuntimeError("Ovis_canadensis vide.")

    return {g: s for g, s in groups.items() if len(s) >= 2}  # ne garde que les groupes avec au moins 2 individus

def compute_region(region, meta, outdir, project):
    """Calcule dXY (Awassi vs chaque groupe cible) pour une région candidate.

    Extrait les SNP bialléliques PASS de la région, calcule dXY par site puis
    moyenne en fenêtres glissantes (WINDOW/STEP), et écrit les tableaux
    résultats (par fenêtre, par site, classement) dans `outdir`.

    Parameters
    ----------
    region : pandas.Series
        Une ligne du tableau de régions (region_id, chr, start, end, ...).
    meta : pandas.DataFrame
        Metadata des échantillons (issue de load_metadata).
    outdir : pathlib.Path
        Dossier de sortie.
    project : pathlib.Path
        Dossier racine des données (pour localiser le VCF).
    """
    region_id = region["region_id"]
    chrom = str(region["chr"])
    start = int(region["start"])
    end = int(region["end"])

    print()
    print("============================================================")
    print(f"Région : {region_id}")
    print(f"chr{chrom}:{start}-{end}")
    print("============================================================")

    vcf = get_vcf(chrom, project)
    vcf_samples = run(["bcftools", "query", "-l", str(vcf)]).splitlines()

    groups = get_groups(meta, vcf_samples)

    selected = set()
    for samples in groups.values():
        selected.update(samples)

    selected_order = [s for s in vcf_samples if s in selected]  # ordre des échantillons sélectionnés, aligné sur le VCF
    sample_to_idx = {s: i for i, s in enumerate(selected_order)}

    group_idx = {}
    for g, samples in groups.items():
        group_idx[g] = np.array([sample_to_idx[s] for s in samples if s in sample_to_idx], dtype=int)

    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        for s in selected_order:
            tmp.write(s + "\n")
        sample_file = tmp.name

    region_str = f"{chrom}:{start}-{end}"

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
        os.remove(sample_file)

    site_rows = []
    n_total = 0
    n_skip_awassi = 0

    targets = [g for g in groups if g != "Awassi"]

    for line in txt.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue

        n_total += 1
        pos = int(parts[1])
        gts = parts[2:]

        alt_counts = np.array([gt_to_alt_count(gt) for gt in gts], dtype=float)

        pA, nA = allele_freq_alt(alt_counts, group_idx["Awassi"])
        if not np.isfinite(pA):
            n_skip_awassi += 1
            continue

        for g in targets:
            pB, nB = allele_freq_alt(alt_counts, group_idx[g])
            if not np.isfinite(pB):
                continue

            dxy_site = pA * (1 - pB) + (1 - pA) * pB  # divergence nucléotidique au site entre Awassi et le groupe cible

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

    site_df = pd.DataFrame(site_rows)

    print(f"SNPs lus : {n_total}")
    print(f"SNPs ignorés car Awassi insuffisant : {n_skip_awassi}")
    print(f"Lignes site dXY : {len(site_df)}")

    if site_df.empty:
        raise RuntimeError(f"Aucun SNP utilisable pour {region_id}")

    win_rows = []

    for wstart in range(start, end - WINDOW + 2, STEP):  # parcourt la région par fenêtres glissantes (pas = STEP)
        wend = wstart + WINDOW - 1
        wmid = (wstart + wend) // 2

        subw = site_df[(site_df["pos"] >= wstart) & (site_df["pos"] <= wend)]

        if subw.empty:
            continue

        for comp, subc in subw.groupby("comparison"):
            if len(subc) < MIN_SNPS:
                continue

            g = subc["group2"].iloc[0]

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

    win_df = pd.DataFrame(win_rows)

    if win_df.empty:
        raise RuntimeError(f"Aucune fenêtre dXY utilisable pour {region_id}")

    out_table = outdir / f"{region_id}_dxy_20kb_step5kb.tsv"
    out_sites = outdir / f"{region_id}_dxy_by_site.tsv.gz"
    ranked_path = outdir / f"{region_id}_dxy_ranked_by_lowest.tsv"

    win_df.to_csv(out_table, sep="\t", index=False)
    site_df.to_csv(out_sites, sep="\t", index=False, compression="gzip")
    win_df.sort_values(["dxy", "group2"]).to_csv(ranked_path, sep="\t", index=False)

    print(f"Table fenêtre : {out_table}")
    print(f"Table site    : {out_sites}")
    print(f"Ranking bas   : {ranked_path}")

def main():
    """Point d'entrée CLI : calcule dXY pour toutes les régions d'un tableau donné."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--project",
        default=os.environ.get("AWASSI_PROJECT_DIR", str(Path.cwd())),
        help="Dossier racine des données (contient data/ et analyses/). "
             "Par défaut : variable d'environnement AWASSI_PROJECT_DIR, sinon le répertoire courant.",
    )
    args = ap.parse_args()

    regions = pd.read_csv(args.regions, sep="\t")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    project = Path(args.project)

    required = {"region_id", "chr", "start", "end", "highlight", "title"}
    missing = required - set(regions.columns)
    if missing:
        raise RuntimeError(f"Colonnes manquantes dans régions : {missing}. Colonnes vues : {regions.columns.tolist()}")

    meta = load_metadata(project)

    for _, row in regions.iterrows():
        compute_region(row, meta, outdir, project)

if __name__ == "__main__":
    main()
