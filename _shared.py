#!/usr/bin/env python3
"""Fonctions communes réutilisées par plusieurs scripts du dépôt (résolution du
dossier de données, extraction de génotypes VCF, metadata, pi), pour éviter les
copies quasi identiques qui existaient auparavant dans fd_genomewide.py, dxy.py,
pi_regions_qc.py, pi_genomewide_baseline.py, heterozygotie_consanguinite.py,
heatmap_haplotypes_arbre_ld.py et dstat_blockjackknife.py.

`default_project()` est importée par tous les scripts Python du dépôt : c'est elle
qui garantit qu'aucun chemin de données n'est codé en dur (voir sa docstring).

Import depuis un script du dépôt (quel que soit son sous-dossier) :

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from _shared import gt_to_alt_count, ...

`Path(__file__).resolve().parents[1]` pointe toujours vers la racine du
dépôt (le script est à un niveau de profondeur), indépendamment du
répertoire courant depuis lequel le script est lancé.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MIN_CALLED_FRAC_GROUP = 0.50  # fraction minimale d'individus génotypés requise par groupe (voir allele_freq_alt)


def default_project():
    """Dossier racine des données (contient data/ et analyses/, hors dépôt).

    Résolution, dans l'ordre : argument `--project <dossier>` de la ligne de
    commande, sinon la variable d'environnement AWASSI_PROJECT_DIR, sinon le
    répertoire courant. Aucun chemin n'est donc codé en dur dans les scripts :
    seuls les sous-chemins relatifs à cette racine le sont.

    La lecture de `--project` se fait directement dans sys.argv (et non via
    argparse) pour rester utilisable au niveau module, y compris quand le
    script est importé par un autre script du dépôt.
    """
    argv = sys.argv[1:]
    if "--project" in argv:
        i = argv.index("--project")
        if i + 1 < len(argv):
            return Path(argv[i + 1])
    return Path(os.environ.get("AWASSI_PROJECT_DIR", str(Path.cwd())))


def run(cmd):
    """Exécute `cmd` et retourne sa sortie standard (texte)."""
    return subprocess.check_output(cmd, text=True)


def parse_bool(x):
    """Convertit une valeur texte de metadata (ex: "oui"/"1") en booléen."""
    return str(x).strip().lower() in {"true", "1", "yes", "oui"}


def read_list(path):
    """Lit un fichier liste (un échantillon par ligne, lignes vides ignorées)."""
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def gt_to_alt_count(gt):
    """Convertit un génotype VCF (ex: "0/1") en nombre de copies de l'allèle alternatif (0, 1 ou 2, NaN si manquant/ambigu)."""
    gt = str(gt).split(":")[0]  # garde seulement le champ GT (avant les ":")
    if gt in {"./.", ".|.", "."}:
        return np.nan

    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)

    if len(parts) != 2 or "." in parts:
        return np.nan

    try:
        return int(parts[0]) + int(parts[1])
    except Exception:
        return np.nan


def allele_freq_alt(counts, idx, min_called_frac=MIN_CALLED_FRAC_GROUP):
    """Fréquence de l'allèle alternatif pour un sous-ensemble d'échantillons.

    Parameters
    ----------
    counts : numpy.ndarray
        Nombre de copies de l'allèle alt par échantillon (voir gt_to_alt_count),
        NaN si non génotypé.
    idx : numpy.ndarray
        Indices des échantillons du groupe dans `counts`.
    min_called_frac : float
        Fraction minimale d'échantillons génotypés requise (sinon NaN).

    Returns
    -------
    tuple[float, int]
        Fréquence alt (NaN si le groupe est vide ou insuffisamment génotypé)
        et nombre d'échantillons effectivement génotypés.
    """
    vals = counts[idx]
    vals = vals[~np.isnan(vals)]

    n_total = len(idx)
    n_called = len(vals)

    if n_total == 0:
        return np.nan, n_called

    if n_called / n_total < min_called_frac:
        return np.nan, n_called

    return vals.sum() / (2 * n_called), n_called  # fréquence alt = somme des comptes / (2 x nb génotypés)


def pi_from_values(vals):
    """Estimateur non biaisé de pi pour un site, à partir de comptes alt (0/1/2) déjà filtrés non-NaN d'un groupe."""
    n_called = len(vals)
    if n_called < 2:
        return np.nan
    two_n = 2 * n_called
    p = vals.sum() / two_n
    return (two_n / (two_n - 1)) * 2 * p * (1 - p)  # correction petit-échantillon (2n/(2n-1))


def pi_site(alt_counts, idx):
    """pi_from_values restreint aux individus d'un groupe (idx). NaN si < 2 génotypés."""
    vals = alt_counts[idx]
    vals = vals[~np.isnan(vals)]
    return pi_from_values(vals)


def load_metadata(project):
    """Charge la table de metadata des échantillons.

    Cherche parmi les emplacements candidats sous `project` et utilise le
    premier fichier trouvé.

    Parameters
    ----------
    project : pathlib.Path
        Dossier racine des données (contient le sous-dossier analyses/).

    Returns
    -------
    pandas.DataFrame
        Metadata dédupliquée par échantillon, avec au moins les colonnes
        sample, group et is_awassi.
    """
    metadata_candidates = [
        project / "analyses/haplotype_heatmap/Awassi_haplo/data/metadata/sample_metadata_387_FST_groups.tsv",
        project / "analyses/synthese_resultats/nnt/06_NNT_indiv_pairwise_base_metadata/results/metadata_used.tsv",
    ]
    for p in metadata_candidates:
        if p.exists():
            meta = pd.read_csv(p, sep="\t")
            print(f"Metadata utilisée : {p}")
            break
    else:
        raise FileNotFoundError("Aucune metadata trouvée.")

    if "sample_id" in meta.columns:
        meta = meta.rename(columns={"sample_id": "sample"})
    if "fst_group" in meta.columns:
        meta = meta.rename(columns={"fst_group": "group"})

    if "sample" not in meta.columns or "group" not in meta.columns:
        raise ValueError("Il faut au minimum les colonnes sample et group dans la metadata.")

    meta["sample"] = meta["sample"].astype(str)
    meta["group"] = meta["group"].astype(str)

    if "is_awassi" in meta.columns:
        meta["is_awassi"] = meta["is_awassi"].apply(parse_bool)
    else:
        meta["is_awassi"] = meta["group"].eq("Awassi")

    return meta.drop_duplicates("sample")


def get_vcf(chrom, project):
    """Trouve le VCF du chromosome demandé sous `project`.

    Parameters
    ----------
    chrom : str
        Nom du chromosome (ex: "16").
    project : pathlib.Path
        Dossier racine des données.

    Returns
    -------
    pathlib.Path
        Chemin du premier VCF candidat qui existe.
    """
    candidates = [
        project / f"data/raw data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
        project / f"data/raw_data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"VCF introuvable pour chr{chrom} : {candidates}")
