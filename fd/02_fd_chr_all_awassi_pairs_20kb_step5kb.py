#!/usr/bin/env python3
#
# fd/D genome-wide par chromosome, fenêtre fine 20kb/step5kb, pour toutes les
# paires de groupes non-Awassi. Version courante utilisée par le pipeline de
# production (03_run_fd_genomewide_all_awassi_pairs_20kb_step5kb.sh et son
# complément 03_resume_..._skip_done.sh).
# Entrée : VCF chr (data/raw data_08_06/ ou raw_data_08_06/), metadata
# Sortie : <outdir>/fd_chr<N>_all_awassi_pairs_20kb_step5kb.tsv.gz
# Usage  : python3 02_fd_chr_all_awassi_pairs_20kb_step5kb.py --chrom <N> --outdir <outdir> [--groups g1,g2,...]
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
import gzip
import math
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
MIN_SNPS = 50  # nombre minimal de SNP requis dans une fenêtre pour la garder

OUTGROUP_FIXED_THRESHOLD = 0.05  # seuil de fréquence pour considérer l'outgroup comme fixé (proche de 0 ou 1)
MIN_CALLED_FRAC_GROUP = 0.50  # fraction minimale d'individus génotypés requise par groupe

# motifs (regex) de noms de groupes à exclure des partenaires (Awassi, outgroups sauvages, etc.)
EXCLUDE_GROUP_PATTERNS = [
    r"Awassi",
    r"Ovis",
    r"OCAN",
    r"canad",
    r"vignei",
    r"orientalis",
    r"mouflon",
    r"wild",
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

    if "sample" not in meta.columns or "group" not in meta.columns:  # vérifie que les colonnes minimales existent
        raise ValueError("Il faut au minimum les colonnes sample et group dans la metadata.")

    meta["sample"] = meta["sample"].astype(str)  # force le type texte pour les ID échantillons
    meta["group"] = meta["group"].astype(str)  # force le type texte pour les groupes

    if "is_awassi" in meta.columns:  # si la colonne existe déjà, la convertit en booléen
        meta["is_awassi"] = meta["is_awassi"].apply(parse_bool)
    else:
        meta["is_awassi"] = meta["group"].eq("Awassi")  # sinon la déduit du nom de groupe "Awassi"

    return meta.drop_duplicates("sample")  # renvoie la metadata sans doublons d'échantillon

def is_excluded_group(g):  # teste si un nom de groupe correspond à un motif exclu
    g = str(g)
    return any(re.search(pat, g, flags=re.IGNORECASE) for pat in EXCLUDE_GROUP_PATTERNS)  # vrai si un motif matche (insensible à la casse)

def get_vcf(chrom):  # trouve le chemin du VCF pour un chromosome donné
    # chemins possibles du VCF pour ce chromosome
    candidates = [
        PROJECT / f"data/raw data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
        PROJECT / f"data/raw_data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
    ]
    for p in candidates:  # teste chaque chemin
        if p.exists():
            return p  # renvoie le premier trouvé
    raise FileNotFoundError(f"VCF introuvable pour chr{chrom} : {candidates}")  # aucun VCF trouvé

def choose_groups(meta, vcf_samples, groups_arg):
    """
    Choix robuste des groupes pour fd.

    P2 = Awassi
    O  = Ovis_canadensis
    P3 = chaque groupe partenaire
    P1 = chaque autre groupe partenaire

    Important :
    - On ignore groups_arg s'il vaut un nombre type 1011, car GROUPS est une variable spéciale de bash.
    - On garde tous les groupes non-Awassi présents dans la metadata et dans le VCF.
    """
    vcf_set = set(vcf_samples)  # ensemble des échantillons présents dans le VCF (recherche rapide)

    # échantillons Awassi (metadata) présents dans le VCF
    awassi = [
        s for s in meta.loc[meta["is_awassi"], "sample"].tolist()
        if s in vcf_set
    ]

    # échantillons outgroup (Ovis canadensis) identifiés par motif dans leur nom
    outgroup = [
        s for s in vcf_samples
        if re.search(r"OCAN|canad|bighorn|ovis_canadensis", s, flags=re.IGNORECASE)
    ]

    # dictionnaire groupe -> liste d'échantillons, initialisé avec Awassi et l'outgroup
    groups = {
        "Awassi": awassi,
        "Ovis_canadensis": outgroup,
    }

    all_meta_groups = sorted([str(g) for g in meta["group"].dropna().unique()])  # tous les groupes présents dans la metadata

    # Si --groups est donné et n'est pas un nombre parasite type 1011, on le respecte.
    if groups_arg and not re.fullmatch(r"[0-9, ]+", str(groups_arg).strip()):  # --groups fourni et pas juste des chiffres (piège variable bash GROUPS)
        candidate_groups = [g.strip() for g in str(groups_arg).split(",") if g.strip()]  # utilise la liste de groupes fournie en argument
    else:
        # ordre préféré des groupes candidats
        preferred = [
            "Asia",
            "Africa",
            "MiddleEastNonAwassi",
            "Europe",
            "America",
            "Australia",
            "EuropeAmericaAustralia",
            "Other",
        ]

        candidate_groups = []  # liste finale des groupes candidats à construire

        for g in preferred:  # ajoute d'abord les groupes préférés présents dans la metadata
            if g in all_meta_groups and g not in candidate_groups:
                candidate_groups.append(g)

        for g in all_meta_groups:  # ajoute ensuite les groupes restants (hors Awassi/outgroup/motifs sauvages)
            if g not in candidate_groups and g not in {"Awassi", "Ovis_canadensis"}:
                if not re.search(r"OCAN|canad|bighorn|vignei|orientalis|mouflon|wild", g, flags=re.IGNORECASE):
                    candidate_groups.append(g)

    print()
    print("Groupes candidats testés comme partenaires :")
    for g in candidate_groups:
        print(" -", g)

    for g in candidate_groups:  # pour chaque groupe candidat, récupère ses échantillons présents dans le VCF
        samples = [
            s for s in meta.loc[meta["group"].astype(str).eq(g), "sample"].tolist()
            if s in vcf_set
        ]

        # On garde les groupes avec au moins 2 individus.
        if len(samples) >= 2:
            groups[g] = samples

    print()
    print("Groupes retenus :")
    for g, s in groups.items():
        print(f"{g}: {len(s)} individus")

    if len(groups["Awassi"]) == 0:  # vérifie qu'il y a au moins un Awassi
        raise RuntimeError("Awassi vide.")

    if len(groups["Ovis_canadensis"]) == 0:  # vérifie qu'il y a au moins un outgroup
        print("Noms wild/outgroups possibles dans le VCF :")
        for s in vcf_samples:  # affiche les noms d'échantillons ressemblant à un outgroup (aide au diagnostic)
            if re.search(r"OCAN|canad|bighorn|ovis|vignei|orientalis|mouflon|wild", s, flags=re.IGNORECASE):
                print(" -", s)
        raise RuntimeError("Ovis_canadensis vide.")

    partners = [g for g in groups if g not in {"Awassi", "Ovis_canadensis"}]  # groupes partenaires (hors Awassi et outgroup)

    if len(partners) < 2:  # il faut au moins 2 groupes partenaires pour comparer
        raise RuntimeError(f"Pas assez de groupes partenaires : {partners}")

    return groups, partners  # renvoie les groupes et la liste des partenaires

def build_comparisons(partners):
    """
    P2 = Awassi fixé.
    Pour chaque groupe partenaire P3, on teste tous les autres groupes comme P1.
    Cela donne tous les signaux Awassi × groupe, mais avec référence explicite.
    """
    comps = []  # liste des comparaisons P1/P2/P3/O à calculer
    for P3 in partners:  # chaque groupe partenaire tour à tour comme P3
        for P1 in partners:  # tous les autres comme P1 potentiel
            if P1 == P3:  # on ignore le cas P1 == P3
                continue
            # enregistre la comparaison (Awassi = P2 fixe, outgroup = Ovis_canadensis)
            comps.append({
                "comparison": f"{P1}__Awassi__{P3}__OCAN",
                "P1": P1,
                "P2": "Awassi",
                "P3": P3,
                "O": "Ovis_canadensis",
            })
    return comps  # renvoie la liste des comparaisons

def gt_to_alt_count(gt):  # convertit un génotype VCF (ex "0/1") en nombre d'allèles alternatifs (0,1,2)
    gt = str(gt).split(":")[0]  # garde seulement le champ GT (avant les ":")
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
    vals = counts[idx]  # comptes d'allèles alternatifs pour les individus du groupe
    vals = vals[~np.isnan(vals)]  # retire les valeurs manquantes

    n_total = len(idx)  # nombre total d'individus dans le groupe
    n_called = len(vals)  # nombre d'individus effectivement génotypés

    if n_total == 0:  # groupe vide
        return np.nan, n_called

    if n_called / n_total < MIN_CALLED_FRAC_GROUP:  # pas assez d'individus génotypés dans ce groupe
        return np.nan, n_called

    return vals.sum() / (2 * n_called), n_called  # fréquence alt = somme des comptes / (2 x nb génotypés)

def extract_and_polarize(chrom, vcf, groups, selected_order, outdir):  # extraction des génotypes VCF et polarisation par l'outgroup
    sample_to_idx = {s: i for i, s in enumerate(selected_order)}  # index de chaque échantillon dans l'ordre sélectionné

    group_idx = {}  # indices (dans selected_order) des échantillons de chaque groupe
    for g, samples in groups.items():  # parcourt chaque groupe
        group_idx[g] = np.array([sample_to_idx[s] for s in samples if s in sample_to_idx], dtype=int)  # indices pour ce groupe

    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:  # fichier temporaire listant les échantillons à extraire
        for s in selected_order:
            tmp.write(s + "\n")
        sample_file = tmp.name  # chemin du fichier temporaire

    # commande bcftools : filtre SNP bialléliques PASS, restreint aux échantillons sélectionnés,
    # extrait CHROM/POS/GT par site (un site par ligne)
    cmd = (
        f"bcftools view "
        f"-f PASS -m2 -M2 -v snps "
        f"-S {sample_file} "
        f"'{vcf}' -Ou | "
        f"bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'"
    )

    print()
    print("Extraction VCF :")
    print(cmd)

    # lance la commande bcftools en flux (stdout lu ligne par ligne, plus économe en mémoire)
    proc = subprocess.Popen(
        ["bash", "-lc", cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    rows = []  # accumulateur des fréquences dérivées par site
    n_total = 0  # compteur de sites lus
    n_outgroup_ambiguous = 0  # compteur de sites où l'outgroup n'est pas assez fixé
    n_missing_group = 0  # compteur de sites avec données insuffisantes dans un groupe

    domestic_groups = [g for g in groups if g not in {"Ovis_canadensis"}]  # tous les groupes sauf l'outgroup

    for line in proc.stdout:  # lit la sortie de bcftools ligne par ligne (un site par ligne)
        n_total += 1
        parts = line.rstrip("\n").split("\t")  # sépare CHROM, POS et les génotypes
        if len(parts) < 3:  # ligne incomplète, on ignore
            continue

        pos = int(parts[1])  # position du site
        gts = parts[2:]  # liste des génotypes bruts pour ce site

        alt_counts = np.array([gt_to_alt_count(gt) for gt in gts], dtype=float)  # convertit chaque génotype en nb d'allèles alternatifs

        pO_alt, _ = allele_freq_alt(alt_counts, group_idx["Ovis_canadensis"])  # fréquence alt chez l'outgroup
        if not np.isfinite(pO_alt):  # outgroup non informatif à ce site
            n_missing_group += 1
            continue

        # Polarisation avec OCAN.
        if pO_alt <= OUTGROUP_FIXED_THRESHOLD:  # outgroup fixé pour l'allèle référence : pas d'inversion
            flip = False
        elif pO_alt >= 1 - OUTGROUP_FIXED_THRESHOLD:  # outgroup fixé pour l'allèle alternatif : on inverse la polarité
            flip = True
        else:
            n_outgroup_ambiguous += 1  # outgroup ni fixé ref ni alt : site ambigu, ignoré
            continue

        rec = {"pos": pos}  # enregistrement des fréquences dérivées pour ce site
        ok = True  # indique si tous les groupes ont une fréquence valide

        for g in domestic_groups:  # calcule la fréquence dérivée pour chaque groupe non-outgroup
            p_alt, _ = allele_freq_alt(alt_counts, group_idx[g])
            if not np.isfinite(p_alt):  # groupe non informatif à ce site
                ok = False
                break

            rec[g] = 1 - p_alt if flip else p_alt  # applique l'inversion de polarité si besoin (fréquence de l'allèle dérivé)

        if not ok:  # un groupe au moins est non informatif : site ignoré
            n_missing_group += 1
            continue

        rows.append(rec)  # garde le site polarisé

    stderr = proc.stderr.read()  # récupère les messages d'erreur du sous-processus
    ret = proc.wait()  # attend la fin du sous-processus et récupère son code retour

    os.remove(sample_file)  # supprime le fichier temporaire d'échantillons

    if ret != 0:  # bcftools a échoué
        print(stderr)
        raise RuntimeError(f"bcftools a échoué pour chr{chrom}")

    freqs = pd.DataFrame(rows)  # assemble tous les sites polarisés en table

    print()
    print(f"chr{chrom} SNPs PASS bialléliques lus : {n_total}")
    print(f"chr{chrom} SNPs polarisés utilisables : {len(freqs)}")
    print(f"chr{chrom} SNPs ignorés outgroup ambigu : {n_outgroup_ambiguous}")
    print(f"chr{chrom} SNPs ignorés données groupe insuffisantes : {n_missing_group}")

    freqs_path = outdir / f"chr{chrom}_derived_freqs.tsv.gz"  # chemin de sauvegarde des fréquences dérivées
    freqs.to_csv(freqs_path, sep="\t", index=False, compression="gzip")  # écrit les fréquences en TSV compressé

    return freqs  # renvoie la table pour la suite du script

def compute_D_fd(p1, p2, p3):  # calcule les statistiques D (ABBA-BABA) et fd pour un jeu de fréquences P1/P2/P3
    ABBA = (1 - p1) * p2 * p3  # poids du motif ABBA (site dérivé chez P2 et P3, pas chez P1)
    BABA = p1 * (1 - p2) * p3  # poids du motif BABA (site dérivé chez P1 et P3, pas chez P2)

    num = np.nansum(ABBA - BABA)  # numérateur commun à D et fd (somme des différences ABBA-BABA)
    den_D = np.nansum(ABBA + BABA)  # dénominateur de la statistique D

    D = num / den_D if den_D > 0 else np.nan  # statistique D (proportion de sites ABBA vs BABA)

    # Martin et al. fd : denominator with dynamic donor PD.
    PD = np.maximum(p2, p3)  # fréquence du "donneur" dynamique (Martin et al. 2015)
    den_fd = np.nansum(PD * (PD - p1))  # dénominateur de fd, calculé avec PD au lieu de p3 fixe

    fd = num / den_fd if den_fd > 0 else np.nan  # statistique fd (estimateur robuste de la fraction introgressée)

    return D, fd, float(np.nansum(ABBA)), float(np.nansum(BABA)), float(num), float(den_D), float(den_fd)  # D, fd + quantités intermédiaires

def scan_fd(chrom, freqs, comparisons, outdir):  # parcourt le chromosome en fenêtres glissantes et calcule D/fd pour chaque comparaison
    if freqs.empty:  # aucun site polarisé disponible
        raise RuntimeError(f"Pas de SNPs polarisés pour chr{chrom}")

    freqs = freqs.sort_values("pos").reset_index(drop=True)  # trie les sites par position
    positions = freqs["pos"].to_numpy(dtype=int)  # tableau des positions triées

    # tableau numpy des fréquences par groupe (toutes les colonnes sauf "pos")
    group_arrays = {
        c: freqs[c].to_numpy(dtype=float)
        for c in freqs.columns
        if c != "pos"
    }

    max_pos = int(positions.max())  # dernière position couverte par les données
    rows = []  # accumulateur des résultats par fenêtre/comparaison

    for start in range(1, max_pos - WINDOW + 2, STEP):  # parcourt le chromosome par fenêtres glissantes (pas = STEP)
        end = start + WINDOW - 1  # fin de la fenêtre
        left = np.searchsorted(positions, start, side="left")  # index du premier SNP dans la fenêtre
        right = np.searchsorted(positions, end, side="right")  # index après le dernier SNP dans la fenêtre

        n_snps = right - left  # nombre de SNP dans la fenêtre
        if n_snps < MIN_SNPS:  # fenêtre avec trop peu de SNP : on l'ignore
            continue

        for comp in comparisons:  # calcule D/fd pour chaque comparaison de groupes dans cette fenêtre
            # récupère les groupes P1/P2/P3 de cette comparaison
            P1 = comp["P1"]
            P2 = comp["P2"]
            P3 = comp["P3"]

            p1 = group_arrays[P1][left:right]  # fréquences P1 dans la fenêtre
            p2 = group_arrays[P2][left:right]  # fréquences P2 (Awassi) dans la fenêtre
            p3 = group_arrays[P3][left:right]  # fréquences P3 dans la fenêtre

            D, fd, ABBA, BABA, num, den_D, den_fd = compute_D_fd(p1, p2, p3)  # calcule D et fd pour cette fenêtre/comparaison

            # enregistre les résultats de la fenêtre (position, comparaison, D, fd, nb SNP, etc.)
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "mid": (start + end) // 2,
                "comparison": comp["comparison"],
                "P1": P1,
                "P2": P2,
                "P3": P3,
                "O": comp["O"],
                "D": D,
                "fd": fd,
                "ABBA": ABBA,
                "BABA": BABA,
                "num": num,
                "den_D": den_D,
                "den_fd": den_fd,
                "n_snps": n_snps,
            })

    out = pd.DataFrame(rows)  # assemble tous les résultats en une table
    out_path = outdir / f"fd_chr{chrom}_all_awassi_pairs_20kb_step5kb.tsv.gz"  # chemin de sortie du fichier fd pour ce chromosome
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")  # écrit les résultats en TSV compressé

    print()
    print(f"chr{chrom} fenêtres/comparaisons écrites : {len(out)}")
    print(f"Résultat : {out_path}")

    return out_path  # renvoie le chemin du fichier écrit

def main():  # point d'entrée du script (gère les arguments et orchestre les étapes)
    ap = argparse.ArgumentParser()  # parseur d'arguments en ligne de commande
    ap.add_argument("--chrom", required=True, help="Chromosome, ex: 16")  # numéro de chromosome à traiter
    ap.add_argument("--outdir", required=True)  # dossier de sortie
    ap.add_argument("--groups", default="", help="Optionnel: groupes séparés par virgule")  # liste optionnelle de groupes à comparer
    args = ap.parse_args()  # parse les arguments fournis

    chrom = str(args.chrom)  # numéro de chromosome en texte
    outdir = Path(args.outdir)  # dossier de sortie en objet Path
    outdir.mkdir(parents=True, exist_ok=True)  # crée le dossier de sortie si besoin

    print("============================================================")
    print(f"fd genome-wide par chromosome — chr{chrom}")
    print("============================================================")
    print(f"WINDOW={WINDOW} STEP={STEP} MIN_SNPS={MIN_SNPS}")

    vcf = get_vcf(chrom)  # localise le VCF du chromosome
    print(f"VCF : {vcf}")

    meta = load_metadata()  # charge la metadata des échantillons

    vcf_samples = run(["bcftools", "query", "-l", str(vcf)]).splitlines()  # liste des échantillons présents dans le VCF
    groups, partners = choose_groups(meta, vcf_samples, args.groups)  # détermine les groupes et les partenaires à comparer
    comparisons = build_comparisons(partners)  # construit la liste des comparaisons P1/P2/P3

    print()
    print(f"Nombre de groupes partenaires : {len(partners)}")
    print(f"Nombre de comparaisons fd : {len(comparisons)}")

    selected = set()  # ensemble de tous les échantillons utilisés (tous groupes confondus)
    for s in groups.values():  # parcourt les échantillons de chaque groupe
        selected.update(s)

    selected_order = [s for s in vcf_samples if s in selected]  # ordre des échantillons sélectionnés, aligné sur l'ordre du VCF

    # Sauvegardes utiles
    # sauvegarde la composition des groupes utilisés (traçabilité)
    pd.DataFrame(
        [{"group": g, "n": len(s), "samples": ",".join(s)} for g, s in groups.items()]
    ).to_csv(outdir / f"chr{chrom}_groups_used.tsv", sep="\t", index=False)

    pd.DataFrame(comparisons).to_csv(outdir / f"chr{chrom}_comparisons_used.tsv", sep="\t", index=False)  # sauvegarde la liste des comparaisons utilisées

    freqs = extract_and_polarize(chrom, vcf, groups, selected_order, outdir)  # extrait les génotypes et calcule les fréquences dérivées polarisées
    scan_fd(chrom, freqs, comparisons, outdir)  # calcule D/fd par fenêtre glissante et écrit les résultats

if __name__ == "__main__":  # point d'entrée du script en exécution directe
    main()
