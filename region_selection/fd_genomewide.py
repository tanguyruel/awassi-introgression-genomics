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

def run(cmd):
    return subprocess.check_output(cmd, text=True)

def parse_bool(x):
    return str(x).strip().lower() in {"true", "1", "yes", "oui"}

def load_metadata():
    for p in METADATA_CANDIDATES:
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

def is_excluded_group(g):
    g = str(g)
    return any(re.search(pat, g, flags=re.IGNORECASE) for pat in EXCLUDE_GROUP_PATTERNS)

def get_vcf(chrom):
    candidates = [
        PROJECT / f"data/raw data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
        PROJECT / f"data/raw_data_08_06/awassi_and_basedata_chr{chrom}.vcf.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"VCF introuvable pour chr{chrom} : {candidates}")

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
    vcf_set = set(vcf_samples)

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

    all_meta_groups = sorted([str(g) for g in meta["group"].dropna().unique()])

    if groups_arg and not re.fullmatch(r"[0-9, ]+", str(groups_arg).strip()):  # --groups fourni et pas juste des chiffres (piège variable bash GROUPS)
        candidate_groups = [g.strip() for g in str(groups_arg).split(",") if g.strip()]
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

        candidate_groups = []

        for g in preferred:  # d'abord les groupes préférés présents dans la metadata
            if g in all_meta_groups and g not in candidate_groups:
                candidate_groups.append(g)

        for g in all_meta_groups:  # ensuite les groupes restants (hors Awassi/outgroup/motifs sauvages)
            if g not in candidate_groups and g not in {"Awassi", "Ovis_canadensis"}:
                if not re.search(r"OCAN|canad|bighorn|vignei|orientalis|mouflon|wild", g, flags=re.IGNORECASE):
                    candidate_groups.append(g)

    print()
    print("Groupes candidats testés comme partenaires :")
    for g in candidate_groups:
        print(" -", g)

    for g in candidate_groups:
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

    if len(groups["Awassi"]) == 0:
        raise RuntimeError("Awassi vide.")

    if len(groups["Ovis_canadensis"]) == 0:
        print("Noms wild/outgroups possibles dans le VCF :")
        for s in vcf_samples:  # aide au diagnostic
            if re.search(r"OCAN|canad|bighorn|ovis|vignei|orientalis|mouflon|wild", s, flags=re.IGNORECASE):
                print(" -", s)
        raise RuntimeError("Ovis_canadensis vide.")

    partners = [g for g in groups if g not in {"Awassi", "Ovis_canadensis"}]

    if len(partners) < 2:
        raise RuntimeError(f"Pas assez de groupes partenaires : {partners}")

    return groups, partners

def build_comparisons(partners):
    """
    P2 = Awassi fixé.
    Pour chaque groupe partenaire P3, on teste tous les autres groupes comme P1.
    Cela donne tous les signaux Awassi × groupe, mais avec référence explicite.
    """
    comps = []
    for P3 in partners:
        for P1 in partners:
            if P1 == P3:
                continue
            comps.append({
                "comparison": f"{P1}__Awassi__{P3}__OCAN",
                "P1": P1,
                "P2": "Awassi",
                "P3": P3,
                "O": "Ovis_canadensis",
            })
    return comps

def gt_to_alt_count(gt):
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

def allele_freq_alt(counts, idx):
    vals = counts[idx]
    vals = vals[~np.isnan(vals)]

    n_total = len(idx)
    n_called = len(vals)

    if n_total == 0:
        return np.nan, n_called

    if n_called / n_total < MIN_CALLED_FRAC_GROUP:
        return np.nan, n_called

    return vals.sum() / (2 * n_called), n_called  # fréquence alt = somme des comptes / (2 x nb génotypés)

def extract_and_polarize(chrom, vcf, groups, selected_order, outdir):  # polarisation par l'outgroup
    sample_to_idx = {s: i for i, s in enumerate(selected_order)}

    group_idx = {}
    for g, samples in groups.items():
        group_idx[g] = np.array([sample_to_idx[s] for s in samples if s in sample_to_idx], dtype=int)

    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:  # fichier temporaire listant les échantillons à extraire
        for s in selected_order:
            tmp.write(s + "\n")
        sample_file = tmp.name

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

    rows = []
    n_total = 0
    n_outgroup_ambiguous = 0
    n_missing_group = 0

    domestic_groups = [g for g in groups if g not in {"Ovis_canadensis"}]

    for line in proc.stdout:
        n_total += 1
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue

        pos = int(parts[1])
        gts = parts[2:]

        alt_counts = np.array([gt_to_alt_count(gt) for gt in gts], dtype=float)

        pO_alt, _ = allele_freq_alt(alt_counts, group_idx["Ovis_canadensis"])
        if not np.isfinite(pO_alt):
            n_missing_group += 1
            continue

        if pO_alt <= OUTGROUP_FIXED_THRESHOLD:  # outgroup fixé pour l'allèle référence : pas d'inversion
            flip = False
        elif pO_alt >= 1 - OUTGROUP_FIXED_THRESHOLD:  # outgroup fixé pour l'allèle alternatif : on inverse la polarité
            flip = True
        else:
            n_outgroup_ambiguous += 1  # outgroup ni fixé ref ni alt : site ambigu, ignoré
            continue

        rec = {"pos": pos}
        ok = True

        for g in domestic_groups:
            p_alt, _ = allele_freq_alt(alt_counts, group_idx[g])
            if not np.isfinite(p_alt):
                ok = False
                break

            rec[g] = 1 - p_alt if flip else p_alt  # applique l'inversion de polarité si besoin (fréquence de l'allèle dérivé)

        if not ok:
            n_missing_group += 1
            continue

        rows.append(rec)

    stderr = proc.stderr.read()
    ret = proc.wait()

    os.remove(sample_file)

    if ret != 0:
        print(stderr)
        raise RuntimeError(f"bcftools a échoué pour chr{chrom}")

    freqs = pd.DataFrame(rows)

    print()
    print(f"chr{chrom} SNPs PASS bialléliques lus : {n_total}")
    print(f"chr{chrom} SNPs polarisés utilisables : {len(freqs)}")
    print(f"chr{chrom} SNPs ignorés outgroup ambigu : {n_outgroup_ambiguous}")
    print(f"chr{chrom} SNPs ignorés données groupe insuffisantes : {n_missing_group}")

    freqs_path = outdir / f"chr{chrom}_derived_freqs.tsv.gz"
    freqs.to_csv(freqs_path, sep="\t", index=False, compression="gzip")

    return freqs

def compute_D_fd(p1, p2, p3):  # statistiques D (ABBA-BABA) et fd pour un jeu de fréquences P1/P2/P3
    ABBA = (1 - p1) * p2 * p3  # poids du motif ABBA (site dérivé chez P2 et P3, pas chez P1)
    BABA = p1 * (1 - p2) * p3  # poids du motif BABA (site dérivé chez P1 et P3, pas chez P2)

    num = np.nansum(ABBA - BABA)  # numérateur commun à D et fd (somme des différences ABBA-BABA)
    den_D = np.nansum(ABBA + BABA)  # dénominateur de la statistique D

    D = num / den_D if den_D > 0 else np.nan  # statistique D (proportion de sites ABBA vs BABA)

    # Martin et al. fd : denominator with dynamic donor PD.
    PD = np.maximum(p2, p3)  # fréquence du "donneur" dynamique (Martin et al. 2015)
    den_fd = np.nansum(PD * (PD - p1))  # dénominateur de fd, calculé avec PD au lieu de p3 fixe

    fd = num / den_fd if den_fd > 0 else np.nan  # statistique fd (estimateur robuste de la fraction introgressée)

    return D, fd, float(np.nansum(ABBA)), float(np.nansum(BABA)), float(num), float(den_D), float(den_fd)

def scan_fd(chrom, freqs, comparisons, outdir):  # parcourt le chromosome en fenêtres glissantes et calcule D/fd pour chaque comparaison
    if freqs.empty:
        raise RuntimeError(f"Pas de SNPs polarisés pour chr{chrom}")

    freqs = freqs.sort_values("pos").reset_index(drop=True)
    positions = freqs["pos"].to_numpy(dtype=int)

    group_arrays = {
        c: freqs[c].to_numpy(dtype=float)
        for c in freqs.columns
        if c != "pos"
    }

    max_pos = int(positions.max())
    rows = []

    for start in range(1, max_pos - WINDOW + 2, STEP):  # fenêtres glissantes ; +2 pour inclure la dernière fenêtre entière
        end = start + WINDOW - 1
        left = np.searchsorted(positions, start, side="left")  # index du premier SNP dans la fenêtre
        right = np.searchsorted(positions, end, side="right")  # index après le dernier SNP dans la fenêtre

        n_snps = right - left
        if n_snps < MIN_SNPS:
            continue

        for comp in comparisons:
            P1 = comp["P1"]
            P2 = comp["P2"]
            P3 = comp["P3"]

            p1 = group_arrays[P1][left:right]
            p2 = group_arrays[P2][left:right]
            p3 = group_arrays[P3][left:right]

            D, fd, ABBA, BABA, num, den_D, den_fd = compute_D_fd(p1, p2, p3)

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

    out = pd.DataFrame(rows)
    out_path = outdir / f"fd_chr{chrom}_all_awassi_pairs_20kb_step5kb.tsv.gz"
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")

    print()
    print(f"chr{chrom} fenêtres/comparaisons écrites : {len(out)}")
    print(f"Résultat : {out_path}")

    return out_path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, help="Chromosome, ex: 16")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--groups", default="", help="Optionnel: groupes séparés par virgule")
    args = ap.parse_args()

    chrom = str(args.chrom)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("============================================================")
    print(f"fd genome-wide par chromosome — chr{chrom}")
    print("============================================================")
    print(f"WINDOW={WINDOW} STEP={STEP} MIN_SNPS={MIN_SNPS}")

    vcf = get_vcf(chrom)
    print(f"VCF : {vcf}")

    meta = load_metadata()

    vcf_samples = run(["bcftools", "query", "-l", str(vcf)]).splitlines()
    groups, partners = choose_groups(meta, vcf_samples, args.groups)
    comparisons = build_comparisons(partners)

    print()
    print(f"Nombre de groupes partenaires : {len(partners)}")
    print(f"Nombre de comparaisons fd : {len(comparisons)}")

    selected = set()
    for s in groups.values():
        selected.update(s)

    selected_order = [s for s in vcf_samples if s in selected]  # ordre aligné sur le VCF (important pour les index de groupe)

    pd.DataFrame(
        [{"group": g, "n": len(s), "samples": ",".join(s)} for g, s in groups.items()]
    ).to_csv(outdir / f"chr{chrom}_groups_used.tsv", sep="\t", index=False)  # traçabilité

    pd.DataFrame(comparisons).to_csv(outdir / f"chr{chrom}_comparisons_used.tsv", sep="\t", index=False)

    freqs = extract_and_polarize(chrom, vcf, groups, selected_order, outdir)
    scan_fd(chrom, freqs, comparisons, outdir)

if __name__ == "__main__":
    main()
