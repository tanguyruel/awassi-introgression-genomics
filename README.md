# Awassi Introgression Genomics

Scripts d'analyse accompagnant le rapport de stage *« Introgression génomique chez le
mouton Awassi — signaux locaux d'introgression et validation multi-critères »* (Tanguy Ruel,
LECA). Recherche de régions du génome où la race Awassi (Moyen-Orient) présente une
affinité locale anormale avec un autre groupe géographique (Afrique, Asie, Europe, Amérique,
Océanie), signal compatible avec une introgression ancienne.

## Contenu de ce dépôt

20 scripts : un par méthode utilisée dans le rapport, exactement ceux cités dans son
**Annexe A — Correspondance méthode / script / référence**. Le dépôt de travail complet
(~340 scripts, itérations de mise au point incluses) reste local ; seule la version de
référence de chaque méthode est publiée ici.

Chaque script commence par un commentaire d'en-tête (rôle, entrée, sortie, usage) et
comporte des commentaires dans le corps du code pour suivre chaque étape du calcul.

### Structure de population

| Méthode | Script |
|---|---|
| ACP (PLINK --pca 20) | `01_pca/01_run_pca_25chr_filtered.sh` |
| ADMIXTURE | `admixture/03_run_admixture_K2_K8_v1.sh` |

### Détection et sélection des régions candidates

| Méthode | Script |
|---|---|
| FST local 20 kb / pas 5 kb | `fst/21_run_FST_local_20kb_step5kb_separate_EU_AM_AUS.sh` |
| FST par SNP | `synthese/13b_plot_fst_fd_ld_aligned_v9_perSNPfst.py` |
| fd | `fd/02_fd_chr_all_awassi_pairs_20kb_step5kb.py` |
| Sélection — score combiné | `synthese/15_rebuild_v10b.py` |
| Sélection — intersection (Manhattan) | `synthese/23_manhattan_fst_fd_genomewide.py` |

### Tests de validation

| Méthode | Script |
|---|---|
| D-stat 150 kb + jackknife | `synthese/27_dstat_9regions_blockjackknife_v1.py` |
| dXY | `dxy/02_compute_dxy_local_candidates.py` |
| π avec contrôle qualité | `synthese/30_pi_9regions_20kb_step5kb_QC_v1.py` |
| π — fond génomique | `synthese/50_pi_genomewide_baseline_v1.py` |
| Phasage (Beagle) | `synthese/10_phase_5regions_beagle_v1.sh` |
| LD par sous-groupe | `synthese/35_ld_per_subgroup_9regions_v1.py` |
| Heatmap + arbre + LD | `synthese/26_heatmap_LD_arbre_9regions_v2.py` |
| Raréfaction (partage d'haplotypes) | `synthese/34_rarefaction_haplotype_sharing_9regions_v1.py` |
| Bootstrap de spécificité | `synthese/65_haplotype_specificity_bootstrap_v1.py` (dépend de `64_haplotype_sharing_all_groups_v1.py`, inclus) |
| F et Ho | `synthese/56_het_inbreeding_v1.py` |
| CCDC91 — reprise sur le cœur du pic | `synthese/69_ccdc91_peak_core_retest_v1.py` (dépend aussi de `64_...py`) |

### Annotation fonctionnelle

| Méthode | Script |
|---|---|
| Annotation fonctionnelle (GFF Oar_v4.0) | `synthese/40_annotate_variants_9regions_v1.py` |

## Logiciels externes utilisés (non fournis dans ce dépôt)

- [PLINK](https://www.cog-genomics.org/plink/) / [PLINK2](https://www.cog-genomics.org/plink/2.0/)
- [bcftools / vcftools](https://samtools.github.io/bcftools/)
- [ADMIXTURE](https://dalexander.github.io/admixture/) (Alexander et al. 2009)
- [Beagle 5.5](https://faculty.washington.edu/browning/beagle/beagle.html) (Browning et al. 2018)
- Python ≥ 3.8 (pandas, numpy, matplotlib, scipy) ; R (ggplot2)

## À savoir avant de relancer un script

- Tous les scripts supposent d'être lancés depuis la racine de ce dépôt.
- Plusieurs scripts ont un chemin absolu en dur vers la machine d'origine (`/home/tanguyruel/...`),
  repéré par un commentaire `# chemin en dur à adapter` — à corriger avant toute exécution.
- `65_...py` et `69_...py` importent `64_haplotype_sharing_all_groups_v1.py` par nom de module :
  les trois doivent rester dans le même dossier `synthese/`.
- `35_ld_per_subgroup_9regions_v1.py` importe directement `26_heatmap_LD_arbre_9regions_v2.py`
  (même dossier requis).

## Données

Les VCF bruts et les résultats intermédiaires volumineux ne sont pas inclus (taille, et
données génétiques dont la redistribution n'est pas garantie). Les tables de résultats
filtrées citées dans le rapport (régions candidates, métriques par test) sont disponibles
sur demande auprès de l'auteur.

## Citation

Code disponible pour accompagner : Ruel T., *Introgression génomique chez le mouton
Awassi*, rapport de stage LECA, 2026. *(mettre à jour avec la référence définitive de
l'article une fois publié.)*

## Licence

Code sous licence MIT (voir `LICENSE`). Les données génomiques associées ne sont pas
couvertes par cette licence.
