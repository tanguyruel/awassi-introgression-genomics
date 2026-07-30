# Jeu de test synthétique

Petit jeu de données **entièrement synthétique** (aucune donnée biologique réelle) contenant
un signal d'introgression connu, pour pouvoir essayer les scripts du dépôt sans les données
réelles du stage (non publiées, cf. [Données](../README.md#données) dans le README principal).

## Contenu

- `generate_test_dataset.py` : script générateur (déterministe, graine fixée) — relire ce
  fichier pour le détail du modèle de simulation.
- `data/raw_data_08_06/awassi_and_basedata_chrTEST1.vcf.gz(.tbi)` : VCF simulé, 1 chromosome
  fictif `TEST1` (300 kb), 1800 SNP, 54 échantillons.
- `analyses/haplotype_heatmap/Awassi_haplo/data/metadata/sample_metadata_387_FST_groups.tsv` :
  metadata (échantillon, groupe, is_awassi).
- `regions_test.tsv` : deux régions pour `dxy.py` — la région introgressée et une région témoin.

L'arborescence `data/…` et `analyses/…` reproduit volontairement les chemins que les scripts
cherchent sous `--project` (voir `load_metadata`/`get_vcf` dans `fd_genomewide.py` et
`dxy.py`), pour que `--project example_data` fonctionne directement sans configuration.

## Le signal simulé

5 groupes domestiques (Awassi, Asia, Africa, Europe, MiddleEastNonAwassi, 10 individus
chacun) + un outgroup sauvage Ovis_canadensis (4 individus, fixé par construction pour
polariser ancestral/dérivé). Les fréquences alléliques par groupe sont tirées par site selon
un modèle de dérive (Balding-Nichols, Fst ≈ 0.06), **sauf** dans la fenêtre
`TEST1:140 000-175 000`, où la fréquence du groupe **Awassi est mélangée à ~75 % avec celle
du groupe Asia** : cela simule un bloc introgressé où Awassi a reçu du flux génique d'Asia.

Attention : chaque SNP est simulé indépendamment (pas de déséquilibre de liaison / structure
haplotypique). Ce jeu convient donc aux méthodes basées sur des fréquences alléliques par
fenêtre (`fd_genomewide.py`, `dxy.py`, et par extension FST/π), mais **pas** aux méthodes qui
ont besoin d'haplotypes réalistes (LD, phasage, partage d'haplotypes) — celles-ci tourneraient
sans erreur mais sans signal biologiquement cohérent.

## Reproduire le signal

```bash
# depuis la racine du dépôt
python3 region_selection/fd_genomewide.py --chrom TEST1 --outdir /tmp/fd_test --project example_data
python3 validation_tests/dxy.py --regions example_data/regions_test.tsv --outdir /tmp/dxy_test --project example_data
```

Résultats obtenus en générant ce jeu (graine fixée, donc reproductibles à l'identique) :

- **fd genome-wide** (`fd_genomewide.py`) : dans la fenêtre introgressée, la comparaison
  `P3=Asia` est la mieux classée (fd le plus élevé) dans 6 des 8 fenêtres 20 kb qui la
  recouvrent (75 %), contre 31 % en dehors — un signal net mais bruité, comme attendu sur un
  aussi petit jeu de SNPs (c'est justement pourquoi le pipeline principal ne conclut jamais
  sur fd seul, cf. validation multi-critères).
- **dXY** (`dxy.py`) : Awassi-Asia a le dXY le plus bas des groupes partenaires dans la région
  test (0.230, contre 0.236-0.246 pour les autres groupes) — écart d'environ 0.006 avec le
  deuxième plus proche. Dans la région témoin (`control_region`, hors fenêtre introgressée),
  Asia n'est plus significativement la plus proche (écart de 0.0006 avec le deuxième groupe,
  dans le bruit) : la baisse de dXY est donc spécifique à la fenêtre introgressée, pas une
  proximité générale entre Awassi et Asia.

## Régénérer / modifier le jeu

```bash
python3 example_data/generate_test_dataset.py
```

Modifier les constantes en tête de `generate_test_dataset.py` (position/taille de la fenêtre
introgressée, groupes, effectifs, Fst, fraction d'introgression) pour explorer d'autres
scénarios.
