Article1
================
Marion MENEC
2023-11-29

# Partie DADA2

## 1.Préparation du répertoire de travail

``` r
refdb_folder <- here::here("article", "refdb")
refdb_folder
```

    ## [1] "/home/rstudio/DADA2/article/refdb"

``` r
if (!dir.exists(refdb_folder)) dir.create(refdb_folder,recursive = TRUE)
```

Lors téléchargement de données, R s’arrête au bout de 60 secondes.
Changement de ce temps à 20 min.

``` r
getOption("timeout")
options(timeout = 1200)
```

Téléchargement des outils.

``` r
devtools::load_all(path = "/home/rstudio/DADA2/course-material-main/R/")
```

``` r
library(dada2); packageVersion("dada2")
```

    ## Loading required package: Rcpp

    ## [1] '1.28.0'

``` r
# Lire les liens FTP à partir du fichier sequences.txt
ftp_links <- readLines("sequences.txt")

# Dossier de sortie
sequences <- "sequences"

# Créer le dossier de sortie s'il n'existe pas
if (!file.exists(sequences)) {
  dir.create(sequences, recursive = TRUE)
}

# Fonction pour télécharger une séquence (forward ou reverse)
download_sra_sequence <- function(ftp_link, output_directory) {
  # Extraire l'identifiant SRA à partir du lien FTP
  sra_id <- sub(".*/(SRR\\d+)/(SRR\\d+_\\d)\\.fastq\\.gz", "\\1_\\2", ftp_link)
  
  # Télécharger le fichier FASTQ
  dest_file <- file.path(output_directory, paste0(sra_id, ".fastq.gz"))
  download.file(ftp_link, dest_file, method = "auto")
}

# Télécharger les séquences pour chaque lien FTP
for (ftp_link in ftp_links) {
  download_sra_sequence(ftp_link, output_directory)
}

# Afficher le contenu du dossier de sortie
dir(sequences)
```

``` r
path_to_fastqs <- here::here("article", "sequences")
```

Sauvegarder le chemin de refdb dans l’espace de travail DANS un objet.

``` r
silva_train_set <- file.path(refdb_folder, 
                             "silva_nr99_v138.1_train_set.fa.gz")
silva_species_assignment <- file.path(refdb_folder, 
                                      "silva_species_assignment_v138.1.fa.gz")
```

Téléchargement des fichiers s’ils n’existent pas.

``` r
if (!file.exists(silva_train_set)){
  download.file(
    "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
    silva_train_set,
    quiet = TRUE
  )
}

if (!file.exists(silva_species_assignment)) {
  download.file(
    "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz",
    silva_species_assignment,
    quiet = TRUE
  )
}
```

## 2.Fichiers d’entrées

Faire une liste avec les fichiers forward.

``` r
fnFs <- sort(list.files(path_to_fastqs,
                        pattern = "_R1.fastq.gz",
                        full.names=TRUE))
```

Faire la même chose avec les fichiers reverse.

``` r
fnRs <- sort(list.files(path_to_fastqs,
                        pattern = "_R2.fastq.gz",
                        full.names = TRUE))
```

``` r
basename(fnFs)
```

    ## [1] "SRR12576537_SRR12576537_R1.fastq.gz" "SRR12576538_SRR12576538_R1.fastq.gz"
    ## [3] "SRR12576539_SRR12576539_R1.fastq.gz" "SRR12576540_SRR12576540_R1.fastq.gz"
    ## [5] "SRR12576541_SRR12576541_R1.fastq.gz" "SRR12576542_SRR12576542_R1.fastq.gz"

basename() : permet de garder le nom du fichier ;

strsplit() : sépare un chaîne de caractère au niveau d’un caractère
donné ; \|\> : permet de mettre plusieurs fonctions à la suite sans
avoir à utiliser de multiples objets ;

sapply() : ici extrait le premier élément de la liste.

``` r
sample_names <- basename(fnFs) |>
  strsplit(split = "_") |>
  sapply(head, 1)
head(sample_names)
```

    ## [1] "SRR12576537" "SRR12576538" "SRR12576539" "SRR12576540" "SRR12576541"
    ## [6] "SRR12576542"

## 3.Vérification de la qualité des séquences

Création d’un fichier où stocker les données de sortie du
qualityprofile. Vérification de la qualité des profils.

``` r
quality_folder <- here::here("outputs",
                             "article",
                             "quality_plots")

if (!dir.create(quality_folder)) {
  dir.create(quality_folder, recursive = TRUE)
}

qualityprofile(fnFs,
               fnRs,
               file.path(quality_folder, "quality_plots.pdf"))
```

## 4.Enlèvement des séquences des amorces

Création d’un dossier dans lequel sera enregistré les reads une fois les
séquences des amorces retirées.

``` r
path_to_trimmed_reads <- here::here(
  "outputs",
  "article",
  "trimmed"
)

if (!dir.exists(path_to_trimmed_reads)) dir.create(path_to_trimmed_reads, recursive = TRUE)
```

Amorces forward (TCGTC-GGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG) et
reverse (GTCTCGTGGG-CTCGGAGATGTGTATAAGAGACAGGACTACHV-GGGTATCTAATCC).

``` r
primer_fwd  <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
primer_rev  <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
# J'ai du enlever les - dans la séquence des amorces sinon la fonction "primer_log" ne fonctionnait pas.
```

Trouver la séquence de l’amorce forward dans les reads forward.

``` r
Biostrings::readDNAStringSet(
  fnFs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]    33 CCTACGGGTGGCAGCAGTGGGGAATATTGCACA                SRR12576537.1 1/1
    ##  [2]    39 CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGGC          SRR12576537.2 2/1
    ##  [3]    40 CCTACGGGAGGCTGCAGTCGAGAATCTTCGTCAATGGACG         SRR12576537.3 3/1
    ##  [4]    73 CCTACGGGGGGCTGCAGTGAGGA...CAGCCATGTCGCGTGCGTGATG SRR12576537.4 4/1
    ##  [5]    33 CCTACGGGTGGCTGCAGTGGGGAATTTTGGACA                SRR12576537.5 5/1
    ##  [6]    32 CCTACGGGGGGCTGCAGTGGGGAATTTTGGTC                 SRR12576537.6 6/1
    ##  [7]    71 CCTACGGGAGGCTGCAGTGAGGA...ACCAGCCATGCCGCGTGCGGGA SRR12576537.7 7/1
    ##  [8]   113 CCTACGGGTGGCAGCAGTGAGGA...AACTGCTTTTATAAGGGAATAA SRR12576537.8 8/1
    ##  [9]    97 CCTACGGGAGGCTGCAGTCGAGA...GGTCTTCTGATTGTTAACCTCT SRR12576537.9 9/1
    ## [10]    32 CCTACGGGGGGCAGCAGTGGGGAATATTGCAC                 SRR12576537.10 10/1

Trouver la séquence de l’amorce reverse dans les reads reverse.

``` r
Biostrings::readDNAStringSet(
  fnRs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]     4 GACT                                             SRR12576537.1 1/2
    ##  [2]     4 GACT                                             SRR12576537.2 2/2
    ##  [3]    43 GACTACTCGGGTATCTAATCCCGTTCGCTCCCCTGGCTTTCGC      SRR12576537.3 3/2
    ##  [4]     4 GACT                                             SRR12576537.4 4/2
    ##  [5]    43 GACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGT      SRR12576537.5 5/2
    ##  [6]    45 GACTACCCGGGTATCTAATCCCGTTTGCTCCCCTCGCTTTCGCAC    SRR12576537.6 6/2
    ##  [7]     9 GACTACAGG                                        SRR12576537.7 7/2
    ##  [8]    27 GACTACCGGGGTATCAGCTAATCTTTT                      SRR12576537.8 8/2
    ##  [9]    67 GACTACAAGGGTATCTAATCCTG...CTGAGTGTCAGTTATCGTCCAG SRR12576537.9 9/2
    ## [10]    42 GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCG       SRR12576537.10 10/2

On enlève les séquences des amorces des reads.

``` r
(primer_log <- primer_trim(
  forward_files = fnFs,
  reverse_files = fnRs,
  primer_fwd = primer_fwd,
  primer_rev = primer_rev,
  output_dir = path_to_trimmed_reads,
  min_size = 200
))
# J'ai du enlever les - dans la séquence des amorces sinon le chunck ne fonctionnait pas
```

## 5.Découpage et filtrage de qualité

Création d’un dossier et des listes de chemins d’accès pour les reads
qui seront filtrés.

``` r
path_to_filtered_reads <- here::here("outputs", "article","filtered")
if (!dir.exists(path_to_filtered_reads)) dir.create(path_to_filtered_reads, recursive = TRUE)

filtFs <- file.path(path_to_filtered_reads, basename(fnFs))
filtRs <- file.path(path_to_filtered_reads, basename(fnRs))
```

Pour faire le lien entre les fichiers et les noms d’échantillons, nommer
simplement le vecteur de noms de fichiers en utilisant les noms
d’échantillons.

``` r
names(filtFs) <- sample_names
names(filtRs) <- sample_names
```

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                                     reads.in reads.out
    ## SRR12576537_SRR12576537_R1.fastq.gz   161823     93851
    ## SRR12576538_SRR12576538_R1.fastq.gz   150276     87167
    ## SRR12576539_SRR12576539_R1.fastq.gz   109001     62603
    ## SRR12576540_SRR12576540_R1.fastq.gz   110147     66861
    ## SRR12576541_SRR12576541_R1.fastq.gz   119114     70356
    ## SRR12576542_SRR12576542_R1.fastq.gz   113636     66636

## 6.Débruitage

Enlever les erreurs de séquençage : regarder la séquence majoritaire et
corriger.

Création d’un modèle d’erreur.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
# temps d'exécution extrêmement long (plus d'une heure. Obliger de passer le chunck en eval=FALSE pour knit, mais le code fonctionne
```

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
# Même chose que le chunck précédent.
```

Visualisation des modèles d’erreur.

``` r
dada2::plotErrors(errF, nominalQ = TRUE)
```

``` r
dada2::plotErrors(errR, nominalQ = TRUE)
```

Les modèles d’erreur sont sur le github dans le fichier
“Article1_files”.

Déreplication :

Ne garder qu’un seul exemplaire des séquences en multiples exemplaires
et mettre dans un fichier combien de fois on a trouvé cette séquence.

``` r
derepFs <- dada2::derepFastq(filtFs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576537_SRR12576537_R1.fastq.gz

    ## Encountered 33557 unique sequences from 93851 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576538_SRR12576538_R1.fastq.gz

    ## Encountered 29598 unique sequences from 87167 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576539_SRR12576539_R1.fastq.gz

    ## Encountered 22362 unique sequences from 62603 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576540_SRR12576540_R1.fastq.gz

    ## Encountered 20537 unique sequences from 66861 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576541_SRR12576541_R1.fastq.gz

    ## Encountered 21462 unique sequences from 70356 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576542_SRR12576542_R1.fastq.gz

    ## Encountered 19607 unique sequences from 66636 total sequences read.

``` r
derepRs <- dada2::derepFastq(filtRs, verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576537_SRR12576537_R2.fastq.gz

    ## Encountered 25645 unique sequences from 93851 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576538_SRR12576538_R2.fastq.gz

    ## Encountered 24507 unique sequences from 87167 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576539_SRR12576539_R2.fastq.gz

    ## Encountered 18475 unique sequences from 62603 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576540_SRR12576540_R2.fastq.gz

    ## Encountered 14943 unique sequences from 66861 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576541_SRR12576541_R2.fastq.gz

    ## Encountered 14989 unique sequences from 70356 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/DADA2/outputs/article/filtered/SRR12576542_SRR12576542_R2.fastq.gz

    ## Encountered 15319 unique sequences from 66636 total sequences read.

Enlèvement des erreurs de séquençage. Obliger de mettre en eval=FALSE
mais fonctionne sinon.

``` r
dadaFs <- dada2::dada(derepFs, err = errF, multithread = TRUE)
```

``` r
dadaRs <- dada2::dada(derepRs, err = errR, multithread = TRUE)
```

## 7.Fusions des reads appariés

Une fois que les erreurs de séquençage enlevés, faire l’alignement des
forward et des reverse, il y a création de contigues. Il n’y a plus de
R1 et R2.

``` r
mergers <- dada2::mergePairs(
  dadaF = dadaFs,
  derepF = derepFs,
  dadaR = dadaRs,
  derepR = derepRs,
  maxMismatch = 0,
  verbose = TRUE
)
```

## 8.Construction de la matrice d’observation

``` r
seqtab <- dada2::makeSequenceTable(mergers)
```

## 9.Enlevement des chimères (artefact)

Lors de la réplication, l’ADN polymérase peut se détacher et arrêter la
polymérisation. A la fin du cycle PCR, il y a 1 brin mère et un petit
brin fille pour cette bactérie. Ce petit brin fille peut venir
s’amorcer/se fixer au brin mère d’une autre bactérie et l’ADN polymérase
va reprendre élongation. Il y aura alors création de chimère : brin
mixte d’une bactérie 1 et d’une bactérie 2.

``` r
seqtab_nochim <- dada2::removeBimeraDenovo(seqtab,
                                           method = "consensus",
                                           multithread = TRUE,
                                           verbose = TRUE)
```

## 10.Assignation taxonomique

Calcul par pourcentage de chance, en partant du haut de l’arbre et
descend peu à peu. Il s’arrête un embranchement quand pour celui d’après
le pourcentage est inférieur au pourcentage défini. L’assignation à une
espèce ne dépend que de la base de données.

``` r
taxonomy <- dada2::assignTaxonomy(
  seqs = seqtab_nochim,
  refFasta = silva_train_set,
  taxLevels = c("Kingdom", "Phylum", "Class",
                "Order", "Family", "Genus",
                "Species"),
  multithread = TRUE,
  minBoot = 60
)
```

Permet l’assignation taxonomique au rang d’espèce, ce que ne permet pas
les lignes de code précédentes.

``` r
taxonomy <- dada2::addSpecies(
  taxonomy,
  silva_species_assignment,
  allowMultiple = FALSE
)
```

## 11.Export des données

Cette section correspond au téléchargement de l’environnement R sur
notre ordinateur.

Les résultats peuvent être exportés comme un objet R, 1 pour le tableau
ASV et un autre pour la taxonomie.

``` r
export_folder <- here::here("outputs", "article", "asv_table")

if (!dir.exists(export_folder)) dir.create(export_folder, recursive = TRUE)

saveRDS(object = seqtab_nochim,
        file = file.path(export_folder, "seqtab_nochim.rds"))

saveRDS(object = taxonomy,
        file = file.path(export_folder, "taxonomy.rds"))
```

Le but des lignes suivantes est de transformer nos résultats en fichier
texte. Création d’une variable pour collecter les séquences ASV.

``` r
asv_seq <- colnames(seqtab_nochim)
```

Création d’identifiant (id) unique pour chaque ASV (plus court que la
séquence en elle-même).

``` r
ndigits <- nchar(length(asv_seq))
asv_id <- sprintf(paste0("ASV_%0", ndigits, "d"), seq_along(asv_seq))
```

Renomer les différentes variables avec les nouveaux id.

``` r
row.names(taxonomy) <- colnames(seqtab_nochim) <- names(asv_seq) <- asv_id
```

Convertir les ASV ids (nom des lignes)dans une nouvelle colonne nommée
“asv”.

``` r
taxonomy_export <- df_export(taxonomy, new_rn = "asv")

seqtab_nochim_export <- t(seqtab_nochim)
seqtab_nochim_export <- df_export(seqtab_nochim_export, new_rn = "asv")
```

Finalement on peut exporter la taxonomie, le tableau ASV et les
séquences (celles là en fichier fasta).

``` r
write.table(taxonomy_export,
            file = file.path(export_folder, "taxonomy.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

write.table(seqtab_nochim_export,
            file = file.path(export_folder, "asv_table.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

cat(paste0(">", names(asv_seq), "\n", asv_seq),
    sep = "\n",
    file = file.path(export_folder, "asv.fasta"))
```

Les statistiques de chaque étapes de prétraitement peuvent également
être exportées. Il faut d’abord assembler le tableau.

``` r
getN <- function(x) sum(dada2::getUniques(x))

log_table <- data.frame(
  input = primer_log$in_reads,
  with_fwd_primer = primer_log$`w/adapters`,
  with_rev_primer = primer_log$`w/adapters2` ,
  with_both_primers = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab_nochim),
  perc_retained = rowSums(seqtab_nochim) / out[, 1] * 100
)

rownames(log_table) <- sample_names


# Ne fonctionne pas puisque que l'objet primer_log ne peut pas être générés car la fonction ne marche pas avec les noms de fichiers.
```

Puis l’exporter.

``` r
df_export(log_table, new_rn = "sample") |>
  write.table(file = file.path(export_folder, "log_table.tsv"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
```

Jusqu’ici le code fonctionne. Je n’ai pas eu le temps de poursuivre plus
à cause du temps de run des chuncks de création des modèles d’erreur.
Tous la fin du code (qui est au dessus de ce texte) à du être mis en
eval=FALSE, sinon je ne pouvais pas knit mais il fonctionne.

Ci-après, les prémices de code de bétâ-diversité théorique qui aurait
permis de continuer. Etant impossible de créer la table asv, il est donc
impossible pour moi de continuer plus loin.

# Partie bétâ-diversité

## Partie 1

Charger les librairies utiles.

``` r
library(phyloseq)
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

Charger les données.

``` r
physeq <- readRDS(here::here("outputs",
                             "article",
                             "asv_table",
                             ""))
```
