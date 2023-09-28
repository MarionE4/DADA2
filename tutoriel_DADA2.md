R Notebook
================

``` bash
wget https://github.com/ANF-MetaBioDiv/course-material/archive/refs/heads/main.zip
unzip main.zip
```

(re)Défini la racine de notre projet sur notre ordinateur.

``` r
refdb_folder <- here::here("data", "refdb")
refdb_folder
```

Création du vrai dossier refdb_folder (si rfdb_folder n’existe pas créer
le, et possibilité créer sous dossier).

``` r
if (!dir.exists(refdb_folder)) dir.create(refdb_folder,recursive = TRUE)
```

``` bash
cp -R course-material-main/data/raw ./data
```

Lors téléchargement de données R s’arrête au bout de 60 secondes,
changement de ce temps à 20 min.

``` r
getOption("timeout")
options(timeout = 1200)
```

Sauver le chemin de refdb dans l’espace de travail DANS un objet.

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

``` r
devtools::load_all(path = "/home/rstudio/DADA2/course-material-main/R/")
```

    ## ℹ Loading ANF_metaB

Sauvegarde du chemin d’accès au répertoire contenant les données brutes
dans un objet.

``` r
path_to_fastqs <- here::here("data", "raw")
```

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

basename() : permet de garder le nom du fichier ;

strsplit() : sépare un chaîne de caractère au niveau d’un caractère
donné ; \|\> : permet de mettre plusieurs fonctions à la suite sans
avoir à utiliser de multiples objets ;

sapply() : ici extrait le premier élément de la liste.

``` r
sample_names <- basename(fnFs) |>
  strsplit(split = "_") |>
  sapply(head, 1) |>
  head()
```

1.  Création d’un fichier où stocker les données de sortie du
    qualityprofile
2.  Vérification de la qualité des profils

``` r
quality_folder <- here::here("outputs",
                             "dada2",
                             "quality_plots")

if (!dir.create(quality_folder)) {
  dir.create(quality_folder, recursive = TRUE)
}

qualityprofile(fnFs,
               fnRs,
               file.path(quality_folder, "quality_plots.pdf"))
```

Création d’un dossier dans lequel sera enregistré les reads une fois les
nucléotides avec une qualité médiocre coupée.

``` r
path_to_trimmed_reads <- here::here(
  "outputs",
  "dada2",
  "trimmed"
)

if (!dir.exists(path_to_trimmed_reads)) dir.create(path_to_trimmed_reads, recursive = TRUE)
```

Les données correspondent à la région V3-V4 (460pb), en utilisant les
amorces Pro341F (CCTACGGGNBGCASCASCAG) et Pro805R
(GACTACNVGGGTATTCTAAT).

``` r
primer_fwd  <- "CCTACGGGNBGCASCAG"
primer_rev  <- "GACTACNVGGGTATCTAAT"
```

Trouver les amorces des forward dans les reads forward.

``` r
Biostrings::readDNAStringSet(
  fnFs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]   293 CCTACGGGGGGCAGCAGTAGGGA...ACATCGGCTTAACCGATGAAGT M01522:260:000000...
    ##  [2]   293 CCTACGGGTGGCACCAGTAGGGA...CGGGGCTTAACCTCGGAACTGC M01522:260:000000...
    ##  [3]   292 CCTACGGGGCGCAGCAGGCGCGA...GGGACCGGGAGAGGTGTGAGGT M01522:260:000000...
    ##  [4]   293 CCTACGGGGTGCAGCAGTAGGGA...TCAAAACTCCCAGTCTAGAGTT M01522:260:000000...
    ##  [5]   291 CCTACGGGTGGCAGCAGTGGGGA...GCAGTGGAAACTGTTGGGCTTG M01522:260:000000...
    ##  [6]   293 CCTACGGGATGCAGCAGGCGCGA...GGGACCGGGAGAGGTGTGGGGG M01522:260:000000...
    ##  [7]   292 CCTACGGGATGCAGCAGTGGGGA...TTTAATCCTGATGAGCTAGAAA M01522:260:000000...
    ##  [8]   293 CCTACGGGGCGCAGCAGTAGGGA...TTAAAACTTTTGTTCTGGAATT M01522:260:000000...
    ##  [9]   292 CCTACGGGTTGCAGCAGTGGGGA...ATTAAAACTTTTCAGCTAGAGT M01522:260:000000...
    ## [10]   293 CCTACGGGAGGCAGCAGTGGGGA...CCCGGGCTCAACCTGGGAACGG M01522:260:000000...

Trouver les amorces reverse dans les reads reverse.

``` r
Biostrings::readDNAStringSet(
  fnRs[1],
  format = "fastq",
  nrec = 10
)
```

    ## DNAStringSet object of length 10:
    ##      width seq                                              names               
    ##  [1]   301 GACTACCAGGGTATCTAATCCTG...GGCTGCTGGCACGAAGTTCGCC M01522:260:000000...
    ##  [2]   301 GACTACCGGGGTATCTAATCCTG...GGCTGCTGGCACGGAGTTAGCC M01522:260:000000...
    ##  [3]   300 AATCCGGTTCGTGCCCCTAGGCT...TCTTTCCCAGCCCTTATTCCAA M01522:260:000000...
    ##  [4]   301 GACTACCGGGGTATCTAATCCTG...GGCTGCTGGCACGGAGTTAGCC M01522:260:000000...
    ##  [5]   301 GACTACCGGGGTATCTAATCCCT...GGCTGCTGGCCCGGAATTAGCC M01522:260:000000...
    ##  [6]   301 GGTATCTAATCCGGTTCGTGCCC...CACCGTCCTTACCCCCCCCTTT M01522:260:000000...
    ##  [7]   301 GGTATCTAATCTTGTTTGCTCCC...CCCGACGTTAGCCGGGGCTTCT M01522:260:000000...
    ##  [8]   301 GACTACGAGGGTATCTAATCCCG...GGCTGCTGGCACGGAATTAGCC M01522:260:000000...
    ##  [9]   301 GGTATCTAATCCTCTTCGCTACC...CACGAAGTTAGCCGGACCTTCT M01522:260:000000...
    ## [10]   301 GACTACGGGGGTATCTAATCCTG...GGCTGCCGGCACGGGGTTAGCC M01522:260:000000...

Copie du dossier bash (comme R avant) dans le dossier DADA2

``` bash
pwd
cp -R /home/rstudio/DADA2/course-material-main/bash .
```

    ## /home/rstudio/DADA2

On enlève les séquences des amorces des reads

``` r
(primer_log <- primer_trim(
  forward_files = fnFs,
  reverse_files = fnRs,
  primer_fwd = primer_fwd,
  primer_rev = primer_rev,
  output_dir = path_to_trimmed_reads,
  min_size = 200
))
```

    ##    sample status in_reads   in_bp too_short too_long too_many_n out_reads
    ## 1    S11B     OK     2000 1186767         0        0          0      1863
    ## 2     S1B     OK     2000 1186613         1        0          0      1855
    ## 3     S2B     OK     2000 1186942         0        0          0      1839
    ## 4     S2S     OK     2000 1186868         0        0          0      1833
    ## 5     S3B     OK     2000 1186650         0        0          0      1860
    ## 6     S3S     OK     2000 1186475         1        0          0      1880
    ## 7     S4B     OK     2000 1186331         2        0          0      1867
    ## 8     S4S     OK     2000 1186681         0        0          0      1872
    ## 9     S5B     OK     2000 1186386         1        0          0      1841
    ## 10    S5S     OK     2000 1186501         1        0          0      1861
    ## 11    S6B     OK     2000 1186261         2        0          0      1839
    ## 12    S6S     OK     2000 1187078         1        0          0      1835
    ## 13    S7B     OK     2000 1186888         0        0          0      1825
    ## 14    S7S     OK     2000 1186299         3        0          0      1845
    ## 15    S8B     OK     2000 1186354         3        0          0      1840
    ## 16    S8S     OK     2000 1186610         1        0          0      1848
    ## 17    S9B     OK     2000 1187038         0        0          0      1834
    ## 18    S9S     OK     2000 1186867         0        0          0      1835
    ##    w/adapters qualtrim_bp out_bp w/adapters2 qualtrim2_bp out2_bp
    ## 1        1986           0 513149        1876            0  528595
    ## 2        1975           0 511096        1877            0  525893
    ## 3        1987           0 506659        1850            0  521371
    ## 4        1989           0 504998        1843            0  519979
    ## 5        1989           0 512326        1870            0  527518
    ## 6        1989           0 517598        1891            0  532758
    ## 7        1980           0 514342        1884            0  529379
    ## 8        1987           0 515511        1884            0  530555
    ## 9        1984           0 506972        1856            0  522013
    ## 10       1991           0 512539        1869            0  527592
    ## 11       1981           0 506577        1857            0  521787
    ## 12       1982           0 505929        1851            0  520562
    ## 13       1987           0 503033        1836            0  517931
    ## 14       1987           0 508524        1857            0  523039
    ## 15       1993           0 507178        1847            0  522137
    ## 16       1982           0 509177        1865            0  524085
    ## 17       1983           0 505424        1851            0  520706
    ## 18       1979           0 505519        1853            0  520103

frgrhr
