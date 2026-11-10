# Galada

Analysis code a paper on plant extracts against *Colletotrichum coccodes* and *Helminthosporium solani* (potato blemish diseases).

## Datasets

The analysis uses two datasets, both on Zenodo (CC BY 4.0) and downloaded automatically by the script.

### 1. Ergosterol, chitin and glucan (Figures 6 and 7)

Zenodo: <https://doi.org/10.5281/zenodo.24052059>
**Experimental dataset of quantification of ergosterol, chitin and glucan content in *C. coccodes* and *H. solani* treated with several plant extracts**

One table (`2026_chitin_glucan_ergosterol.csv`, `;`-separated) with the chitin, glucan and ergosterol content of fungal spores treated with plant extracts. It comes from the project "Integrierte Bekämpfung von *C. coccodes* und *H. solani* in der Kartoffelwirtschaft" (Innosuisse/CTI grant 18536.1 PFLS-LS).

| Column | Meaning |
|---|---|
| `Number` | internal sample code |
| `Treatment` | product applied to the spores (PDB, PDB+MeOH, Frangula, Rheum, Curcuma) |
| `Pathogen` | *C. coccodes* or *H. solani* |
| `chitin ug/mg` | chitin, in glucosamine equivalents |
| `glucan ug/mg` | glucan, in glucan equivalents |
| `Ergosterol (ug/mg DW)` | ergosterol |

Empty cells mean the variable was not measured for that sample.

### 2. Post-harvest disease severity (Figure 8)

Zenodo: <https://doi.org/10.5281/zenodo.24111700>, file `postharvest_data_analysis.csv` (`;`-separated, 134 rows)

Post-harvest assessment of black dot (*C. coccodes*, `Cc`) and silver scurf (*H. solani*, `Hs`) on tubers treated with different products.

| Column | Meaning |
|---|---|
| `Unit number` | experimental unit |
| `Treatment`, `Product` | treatment group and product applied (Water, Amistar, Curcuma, Frangula, Rheum) |
| `Year` | year of the trial |
| `Observator` | person who scored the tubers |
| `Cc 0`–`Cc 4`, `Hs 0`–`Hs 4` | number of tubers in each disease severity class (0–4) |
| `Total` | number of tubers scored |
| `Cc incidence`, `Hs incidence` | % of tubers with symptoms |
| `Cc Tax`, `Hs Tax` | disease severity index (%) |

## Script: `Massana_JoF2026.R`

Run from the project folder (`Rscript Massana_JoF2026.R`). It downloads the datasets from Zenodo (Figures 6–7: record 24052059; Figure 8: record 24111700), so no local data file is needed. Required R packages: dplyr, ggplot2, readr, tidyr, tidyverse, gridExtra, RColorBrewer, emmeans.

| Section | What it does |
|---|---|
| **Figure 6**: ergosterol (dataset 1) | Summary statistics (n, mean, SD, median, quartiles, min, max) per pathogen × treatment. Four prespecified Wilcoxon rank-sum tests per pathogen (PDB vs PDB+MeOH, PDB vs Frangula, PDB+MeOH vs Rheum, PDB+MeOH vs Curcuma), Holm-adjusted within each pathogen. Boxplots with individual points, one panel per pathogen. |
| **Figure 7**: chitin and glucan (dataset 1) | Same summaries and the same four planned Wilcoxon comparisons, for glucan and chitin in each pathogen. Four-panel figure (glucan and chitin × pathogen). |
| **Figure 8**: post-harvest disease severity | Black dot and silver scurf severity per product and year (dataset 2). Linear models with year adjustment, then `emmeans` contrasts of each product against Water (Holm). |

Results (summary tables and test results) are printed to the console.

### Outputs (folder `Figure/`, created automatically)

- `figure_6.png`: ergosterol boxplots 
- `figure_7.png`: chitin and glucan panels 
- Figure 8 (`figure_8.png`) and the supplementary tables S3a/S3b are written by lines currently commented out at the end of the script.
