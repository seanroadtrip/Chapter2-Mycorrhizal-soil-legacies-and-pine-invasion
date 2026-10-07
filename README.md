# Soil legacy and shade effects on wilding conifer seedlings (PhD thesis, Chapter 2)

Data and R code for Chapter 2 of the PhD thesis of Sean Davis, School of Biological Sciences, University of Canterbury, Christchurch, New Zealand, 2026.

The chapter tests whether pine legacy soils (from sites previously invaded by wilding conifers) change how wilding conifer seedlings respond to shade, compared with uninvaded grassland soils. Seedlings of *Pinus contorta*, *Pinus radiata*, and *Pseudotsuga menziesii* were grown in mesocosms across a shade gradient of 45 shade boxes (15 nominal shade levels × 3 replicates), in soils from three South Island regions (Craigieburn/Cass, Lake Pukaki, Molesworth). Each soil source was used both as a sterilised bulk substrate and as a live inoculum, giving four soil treatments. We measured final biomass and root tip ectomycorrhization.

## Repository contents

| File | Description |
|---|---|
| `chapter_2_script.R` | Analysis script: fits all models and produces most tables and figures reported in the chapter |
| `potdat_sup.csv` | Seedling-level dataset (one row per pot, 540 pots), described below |
| `potdat_sup.rds` | The same dataset saved as an R object, with factor levels set; this is the file the analysis script reads |

## How to reproduce the analysis

1. Install R (version 4.6.1 was used) and the packages below.
2. Download or clone this repository and set it as your working directory.
3. Run `chapter_2_script.R` from top to bottom.


Fitting requires a working C++ toolchain for Stan (on Windows, the version of Rtools matching your R version).

### Software

| Software | Version |
|---|---|
| R | 4.6.1 |
| brms | 2.23.0 |
| loo | 2.10.1 |

`sessionInfo()` at the end of the script records the full environment.

## Analysis overview

- **Biomass:** a log-scale Ricker function of % PPFD transmitted, with amplitude and decay estimated for each species × inoculum group, lognormal errors, and random intercepts for source region and region × shade box. Derived quantities (peak shade, peak biomass, post-peak decline, probability of the soil curves crossing) are calculated from the posterior draws.
- **Bulk soil check:** the biomass model refitted with an additional effect of pine legacy bulk soil.
- **Ectomycorrhization:** zero-inflated beta regression of the proportion of mycorrhized root tips on light × inoculum × species.
- **Ectomycorrhization and biomass:** biomass models with soil treatment, mycorrhization, or both, compared by approximate leave-one-out cross-validation, plus an estimate of the share of the soil effect accounted for by ectomycorrhization.

All models were fitted in brms with 4 chains × 10,000 iterations (2,500 warmup) and seed 67, pronounced "six-seven".

## Data dictionary: `potdat_sup.csv`

One row per pot. `NA` indicates a missing value (seedling died or a measurement was not collected).

### Identifiers and design

| Column | Description |
|---|---|
| `pot` | Pot number |
| `tent` | Shade box number (1–45) |
| `label` | Treatment label as written on the pot, e.g. `MWG-P PIRA` |
| `region` | Source region of soil: `CB` Craigieburn/Cass, `LP` Lake Pukaki, `MW` Molesworth (derived from label)|
| `bulk` | Sterilised bulk substrate: `G` uninvaded grassland, `P` pine legacy (derived from label)|
| `inoculate` | Live soil inoculum: `G` uninvaded grassland, `P` pine legacy (derived from label)|
| `species` | `PICO` *Pinus contorta*, `PIRA` *Pinus radiata*, `PSME` *Pseudotsuga menziesii* (derived from label)|
| `analysed` | `TRUE` for the 506 seedlings with complete data, used in all models |

### Light

| Column | Description |
|---|---|
| `ppfd_trans` | Photosynthetic photon flux density transmitted to pot height, % of above-box PPFD (one value per shade box) |
| `ppfd_red` | PPFD reduced (shade), % (= 100 − `ppfd_trans`) |
| `ppfd_prop` | `ppfd_trans` as a proportion (0–1), used in the ectomycorrhization model |

### Growth and biomass

| Column | Description |
|---|---|
| `height` | Seedling height at harvest [mm] |
| `diameter_mm` | Stem diameter at harvest (mm) |
| `shoot_dry_g` | Shoot dry mass (g) |
| `root_wet_g` | Whole root system wet mass (g) |
| `root_sub_wet_g` | Wet mass of the root subsample (g) |
| `root_sub_dry_g` | Dry mass of the root subsample (g) |
| `root_dry_g` | Whole root dry mass (g), scaled from the subsample: `root_sub_dry_g × root_wet_g / root_sub_wet_g` |
| `final_mass` | Total seedling dry mass (g): `shoot_dry_g + root_dry_g`. Response variable for the biomass models |
| `root_shoot` | Root:shoot dry mass ratio |

### Roots and ectomycorrhization

Root length was estimated by the line intersect method (Newman 1966; Tennant 1975). Root tips were classified as mycorrhized or non-mycorrhized.

| Column | Description |
|---|---|
| `nonecm_tips` | Number of non-mycorrhized root tips |
| `ecm_tips` | Number of ectomycorrhized root tips |
| `total_tips` | `nonecm_tips + ecm_tips` |
| `ecm_colonization` | Root tips ectomycorrhized (%): `100 × ecm_tips / total_tips`. Response variable for the ectomycorrhization model |
| `ecm_prop` | `ecm_colonization` as a proportion (0–1) |
| `nonecm_intersects` | Gridline intersections with non-mycorrhized roots |
| `ecm_intersects` | Gridline intersections with ectomycorrhized roots |
| `total_intersects` | `nonecm_intersects + ecm_intersects` |
| `grid_area` | Area of the counting grid [mm]; 236.8 for the small grid used in shade boxes 8, 13, 22, 37, 40, and 43, otherwise 259.4 |
| `grid_line_len` | Total length of gridlines [mm]; 113.1 (small grid) or 132.7 |
| `root_length` | Root length [mm]: `π × total_intersects × grid_area / (2 × grid_line_len)` |
| `ecm_root_length` | Ectomycorrhized root length [mm], same formula using `ecm_intersects` |
| `ecm_length_pct` | Root length ectomycorrhized (%) |

## Citation

If you use these data or code, please cite:

Davis, S. (2026). Restoration of indigenous shade canopies to resist wilding conifer reinvasion. PhD thesis, University of Canterbury, Christchurch, New Zealand.

and the archived version of this repository.

## Licence

 CC BY 4.0 for the data and MIT for the code.

## Contact

Sean Davis, School of Biological Sciences, University of Canterbury. [seanroadtrip@gmail.com]
