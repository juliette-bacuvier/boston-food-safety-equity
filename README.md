# boston-food-safety-equity

**Analysis of food safety inspection outcomes in Boston by neighbourhood socioeconomic and demographic characteristics** — Wellesley College, DS340H Data Science Capstone, Spring 2026.

Author: Juliette Bacuvier

---

## Project summary

This project asks whether food safety inspection outcomes in Boston systematically differ across neighbourhoods of varying socioeconomic status, and whether any such gap has widened between 2015 and 2025. Boston food establishment inspection records are linked to tract-level demographics from the American Community Survey (ACS) and modelled using mixed-effects logistic regression.

**Main finding:** Neighbourhood income has a real but modest association with inspection pass rates. Establishments in low-income tracts have roughly a 3.4 percentage-point lower predicted pass probability than those in high-income tracts. The gap is statistically stable across the ten-year window and is dwarfed by within-establishment variability and by differences across business types.

The full write-up lives in the compiled poster (`DS_340_H___Project_Poster (2).pdf`).

---

## ⚠️ Dataset not included ⚠️

**The primary inspection dataset is not in this repository because the file is too large to commit to GitHub.** The raw Boston Food Establishment Inspections file contains 872,111 violation records and exceeds GitHub's 100 MB file size limit.

To reproduce the analysis, download the data manually:

- **Boston Food Establishment Inspections (2015–2025):** available from [Analyse Boston (City of Boston Open Data)](https://data.boston.gov/dataset/food-establishment-inspections). Download the CSV and place it in the repository root (or update the file path at the top of `BacuvierJ.DS340H.Project.R`).
- **ACS 2019 5-year estimates:** pulled via the `tidycensus` R package inside the script. Requires a free Census API key: set it once with `tidycensus::census_api_key("YOUR_KEY", install = TRUE)`.
- **Census tract shapefiles (Massachusetts, 2019):**  fetched via the `tigris` R package inside the script; cached locally on first run.

---

## Repository contents

| File | Description |
|---|---|
| `BacuvierJ.DS340H.Project.R` | Main analysis script — data cleaning, spatial join to ACS tracts, EDA, mixed-effects modelling, diagnostics, and figure export. |
| `DS_340_H___Project_Poster (2).pdf` | Final version of the capstone poster. |
| `README.md` | This file. |
| `LICENSE` | MIT license. |
| `.gitignore` | Excludes the raw inspection CSV and cached shapefiles from version control. |

---

## Data

| Source | Description | Records |
|---|---|---|
| Boston Food Establishment Inspections (2015–2025) | Inspection results with violations, dates, business types, addresses | 872,111 raw $\rightarrow$ 100,211 unique inspections after deduplication |
| ACS 2019 5-year estimates | Tract-level median household income, poverty rate, racial composition, foreign-born share | 174 census tracts |
| TIGER/Line shapefiles | Census tract polygons for spatial join | 174 polygons |

**Final analysis sample:** ~91,600 inspections across 6,912 establishments in 174 census tracts. Overall fail rate: 46% (Boston uses an inclusive definition: any inspection with at least one violation of any severity counts as a fail).

---

## Methods

- **Mixed-effects logistic regression** via `glmmTMB`, with random intercepts by establishment to account for repeated inspections.
- **Sequential model building:** null $\rightarrow$ income-only $\rightarrow$ SES composite $\rightarrow$ fully adjusted $\rightarrow$ income × time interaction.
- **PCA-based SES index** to address multicollinearity between income and poverty (r = −0.81). PC1 captures 79% of variance across income, poverty, % white, and foreign-born share.
- **Train/test split by establishment** (80/20) to avoid information leakage from repeat visits.
- **Robustness checks:** COVID-period exclusion, nonwhite-composition sensitivity, DHARMa residual diagnostics, held-out AUC.

---

## Key results

- **Null model ICC = 0.046** $\rightarrow$ only 4.6% of pass/fail variation lies between establishments; the remaining 95% is within-establishment (visit-level) noise.
- **Income effect:** +1 SD in tract income (~$40K) $\rightarrow$ odds of passing ↑ 6.4%. Low-income Food Service: 50.1% predicted pass; high-income: 53.5%.
- **Time interaction:** income × year is statistically significant (p = 0.006) but practically negligible (OR = 1.022) $\rightarrow$ the gap has not meaningfully widened.
- **Business type matters more than SES:** Mobile food OR ≈ 1.42; Retail food OR ≈ 1.18.

---

## How to reproduce

### Prerequisites

- R ≥ 4.2
- R packages: `glmmTMB`, `tidycensus`, `tigris`, `sf`, `tidyverse`, `DHARMa`, `pROC`, `broom.mixed`

```r
install.packages(c("glmmTMB", "tidycensus", "tigris", "sf", "tidyverse", "DHARMa", "pROC", "broom.mixed"))
```

### Steps

1. Download the raw Boston inspections CSV from [Analyse Boston](https://data.boston.gov/dataset/food-establishment-inspections) and place it in the repository root.
2. Set your Census API key (one-time): `tidycensus::census_api_key("YOUR_KEY", install = TRUE)`
3. Open `BacuvierJ.DS340H.Project.R` in RStudio and run it end-to-end, or from the terminal:
   ```bash
   Rscript BacuvierJ.DS340H.Project.R
   ```

---

## Limitations

- Single-city analysis $\rightarrow$ findings may not generalise.
- ACS 2019 demographics applied to the full 2015–2025 window.
- Pass/fail is coarse $\rightarrow$ doesn't capture violation severity.
- Inspection frequency is endogenous to outcomes (failures trigger re-inspections), which risks collider bias on the SES coefficient.
- Spatial autocorrelation is not modelled.

See the poster for the full caveat discussion.

---

## Citation

> Bacuvier, J. (2026). *Mapping Inequality on the Boston Plate: Neighbourhood Socioeconomics and Food Inspection Outcomes in Boston, 2015–2025.* Wellesley College, DS340H Data Science Capstone.

---

## License

MIT $\rightarrow$ see `LICENSE`.
