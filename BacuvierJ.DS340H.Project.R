library(tidyverse)
library(lubridate)
library(sf)
library(tigris)
library(leaflet)
library(tidycensus)
library(glmmTMB)
library(ggcorrplot)
library(pROC)
library(DHARMa)
library(car)
library(tmap)

# Resolve namespace conflicts: car masks dplyr::recode and can interfere
# with dplyr::filter in some edge cases
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("recode", "dplyr")

options(tigris_use_cache = TRUE)

# ==============================================================================
# 1. LOADING DATA
# ==============================================================================

raw_data = read_csv('/Users/juliettebacuvier/desktop/tmpp7jeda9f.csv')

# ==============================================================================
# 2. DATE PARSING & FILTERING
# ==============================================================================

raw_data = raw_data %>%
  mutate(insp_date = ymd_hms(resultdttm), year = year(insp_date)) %>%
  filter(year >= 2015, year <= 2025, !is.na(insp_date))

# ==============================================================================
# 3. OUTCOME VARIABLE
# ==============================================================================

cat("=== viol_status distribution ===\n")
raw_data %>% count(viol_status, sort = TRUE) %>% print(n = 20)

cat("\n=== result distribution ===\n")
raw_data %>% count(result, sort = TRUE) %>% print(n = 20)

# If viol_status varies within a visit, it's per-violation not per-inspection
conflicting = raw_data %>%
  group_by(licenseno, insp_date) %>%
  summarise(n_statuses = n_distinct(viol_status), .groups = "drop") %>%
  filter(n_statuses > 1)
cat("\nVisits with conflicting viol_status:", nrow(conflicting), "\n")

# pass = 1, fail = 0
raw_data = raw_data %>%
  mutate(pass_fail = case_when(
    str_detect(result, regex("pass", ignore_case = TRUE)) ~ 1L,
    str_detect(result, regex("fail", ignore_case = TRUE)) ~ 0L,
    TRUE ~ NA_integer_))

# ==============================================================================
# 4. DEDUPLICATION (moved before spatial join for efficiency)
# ==============================================================================

'
Raw data is one row per violation per visit.
We want one row per visit, keeping the worst outcome (fail = 0).
Moving this before the spatial join avoids expensive geometry operations
on ~2-3x more rows than necessary.
'
raw_data = raw_data %>%
  filter(!is.na(pass_fail)) %>%
  group_by(licenseno, insp_date) %>%
  arrange(pass_fail) %>%
  slice(1) %>%
  ungroup()

cat("\nAfter deduplication:", nrow(raw_data), "unique inspections\n")

# ==============================================================================
# 5. GEOGRAPHY PARSING
# ==============================================================================

raw_data = raw_data %>%
  mutate(
    location = str_remove_all(location, "[()]"),
    lat = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 1])),
    lon = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 2])))

# ==============================================================================
# 6. ACS CENSUS DATA
# ==============================================================================

'
Note: Using ACS 2019 5-year estimates for inspections spanning 2015-2025.
This is a deliberate choice -- 2019 is the centre of the time range and
pre-dates COVID disruptions to Census data collection.
Tract demographics do shift over a decade; this is acknowledged as a limitation.
'

acs_data = get_acs(
  geography = "tract", state = "MA",
  variables = c(
    median_income = "B19013_001",
    total_pop = "B01003_001",
    below_poverty = "B17001_002",
    total_poverty = "B17001_001",
    white_alone = "B02001_002",
    black_alone = "B02001_003",
    indigenous_alone = "B02001_004",
    asian_alone = "B02001_005",
    hawaiian_alone = "B02001_006",
    other_race_alone = "B02001_007",
    total_race = "B02001_001",
    hispanic = "B03002_012",
    total_hispanic = "B03002_001",
    foreign_born = "B05002_013",
    total_nativity = "B05002_001"),
  year = 2019, survey = "acs5")

# Pivot wide and compute derived rates
acs_clean = acs_data %>%
  select(GEOID, variable, estimate) %>%
  pivot_wider(names_from = variable, values_from = estimate) %>%
  mutate(
    poverty_rate = ifelse(total_poverty == 0, NA_real_, below_poverty / total_poverty),
    pct_white = ifelse(total_race == 0, NA_real_, white_alone / total_race),
    pct_black = ifelse(total_race == 0, NA_real_, black_alone / total_race),
    pct_indigenous = ifelse(total_race == 0, NA_real_, indigenous_alone / total_race),
    pct_asian = ifelse(total_race == 0, NA_real_, asian_alone / total_race),
    pct_hawaiian = ifelse(total_race == 0, NA_real_, hawaiian_alone / total_race),
    pct_other_race = ifelse(total_race == 0, NA_real_, other_race_alone / total_race),
    pct_hispanic = ifelse(total_hispanic == 0, NA_real_, hispanic / total_hispanic),
    pct_foreign_born = ifelse(total_nativity == 0, NA_real_, foreign_born / total_nativity),
    pct_nonwhite = 1 - pct_white)

# ==============================================================================
# 7. SPATIAL JOIN
# ==============================================================================

ma_tracts = tracts(state = "MA", year = 2019, cb = TRUE) %>%
  st_transform(crs = 4326)

analysis_data = raw_data %>%
  filter(!is.na(lat), !is.na(lon))

inspections_sf = analysis_data %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

inspections_sf = st_join(inspections_sf, ma_tracts %>% select(GEOID), join = st_within)

analysis_data$tract_geoid = inspections_sf$GEOID

analysis_data = analysis_data %>%
  left_join(acs_clean, by = c("tract_geoid" = "GEOID"))

# ==============================================================================
# 8. CONTROL VARIABLES
# ==============================================================================

# inspection frequency per establishment per year
insp_freq = analysis_data %>%
  group_by(licenseno, year) %>%
  summarise(insp_per_year = n(), .groups = "drop")

analysis_data = analysis_data %>%
  left_join(insp_freq, by = c("licenseno", "year"))

# license age in years at time of inspection
analysis_data = analysis_data %>%
  mutate(iss_date = ymd_hms(issdttm),
         license_age = as.numeric(difftime(insp_date, iss_date, units = "days")) / 365.25,
         license_age = pmax(license_age, 0, na.rm = TRUE))

# business type
analysis_data = analysis_data %>%
  mutate(business_type = case_when(
    licensecat == "FS" ~ "Food Service",
    licensecat == "FT" ~ "Food Truck/Mobile",
    licensecat == "RF" ~ "Retail Food",
    licensecat == "MFW" ~ "Manufacturing/Wholesale",
    TRUE ~ "Other"),
    business_type = factor(business_type, levels = c("Food Service", "Retail Food",
                                                     "Food Truck/Mobile", "Manufacturing/Wholesale", "Other")))

# continuous time variable for regression
analysis_data = analysis_data %>%
  mutate(time_scaled = scale(as.numeric(insp_date))[, 1])

# income groupings
analysis_data = analysis_data %>%
  mutate(
    income_tercile = ntile(median_income, 3),
    income_tercile = factor(income_tercile, labels = c("Low Income", "Middle Income", "High Income")),
    income_quartile = ntile(median_income, 4),
    income_quartile = factor(income_quartile, labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)")))

analysis_data = analysis_data %>%
  mutate(acs_missing = is.na(median_income), zip = as.character(zip), year = as.integer(year))

# ==============================================================================
# 9. DIAGNOSTICS
# ==============================================================================

cat("\n=== Final dataset ===\n")
cat("Rows:", nrow(analysis_data), "\n")
cat("Unique establishments:", n_distinct(analysis_data$licenseno), "\n")
cat("Unique tracts:", n_distinct(analysis_data$tract_geoid, na.rm = TRUE), "\n")
cat("Pass rate:", round(mean(analysis_data$pass_fail == 1, na.rm = TRUE) * 100, 1), "%\n")
cat("Fail rate:", round(mean(analysis_data$pass_fail == 0, na.rm = TRUE) * 100, 1), "%\n")
cat("Year range:", min(analysis_data$year), "-", max(analysis_data$year), "\n")

cat("\n=== Missingness ===\n")
analysis_data %>%
  summarise(across(c(pass_fail, lat, lon, tract_geoid, median_income, poverty_rate,
                     pct_white, pct_black, pct_hispanic, pct_foreign_born, license_age,
                     insp_per_year, business_type),
                   ~ sum(is.na(.)) )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(analysis_data) * 100, 1)) %>%
  arrange(desc(n_missing)) %>%
  print(n = 20)

# ==============================================================================
# 10. EDA VISUALISATIONS
# ==============================================================================

# Naming colours
red = "#b90d0d"
yellow = "#f3d357"
blue = "#2a7b9b"
white = "white"
gold = "#f0a500"
colour_palette = c(red, gold, blue)
heat_palette = c(blue, white, red)

# --- Overall outcome distribution ---
'
This is a nearly 50/50 outcome, making prediction inherently difficult and meaning 
even small shifts in probability are hard to detect visually.
'
yearly_outcomes = analysis_data %>%
  filter(!is.na(pass_fail)) %>%
  group_by(year) %>%
  summarise(n = n(), fail_rate = mean(1 - pass_fail, na.rm = TRUE), .groups = "drop")

ggplot(yearly_outcomes, aes(x = year, y = fail_rate)) +
  geom_col(fill = red, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = paste0(round(fail_rate * 100, 1), "%")), vjust = -0.5, size = 3.5) +
  geom_text(aes(label = paste0("n=", scales::comma(n))), y = 0.02, size = 2.8, color = "white") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.6),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title = "Overall Inspection Fail Rate by Year",
    subtitle = "Fail rates hover around 45-50% across the study period",
    x = "Year", y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot0_overall_fail_rate_by_year.png", width = 9, height = 5, dpi = 300)

# --- Choropleth map (interactive) ---
neighbourhood_summary = analysis_data %>%
  dplyr::filter(!is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(n_inspections = n(), fail_rate = mean(1 - pass_fail, na.rm = TRUE), .groups = "drop")

plot_data = neighbourhood_summary %>% filter(n_inspections >= 10)

map_tracts = ma_tracts %>%
  filter(GEOID %in% plot_data$tract_geoid) %>%
  left_join(plot_data, by = c("GEOID" = "tract_geoid"))

pal = colorNumeric(palette = c(yellow, red), domain = map_tracts$fail_rate)

leaflet(map_tracts) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(fillColor = ~pal(fail_rate), fillOpacity = 0.7, color = "white", weight = 1,
              popup = ~paste0(
                "<b>Tract ", GEOID, "</b><br>",
                "Fail rate: ", sprintf("%.1f%%", fail_rate * 100), "<br>",
                "Inspections: ", n_inspections)
  ) %>%
  addLegend(pal = pal, values = ~fail_rate, title = "Proportion Failed",
            labFormat = labelFormat(suffix = "%", transform = function(x) x * 100),
            position = "bottomright")

# --- Static choropleth (tmap) for portfolio/report ---
tmap_mode("plot")

static_map = tm_shape(map_tracts) +
  tm_fill("fail_rate", palette = c(yellow, red), style = "cont",
          title = "Fail Rate", legend.format = list(fun = function(x) paste0(round(x * 100), "%"))) +
  tm_borders(col = "white", lwd = 0.5) +
  tm_layout(
    title = "Food Inspection Fail Rate by Census Tract",
    title.size = 1.2, title.fontface = "bold",
    legend.position = c("right", "bottom"),
    frame = FALSE)

tmap_save(static_map, "plot_choropleth_static.png", width = 9, height = 7, dpi = 300)

# --- Tract-level summary for scatterplots ---
tract_summary = analysis_data %>%
  filter(!is.na(median_income), !is.na(pass_fail), !is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(
    n_inspections = n(),
    fail_rate = mean(1 - pass_fail, na.rm = TRUE),
    median_income = first(median_income),
    poverty_rate = first(poverty_rate),
    pct_white = first(pct_white),
    pct_black = first(pct_black),
    pct_hispanic = first(pct_hispanic),
    pct_foreign_born = first(pct_foreign_born),
    pct_nonwhite = first(pct_nonwhite),
    .groups = "drop") %>%
  filter(n_inspections >= 10)

# --- Fail rate vs. median income ---
ggplot(tract_summary, aes(x = median_income / 1000, y = fail_rate)) +
  geom_point(aes(size = n_inspections), alpha = 0.5, color = red) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.8) +
  scale_x_continuous(labels = scales::dollar_format(suffix = "K")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_size_continuous(range = c(1.5, 8), name = "Inspections") +
  labs(
    title = "Food Inspection Fail Rate vs. Median Household Income",
    subtitle = "Each point is a census tract (min. 10 inspections)",
    x = "Median Household Income (ACS 2019, tract level)",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot1_failrate_vs_income_tract.png", width = 9, height = 6, dpi = 300)

# --- Fail rate over time by income tercile ---

# one income value per tract to avoid duplicate rows in the join
tract_income = analysis_data %>%
  filter(!is.na(median_income), !is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(median_income = first(median_income), .groups = "drop") %>%
  mutate(income_tercile = ntile(median_income, 3),
         income_tercile = factor(income_tercile, labels = c("Low Income", "Middle Income", "High Income")))

time_trends = analysis_data %>%
  filter(!is.na(median_income), !is.na(pass_fail), !is.na(tract_geoid)) %>%
  select(-income_tercile) %>%
  left_join(tract_income %>% select(tract_geoid, income_tercile), by = "tract_geoid") %>%
  group_by(year, income_tercile) %>%
  summarise(fail_rate = mean(1 - pass_fail, na.rm = TRUE), n = n(), .groups = "drop")

ggplot(time_trends, aes(x = year, y = fail_rate, color = income_tercile)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_color_manual(values = colour_palette, name = "Neighbourhood\nIncome Level") +
  labs(
    title = "Inspection Fail Rate Over Time by Neighbourhood Income Level",
    subtitle = "Census tracts grouped into income terciles (ACS 2019 median household income)\nNote: y-axis is compressed; the absolute gap between groups is ~3-6 pp",
    x = "Year",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot2_failrate_time_trends_tract.png", width = 9, height = 6, dpi = 300)

# --- Income gap over time ---
# Makes the magnitude of the gap explicit rather than relying on
# visual comparison across compressed y-axis
income_gap = time_trends %>%
  select(year, income_tercile, fail_rate) %>%
  pivot_wider(names_from = income_tercile, values_from = fail_rate) %>%
  mutate(gap_pp = (`Low Income` - `High Income`) * 100)

ggplot(income_gap, aes(x = year, y = gap_pp)) +
  geom_line(linewidth = 1, color = red) +
  geom_point(size = 2.5, color = red) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title = "Fail Rate Gap: Low-Income vs. High-Income Neighbourhoods",
    subtitle = "Positive values = low-income tracts fail more often; gap is modest but persistent",
    x = "Year",
    y = "Fail Rate Gap (percentage points)",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot2b_income_gap_over_time.png", width = 9, height = 5, dpi = 300)

# --- Fail rate by business type ---
business_summary = analysis_data %>%
  filter(!is.na(pass_fail)) %>%
  group_by(business_type) %>%
  summarise(n_inspections = n(), fail_rate = mean(1 - pass_fail, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(fail_rate))

ggplot(business_summary, aes(x = reorder(business_type, fail_rate), y = fail_rate)) +
  geom_col(fill = red, alpha = 0.8, width = 0.6) +
  geom_text(aes(label = paste0(round(fail_rate * 100, 1), "%")), hjust = -0.15, size = 3.5) +
  geom_text(aes(label = paste0("n=", scales::comma(n_inspections))), y = 0.002, hjust = 0, size = 3, color = "white") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = "Inspection Fail Rate by Business Type",
    subtitle = "Justifies including business type as a control variable",
    x = NULL, y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data"
  ) + theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), panel.grid.major.y = element_blank())

ggsave("plot3_failrate_by_business.png", width = 9, height = 5, dpi = 300)

# --- Correlation heatmap ---
cor_vars = analysis_data %>%
  filter(!acs_missing, !is.na(tract_geoid)) %>%
  distinct(tract_geoid, .keep_all = TRUE) %>%
  select(median_income, poverty_rate, pct_white, pct_black, pct_asian, pct_hispanic, pct_foreign_born) %>%
  drop_na()

cor_matrix = cor(cor_vars, use = "complete.obs")
colnames(cor_matrix) = rownames(cor_matrix) =
  c("Median Income", "Poverty Rate", "% White", "% Black", "% Asian", "% Hispanic", "% Foreign Born")

ggcorrplot(cor_matrix, type = "lower", lab = TRUE, lab_size = 3.5, colors = heat_palette,
           title = "Correlation Between Tract-Level Demographic Variables",
           ggtheme = theme_minimal(base_size = 13)) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plot4_correlation_heatmap_tract.png", width = 8, height = 7, dpi = 300)

# --- Fail rate vs. poverty rate ---
ggplot(tract_summary, aes(x = poverty_rate, y = fail_rate)) +
  geom_point(aes(size = n_inspections), alpha = 0.5, color = red) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.8) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_size_continuous(range = c(1.5, 8), name = "Inspections") +
  labs(
    title = "Food Inspection Fail Rate vs. Poverty Rate",
    subtitle = "Each point is a census tract (min. 10 inspections); relationship is near-zero",
    x = "Poverty Rate (ACS 2019, tract level)",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot5_failrate_vs_poverty_tract.png", width = 9, height = 6, dpi = 300)

# --- Inspection frequency by income tercile ---
# Result: minimal difference (3.64 vs 3.45) across income groups.
# This is a useful null finding -- argues against systematic targeting bias.
insp_freq_by_income = analysis_data %>%
  filter(!is.na(income_tercile), !is.na(insp_per_year)) %>%
  group_by(income_tercile) %>%
  summarise(mean_insp_freq = mean(insp_per_year, na.rm = TRUE), median_insp_freq = median(insp_per_year, na.rm = TRUE),
            n = n(), .groups = "drop")

ggplot(insp_freq_by_income, aes(x = income_tercile, y = mean_insp_freq, fill = income_tercile)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = round(mean_insp_freq, 2)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = colour_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Mean Inspection Frequency by Neighbourhood Income Level",
    subtitle = "Minimal difference across groups -- little evidence of systematic targeting bias",
    x = "Neighbourhood Income Level",
    y = "Mean Inspections per Establishment per Year",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot6_insp_freq_by_income.png", width = 8, height = 5, dpi = 300)

# ==============================================================================
# 11. CROSS-VALIDATION SPLIT
# Holds back 20% of establishments before modelling
# Splits by establishment, not by row, to avoid leakage
# ==============================================================================

set.seed(42)

all_establishments = analysis_data %>%
  filter(!is.na(median_income), !is.na(pass_fail), !is.na(license_age), !is.na(insp_per_year)) %>%
  distinct(licenseno) %>%
  pull(licenseno)

test_establishments = sample(all_establishments, size = round(0.2 * length(all_establishments)))

train_data = analysis_data %>%
  filter(licenseno %in% all_establishments, !licenseno %in% test_establishments)

test_data = analysis_data %>%
  filter(licenseno %in% test_establishments)

cat("\n=== Train/test split ===\n")
cat("Training:", n_distinct(train_data$licenseno), "establishments,", nrow(train_data), "inspections\n")
cat("Test:", n_distinct(test_data$licenseno), "establishments,", nrow(test_data), "inspections\n")

# ==============================================================================
# 12. MODEL DATA PREPARATION
# ==============================================================================

'
Keeps only establishments with 2+ inspections.
Rationale: establishments with a single inspection cannot contribute to
the within-establishment variance estimate needed for the random intercept.
'
model_data = train_data %>%
  group_by(licenseno) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  filter(!is.na(median_income), !is.na(poverty_rate), !is.na(pct_white), !is.na(license_age),
         !is.na(insp_per_year), !is.na(business_type), !is.na(time_scaled), !is.na(pct_nonwhite))

# Store scaling parameters explicitly so train and test use the same values
scaling_params = list(
  income_mean = mean(model_data$median_income, na.rm = TRUE),
  income_sd = sd(model_data$median_income, na.rm = TRUE),
  poverty_mean = mean(model_data$poverty_rate, na.rm = TRUE),
  poverty_sd = sd(model_data$poverty_rate, na.rm = TRUE),
  pct_white_mean = mean(model_data$pct_white, na.rm = TRUE),
  pct_white_sd = sd(model_data$pct_white, na.rm = TRUE),
  pct_nonwhite_mean = mean(model_data$pct_nonwhite, na.rm = TRUE),
  pct_nonwhite_sd = sd(model_data$pct_nonwhite, na.rm = TRUE),
  license_age_mean = mean(model_data$license_age, na.rm = TRUE),
  license_age_sd = sd(model_data$license_age, na.rm = TRUE),
  insp_freq_mean = mean(model_data$insp_per_year, na.rm = TRUE),
  insp_freq_sd = sd(model_data$insp_per_year, na.rm = TRUE))

model_data = model_data %>%
  mutate(
    income_scaled = (median_income - scaling_params$income_mean) / scaling_params$income_sd,
    poverty_scaled = (poverty_rate - scaling_params$poverty_mean) / scaling_params$poverty_sd,
    pct_white_scaled = (pct_white - scaling_params$pct_white_mean) / scaling_params$pct_white_sd,
    pct_nonwhite_scaled = (pct_nonwhite - scaling_params$pct_nonwhite_mean) / scaling_params$pct_nonwhite_sd,
    license_age_scaled = (license_age - scaling_params$license_age_mean) / scaling_params$license_age_sd,
    insp_freq_scaled = (insp_per_year - scaling_params$insp_freq_mean) / scaling_params$insp_freq_sd)

cat("\nModel dataset:", nrow(model_data), "inspections,", n_distinct(model_data$licenseno), "establishments\n")

# ==============================================================================
# 13. COLLINEARITY DIAGNOSTICS
# ==============================================================================

'
The correlation heatmap flagged high correlations between income, poverty,
and pct_white (income-poverty r = -0.81). VIFs from a fixed-effects GLM
check whether this is problematic for the regression.

RESULT: All VIFs < 5 (income = 4.13, poverty = 2.84, pct_white = 2.19),
within conventional thresholds. However, the poverty coefficient flips sign
in the full model (see Section 14 notes), confirming that individual
coefficients on correlated SES variables should be interpreted with caution.
The PCA composite addresses this.
'

vif_model = glm(pass_fail ~ income_scaled + poverty_scaled + pct_white_scaled +
                  license_age_scaled + insp_freq_scaled + business_type + time_scaled,
                data = model_data, family = binomial)

vif_values = vif(vif_model)
print(vif_values)

'
Composite SES index via PCA.
PC1 captures the shared socioeconomic dimension across income, poverty, and
racial composition without the instability of including all three separately.
PC1 explains ~79% of variance with loadings:
  income (+0.62), poverty (-0.57), pct_white (+0.55)
This aligns with expectations: higher SES = higher income, lower poverty, whiter.
'
ses_vars = model_data %>%
  select(income_scaled, poverty_scaled, pct_white_scaled) %>%
  drop_na()

ses_pca = prcomp(ses_vars, center = FALSE, scale. = FALSE)
model_data$ses_index = ses_pca$x[, 1]

# Store PCA rotation for applying to test data later
ses_rotation = ses_pca$rotation[, 1]

# ==============================================================================
# 14. MIXED-EFFECTS LOGISTIC REGRESSION (using glmmTMB)
# ==============================================================================

'
Model 0: Null model.
Random intercept only. The ICC tells us what share of pass/fail variation
is between vs. within establishments.

RESULT: ICC = 0.046 -- only 4.6% of variation is between establishments.
This is low: most variation is WITHIN establishments over time (the same
restaurant sometimes passes, sometimes fails). Pass/fail outcomes are
noisy and variable for everyone, not clustered by establishment identity.
'

m0 = glmmTMB(pass_fail ~ 1 + (1 | licenseno), data = model_data, family = binomial)

re_var_m0 = VarCorr(m0)$cond$licenseno[1, 1]
icc = re_var_m0 / (re_var_m0 + (pi^2 / 3))
cat("\nNull model ICC:", round(icc, 3), "(", round(icc * 100, 1), "% between-establishment)\n")

# --- Model 1a: Income only ---
m1a = glmmTMB(pass_fail ~ income_scaled + (1 | licenseno), data = model_data, family = binomial)
income_or_unadjusted = exp(fixef(m1a)$cond["income_scaled"])

# --- Model 1b: Poverty only (robustness) ---
# Unadjusted poverty OR = 0.979 (higher poverty -> slightly worse outcomes).
# This is the expected direction but barely significant (p = 0.035).
m1b = glmmTMB(pass_fail ~ poverty_scaled + (1 | licenseno), data = model_data, family = binomial)

# --- Model 1c: Composite SES index ---
# Cleanest neighbourhood-SES predictor: avoids collinearity entirely.
# OR = 1.065 per unit of PC1 (higher SES -> higher odds of passing).
m1c = glmmTMB(pass_fail ~ ses_index + (1 | licenseno), data = model_data, family = binomial)

'
Model 2: Fully adjusted.

KEY INTERPRETATION NOTE: The poverty coefficient FLIPS SIGN here.
Unadjusted, poverty OR = 0.979 (more poverty -> worse outcomes, as expected).
Adjusted, poverty OR = 1.114 (more poverty -> BETTER outcomes?!).

This is a classic suppression effect from collinearity. When you condition
on income AND pct_white simultaneously, the residual variation in poverty
no longer represents "disadvantage" -- it captures whatever is left over
after income and racial composition are accounted for, which could be
urbanicity, density, or other confounded characteristics.

This is exactly why the SES composite (m1c) is the more reliable measure
of the neighbourhood SES effect. The individual coefficients on income,
poverty, and pct_white in this model should NOT be interpreted as
independent causal effects.

The strongest predictor by far is INSPECTION FREQUENCY (OR = 0.665):
establishments inspected more often are much less likely to pass on any
given visit. IMPORTANT CAVEAT: this likely reflects a mechanical artifact
rather than a causal effect. When an establishment fails, it gets
re-inspected, which simultaneously inflates its inspection count AND its
fail count. The negative correlation between frequency and passing is at
least partly an artifact of this enforcement feedback loop.
'

m2 = glmmTMB(pass_fail ~ income_scaled + poverty_scaled + pct_white_scaled + license_age_scaled +
               insp_freq_scaled + business_type + time_scaled + (1 | licenseno),
             data = model_data, family = binomial)

or_adjusted = exp(fixef(m2)$cond)
ci_adjusted = exp(confint(m2, parm = "beta_", method = "Wald"))

or_table = data.frame(
  term = names(or_adjusted),
  OR = round(or_adjusted, 3),
  CI_lower = round(ci_adjusted[, 1], 3),
  CI_upper = round(ci_adjusted[, 2], 3)
) %>%
  filter(term != "(Intercept)")

cat("\n=== Odds Ratios: fully adjusted model ===\n")
print(or_table)

income_or_adjusted = exp(fixef(m2)$cond["income_scaled"])

# --- Model 2b: pct_nonwhite robustness check ---
# Mirrors pct_white coefficient in opposite direction, as expected.
# pct_nonwhite OR = 0.872 vs pct_white OR = 1.147 -- reciprocal relationship holds.
m2b = glmmTMB(pass_fail ~ income_scaled + poverty_scaled + pct_nonwhite_scaled + license_age_scaled +
                insp_freq_scaled + business_type + time_scaled + (1 | licenseno),
              data = model_data, family = binomial)

'
Model 3: Income x time interaction.

RESULT: The interaction is statistically significant (OR = 1.022, p = 0.006)
but the effect size is practically negligible. A 2.2% change in the income
OR per SD of time means the income gradient is shifting by a trivial amount
per decade. The predicted probability plot confirms this: the lines for
different income groups are visually parallel.

Interpretation: the income-outcome gap is widening very slightly over time,
but the magnitude is too small to be substantively meaningful. For practical
purposes, the gap is stable.
'

m3 = glmmTMB(pass_fail ~ income_scaled * time_scaled + poverty_scaled + pct_white_scaled +
               license_age_scaled + insp_freq_scaled + business_type + (1 | licenseno),
             data = model_data, family = binomial)

interaction_coef = fixef(m3)$cond["income_scaled:time_scaled"]
interaction_or = exp(interaction_coef)
interaction_p = summary(m3)$coefficients$cond["income_scaled:time_scaled", "Pr(>|z|)"]

# ==============================================================================
# 15. RANDOM EFFECT VARIANCE DECOMPOSITION
# ==============================================================================

'
The RE variance drops dramatically from null (0.160) to fully adjusted (0.020),
an 87.7% reduction. This means the fixed effects -- primarily inspection
frequency and time -- absorb nearly all of the between-establishment variance
that did exist.

Note that the ICC was already low (4.6%), so the between-establishment
variance was small to begin with. The fixed effects explain that small
between-establishment component almost completely, but the vast majority
of variation was always within-establishment noise.
'

re_variances = data.frame(
  Model = c("Null", "Income only", "SES index", "Fully adjusted", "Income x time"),
  RE_var = round(c(
    VarCorr(m0)$cond$licenseno[1, 1],
    VarCorr(m1a)$cond$licenseno[1, 1],
    VarCorr(m1c)$cond$licenseno[1, 1],
    VarCorr(m2)$cond$licenseno[1, 1],
    VarCorr(m3)$cond$licenseno[1, 1]), 4))
re_variances$pct_reduction = round((1 - re_variances$RE_var / re_variances$RE_var[1]) * 100, 1)

cat("\n=== Random effect variance across models ===\n")
print(re_variances)

# ==============================================================================
# 16. COEFFICIENT PLOT (redesigned)
# ==============================================================================

coef_plot_data = data.frame(
  term = names(or_adjusted),
  OR = or_adjusted,
  CI_lower = ci_adjusted[, 1],
  CI_upper = ci_adjusted[, 2]
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = case_match(term,
                      "income_scaled" ~ "Median Income",
                      "poverty_scaled" ~ "Poverty Rate",
                      "pct_white_scaled" ~ "% White",
                      "license_age_scaled" ~ "License Age",
                      "insp_freq_scaled" ~ "Inspection Frequency",
                      "time_scaled" ~ "Time",
                      "business_typeRetail Food" ~ "Retail Food",
                      "business_typeFood Truck/Mobile" ~ "Food Truck / Mobile",
                      "business_typeManufacturing/Wholesale" ~ "Manufacturing / Wholesale",
                      "business_typeOther" ~ "Other Business",
                      .default = term),
    # group variables for visual separation
    category = case_when(
      term %in% c("Median Income", "Poverty Rate", "% White") ~ "Neighbourhood SES",
      term %in% c("License Age", "Inspection Frequency", "Time") ~ "Establishment & Time",
      TRUE ~ "Business Type"),
    category = factor(category, levels = c("Neighbourhood SES", "Establishment & Time", "Business Type")),
    # colour-code by whether CI crosses 1
    significant = CI_lower > 1 | CI_upper < 1,
    direction = case_when(
      !significant ~ "Not significant",
      OR > 1 ~ "Higher odds of passing",
      TRUE ~ "Lower odds of passing"))

ggplot(coef_plot_data, aes(x = OR, y = reorder(term, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = direction),
                 height = 0, linewidth = 0.9) +
  geom_point(aes(color = direction, fill = direction),
             size = 2.2, shape = 21, stroke = 0.8) +
  geom_text(aes(label = sprintf("%.2f", OR)),
            hjust = -0.5, vjust = 0.5, size = 3, color = "grey30", fontface = "bold") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.5, 0.7, 1, 1.5, 2),
                labels = c("0.5", "0.7", "1", "1.5", "2"),
                expand = expansion(mult = c(0.08, 0.18))) +
  scale_color_manual(values = c(
    "Higher odds of passing" = "#2a7b9b",
    "Lower odds of passing" = "#b90d0d",
    "Not significant" = "grey60")) +
  scale_fill_manual(values = c(
    "Higher odds of passing" = "#2a7b9b",
    "Lower odds of passing" = "#b90d0d",
    "Not significant" = "grey80")) +
  labs(
    title = "What Predicts Passing a Food Inspection?",
    subtitle = "Odds ratios from mixed-effects logistic regression (95% CI)",
    x = "Odds Ratio (log scale)",
    y = NULL,
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates (2019)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10, margin = margin(b = 8)),
    plot.caption = element_text(color = "grey50", size = 8, hjust = 0),
    strip.text.y = element_text(angle = 0, face = "bold", size = 10,
                                color = "grey30", hjust = 0),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.placement = "outside",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.spacing.y = unit(0.4, "lines"),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 9, color = "grey40", margin = margin(t = 6)),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.margin = margin(b = -4),
    legend.key.size = unit(0.4, "cm"))

ggsave("plot7_coefficient_plot.png", width = 8, height = 5.5, dpi = 300)

# ==============================================================================
# 17. PREDICTED PROBABILITY PLOT
# ==============================================================================

pred_grid = expand.grid(
  income_scaled = c(
    quantile(model_data$income_scaled, 0.17),
    quantile(model_data$income_scaled, 0.50),
    quantile(model_data$income_scaled, 0.83)),
  time_scaled = seq(min(model_data$time_scaled), max(model_data$time_scaled), length.out = 50),
  poverty_scaled = 0, pct_white_scaled = 0, license_age_scaled = 0, insp_freq_scaled = 0,
  business_type = "Food Service"
) %>%
  mutate(business_type = factor(business_type, levels = levels(model_data$business_type)),
         income_group = case_when(
           income_scaled == quantile(model_data$income_scaled, 0.17) ~ "Low Income",
           income_scaled == quantile(model_data$income_scaled, 0.50) ~ "Middle Income",
           TRUE ~ "High Income"),
         income_group = factor(income_group, levels = c("Low Income", "Middle Income", "High Income")),
         year_approx = min(model_data$year) + (time_scaled - min(model_data$time_scaled)) /
           (max(model_data$time_scaled) - min(model_data$time_scaled)) * (max(model_data$year) - min(model_data$year)))

pred_grid$pred_prob = predict(m3, newdata = pred_grid, type = "response", re.form = NA)

ggplot(pred_grid, aes(x = year_approx, y = pred_prob, color = income_group)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_color_manual(values = colour_palette, name = "Neighbourhood\nIncome Level") +
  labs(
    title = "Predicted Probability of Passing Inspection by Income Level Over Time",
    subtitle = "Mixed-effects logistic regression with income x time interaction\nOther variables held at mean; Food Service reference category\nLines are nearly parallel: the income effect is small and essentially stable",
    x = "Year",
    y = "Predicted Probability of Passing",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates (2019)"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot8_predicted_probs_income_time.png", width = 9, height = 6, dpi = 300)

# ==============================================================================
# 18. MARGINAL EFFECTS TABLE
# ==============================================================================

'
Concrete predicted probabilities at representative values.
This is the clearest way to communicate the income effect:
  Low income tract (~$41K):  50.1% pass probability
  High income tract (~$128K): 53.5% pass probability
  Gap: ~3.4 percentage points

For context, the difference between Food Service and Retail Food fail rates
(47.7% vs 40.1%) is much larger than the income effect.
'

marginal_grid = expand.grid(
  income_scaled = c(
    quantile(model_data$income_scaled, 0.17),
    quantile(model_data$income_scaled, 0.50),
    quantile(model_data$income_scaled, 0.83)),
  poverty_scaled = 0, pct_white_scaled = 0,
  license_age_scaled = 0, insp_freq_scaled = 0,
  time_scaled = 0,
  business_type = "Food Service"
) %>%
  mutate(business_type = factor(business_type, levels = levels(model_data$business_type)),
         income_label = c("Low Income (17th pctile)", "Middle Income (median)", "High Income (83rd pctile)"))

marginal_grid$pred_pass_prob = predict(m2, newdata = marginal_grid, type = "response", re.form = NA)
marginal_grid$approx_income = marginal_grid$income_scaled * scaling_params$income_sd + scaling_params$income_mean

marginal_table = marginal_grid %>%
  select(income_label, approx_income, pred_pass_prob) %>%
  mutate(
    approx_income = paste0("$", scales::comma(round(approx_income))),
    pred_pass_prob = paste0(round(pred_pass_prob * 100, 1), "%"))

cat("\n=== Predicted pass probability by income level ===\n")
print(marginal_table)

# ==============================================================================
# 19. RANDOM EFFECTS DISTRIBUTION
# ==============================================================================

'
The random intercepts are tightly clustered around zero (SD = 0.14 after
adjustment). This confirms that once we account for observable characteristics,
establishments do not differ much from each other in baseline pass probability.
The raw ICC of 4.6% was already small, and the fixed effects explain most of it.
'

re_estimates = ranef(m2)$cond$licenseno
re_df = data.frame(random_intercept = re_estimates[, 1])

ggplot(re_df, aes(x = random_intercept)) +
  geom_histogram(bins = 60, fill = red, alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  labs(
    title = "Distribution of Establishment-Level Random Intercepts",
    subtitle = paste0("From fully adjusted model | SD = ", round(sd(re_df$random_intercept), 3),
                      "\nTightly clustered: observable covariates explain most between-establishment variation"),
    x = "Random Intercept (log-odds scale)",
    y = "Number of Establishments",
    caption = "Source: Boston Food Inspection data"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot9_random_effects_distribution.png", width = 9, height = 5, dpi = 300)

# ==============================================================================
# 20. MODEL DIAGNOSTICS
# ==============================================================================

'
DHARMa simulation-based residual diagnostics.
RESULTS:
  KS uniformity test: p = 0.054 (borderline, but expected at n = 72k;
  tiny deviations from uniformity are normal with large samples).
  Overdispersion test: dispersion = 1.0003, p = 0.80.
  No evidence of misspecification. The model fits well.
'

dharma_res = simulateResiduals(m2, n = 250, plot = FALSE)

png("plot10_dharma_diagnostics.png", width = 10, height = 5, units = "in", res = 300)
plot(dharma_res, main = "DHARMa Residual Diagnostics")
dev.off()

# ==============================================================================
# 21. TEST SET EVALUATION
# ==============================================================================

# Scale test data using stored training parameters
test_model_data = test_data %>%
  filter(!is.na(median_income), !is.na(poverty_rate), !is.na(pct_white), !is.na(license_age),
         !is.na(insp_per_year), !is.na(business_type), !is.na(time_scaled)) %>%
  mutate(
    income_scaled = (median_income - scaling_params$income_mean) / scaling_params$income_sd,
    poverty_scaled = (poverty_rate - scaling_params$poverty_mean) / scaling_params$poverty_sd,
    pct_white_scaled = (pct_white - scaling_params$pct_white_mean) / scaling_params$pct_white_sd,
    license_age_scaled = (license_age - scaling_params$license_age_mean) / scaling_params$license_age_sd,
    insp_freq_scaled = (insp_per_year - scaling_params$insp_freq_mean) / scaling_params$insp_freq_sd)

test_model_data$pred_prob = predict(m2, newdata = test_model_data, type = "response", re.form = NA)
roc_obj = roc(test_model_data$pass_fail, test_model_data$pred_prob)

# Null model AUC for comparison
test_model_data$pred_prob_null = predict(m0, newdata = test_model_data, type = "response", re.form = NA)
roc_null = roc(test_model_data$pass_fail, test_model_data$pred_prob_null)

'
RESULTS:
  Full model AUC:  0.624
  Null model AUC:  0.500
  AUC gain:        0.124

AUC of 0.624 is modest but expected for a ~50/50 outcome with weak
neighbourhood-level predictors. The 0.124 AUC gain shows fixed effects
do add discriminative power vs. a constant prediction, but the practical
ceiling is low. The model value is in identifying population-level
associations, not predicting individual inspection outcomes.
'

cat("\n=== Test set evaluation ===\n")
cat("Full model AUC:", round(auc(roc_obj), 3), "\n")
cat("Null model AUC:", round(auc(roc_null), 3), "\n")
cat("AUC gain from fixed effects:", round(auc(roc_obj) - auc(roc_null), 3), "\n")

# ==============================================================================
# 22. MODEL COMPARISON
# ==============================================================================

'
INTERPRETATION:
Poverty alone barely improves on null (Delta AIC = -2.4, marginal).
Income alone is clearly better (Delta AIC = -79).
SES composite is best single predictor (Delta AIC = -92).
Full model is a massive improvement (Delta AIC = -2730), driven primarily
by inspection frequency and time, not by demographics.
'

model_comparison = data.frame(
  Model = c("Null (intercept only)", "Income only", "Poverty only", "SES composite index",
            "Fully adjusted", "Income x time interaction"),
  AIC = round(c(AIC(m0), AIC(m1a), AIC(m1b), AIC(m1c), AIC(m2), AIC(m3)), 1),
  Delta_AIC = round(c(AIC(m0), AIC(m1a), AIC(m1b), AIC(m1c), AIC(m2), AIC(m3)) - AIC(m0), 1),
  logLik = round(c(logLik(m0), logLik(m1a), logLik(m1b), logLik(m1c), logLik(m2), logLik(m3)), 1))

cat("\n=== Model comparison ===\n")
print(model_comparison)

# ==============================================================================
# 23. SENSITIVITY ANALYSIS: EXCLUDING COVID YEARS
# ==============================================================================

# RESULT: All ORs shift by < 0.02 when excluding 2020-2021.
# Main findings are not driven by pandemic-era disruptions.

model_data_nocovid = model_data %>%
  filter(!year %in% c(2020, 2021))

m2_nocovid = glmmTMB(pass_fail ~ income_scaled + poverty_scaled + pct_white_scaled +
                       license_age_scaled + insp_freq_scaled + business_type + time_scaled + (1 | licenseno),
                     data = model_data_nocovid, family = binomial)

or_nocovid = exp(fixef(m2_nocovid)$cond)

sensitivity_compare = data.frame(
  term = names(or_adjusted),
  OR_full = round(or_adjusted, 3),
  OR_no_covid = round(or_nocovid, 3)
) %>%
  filter(term != "(Intercept)") %>%
  mutate(difference = round(OR_no_covid - OR_full, 4))

cat("\n=== COVID sensitivity check ===\n")
print(sensitivity_compare)

# ==============================================================================
# 24. SUMMARY OF KEY FINDINGS
# ==============================================================================
'
SUMMARY OF KEY FINDINGS

1. WITHIN-ESTABLISHMENT VARIATION DOMINATES
   ICC = 0.046 --> only 4.6% of pass/fail variation is between establishments.
   Most variation occurs within the same establishment over time: the same
   restaurant sometimes passes, sometimes fails. Pass/fail is a noisy outcome.
   The modest between-establishment variance that exists is almost entirely
   explained by observable characteristics (87.7% RE variance reduction),
   especially inspection frequency and time trends.

2. NEIGHBOURHOOD INCOME HAS A SIGNIFICANT BUT SMALL EFFECT
   Adjusted OR = 1.064 (6.4% increase in odds of passing per 1 SD income).
   In concrete terms: a Food Service establishment in a low-income tract
   (~$41K median income) has a 50.1% predicted pass probability, vs 53.5%
   in a high-income tract (~$128K). The 3.4 pp gap is real but modest -->
   much smaller than the difference between business types (~8 pp between
   Food Service and Retail Food).

3. THE INCOME GRADIENT IS ESSENTIALLY STABLE OVER TIME
   Income x time interaction OR = 1.022, p = 0.006. Statistically significant
   but practically negligible: a ~2% shift in the income OR per SD of time.
   The predicted probability lines are visually parallel. The gap is widening
   very slightly, but the magnitude is too small to be substantively meaningful.

4. COLLINEARITY COMPLICATES THE FULL MODEL
   The poverty coefficient flips sign when income and % White are included
   (unadjusted OR = 0.979, adjusted OR = 1.114). This is a suppression
   effect, not evidence that poverty improves outcomes. The PCA composite
   SES index (OR = 1.065) provides a cleaner single-predictor estimate
   of the neighbourhood SES effect.

5. INSPECTION FREQUENCY IS THE DOMINANT PREDICTOR (OR = 0.665)
   But this likely reflects enforcement mechanics rather than food safety:
   failures trigger re-inspections, which simultaneously inflate frequency
   and fail counts, creating an artificial negative correlation. This
   variable should be interpreted as a control for enforcement patterns,
   not as a causal predictor of food safety quality.

6. INSPECTION TARGETING SHOWS MINIMAL INCOME BIAS
   Mean inspections per year: Low income = 3.64, High income = 3.45.
   The difference is negligible, arguing against systematic targeting.

7. COVID SENSITIVITY: ROBUST
   All ORs shift by < 0.02 when excluding 2020-2021.

8. MODEL DIAGNOSTICS: CLEAN
   No overdispersion (1.0003). KS test borderline (p = 0.054) but expected
   at n = 72k. AUC = 0.624 on held-out test set --> modest but appropriate
   for a model identifying population-level associations rather than
   predicting individual outcomes.
'


