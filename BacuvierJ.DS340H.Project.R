library(tidyverse)
library(lubridate)
library(sf)
library(tigris)
library(leaflet)
library(tidycensus)
library(glmmTMB)
library(ggcorrplot)
library(pROC)

options(tigris_use_cache = TRUE)

# ==============================================================================
# 1. LOADING DATA
# ==============================================================================

raw_data = read_csv('/Users/juliettebacuvier/desktop/tmpp7jeda9f.csv')
names(raw_data)

# ==============================================================================
# 2. DATE PARSING & FILTERING
# ==============================================================================

raw_data = raw_data %>%
  mutate(insp_date = ymd_hms(resultdttm), year = year(insp_date)) %>%
  filter(year >= 2015, year <= 2025, !is.na(insp_date))

# ==============================================================================
# 3. OUTCOME VARIABLE
# ==============================================================================

# Check both candidate outcome columns before deciding which to use
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
# 4. GEOGRAPHY PARSING
# ==============================================================================

raw_data = raw_data %>%
  mutate(
    location = str_remove_all(location, "[()]"),
    lat = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 1])),
    lon = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 2])))

# ==============================================================================
# 5. ACS CENSUS DATA
# ==============================================================================

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
    poverty_rate = ifelse(total_poverty == 0, NA_real_, below_poverty/ total_poverty),
    pct_white = ifelse(total_race == 0, NA_real_, white_alone / total_race),
    pct_black = ifelse(total_race == 0, NA_real_, black_alone / total_race),
    pct_indigenous = ifelse(total_race == 0, NA_real_, indigenous_alone / total_race),
    pct_asian = ifelse(total_race == 0, NA_real_, asian_alone / total_race),
    pct_hawaiian = ifelse(total_race == 0, NA_real_, hawaiian_alone / total_race),
    pct_other_race = ifelse(total_race == 0, NA_real_, other_race_alone / total_race),
    pct_hispanic = ifelse(total_hispanic == 0, NA_real_, hispanic / total_hispanic),
    pct_foreign_born = ifelse(total_nativity == 0, NA_real_, foreign_born / total_nativity))

cat("\nACS tract data:", nrow(acs_clean), "tracts\n")

# ==============================================================================
# 6. SPATIAL JOIN
# ==============================================================================

ma_tracts = tracts(state = "MA", year = 2019, cb = TRUE) %>%
  st_transform(crs = 4326)

analysis_data = raw_data %>%
  filter(!is.na(pass_fail), !is.na(lat), !is.na(lon))

inspections_sf = analysis_data %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

inspections_sf = st_join(inspections_sf, ma_tracts %>% select(GEOID), join = st_within)

analysis_data$tract_geoid = inspections_sf$GEOID

cat("\nTract spatial join match rate:", round(mean(!is.na(analysis_data$tract_geoid)) * 100, 1), "%\n")

analysis_data = analysis_data %>%
  left_join(acs_clean, by = c("tract_geoid" = "GEOID"))

# ==============================================================================
# 7. DEDUPLICATION
# ==============================================================================

# Raw data is one row per violation per visit
# We want one row per visit, keeping the worst outcome (fail = 0)
analysis_data = analysis_data %>%
  group_by(licenseno, insp_date) %>%
  arrange(pass_fail) %>%
  slice(1) %>%
  ungroup()

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
cat("\n=== Top licensecat values ===\n")
raw_cats = sort(table(analysis_data$licensecat), decreasing = TRUE)
print(head(raw_cats, 30))

analysis_data = analysis_data %>%
  mutate(business_type = case_when(
      licensecat == "FS"  ~ "Food Service",
      licensecat == "FT"  ~ "Food Truck/Mobile",
      licensecat == "RF"  ~ "Retail Food",
      licensecat == "MFW" ~ "Manufacturing/Wholesale",
      TRUE ~ "Other"),
    business_type = factor(business_type, levels = c("Food Service", "Retail Food",
                            "Food Truck/Mobile", "Manufacturing/Wholesale", "Other")))

# if >15% fall into Other, the taxonomy needs more work
cat("\n=== Business type distribution ===\n")
analysis_data %>%
  count(business_type, sort = TRUE) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  print()

# continuous time variable for regression
# Using scaled numeric date rather than year dummies -- no reason to expect jumps at year boundaries
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
cat("ACS match rate:", round(mean(!analysis_data$acs_missing) * 100, 1), "%\n")
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

cat("\n=== Fail rate range across tracts (min 10 inspections) ===\n")
tract_fail = analysis_data %>%
  filter(!is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(n = n(),,fail_rate = mean(1 - pass_fail, na.rm = TRUE), .groups = "drop") %>%
  filter(n >= 10)

cat("Tracts:", nrow(tract_fail), "\n")
cat("Fail rate range:", round(min(tract_fail$fail_rate) * 100, 1), "% to", round(max(tract_fail$fail_rate) * 100, 1), "%\n")
cat("Fail rate SD:", round(sd(tract_fail$fail_rate) * 100, 1), "pp\n")

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

# --- Choropleth map ---
neighbourhood_summary = analysis_data %>%
  filter(!is.na(tract_geoid)) %>%
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
    subtitle = "Each point is a census tract (min. 10 inspections); shaded band = 95% CI",
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

# dropping income_tercile from analysis_data first to avoid .x/.y clash on join
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
    subtitle = "Census tracts grouped into income terciles (ACS 2019 median household income)",
    x = "Year",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot2_failrate_time_trends_tract.png", width = 9, height = 6, dpi = 300)

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
    subtitle = "Each point is a census tract (min. 10 inspections)",
    x = "Poverty Rate (ACS 2019, tract level)",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot5_failrate_vs_poverty_tract.png", width = 9, height = 6, dpi = 300)

# --- Inspection frequency by income tercile ---
# checks whether lower-income areas are inspected more often (targeting bias)
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
    subtitle = "Checks for inspection targeting bias: are lower-income areas inspected more?",
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
cat("Training establishments:", n_distinct(train_data$licenseno), "\n")
cat("Training inspections:   ", nrow(train_data), "\n")
cat("Test establishments:    ", n_distinct(test_data$licenseno), "\n")
cat("Test inspections:       ", nrow(test_data), "\n")

# ==============================================================================
# 12. MODEL DATA PREPARATION
# ==============================================================================

# Keeps only establishments with 2+ inspections (helps model convergence)
# and drops any remaining NAs on model variables
model_data = train_data %>%
  group_by(licenseno) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  filter(!is.na(median_income), !is.na(poverty_rate), !is.na(pct_white), !is.na(license_age),
    !is.na(insp_per_year), !is.na(business_type), !is.na(time_scaled)
    ) %>%
  mutate(
    # scaled to mean 0, SD 1 so coefficients are comparable
    income_scaled = scale(median_income)[, 1],
    poverty_scaled = scale(poverty_rate)[, 1],
    pct_white_scaled = scale(pct_white)[, 1],
    license_age_scaled = scale(license_age)[, 1],
    insp_freq_scaled = scale(insp_per_year)[, 1] )

cat("\nModel dataset:", nrow(model_data), "inspections,", n_distinct(model_data$licenseno), "establishments\n")

# ==============================================================================
# 13. MIXED-EFFECTS LOGISTIC REGRESSION (using glmmTMB)
# ==============================================================================

# --- Model 0: Null model ---
# Random intercept only, no predictors
# Tells us how much variation is between vs. within establishments
cat("\nFitting null model...\n")
m0 = glmmTMB(pass_fail ~ 1 + (1 | licenseno), data = model_data, family = binomial)
summary(m0)

# ICC = between-establishment variance / total variance
re_var = VarCorr(m0)$cond$licenseno[1, 1]
icc = re_var / (re_var + (pi^2 / 3))
cat("Null model ICC:", round(icc, 3), "\n")
cat(round(icc * 100, 1), "% of outcome variation is between establishments\n")

# --- Model 1: Income only ---
cat("\nFitting income-only model...\n")
m1 = glmmTMB(pass_fail ~ income_scaled + (1 | licenseno), data = model_data, family = binomial)
summary(m1)

income_or_unadjusted = exp(fixef(m1)$cond["income_scaled"])
cat("\nUnadjusted OR for income (per 1 SD increase):", round(income_or_unadjusted, 3), "\n")

# --- Model 2: Fully adjusted ---
cat("\nFitting fully adjusted model...\n")
m2 = glmmTMB(pass_fail ~ income_scaled + poverty_scaled + pct_white_scaled + license_age_scaled +
    insp_freq_scaled + business_type + time_scaled + (1 | licenseno),
  data = model_data, family = binomial)
summary(m2)

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

# How much does the income OR shrink after adding controls?
income_or_adjusted = exp(fixef(m2)$cond["income_scaled"])
cat("\nUnadjusted income OR:", round(income_or_unadjusted, 3))
cat("\nAdjusted income OR:  ", round(income_or_adjusted, 3))
cat("\nAttenuation:         ", round((1 - income_or_adjusted / income_or_unadjusted) * 100, 1), "%\n")

# --- Model 3: Income x time interaction ---
# Tests whether the income gradient is changing over time
cat("\nFitting interaction model...\n")
m3 = glmmTMB(pass_fail ~ income_scaled * time_scaled + poverty_scaled + pct_white_scaled +
    license_age_scaled + insp_freq_scaled + business_type + (1 | licenseno),
  data = model_data, family = binomial)
summary(m3)

interaction_coef = fixef(m3)$cond["income_scaled:time_scaled"]
cat("\nIncome x time interaction OR:", round(exp(interaction_coef), 4), "\n")

# ==============================================================================
# 14. COEFFICIENT PLOT
# ==============================================================================

coef_plot_data = data.frame(
  term = names(or_adjusted),
  OR = or_adjusted,
  CI_lower = ci_adjusted[, 1],
  CI_upper = ci_adjusted[, 2]
) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = recode(term,
                  "income_scaled" = "Median Income (scaled)",
                  "poverty_scaled" = "Poverty Rate (scaled)",
                  "pct_white_scaled" = "% White (scaled)",
                  "license_age_scaled" = "License Age (scaled)",
                  "insp_freq_scaled" = "Inspection Frequency (scaled)",
                  "time_scaled" = "Time (scaled)",
                  "business_typeRetail Food" = "Business: Retail Food",
                  "business_typeFood Truck/Mobile" = "Business: Food Truck/Mobile",
                  "business_typeManufacturing/Wholesale" = "Business: Manufacturing/Wholesale",
                  "business_typeOther" = "Business: Other"))

ggplot(coef_plot_data, aes(x = OR, y = reorder(term, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.25, color = "grey40") +
  geom_point(size = 3, color = red) +
  scale_x_log10() +
  labs(
    title = "Odds Ratios from Mixed-Effects Logistic Regression",
    subtitle = "Outcome: passing a food safety inspection | Random intercept by establishment\nOR > 1 = higher odds of passing; OR < 1 = lower odds of passing",
    x = "Odds Ratio (log scale, 95% CI)",
    y = NULL,
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates (2019)"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot7_coefficient_plot.png", width = 9, height = 6, dpi = 300)

# ==============================================================================
# 15. PREDICTED PROBABILITY PLOT
# ==============================================================================

# Prediction grid across income levels & time
# All other variables held at their mean (= 0 after scaling)
pred_grid = expand.grid(
  income_scaled = c(
    quantile(model_data$income_scaled, 0.17), # low income tercile
    quantile(model_data$income_scaled, 0.50), # middle
    quantile(model_data$income_scaled, 0.83) # high income tercile
    ),
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
    # back-transformed scaled time to approximate year for the x axis
    year_approx = min(model_data$year) + (time_scaled - min(model_data$time_scaled)) /
      (max(model_data$time_scaled) - min(model_data$time_scaled)) * (max(model_data$year) - min(model_data$year)))

# re.form = NA gives population-average predictions (ignores random effects)
pred_grid$pred_prob = predict(m3, newdata = pred_grid, type = "response", re.form = NA)

ggplot(pred_grid, aes(x = year_approx, y = pred_prob, color = income_group)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_color_manual(values = colour_palette, name = "Neighbourhood\nIncome Level") +
  labs(
    title = "Predicted Probability of Passing Inspection by Income Level Over Time",
    subtitle = "Mixed-effects logistic regression with income × time interaction\nOther variables held at mean; Food Service reference category",
    x = "Year",
    y = "Predicted Probability of Passing",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates (2019)"
  ) + theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))

ggsave("plot8_predicted_probs_income_time.png", width = 9, height = 6, dpi = 300)

# ==============================================================================
# 16. TEST SET EVALUATION
# ==============================================================================

# Scale test data using training means/SDs to avoid leakage
test_model_data = test_data %>%
  filter(!is.na(median_income), !is.na(poverty_rate), !is.na(pct_white), !is.na(license_age), 
         !is.na(insp_per_year), !is.na(business_type), !is.na(time_scaled)
  ) %>%
  mutate(
    income_scaled = (median_income - mean(model_data$median_income, na.rm = TRUE)) /
      sd(model_data$median_income, na.rm = TRUE),
    poverty_scaled = (poverty_rate - mean(model_data$poverty_rate, na.rm = TRUE)) /
      sd(model_data$poverty_rate, na.rm = TRUE),
    pct_white_scaled = (pct_white - mean(model_data$pct_white, na.rm = TRUE)) /
      sd(model_data$pct_white, na.rm = TRUE),
    license_age_scaled = (license_age - mean(model_data$license_age, na.rm = TRUE)) /
      sd(model_data$license_age, na.rm = TRUE),
    insp_freq_scaled = (insp_per_year - mean(model_data$insp_per_year, na.rm = TRUE)) /
      sd(model_data$insp_per_year, na.rm = TRUE))

test_model_data$pred_prob = predict(m2, newdata = test_model_data, type = "response", re.form = NA)

roc_obj = roc(test_model_data$pass_fail, test_model_data$pred_prob)
cat("\n=== Test set evaluation ===\n")
cat("AUC:", round(auc(roc_obj), 3), "\n")
cat("Mean predicted pass probability:", round(mean(test_model_data$pred_prob, na.rm = TRUE), 3), "\n")
cat("Observed pass rate:             ", round(mean(test_model_data$pass_fail, na.rm = TRUE), 3), "\n")

# ==============================================================================
# 17. MODEL COMPARISON
# ==============================================================================

model_comparison = data.frame(
  Model = c("Null (intercept only)", "Income only (unadjusted)", "Fully adjusted", "Income x time interaction"),
  AIC = round(c(AIC(m0), AIC(m1), AIC(m2), AIC(m3)), 1),
  Delta_AIC = round(c(AIC(m0), AIC(m1), AIC(m2), AIC(m3)) - AIC(m0), 1),
  logLik = round(c(logLik(m0), logLik(m1), logLik(m2), logLik(m3)), 1))

cat("\n=== Model comparison ===\n")
cat("Delta AIC > 2 = meaningful improvement; Delta AIC > 10 = strong\n")
print(model_comparison)

# LRT: does income improve on null?
cat("\n=== LRT: Null vs Income only ===\n")
lrt_m0_m1 = anova(m0, m1)
print(lrt_m0_m1)

# LRT: does full adjustment improve on income only?
cat("\n=== LRT: Income only vs Fully adjusted ===\n")
lrt_m1_m2 = anova(m1, m2)
print(lrt_m1_m2)

# LRT: does the interaction term improve fit?
# Key test -- is the income gradient actually changing over time?
cat("\n=== LRT: Fully adjusted vs Income x time ===\n")
lrt_m2_m3 = anova(m2, m3)
print(lrt_m2_m3)

# Random effect variance should normally shrink as fixed effects absorb between-unit variation
cat("\n=== Random effect variance across models ===\n")
re_variances = data.frame(
  Model = c("Null", "Income only", "Fully adjusted", "Income x time"),
  RE_var = round(c(VarCorr(m0)$cond$licenseno[1, 1], VarCorr(m1)$cond$licenseno[1, 1],
    VarCorr(m2)$cond$licenseno[1, 1], VarCorr(m3)$cond$licenseno[1, 1]), 4))
print(re_variances)




