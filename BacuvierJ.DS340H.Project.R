library(tidyverse)
library(lubridate)
library(sf)
library(tigris)
library(tidycensus)
library(glmmTMB)
library(ggcorrplot)
library(pROC)
library(car)
library(tmap)

# Resolve namespace conflicts
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

# pass = 1, fail = 0
raw_data = raw_data %>%
  mutate(pass_fail = case_when(
    str_detect(result, regex("pass", ignore_case = TRUE)) ~ 1L,
    str_detect(result, regex("fail", ignore_case = TRUE)) ~ 0L,
    TRUE ~ NA_integer_))

# ==============================================================================
# 4. DEDUPLICATION
# ==============================================================================

# Raw data is one row per violation per visit.
# Collapse to one row per visit, keeping the worst outcome (fail = 0).
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

acs_data = get_acs(
  geography = "tract", state = "MA",
  variables = c(
    median_income = "B19013_001",
    total_pop = "B01003_001",
    below_poverty = "B17001_002",
    total_poverty = "B17001_001",
    white_alone = "B02001_002",
    total_race = "B02001_001"),
  year = 2019, survey = "acs5")

acs_clean = acs_data %>%
  select(GEOID, variable, estimate) %>%
  pivot_wider(names_from = variable, values_from = estimate) %>%
  mutate(
    poverty_rate = ifelse(total_poverty == 0, NA_real_, below_poverty / total_poverty),
    pct_white = ifelse(total_race == 0, NA_real_, white_alone / total_race),
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

# Inspection frequency per establishment per year
insp_freq = analysis_data %>%
  group_by(licenseno, year) %>%
  summarise(insp_per_year = n(), .groups = "drop")

analysis_data = analysis_data %>%
  left_join(insp_freq, by = c("licenseno", "year"))

# License age in years at time of inspection
analysis_data = analysis_data %>%
  mutate(iss_date = ymd_hms(issdttm),
         license_age = as.numeric(difftime(insp_date, iss_date, units = "days")) / 365.25,
         license_age = pmax(license_age, 0, na.rm = TRUE))

# Business type
analysis_data = analysis_data %>%
  mutate(business_type = case_when(
    licensecat == "FS" ~ "Food Service",
    licensecat == "FT" ~ "Food Truck/Mobile",
    licensecat == "RF" ~ "Retail Food",
    licensecat == "MFW" ~ "Manufacturing/Wholesale",
    TRUE ~ "Other"),
    business_type = factor(business_type, levels = c("Food Service", "Retail Food",
                                                     "Food Truck/Mobile", "Manufacturing/Wholesale", "Other")))

# Continuous time variable for regression
analysis_data = analysis_data %>%
  mutate(time_scaled = scale(as.numeric(insp_date))[, 1])

# Income tercile grouping
analysis_data = analysis_data %>%
  mutate(
    income_tercile = ntile(median_income, 3),
    income_tercile = factor(income_tercile, labels = c("Low Income", "Middle Income", "High Income")),
    year = as.integer(year))

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
                     pct_white, license_age, insp_per_year, business_type),
                   ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(analysis_data) * 100, 1)) %>%
  arrange(desc(n_missing)) %>%
  print(n = 20)

# ==============================================================================
# 10. EDA VISUALISATIONS
# ==============================================================================

red = "#b90d0d"
blue = "#2a7b9b"
gold = "#f0a500"
white = "white"
colour_palette = c(red, gold, blue)
heat_palette = c(blue, white, red)

# --- Choropleth maps ---
neighbourhood_summary = analysis_data %>%
  filter(!is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(
    n_inspections = n(),
    pass_rate = mean(pass_fail, na.rm = TRUE),
    median_income = first(median_income),
    .groups = "drop")

plot_data = neighbourhood_summary %>% filter(n_inspections >= 10)

map_tracts = ma_tracts %>%
  filter(GEOID %in% plot_data$tract_geoid) %>%
  left_join(plot_data, by = c("GEOID" = "tract_geoid"))

# Boston Neighbourhood boundaries from Analyze Boston
boston_neigh = st_read(
  "https://bostonopendata-boston.opendata.arcgis.com/datasets/boston::boston-neighborhoods.geojson",
  quiet = TRUE
) %>% st_transform(4326)

tmap_mode("plot")

# --- Pass Rate Map ---
static_map_pass = tm_shape(map_tracts) +
  tm_fill("pass_rate", palette = c(red, gold, blue), style = "cont",
          title = "Pass Rate",
          legend.format = list(fun = function(x) paste0(round(x * 100), "%"))) +
  tm_borders(col = "white", lwd = 0.2) +
  tm_shape(boston_neigh) +
  tm_borders(col = "black", lwd = 1.5) +
  tm_layout(
    legend.position = c("left", "bottom"),
    legend.bg.color = "white",
    legend.bg.alpha = 0.85,
    legend.frame = FALSE,
    frame = FALSE,
    inner.margins = c(0.02, 0.02, 0.05, 0.02))

static_map_pass
#tmap_save(static_map_pass,
          "/Users/juliettebacuvier/Desktop/plot_choropleth_pass_static.png",
          width = 8, height = 9, dpi = 300)

# --- Median Income Map ---
static_map_income = tm_shape(map_tracts) +
  tm_fill("median_income", palette = "brewer.yl_gn", style = "cont",
          title = "Median Income",
          legend.format = list(fun = function(x) scales::dollar(x))) +
  tm_borders(col = "white", lwd = 0.2) +
  tm_shape(boston_neigh) +
  tm_borders(col = "black", lwd = 1.5) +
  tm_layout(
    legend.position = c("right"),
    legend.bg.color = "white",
    legend.bg.alpha = 0.85,
    legend.frame = FALSE,
    frame = FALSE,
    inner.margins = c(0.02, 0.02, 0.05, 0.02))

static_map_income
#tmap_save(static_map_income,
          "/Users/juliettebacuvier/Desktop/plot_choropleth_income_static.png",
          width = 8, height = 9, dpi = 300)

# ==============================================================================
# 11. CROSS-VALIDATION SPLIT
# Holds back 20% of establishments before modelling.
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

# Keep only establishments with 2+ inspections.
# Establishments with a single inspection cannot contribute to the
# within-establishment variance estimate needed for the random intercept.
model_data = train_data %>%
  group_by(licenseno) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  filter(!is.na(median_income), !is.na(poverty_rate), !is.na(pct_white),
         !is.na(license_age), !is.na(insp_per_year), !is.na(business_type),
         !is.na(time_scaled))

# Store scaling parameters explicitly so train and test use the same values
scaling_params = list(
  income_mean = mean(model_data$median_income, na.rm = TRUE),
  income_sd = sd(model_data$median_income, na.rm = TRUE),
  poverty_mean = mean(model_data$poverty_rate, na.rm = TRUE),
  poverty_sd = sd(model_data$poverty_rate, na.rm = TRUE),
  pct_white_mean = mean(model_data$pct_white, na.rm = TRUE),
  pct_white_sd = sd(model_data$pct_white, na.rm = TRUE),
  license_age_mean = mean(model_data$license_age, na.rm = TRUE),
  license_age_sd = sd(model_data$license_age, na.rm = TRUE),
  insp_freq_mean = mean(model_data$insp_per_year, na.rm = TRUE),
  insp_freq_sd = sd(model_data$insp_per_year, na.rm = TRUE))

model_data = model_data %>%
  mutate(
    income_scaled = (median_income - scaling_params$income_mean) / scaling_params$income_sd,
    poverty_scaled = (poverty_rate - scaling_params$poverty_mean) / scaling_params$poverty_sd,
    pct_white_scaled = (pct_white - scaling_params$pct_white_mean) / scaling_params$pct_white_sd,
    license_age_scaled = (license_age - scaling_params$license_age_mean) / scaling_params$license_age_sd,
    insp_freq_scaled = (insp_per_year - scaling_params$insp_freq_mean) / scaling_params$insp_freq_sd)

cat("\nModel dataset:", nrow(model_data), "inspections,", n_distinct(model_data$licenseno), "establishments\n")

# ==============================================================================
# 13. COLLINEARITY DIAGNOSTICS
# ==============================================================================

vif_model = glm(pass_fail ~ income_scaled + poverty_scaled + pct_white_scaled +
                  license_age_scaled + insp_freq_scaled + business_type + time_scaled,
                data = model_data, family = binomial)

vif_values = vif(vif_model)
print(vif_values)

ses_vars = model_data %>%
  select(income_scaled, poverty_scaled, pct_white_scaled) %>%
  drop_na()

ses_pca = prcomp(ses_vars, center = FALSE, scale. = FALSE)
model_data$ses_index = ses_pca$x[, 1]

print(summary(ses_pca))

# ==============================================================================
# 14. MIXED-EFFECTS LOGISTIC REGRESSION (glmmTMB)
# ==============================================================================

# Model 0: Null model (random intercept only).
m0 = glmmTMB(pass_fail ~ 1 + (1 | licenseno), data = model_data, family = binomial)

re_var_m0 = VarCorr(m0)$cond$licenseno[1, 1]
icc = re_var_m0 / (re_var_m0 + (pi^2 / 3))
cat("\nNull model ICC:", round(icc, 3), "(", round(icc * 100, 1), "% between-establishment)\n")

# Model 1a: Income only
m1a = glmmTMB(pass_fail ~ income_scaled + (1 | licenseno),
              data = model_data, family = binomial)

# Model 1b: Poverty only (robustness check)
m1b = glmmTMB(pass_fail ~ poverty_scaled + (1 | licenseno),
              data = model_data, family = binomial)

# Model 2: Fully adjusted.
m2 = glmmTMB(pass_fail ~ income_scaled + poverty_scaled + pct_white_scaled +
               license_age_scaled + insp_freq_scaled + business_type +
               time_scaled + (1 | licenseno),
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

# Model 3: Income x time interaction.
m3 = glmmTMB(pass_fail ~ income_scaled * time_scaled + poverty_scaled +
               pct_white_scaled + license_age_scaled + insp_freq_scaled +
               business_type + (1 | licenseno),
             data = model_data, family = binomial)

# ==============================================================================
# 15. RANDOM EFFECT VARIANCE DECOMPOSITION
# ==============================================================================

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
# 16. COEFFICIENT PLOT
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
    category = case_when(
      term %in% c("Median Income", "Poverty Rate", "% White") ~ "Neighbourhood SES",
      term %in% c("License Age", "Inspection Frequency", "Time") ~ "Establishment & Time",
      TRUE ~ "Business Type"),
    category = factor(category, levels = c("Neighbourhood SES", "Establishment & Time", "Business Type")),
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
# 17. COMBINED PASS RATE OVER TIME (Raw vs Modeled)
# ==============================================================================

# Tract-level income tercile assignment (one income per tract)
tract_income = analysis_data %>%
  filter(!is.na(median_income), !is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(median_income = first(median_income), .groups = "drop") %>%
  mutate(income_tercile = ntile(median_income, 3),
         income_tercile = factor(income_tercile,
                                 labels = c("Low Income", "Middle Income", "High Income")))

time_trends = analysis_data %>%
  filter(!is.na(median_income), !is.na(pass_fail), !is.na(tract_geoid)) %>%
  select(-income_tercile) %>%
  left_join(tract_income %>% select(tract_geoid, income_tercile), by = "tract_geoid") %>%
  group_by(year, income_tercile) %>%
  summarise(pass_rate = mean(pass_fail, na.rm = TRUE), n = n(), .groups = "drop")

pred_grid = expand.grid(
  income_scaled = c(
    quantile(model_data$income_scaled, 0.17),
    quantile(model_data$income_scaled, 0.50),
    quantile(model_data$income_scaled, 0.83)),
  time_scaled = seq(min(model_data$time_scaled), max(model_data$time_scaled), length.out = 50),
  poverty_scaled = 0, pct_white_scaled = 0, license_age_scaled = 0, insp_freq_scaled = 0,
  business_type = "Food Service") %>%
  mutate(business_type = factor(business_type, levels = levels(model_data$business_type)),
         income_tercile = case_when(
           income_scaled == quantile(model_data$income_scaled, 0.17) ~ "Low Income",
           income_scaled == quantile(model_data$income_scaled, 0.50) ~ "Middle Income",
           TRUE ~ "High Income"),
         income_tercile = factor(income_tercile,
                                 levels = c("Low Income", "Middle Income", "High Income")),
         year_approx = min(model_data$year) + (time_scaled - min(model_data$time_scaled)) /
           (max(model_data$time_scaled) - min(model_data$time_scaled)) *
           (max(model_data$year) - min(model_data$year)))

pred_grid$pred_prob = predict(m3, newdata = pred_grid, type = "response", re.form = NA)

df_raw = time_trends %>%
  select(year, pass_rate, income_tercile) %>%
  mutate(data_type = "Raw", legend_key = paste0(income_tercile, " (Raw)"))

df_model = pred_grid %>%
  select(year = year_approx, pass_rate = pred_prob, income_tercile) %>%
  mutate(data_type = "Modeled", legend_key = paste0(income_tercile, " (Modeled)"))

df_combined = bind_rows(df_raw, df_model) %>%
  mutate(legend_key = factor(legend_key, levels = c(
    "High Income (Modeled)", "High Income (Raw)",
    "Middle Income (Modeled)", "Middle Income (Raw)",
    "Low Income (Modeled)", "Low Income (Raw)")))

legend_colors = c(
  "High Income (Modeled)" = blue, "High Income (Raw)" = blue,
  "Middle Income (Modeled)" = gold, "Middle Income (Raw)" = gold,
  "Low Income (Modeled)" = red, "Low Income (Raw)" = red)

legend_linetypes = c(
  "High Income (Modeled)" = "solid", "High Income (Raw)" = "dashed",
  "Middle Income (Modeled)" = "solid", "Middle Income (Raw)" = "dashed",
  "Low Income (Modeled)" = "solid", "Low Income (Raw)" = "dashed")

ggplot(df_combined, aes(x = year, y = pass_rate, color = legend_key, linetype = legend_key)) +
  geom_line(aes(linewidth = data_type, alpha = data_type)) +
  geom_point(data = filter(df_combined, data_type == "Raw"), size = 2.5, alpha = 0.4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_color_manual(values = legend_colors, name = "Income Level & Data Type") +
  scale_linetype_manual(values = legend_linetypes, name = "Income Level & Data Type") +
  scale_linewidth_manual(values = c("Modeled" = 1.5, "Raw" = 0.8), guide = "none") +
  scale_alpha_manual(values = c("Modeled" = 1, "Raw" = 0.4), guide = "none") +
  labs(
    title = "Inspection Pass Rate Over Time by Income Level",
    subtitle = "Solid lines: modelled predicted probability | Faded dashed lines: raw observed pass rate\nModel adjusts for business type, inspection frequency, and license age",
    x = "Year", y = "Pass Rate / Predicted Probability",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.width = unit(1.2, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4))

ggsave("plot8_combined_passrate_time.png", width = 10, height = 6, dpi = 300)

# ==============================================================================
# 18. MARGINAL EFFECTS TABLE
# ==============================================================================

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
         income_label = c("Low Income (17th pctile)", "Middle Income (median)",
                          "High Income (83rd pctile)"))

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
# 19. TEST SET EVALUATION
# ==============================================================================

test_model_data = test_data %>%
  filter(!is.na(median_income), !is.na(poverty_rate), !is.na(pct_white),
         !is.na(license_age), !is.na(insp_per_year), !is.na(business_type),
         !is.na(time_scaled)) %>%
  mutate(
    income_scaled = (median_income - scaling_params$income_mean) / scaling_params$income_sd,
    poverty_scaled = (poverty_rate - scaling_params$poverty_mean) / scaling_params$poverty_sd,
    pct_white_scaled = (pct_white - scaling_params$pct_white_mean) / scaling_params$pct_white_sd,
    license_age_scaled = (license_age - scaling_params$license_age_mean) / scaling_params$license_age_sd,
    insp_freq_scaled = (insp_per_year - scaling_params$insp_freq_mean) / scaling_params$insp_freq_sd)

test_model_data$pred_prob = predict(m2, newdata = test_model_data, type = "response", re.form = NA)
roc_obj = roc(test_model_data$pass_fail, test_model_data$pred_prob)

test_model_data$pred_prob_null = predict(m0, newdata = test_model_data, type = "response", re.form = NA)
roc_null = roc(test_model_data$pass_fail, test_model_data$pred_prob_null)

cat("\n=== Test set evaluation ===\n")
cat("Full model AUC:", round(auc(roc_obj), 3), "\n")
cat("Null model AUC:", round(auc(roc_null), 3), "\n")
cat("AUC gain from fixed effects:", round(auc(roc_obj) - auc(roc_null), 3), "\n")

# ==============================================================================
# 20. MODEL COMPARISON
# ==============================================================================

model_comparison = data.frame(
  Model = c("Null (intercept only)",
            "Income only",
            "Poverty only",
            "SES composite (PCA)",
            "Fully adjusted",
            "Income x time interaction"),
  AIC = c(AIC(m0), AIC(m1a), AIC(m1b), AIC(m2), AIC(m3)),
  logLik = as.numeric(c(logLik(m0), logLik(m1a), logLik(m1b),
                        logLik(m2), logLik(m3))),
  df = c(attr(logLik(m0), "df"), attr(logLik(m1a), "df"), attr(logLik(m1b), "df"),
         attr(logLik(m2), "df"), attr(logLik(m3), "df"))
) %>%
  mutate(
    Delta_AIC = AIC - AIC[1],
    AIC = round(AIC, 1),
    logLik = round(logLik, 1),
    Delta_AIC = round(Delta_AIC, 1)
  ) %>%
  select(Model, df, logLik, AIC, Delta_AIC)

cat("\n=== Model comparison ===\n")
print(model_comparison)

write_csv(model_comparison, "model_comparison.csv")



