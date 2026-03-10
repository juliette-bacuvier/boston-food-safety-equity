library(tidyverse)
library(lubridate)
library(sf)
library(tigris)
library(leaflet)
library(tidycensus)

options(tigris_use_cache = TRUE)

#--------------
# Loading data
#-------------
raw_data = read_csv('/Users/juliettebacuvier/desktop/tmpp7jeda9f.csv')

names(raw_data)
#glimpse(raw_data)

# parsing dates & filtering to last 5 years
raw_data = raw_data %>%
   mutate(
      insp_date = ymd_hms(resultdttm),
      year = year(insp_date)
   ) %>%
   filter(year >= 2015, !is.na(insp_date))

# -------------------------------------------------------------
# Outcome variable —-> investigate both viol_status & result
# -------------------------------------------------------------

cat("=== viol_status distribution ===\n")
raw_data %>% count(viol_status, sort = TRUE) %>% print(n = 20)

cat("\n=== result distribution ===\n")
raw_data %>% count(result, sort = TRUE) %>% print(n = 20)

# checking if viol_status varies within a single inspection visit
# (if it does, it's per-violation, not per-inspection)
conflicting = raw_data %>%
  group_by(licenseno, insp_date) %>%
  summarise(n_statuses = n_distinct(viol_status), .groups = "drop") %>%
  filter(n_statuses > 1)
cat("\nVisits with conflicting viol_status:", nrow(conflicting), "\n")

raw_data = raw_data %>%
  mutate(
    pass_fail = case_when(
      str_detect(result, regex("pass", ignore_case = TRUE)) ~ 1L,
      str_detect(result, regex("fail", ignore_case = TRUE)) ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# ----------------------------
# Geography parsing
# ----------------------------

raw_data = raw_data %>%
  mutate(
    location = str_remove_all(location, "[()]"),
    lat = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 1])),
    lon = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 2]))
  )

# ---------------------------
# Census tract ACS data
# ---------------------------

acs_data = get_acs(
   geography = "tract",
   state = "MA",
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
      total_nativity = "B05002_001"
   ),
   year = 2019,
   survey = "acs5"
)

glimpse(acs_data)

# cleaning acs data: pivot wide and compute derived variables
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
      pct_foreign_born = ifelse(total_nativity == 0, NA_real_, foreign_born / total_nativity)
   )

cat("\nACS tract data:", nrow(acs_clean), "tracts\n")

# -------------------------------------------------------------
# Spatial join --> assigning each inspection to a census tract
# -------------------------------------------------------------

# loading MA census tract geometries
ma_tracts = tracts(state = "MA", year = 2019, cb = TRUE) %>%
  st_transform(crs = 4326)

# filtering to valid inspection records
analysis_data = raw_data %>%
  filter(!is.na(pass_fail), !is.na(lat), !is.na(lon))

# converting inspections to sf points
inspections_sf = analysis_data %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# spatial join: find which tract each inspection point falls in
inspections_sf = st_join(inspections_sf, ma_tracts %>% select(GEOID), join = st_within)

# pulling the tract GEOID back into the regular data frame
analysis_data$tract_geoid = inspections_sf$GEOID

# reporting match rate
cat("\nTract spatial join match rate:",
    round(mean(!is.na(analysis_data$tract_geoid)) * 100, 1), "%\n")

# -------------------------------------------------------------
# Joining ACS to inspections via tract GEOID
# -------------------------------------------------------------

analysis_data = analysis_data %>%
  left_join(acs_clean, by = c("tract_geoid" = "GEOID"))

# -------------------------------------------------------------
# Deduplication --> one row per inspection visit
# -------------------------------------------------------------

# The raw data has one row per violation within a visit, but pass_fail
# reflects the overall visit outcome. Keeping duplicates inflates counts.
# Take the WORST outcome per visit (Fail = 0 beats Pass = 1).
analysis_data = analysis_data %>%
  group_by(licenseno, insp_date) %>%
  arrange(pass_fail) %>% # puts 0 (Fail) before 1 (Pass)
  slice(1) %>% # keeps the worst outcome

ungroup()

# -------------------------------------------------------------
# "Control variables"
# -------------------------------------------------------------

# inspection frequency per establishment per year
insp_freq = analysis_data %>%
  group_by(licenseno, year) %>%
  summarise(insp_per_year = n(), .groups = "drop")

analysis_data = analysis_data %>%
  left_join(insp_freq, by = c("licenseno", "year"))

# license age in years since issuance
analysis_data = analysis_data %>%
  mutate(
    iss_date = ymd_hms(issdttm),
    license_age = as.numeric(difftime(insp_date, iss_date, units = "days")) / 365.25,
    license_age = pmax(license_age, 0, na.rm = TRUE)
  )

# broad business type categories
analysis_data = analysis_data %>%
  mutate(
    business_type = case_when(
      str_detect(licensecat, regex("eating|drinking|food service|restaurant", ignore_case = TRUE)) ~ "Restaurant/Food Service",
      str_detect(licensecat, regex("retail|grocery|market|store|supermarket", ignore_case = TRUE)) ~ "Retail Food",
      str_detect(licensecat, regex("mobile|pushcart|cart|truck", ignore_case = TRUE)) ~ "Mobile Food",
      str_detect(licensecat, regex("cater", ignore_case = TRUE)) ~ "Caterer",
      str_detect(licensecat, regex("school|institution|hospital|care|clinic", ignore_case = TRUE)) ~ "Institutional",
      TRUE ~ "Other"
    ),
    business_type = factor(business_type)
  )

# -------------------------------------------------------------
# Data quality flags & type conversions
# -------------------------------------------------------------

analysis_data = analysis_data %>%
  mutate(
    acs_missing = is.na(median_income),
    zip = as.character(zip),
    year = as.integer(year)
  )

# -------------------------------------------------------------
# Summary diagnositcs
# -------------------------------------------------------------

cat("\n=== Final dataset ===\n")
cat("Rows:", nrow(analysis_data), "\n")
cat("Unique tracts:", n_distinct(analysis_data$tract_geoid, na.rm = TRUE), "\n")
cat("Unique ZIPs:", n_distinct(analysis_data$zip, na.rm = TRUE), "\n")
cat("ACS match rate:", round(mean(!analysis_data$acs_missing) * 100, 1), "%\n")
cat("Year range:", min(analysis_data$year), "-", max(analysis_data$year), "\n\n")

# Missingness summary
cat("=== Missingness ===\n")
analysis_data %>%
  summarise(across(
    c(pass_fail, lat, lon, tract_geoid, median_income, poverty_rate,
      pct_white, pct_black, pct_hispanic, pct_foreign_born,
      license_age, insp_per_year, business_type),
    ~ sum(is.na(.))
  )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(analysis_data) * 100, 1)) %>%
  arrange(desc(n_missing)) %>%
  print(n = 20)

# Fail rate range across tracts
cat("\n=== Fail rate range across tracts (min 10 inspections) ===\n")
tract_fail = analysis_data %>%
  filter(!is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(
    n = n(),
    fail_rate = mean(1 - pass_fail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

cat("Tracts:", nrow(tract_fail), "\n")
cat("Fail rate range:", round(min(tract_fail$fail_rate) * 100, 1), "% to",
    round(max(tract_fail$fail_rate) * 100, 1), "%\n")
cat("Fail rate SD:", round(sd(tract_fail$fail_rate) * 100, 1), "pp\n")

glimpse(analysis_data)



#===============================================================================
# EXPLORATORY DATA ANALYSIS
# ==============================================================================
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("ggcorrplot", quietly = TRUE)) install.packages("ggcorrplot", repos = "https://cloud.r-project.org")
library(ggplot2)
library(ggcorrplot)

#---------------------------------------------------------------------------
# Chloropleth Map
#---------------------------------------------------------------------------
neighbourhood_summary = analysis_data %>%
  filter(!is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarise(
    n_inspections = n(),
    fail_rate = mean(1 - pass_fail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(fail_rate))

print(neighbourhood_summary)

plot_data = neighbourhood_summary %>%
  filter(n_inspections >= 10)

# Filter tract geometries to those in our data
map_tracts = ma_tracts %>%
  filter(GEOID %in% plot_data$tract_geoid) %>%
  left_join(plot_data, by = c("GEOID" = "tract_geoid"))

pal = colorNumeric(
  palette = c("#f3d357", "#b90d0d"),
  domain = map_tracts$fail_rate
)

leaflet(map_tracts) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~pal(fail_rate),
    fillOpacity = 0.7,
    color = "white",
    weight = 1,
    popup = ~paste0(
      "<b>Tract ", GEOID, "</b><br>",
      "Fail rate: ", sprintf("%.1f%%", fail_rate * 100), "<br>",
      "Inspections: ", n_inspections
    )
  ) %>%
  addLegend(
    pal = pal,
    values = ~fail_rate,
    title = "Proportion Failed",
    labFormat = labelFormat(suffix = "%", transform = function(x) x * 100),
    position = "bottomright"
  )

#---------------------------------------------------------------------------
# Fail rate vs. Median income (ZIP-level scatterplot)
#---------------------------------------------------------------------------

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
    .groups = "drop"
  )

ggplot(tract_summary, aes(x = median_income / 1000, y = fail_rate)) +
  geom_point(aes(size = n_inspections), alpha = 0.5, color = "#b90d0d") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.8) +
  scale_x_continuous(labels = scales::dollar_format(suffix = "K")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_size_continuous(range = c(1.5, 8), name = "Inspections") +
  labs(
    title = "Food Inspection Fail Rate vs. Median Household Income",
    subtitle = "Each point is a census tract (min. 10 inspections); line shows linear trend",
    x = "Median Household Income (ACS 2019, tract level)",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates (tract level)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("plot1_failrate_vs_income_tract.png", width = 9, height = 6, dpi = 300)


#---------------------------------------------------------------------------
# Fail rate over time by income tercile --> to show whether disparities are growing, shrinking or stable
#---------------------------------------------------------------------------

tract_income = analysis_data %>%
  filter(!is.na(median_income), !is.na(tract_geoid)) %>%
  distinct(tract_geoid, median_income) %>%
  mutate(
    income_tercile = ntile(median_income, 3),
    income_tercile = factor(income_tercile,
                            labels = c("Low Income", "Middle Income", "High Income"))
  )

time_trends = analysis_data %>%
  filter(!is.na(median_income), !is.na(pass_fail), !is.na(tract_geoid)) %>%
  left_join(tract_income %>% select(tract_geoid, income_tercile), by = "tract_geoid") %>%
  group_by(year, income_tercile) %>%
  summarise(
    fail_rate = mean(1 - pass_fail, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

ggplot(time_trends, aes(x = year, y = fail_rate, color = income_tercile)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_color_manual(
    values = c("Low Income" = "#b90d0d", "Middle Income" = "#f0a500", "High Income" = "#2a7b9b"),
    name = "Neighbourhood\nIncome Level"
  ) +
  labs(
    title = "Inspection Fail Rate Over Time by Neighbourhood Income Level",
    subtitle = "Census tracts grouped into income terciles (ACS median household income)",
    x = "Year",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates (tract level)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("plot2_failrate_time_trends_tract.png", width = 9, height = 6, dpi = 300)

#---------------------------------------------------------------------------
# Fail rate by business type
#---------------------------------------------------------------------------

business_summary = analysis_data %>%
  filter(!is.na(pass_fail)) %>%
  group_by(business_type) %>%
  summarise(
    n_inspections = n(),
    fail_rate = mean(1 - pass_fail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(fail_rate))

table(analysis_data$licensecat)
#will need to modify this in the future bc graph is not what i want :(

ggplot(business_summary, aes(x = reorder(business_type, fail_rate), y = fail_rate)) +
  geom_col(fill = "#b90d0d", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = paste0(round(fail_rate * 100, 1), "%")),
            hjust = -0.15, size = 3.5) +
  geom_text(aes(label = paste0("n=", scales::comma(n_inspections))),
            y = 0.002, hjust = 0, size = 3, color = "white") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = "Inspection Fail Rate by Business Type",
    subtitle = "Justifies including business type as a control variable",
    x = NULL,
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave("plot3_failrate_by_business.png", width = 9, height = 5, dpi = 300)

#---------------------------------------------------------------------------
# Correlation heatmap of ACS variables
#---------------------------------------------------------------------------

cor_vars = analysis_data %>%
  filter(!acs_missing, !is.na(tract_geoid)) %>%
  distinct(tract_geoid, .keep_all = TRUE) %>%
  select(median_income, poverty_rate, pct_white, pct_black, pct_asian,
         pct_hispanic, pct_foreign_born) %>%
  drop_na()

cor_matrix = cor(cor_vars, use = "complete.obs")

colnames(cor_matrix) = c("Median Income", "Poverty Rate", "% White", "% Black",
                          "% Asian", "% Hispanic", "% Foreign Born")
rownames(cor_matrix) = colnames(cor_matrix)

ggcorrplot(
  cor_matrix,
  type = "lower",
  lab = TRUE,
  lab_size = 3.5,
  colors = c("#2a7b9b", "white", "#b90d0d"),
  title = "Correlation Between Tract-Level Demographic Variables",
  ggtheme = theme_minimal(base_size = 13)
) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plot4_correlation_heatmap_tract.png", width = 8, height = 7, dpi = 300)

# -----------------------------------------------------
# Fail rate vs. Poverty rate (tract-lvl)
# -----------------------------------------------------

ggplot(tract_summary, aes(x = poverty_rate, y = fail_rate)) +
  geom_point(aes(size = n_inspections), alpha = 0.5, color = "#b90d0d") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.8) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_size_continuous(range = c(1.5, 8), name = "Inspections") +
  labs(
    title = "Food Inspection Fail Rate vs. Poverty Rate",
    subtitle = "Each point is a census tract (min. 10 inspections)",
    x = "Poverty Rate (ACS 2019, tract level)",
    y = "Inspection Fail Rate",
    caption = "Source: Boston Food Inspection data & ACS 5-year estimates (tract level)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("plot5_failrate_vs_poverty_tract.png", width = 9, height = 6, dpi = 300)

