library(tidyverse)
library(lubridate)
library(sf)
library(tigris)
library(leaflet)

# =============================================================================
# DATA CLEANING OVERVIEW
# =============================================================================
# This script prepares the analytical dataset for studying whether food safety
# outcomes in Boston vary with neighbourhood socioeconomic and demographic
# characteristics, after controlling for business type, license age, and
# inspection frequency.
#
# The cleaning pipeline proceeds in the following steps:
#
# 1. LOADING & DATE PARSING (raw_data)
#    - Read the Boston Food Inspection CSV.
#    - Parse resultdttm (inspection datetime) using ymd_hms() and extract year.
#    - Filter to inspections from 2015 onwards to focus on recent trends and
#      ensure sufficient sample size, dropping records with missing dates.
#
# 2. OUTCOME VARIABLE (pass_fail)
#    - Recode viol_status into a binary numeric variable: Pass = 1, Fail = 0.
#    - Records with any other value (e.g. NA, "HE_Fail") are set to NA and
#      excluded downstream via filter(!is.na(pass_fail)).
#
# 3. GEOLOCATION PARSING (lat, lon)
#    - Extract latitude and longitude from the location column, which stores
#      coordinates as "(lat, lon)" strings.
#    - Records with missing coordinates are excluded downstream.
#
# 4. CENSUS DATA (acs_data -> acs_clean)
#    - Pull 2019 ACS 5-year estimates at the ZCTA (ZIP code tabulation area)
#      level via tidycensus. The 5-year survey is used because 1-year estimates
#      are not published for small geographies like ZCTAs.
#    - Variables pulled: median household income, total population, poverty
#      status counts, full racial breakdown (White, Black, Indigenous, Asian,
#      Hawaiian, Other), Hispanic/Latino ethnicity, and foreign-born status.
#    - Data arrives in long format (one row per ZCTA per variable); pivoted
#      wide and derived rates computed (poverty_rate, pct_* for each group,
#      pct_foreign_born).
#
# 5. RECORD-LEVEL ANALYTICAL DATASET (analysis_data)
#    - Filter raw_data to valid pass_fail, lat, and lon values.
#    - Left-join ACS neighbourhood characteristics to each inspection record
#      via ZIP code (padded to 5 digits to match ZCTA GEOIDs).
#    - Deduplicate to one row per inspection visit using distinct(licenseno,
#      insp_date): the raw data has one row per violation within a visit, but
#      viol_status (and hence pass_fail) reflects the overall visit outcome,
#      so keeping multiple rows per visit would inflate inspection counts.
#
# 6. CONTROL VARIABLES
#    - insp_per_year: number of inspection visits per establishment per year,
#      computed after deduplication so it counts visits, not violations.
#      Used as a control for inspection targeting bias.
#    - license_age: years elapsed between license issuance and inspection date,
#      floored at 0 to handle any data entry errors. Proxy for establishment
#      experience/maturity.
#    - business_type: licensecat collapsed into 6 broad categories
#      (Restaurant/Food Service, Retail Food, Mobile Food, Caterer,
#      Institutional, Other) via regex matching. Stored as a factor.
#
# 7. DATA QUALITY FLAGS & TYPE CONVERSIONS
#    - acs_missing: TRUE if the inspection ZIP did not match any ZCTA in the
#      ACS data, indicating the neighbourhood-level covariates are unavailable
#      for that record.
#    - zip converted to character; year converted to integer (suitable for use
#      as a trend term or factor in regression models).
#
# OUTPUT DATASETS:
#    - analysis_data : inspection-level dataset ready for regression modelling
#                      and trend analysis (one row per inspection visit)
#    - neighbourhood_summary : ZIP-level aggregation of fail rates, used for
#                              the EDA choropleth map below
# =============================================================================

# Loading data
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

# creating binary pass/fail variable from viol_status
raw_data = raw_data %>%
  mutate(
    pass_fail = case_when(
      viol_status == "Pass" ~ 1L,
      viol_status == "Fail" ~ 0L,
      TRUE ~ NA_integer_)
  )

# parsing location into lat/lon
raw_data = raw_data %>%
  mutate(
    location = str_remove_all(location, "[()]"),
    lat = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 1])),
    lon = as.numeric(str_trim(str_split_fixed(location, ",", 2)[, 2])))

# extracting acs data
library(tidycensus)

acs_data = get_acs(
  geography = "zcta",
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

# building record-level analytical dataset: join ACS to each inspection record
analysis_data = raw_data %>%
  filter(!is.na(pass_fail), !is.na(lat), !is.na(lon)) %>%
  mutate(zip_pad = str_pad(as.character(zip), width = 5, pad = "0")) %>%
  left_join(acs_clean, by = c("zip_pad" = "GEOID")) %>%
  select(-zip_pad)

# deduplicating to one row per inspection visit
# (raw data has one row per violation, & viol_status reflects the overall inspection outcome)
analysis_data = analysis_data %>%
  distinct(licenseno, insp_date, .keep_all = TRUE)

# computing inspection frequency per establishment per year (control variable)
insp_freq = analysis_data %>%
  group_by(licenseno, year) %>%
  summarise(insp_per_year = n(), .groups = "drop")

analysis_data = analysis_data %>%
  left_join(insp_freq, by = c("licenseno", "year"))

# computing license age in years since issuance
analysis_data = analysis_data %>%
  mutate(
    iss_date = ymd_hms(issdttm),
    license_age = as.numeric(difftime(insp_date, iss_date, units = "days")) / 365.25,
    license_age = pmax(license_age, 0, na.rm = TRUE)
  )

# collapsing licensecat into broad business type categories
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

# flagging records with missing ACS data (zip did not match any ZCTA)
analysis_data = analysis_data %>%
  mutate(acs_missing = is.na(median_income))

# converting types for modelling
analysis_data = analysis_data %>%
  mutate(
    zip = as.character(zip),
    year = as.integer(year))

glimpse(analysis_data)

# neighbourhood-level summary (derived from analysis_data, used for EDA map)
neighbourhood_summary = analysis_data %>%
  group_by(zip) %>%
  summarise(
    n_inspections = n(),
    fail_rate = mean(1 - pass_fail, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(fail_rate))

print(neighbourhood_summary)

# visualising fail rate by neighbourhood
plot_data = neighbourhood_summary %>%
  filter(n_inspections >= 10) %>%
  arrange(desc(fail_rate))

# ZIP to neighbourhood crosswalk for Boston and surrounding areas
boston_neighbourhoods = tibble(
  zip = c("02101", "02108", "02109", "02110", "02111", "02113", "02114", "02115",
          "02116", "02118", "02119", "02120", "02121", "02122", "02124", "02125",
          "02126", "02127", "02128", "02129", "02130", "02131", "02132", "02134",
          "02135", "02136", "02163", "02199", "02210", "02215", "02445", "02446", 
          "02467", "02138", "02139", "02140", "02141", "02142", "02143", "02144", 
          "02145", "02472"),
  neighbourhood = c("Downtown", "Beacon Hill", "North End/Waterfront", "Financial 
   District", "Chinatown", "North End", "Beacon Hill/West End", "Fenway/Longwood",
                    "Back Bay", "South End", "Roxbury", "Mission Hill", "Dorchester (Grove Hall)",
                    "Dorchester (Neponset)", "Dorchester (Codman Square)", "Dorchester (Savin Hill)",
                    "Mattapan", "South Boston", "East Boston", "Charlestown", "Jamaica Plain", "Roslindale",
                    "West Roxbury", "Allston", "Brighton", "Hyde Park", "Allston (HBS)", 
                    "Back Bay (Prudential)", "Seaport/Fort Point", "Fenway/Kenmore", "Brookline", 
                    "Brookline (Longwood/Coolidge Corner)", "Chestnut Hill", "Cambridge (Harvard Square)",
                    "Cambridge (Central/MIT)", "Cambridge (North)", "Cambridge (East)", "Cambridge (Kendall Square)", "Somerville (Davis Square)", "Somerville", "Somerville (East)", "Watertown"
  )
)

options(tigris_use_cache = TRUE)
boston_zctas = zctas(year = 2020, cb = TRUE) %>%
  filter(ZCTA5CE20 %in% as.character(plot_data$zip))

map_data = boston_zctas %>%
  left_join(
    plot_data %>% mutate(zip = as.character(zip)),
    by = c("ZCTA5CE20" = "zip")
  ) %>%
  left_join(boston_neighbourhoods, by = c("ZCTA5CE20" = "zip")) %>%
  mutate(neighbourhood = if_else(is.na(neighbourhood), ZCTA5CE20, neighbourhood)) %>%
  st_transform(crs = 4326)

pal = colorNumeric(palette = c("#f3d357", "#b90d0d"), domain = map_data$fail_rate)

leaflet(map_data) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~pal(fail_rate),
    fillOpacity = 0.7,
    color = "white",
    weight = 1,
    popup = ~paste0("<b>", neighbourhood, "</b><br>Fail rate: ", sprintf("%.1f%%", fail_rate * 100))
  ) %>%
  addLegend(
    pal = pal,
    values = ~fail_rate,
    title = "Proportion Failed",
    labFormat = labelFormat(suffix = "%", transform = function(x) x * 100),
    position = "bottomright"
    #decreasing = TRUE #commented it out because it's not working :(
  )

