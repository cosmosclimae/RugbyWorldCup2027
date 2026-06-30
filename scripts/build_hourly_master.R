#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: Rscript build_hourly_master.R <stadium> <timezone> <era5land_csv> <era5gust_csv> <output_csv>")
}

stadium   <- args[1]
tz_local  <- args[2]
f_land    <- args[3]
f_gust    <- args[4]
f_out     <- args[5]

message("[INFO] stadium = ", stadium)
message("[INFO] timezone = ", tz_local)
message("[INFO] ERA5-Land = ", f_land)
message("[INFO] ERA5 gust = ", f_gust)
message("[INFO] output = ", f_out)

#---------------------------
# Helpers
#---------------------------

parse_datetime_utc <- function(x) {
  x <- as.character(x)

  out <- suppressWarnings(ymd_hms(x, tz = "UTC"))
  if (!all(is.na(out))) return(out)

  out <- suppressWarnings(ymd_hms(gsub("T", " ", x), tz = "UTC"))
  if (!all(is.na(out))) return(out)

  out <- suppressWarnings(ymd_hm(x, tz = "UTC"))
  if (!all(is.na(out))) return(out)

  out <- suppressWarnings(ymd_hm(gsub("T", " ", x), tz = "UTC"))
  if (!all(is.na(out))) return(out)

  out <- suppressWarnings(ymd(x, tz = "UTC"))
  out
}

# Relative humidity from T and Td (both in °C)
calc_rh <- function(t_c, td_c) {
  es <- 6.112 * exp((17.67 * t_c) / (t_c + 243.5))
  e  <- 6.112 * exp((17.67 * td_c) / (td_c + 243.5))
  rh <- 100 * e / es
  pmin(pmax(rh, 0), 100)
}

# Stull wet-bulb approximation (T in °C, RH in %)
calc_tw_stull <- function(t_c, rh) {
  t_c * atan(0.151977 * sqrt(rh + 8.313659)) +
    atan(t_c + rh) -
    atan(rh - 1.676331) +
    0.00391838 * (rh^(3/2)) * atan(0.023101 * rh) -
    4.686035
}

# Simple globe temperature proxy using hourly RSDS (W m-2)
calc_tg_proxy <- function(t_c, rh, rsds_wm2) {
  0.01498 * rsds_wm2 + 1.184 * t_c - 0.0789 * rh - 2.739
}

# WBGT
calc_wbgt <- function(t_c, rh, rsds_wm2) {
  tw <- calc_tw_stull(t_c, rh)
  tg <- calc_tg_proxy(t_c, rh, rsds_wm2)
  0.7 * tw + 0.2 * tg + 0.1 * t_c
}

# Heat Index in °C
calc_hi <- function(t_c, rh) {
  rh <- pmax(pmin(rh, 99), 1)
  t_f <- t_c * 9/5 + 32

  hi_f_raw <- -42.379 +
    2.04901523 * t_f +
    10.14333127 * rh -
    0.22475541 * t_f * rh -
    6.83783e-3 * t_f^2 -
    5.481717e-2 * rh^2 +
    1.22874e-3 * t_f^2 * rh +
    8.5282e-4 * t_f * rh^2 -
    1.99e-6 * t_f^2 * rh^2

  hi_f <- ifelse(t_f >= 80 & rh >= 40, hi_f_raw, t_f)
  (hi_f - 32) * 5/9
}

#---------------------------
# Read raw CSVs
#---------------------------

land_raw <- read_csv(f_land, show_col_types = FALSE)
gust_raw <- read_csv(f_gust, show_col_types = FALSE)

message("[INFO] ERA5-Land raw columns: ", paste(names(land_raw), collapse = ", "))
message("[INFO] ERA5 gust raw columns: ", paste(names(gust_raw), collapse = ", "))

#---------------------------
# Check expected columns
#---------------------------

required_land <- c("valid_time", "latitude", "longitude", "t2m", "d2m", "ssrd", "tp", "u10", "v10")
missing_land <- setdiff(required_land, names(land_raw))
if (length(missing_land) > 0) {
  stop("Missing required ERA5-Land columns: ",
       paste(missing_land, collapse = ", "),
       "\nAvailable columns: ", paste(names(land_raw), collapse = ", "))
}

required_gust <- c("valid_time", "fg10")
missing_gust <- setdiff(required_gust, names(gust_raw))
if (length(missing_gust) > 0) {
  stop("Missing required ERA5 gust columns: ",
       paste(missing_gust, collapse = ", "),
       "\nAvailable columns: ", paste(names(gust_raw), collapse = ", "))
}

#---------------------------
# Parse datetimes
#---------------------------

land <- land_raw %>%
  mutate(datetime_utc = parse_datetime_utc(valid_time))

gust <- gust_raw %>%
  mutate(datetime_utc = parse_datetime_utc(valid_time))

if (all(is.na(land$datetime_utc))) {
  stop("ERA5-Land datetime could not be parsed.")
}
if (all(is.na(gust$datetime_utc))) {
  stop("ERA5 gust datetime could not be parsed.")
}

#---------------------------
# Keep only useful columns and sort
#---------------------------

land <- land %>%
  filter(!is.na(datetime_utc)) %>%
  select(datetime_utc, latitude, longitude, t2m, d2m, ssrd, tp, u10, v10) %>%
  arrange(datetime_utc)

gust <- gust %>%
  filter(!is.na(datetime_utc)) %>%
  select(datetime_utc, gust = fg10) %>%
  arrange(datetime_utc)

message("[INFO] nrow land = ", nrow(land))
message("[INFO] nrow gust = ", nrow(gust))
message("[INFO] duplicated datetimes in land = ", sum(duplicated(land$datetime_utc)))
message("[INFO] duplicated datetimes in gust = ", sum(duplicated(gust$datetime_utc)))

#---------------------------
# Hard checks before bind
#---------------------------

if (sum(duplicated(land$datetime_utc)) > 0) {
  dup_land <- land %>%
    count(datetime_utc) %>%
    filter(n > 1)
  stop("Duplicated datetimes in ERA5-Land. Example:\n",
       paste(utils::capture.output(print(head(dup_land, 10))), collapse = "\n"))
}

if (sum(duplicated(gust$datetime_utc)) > 0) {
  dup_gust <- gust %>%
    count(datetime_utc) %>%
    filter(n > 1)
  stop("Duplicated datetimes in ERA5 gust. Example:\n",
       paste(utils::capture.output(print(head(dup_gust, 10))), collapse = "\n"))
}

if (nrow(land) != nrow(gust)) {
  stop("ERA5-Land and ERA5 gust do not have the same number of rows after parsing.")
}

if (!identical(land$datetime_utc, gust$datetime_utc)) {
  bad_land <- setdiff(land$datetime_utc, gust$datetime_utc)
  bad_gust <- setdiff(gust$datetime_utc, land$datetime_utc)
  stop(
    "Datetime vectors are not identical between ERA5-Land and ERA5 gust.\n",
    "First unmatched in land: ", paste(head(bad_land, 5), collapse = ", "), "\n",
    "First unmatched in gust: ", paste(head(bad_gust, 5), collapse = ", ")
  )
}

#---------------------------
# Safe bind
#---------------------------

df <- land %>%
  mutate(gust = gust$gust)

#---------------------------
# Unit conversions + local time
#---------------------------

df <- df %>%
  mutate(
    stadium = stadium,
    timezone = tz_local,

    t2m_c = t2m - 273.15,
    d2m_c = d2m - 273.15,
    tp_mm = tp * 1000.0,
    rsds_wm2 = ssrd / 3600.0,

    ws10 = sqrt(u10^2 + v10^2),

    # vrai datetime local pour les calculs
    datetime_local = with_tz(datetime_utc, tzone = tz_local),

    # composantes locales
    date_local = as.Date(datetime_local),
    hour_local = hour(datetime_local),
    month_local = month(datetime_local),
    year_local = year(datetime_local),

    # colonnes texte explicites pour export CSV
    datetime_utc_str   = format(datetime_utc, tz = "UTC", format = "%Y-%m-%d %H:%M:%S %Z"),
    datetime_local_str = format(datetime_local, tz = tz_local, format = "%Y-%m-%d %H:%M:%S %Z")
  )

#---------------------------
# Derived metrics
#---------------------------

df <- df %>%
  mutate(
    rh = calc_rh(t2m_c, d2m_c),
    wbgt = calc_wbgt(t2m_c, rh, rsds_wm2),
    hi = calc_hi(t2m_c, rh)
  ) %>%
  arrange(datetime_utc)

#---------------------------
# Final output
#---------------------------

out <- df %>%
  select(
    stadium, timezone,
    latitude, longitude,
    datetime_utc_str, datetime_local_str, date_local, year_local, month_local, hour_local,
    t2m_c, d2m_c, rh,
    rsds_wm2, tp_mm,
    u10, v10, ws10, gust,
    wbgt, hi
  )

dir.create(dirname(f_out), recursive = TRUE, showWarnings = FALSE)
write_csv(out, f_out)

message("[OK] Wrote hourly master: ", f_out)