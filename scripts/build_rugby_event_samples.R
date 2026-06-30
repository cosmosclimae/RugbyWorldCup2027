suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(lubridate)
})

# ---- paths: edit as needed ----
schedule_file <- "f:/climae/RugbyWC2027/data/match_schedule_2027.csv"
stadium_dir   <- "f:/climae/RugbyWC2027/processed/hourly_master"     # contains adelaide_oval.csv, perth_stadium.csv, etc.
output_file    <- "f:/climae/RugbyWC2027/processed/rugbywc2027_event_samples_master.csv"



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
})

# =========================
# PARAMETRES
# =========================

schedule_file <- "f:/climae/RugbyWC2027/data/match_schedule_2027.csv"

# Dossier contenant les 8 CSV stades
stadium_dir   <- "f:/climae/RugbyWC2027/processed/hourly_master" 

# Fichier de sortie
out_file      <- "f:/climae/RugbyWC2027/processed/rugbywc2027_event_samples_master.csv"

# Offsets et périodes
offset_days <- -2:2
p1_years <- 1996:2010
p2_years <- 2011:2025

# Seuils
WBGT24_THRESHOLD <- 24
WBGT26_THRESHOLD <- 26
HI32_THRESHOLD   <- 32
RAIN_MATCH_THRESHOLD <- 1      # mm sur la fenêtre match
RAIN24H_THRESHOLD    <- 5      # mm sur 24h avant match
RAIN5D_THRESHOLD     <- 25     # mm cumul 5 jours avant match (antecedent / saturation surface)
GUST_THRESHOLD       <- 12     # m/s

# Fenêtres
MATCH_WINDOW_BEFORE_HOURS <- 1
MATCH_WINDOW_AFTER_HOURS  <- 2

# =========================
# FONCTIONS
# =========================

safe_na_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_na_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_na_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_na_sum <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
}

find_stadium_file <- function(stadium_slug, stadium_dir) {
  f <- file.path(stadium_dir, paste0(stadium_slug, ".csv"))
  if (!file.exists(f)) {
    stop(sprintf("Fichier stade introuvable: %s", f))
  }
  f
}

prepare_stadium_dt <- function(file) {
  dt <- fread(file)

  required_cols <- c(
    "stadium", "timezone", "latitude", "longitude",
    "datetime_local_str", "date_local", "year_local", "hour_local",
    "tp_mm", "ws10", "gust", "wbgt", "hi"
  )

  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Colonnes manquantes dans %s : %s",
      file, paste(missing_cols, collapse = ", ")
    ))
  }

  # datetime_local_str du type "1996-01-01 11:30:00 ACDT"
  # On enlève l'abréviation de timezone finale pour parser proprement
  dt[, datetime_local_clean := sub(" [A-Z]{3,4}$", "", datetime_local_str)]
  dt[, datetime_local := ymd_hms(datetime_local_clean, tz = "UTC")]

  # IMPORTANT :
  # On garde ces datetimes comme "horloge locale" pour faire les comparaisons
  # entre match schedule local et données locales. On ne reconvertit pas.
  # Ici tz="UTC" sert juste à stocker des datetimes cohérents sans décalage.

  setorder(dt, datetime_local)

  dt
}

make_local_datetime <- function(date_str, time_str) {
  ymd_hm(paste(date_str, time_str), tz = "UTC")
}

nearest_available_time <- function(target_dt, available_dts) {
  idx <- which.min(abs(as.numeric(difftime(available_dts, target_dt, units = "secs"))))
  available_dts[idx]
}

format_dt <- function(x) {
  if (is.na(x)) return(NA_character_)
  format(x, "%Y-%m-%d %H:%M:%S")
}

extract_one_event_line <- function(match_row, stadium_dt, target_year, period_label, offset_day) {
  # Match original 2027
  kickoff_orig <- make_local_datetime(match_row$match_date, match_row$match_time_local)

  # Transposition dans l'année cible avec offset
  month_day <- format(kickoff_orig, "%m-%d")
  hhmmss    <- format(kickoff_orig, "%H:%M:%S")

  sample_kickoff_nominal <- ymd_hms(
    paste0(target_year, "-", month_day, " ", hhmmss),
    tz = "UTC"
  ) + days(offset_day)

  # Nearest timestamp réellement disponible dans le fichier stade
  kickoff_nearest <- nearest_available_time(
    target_dt = sample_kickoff_nominal,
    available_dts = stadium_dt$datetime_local
  )

  # Fenêtre match fondée sur le nearest réel
  window_start <- kickoff_nearest - hours(MATCH_WINDOW_BEFORE_HOURS)
  window_end   <- kickoff_nearest + hours(MATCH_WINDOW_AFTER_HOURS)

  # Fenêtre pluie 24h avant
  rain24h_start <- kickoff_nearest - hours(24)
  rain24h_end   <- kickoff_nearest

  # Fenêtre pluie 5 jours avant (cumul antecedent / saturation du terrain)
  rain5d_start <- kickoff_nearest - hours(120)
  rain5d_end   <- kickoff_nearest

  # Nuit de récupération: 22:00 → 06:00 après le match
  # Si le match se termine tard, on garde quand même la nuit civile suivante du jour de match
  sample_date <- as.Date(kickoff_nearest)
  night_start <- ymd_hms(paste0(sample_date, " 22:00:00"), tz = "UTC")
  night_end   <- night_start + hours(8)

  # Sous-ensembles
  dt_match <- stadium_dt[datetime_local >= window_start & datetime_local <= window_end]
  dt_r24   <- stadium_dt[datetime_local > rain24h_start & datetime_local <= rain24h_end]
  dt_r5d   <- stadium_dt[datetime_local > rain5d_start  & datetime_local <= rain5d_end]
  dt_night <- stadium_dt[datetime_local >= night_start & datetime_local <= night_end]

  # Ligne kickoff exacte = nearest timestamp
  dt_kick <- stadium_dt[datetime_local == kickoff_nearest]

  if (nrow(dt_kick) != 1) {
    stop(sprintf(
      "Kickoff nearest non unique pour match_id=%s, year=%s, offset=%s",
      match_row$match_id, target_year, offset_day
    ))
  }

  # Indicateurs
  WBGT_kickoff <- dt_kick$wbgt[1]
  WBGT_max     <- safe_na_max(dt_match$wbgt)

  HI_kickoff   <- dt_kick$hi[1]
  HI_max       <- safe_na_max(dt_match$hi)

  rain_window  <- safe_na_sum(dt_match$tp_mm)
  wind_mean    <- safe_na_mean(dt_match$ws10)
  gust_max     <- safe_na_max(dt_match$gust)

  rain_24h_before <- safe_na_sum(dt_r24$tp_mm)
  rain_5d_before  <- safe_na_sum(dt_r5d$tp_mm)

  Tmin_night <- safe_na_min(dt_night$t2m_c)
  HI_night   <- safe_na_max(dt_night$hi)

  # Flags
  WBGT24_flag    <- as.integer(!is.na(WBGT_max) && WBGT_max > WBGT24_THRESHOLD)
  WBGT26_flag    <- as.integer(!is.na(WBGT_max) && WBGT_max > WBGT26_THRESHOLD)
  HI32_flag      <- as.integer(!is.na(HI_max)   && HI_max   > HI32_THRESHOLD)

  rain_match_flag <- as.integer(!is.na(rain_window)     && rain_window     > RAIN_MATCH_THRESHOLD)
  rain24h_flag    <- as.integer(!is.na(rain_24h_before) && rain_24h_before > RAIN24H_THRESHOLD)
  rain5d_flag     <- as.integer(!is.na(rain_5d_before)  && rain_5d_before  > RAIN5D_THRESHOLD)
  gust_flag       <- as.integer(!is.na(gust_max)        && gust_max        > GUST_THRESHOLD)

  compound_count <- WBGT24_flag + rain_match_flag + gust_flag

  data.table(
    match_id = match_row$match_id,
    match_date = match_row$match_date,
    match_time_local = match_row$match_time_local,
    kickoff_datetime_local = match_row$kickoff_datetime_local,
    kickoff_hour_local = hour(kickoff_orig) + minute(kickoff_orig) / 60,

    phase = match_row$phase,
    phase_order = match_row$phase_order,
    pool = match_row$pool,
    round_label = match_row$round_label,
    team_A = match_row$team_A,
    team_B = match_row$team_B,
    side_A_type = match_row$side_A_type,
    side_B_type = match_row$side_B_type,

    stadium = match_row$stadium,
    venue_name = match_row$venue_name,
    city = match_row$city,
    venue_label = match_row$venue_label,
    lat = stadium_dt$latitude[1],
    lon = stadium_dt$longitude[1],

    year = target_year,
    period = period_label,
    offset_day = offset_day,

    sample_date_local = as.character(as.Date(sample_kickoff_nominal)),
    sample_kickoff_datetime_local = format_dt(sample_kickoff_nominal),
    kickoff_hour_nearest = format(kickoff_nearest, "%H:%M:%S"),

    window_start_local = format_dt(window_start),
    window_end_local = format_dt(window_end),

    rain24h_start_local = format_dt(rain24h_start),
    rain24h_end_local = format_dt(rain24h_end),

    rain5d_start_local = format_dt(rain5d_start),
    rain5d_end_local = format_dt(rain5d_end),

    night_start_local = format_dt(night_start),
    night_end_local = format_dt(night_end),

    WBGT_kickoff = WBGT_kickoff,
    WBGT_max = WBGT_max,
    HI_kickoff = HI_kickoff,
    HI_max = HI_max,
    rain_window = rain_window,
    wind_mean = wind_mean,
    gust_max = gust_max,
    rain_24h_before = rain_24h_before,
    rain_5d_before = rain_5d_before,
    Tmin_night = Tmin_night,
    HI_night = HI_night,

    WBGT24_flag = WBGT24_flag,
    WBGT26_flag = WBGT26_flag,
    HI32_flag = HI32_flag,
    rain_match_flag = rain_match_flag,
    rain24h_flag = rain24h_flag,
    rain5d_flag = rain5d_flag,
    gust_flag = gust_flag,

    compound_count = compound_count
  )
}

# =========================
# LECTURE CALENDRIER
# =========================

schedule_dt <- fread(schedule_file)

required_schedule_cols <- c(
  "match_id", "match_date", "match_time_local", "kickoff_datetime_local",
  "stadium", "venue_name", "city", "venue_label",
  "phase", "phase_order", "pool", "round_label",
  "team_A", "team_B", "side_A_type", "side_B_type"
)

missing_schedule_cols <- setdiff(required_schedule_cols, names(schedule_dt))
if (length(missing_schedule_cols) > 0) {
  stop(sprintf(
    "Colonnes manquantes dans le calendrier : %s",
    paste(missing_schedule_cols, collapse = ", ")
  ))
}

# =========================
# PREP STADES
# =========================

stadium_ids <- unique(schedule_dt$stadium)

stadium_list <- vector("list", length(stadium_ids))
names(stadium_list) <- stadium_ids

for (s in stadium_ids) {
  f <- find_stadium_file(s, stadium_dir)
  cat(sprintf("Lecture stade: %s\n", f))
  stadium_list[[s]] <- prepare_stadium_dt(f)
}

# =========================
# CONSTRUCTION CSV MAITRE
# =========================

all_rows <- vector("list", nrow(schedule_dt) * (length(p1_years) + length(p2_years)) * length(offset_days))
k <- 1L

for (i in seq_len(nrow(schedule_dt))) {
  match_row <- schedule_dt[i]
  stadium_dt <- stadium_list[[match_row$stadium]]

  # Période 1
  for (yy in p1_years) {
    for (od in offset_days) {
      all_rows[[k]] <- extract_one_event_line(
        match_row = match_row,
        stadium_dt = stadium_dt,
        target_year = yy,
        period_label = "p1",
        offset_day = od
      )
      k <- k + 1L
    }
  }

  # Période 2
  for (yy in p2_years) {
    for (od in offset_days) {
      all_rows[[k]] <- extract_one_event_line(
        match_row = match_row,
        stadium_dt = stadium_dt,
        target_year = yy,
        period_label = "p2",
        offset_day = od
      )
      k <- k + 1L
    }
  }
}

result_dt <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)

# Ordre final
setcolorder(result_dt, c(
  "match_id",
  "match_date",
  "match_time_local",
  "kickoff_datetime_local",
  "kickoff_hour_local",
  "phase",
  "phase_order",
  "pool",
  "round_label",
  "team_A",
  "team_B",
  "side_A_type",
  "side_B_type",
  "stadium",
  "venue_name",
  "city",
  "venue_label",
  "lat",
  "lon",
  "year",
  "period",
  "offset_day",
  "sample_date_local",
  "sample_kickoff_datetime_local",
  "kickoff_hour_nearest",
  "window_start_local",
  "window_end_local",
  "rain24h_start_local",
  "rain24h_end_local",
  "rain5d_start_local",
  "rain5d_end_local",
  "night_start_local",
  "night_end_local",
  "WBGT_kickoff",
  "WBGT_max",
  "HI_kickoff",
  "HI_max",
  "rain_window",
  "wind_mean",
  "gust_max",
  "rain_24h_before",
  "rain_5d_before",
  "Tmin_night",
  "HI_night",
  "WBGT24_flag",
  "WBGT26_flag",
  "HI32_flag",
  "rain_match_flag",
  "rain24h_flag",
  "rain5d_flag",
  "gust_flag",
  "compound_count"
))

fwrite(result_dt, output_file)

cat("\nFichier créé : ", output_file, "\n", sep = "")
cat("Nombre de lignes : ", nrow(result_dt), "\n", sep = "")
cat("Nombre attendu    : ", nrow(schedule_dt) * 30 * 5, "\n", sep = "")