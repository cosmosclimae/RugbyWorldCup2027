shell.executable("/bin/bash")
from pathlib import Path
VENV_ACTIVATE = str(config.get("venv_activate", "~/venv_cds/bin/activate"))
STADIUMS = {
    "adelaide_oval": {
        "lat": -34.915556,
        "lon": 138.596111,
    },
    "brisbane_stadium": {
        "lat": -27.464722,
        "lon": 153.009444,
    },
    "docklands_stadium": {
        "lat": -37.816389,
        "lon": 144.947500,
    },
    "newcastle_stadium": {
        "lat": -32.918889,
        "lon": 151.726667,
    },
    "perth_stadium": {
        "lat": -31.951111,
        "lon": 115.889167,
    },
    "stadium_australia": {
        "lat": -33.847222,
        "lon": 151.063333,
    },
    "sydney_football_stadium": {
        "lat": -33.889167,
        "lon": 151.225278,
    },
    "north_queensland_stadium": {
        "lat": -19.266111,
        "lon": 146.816667,
    },
}
TIMEZONES = {
    "adelaide_oval": "Australia/Adelaide",
    "brisbane_stadium": "Australia/Brisbane",
    "docklands_stadium": "Australia/Melbourne",
    "newcastle_stadium": "Australia/Sydney",
    "perth_stadium": "Australia/Perth",
    "stadium_australia": "Australia/Sydney",
    "sydney_football_stadium": "Australia/Sydney",
    "north_queensland_stadium": "Australia/Brisbane",
}
START_DATE = "1996-01-01"
END_DATE = "2025-12-31"

RAW_ERA5LAND_DIR = Path("/mnt/f/climae/RugbyWC2027/raw/era5land_timeseries")
RAW_ERA5_DIR = Path("/mnt/f/climae/RugbyWC2027/raw/era5_timeseries")
HOURLY_MASTER_DIR = Path("/mnt/f/climae/RugbyWC2027/processed/hourly_master")
SCRIPT_ERA5LAND = Path(workflow.basedir) / "scripts" / "download_era5land_fifa_month.py"
SCRIPT_ERA5 = Path(workflow.basedir) / "scripts" / "download_era5_fifa_month.py"
SCRIPT_HOURLY_MASTER = Path(workflow.basedir) / "scripts" / "build_hourly_master.R"

rule all:
    input:
        expand(str(RAW_ERA5LAND_DIR / "{stadium}.csv"), stadium=STADIUMS.keys()),
        expand(str(RAW_ERA5_DIR / "{stadium}.csv"), stadium=STADIUMS.keys()),
        expand(str(HOURLY_MASTER_DIR / "{stadium}.csv"), stadium=STADIUMS.keys()),

rule download_era5land_timeseries:
    output:
        str(RAW_ERA5LAND_DIR / "{stadium}.csv")
    params:
        lat=lambda wc: STADIUMS[wc.stadium]["lat"],
        lon=lambda wc: STADIUMS[wc.stadium]["lon"],
        start=START_DATE,
        end=END_DATE,
    log:
        "logs/era5land/{stadium}.log"
    shell:
        r"""
        set -euo pipefail           
        mkdir -p {RAW_ERA5LAND_DIR} logs/era5land
            source {VENV_ACTIVATE} && \
            python '{SCRIPT_ERA5LAND}' \
                --stadium {wildcards.stadium} \
                --lat {params.lat} \
                --lon {params.lon} \
                --start-date {params.start} \
                --end-date {params.end} \
                --output {output} \
                > {log} 2>&1
        """

rule download_era5_gust_timeseries:
    output:
        str(RAW_ERA5_DIR / "{stadium}.csv")
    params:
        lat=lambda wc: STADIUMS[wc.stadium]["lat"],
        lon=lambda wc: STADIUMS[wc.stadium]["lon"],
        start=START_DATE,
        end=END_DATE,
    log:
        "logs/era5/{stadium}.log"
    shell:
        r"""
        set -euo pipefail           
        mkdir -p {RAW_ERA5LAND_DIR} logs/era5land
            mkdir -p {RAW_ERA5_DIR} logs/era5
            source {VENV_ACTIVATE} && \
            python '{SCRIPT_ERA5}' \
                --stadium {wildcards.stadium} \
                --lat {params.lat} \
                --lon {params.lon} \
                --start-date {params.start} \
                --end-date {params.end} \
                --output {output} \
                > {log} 2>&1
            """

rule build_hourly_master:
    input:
        land=str(RAW_ERA5LAND_DIR / "{stadium}.csv"),
        gust=str(RAW_ERA5_DIR / "{stadium}.csv"),
    output:
        str(HOURLY_MASTER_DIR / "{stadium}.csv")
    params:
        timezone=lambda wc: TIMEZONES[wc.stadium],
    log:
        "logs/hourly_master/{stadium}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {HOURLY_MASTER_DIR} logs/hourly_master
        Rscript {SCRIPT_HOURLY_MASTER} \
            {wildcards.stadium} \
            {params.timezone} \
            {input.land} \
            {input.gust} \
            {output} \
            > {log} 2>&1
        """