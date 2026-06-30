# Rugby World Cup 2027 climate-mode exposure analysis

This repository contains the code used to analyse the climate-mode-conditioned
multi-hazard exposure of the Rugby World Cup 2027 in Australia.

The analysis combines the official tournament schedule with hourly ERA5-Land
and ERA5 reanalysis data to build historical analogue realisations of each
scheduled match over 1996--2025. Each match realisation is classified into one
of five mutually exclusive match-weather regimes: **Fair**, **Hot**, **Wet**,
**Windy** or **Compound**. The resulting type distribution is then estimated
conditionally on ENSO and Indian Ocean Dipole (IOD) phase.

## Research aim

The project asks whether the phase of large-scale climate modes, especially
ENSO and the IOD, changes the distribution of playing-condition regimes faced
by a fixed sporting mega-event schedule.

The analysis does **not** forecast the weather of the 2027 tournament itself.
Instead, it asks how the fixed 2027 schedule would have sampled historical
October--November weather conditions under different ENSO and IOD phases.

## Data

The workflow uses:

* ERA5-Land hourly data for:

  * 2-m air temperature
  * 2-m dew-point temperature
  * surface solar radiation downwards
  * total precipitation
  * 10-m wind components
* ERA5 hourly data for:

  * 10-m wind gusts
* Oceanic Niño Index (ONI) for ENSO phase classification
* Dipole Mode Index (DMI) for IOD phase classification
* Rugby World Cup 2027 schedule information:

  * venue
  * stadium
  * match date
  * local kick-off time

Raw ERA5, ERA5-Land, ONI and DMI data are not redistributed in this repository.
They must be downloaded from their original providers.

## Event-based sampling

For each scheduled match, the 2027 calendar date and local kick-off hour are
repeated in every year from 1996 to 2025. Each match is also sampled over a
±2-day window around the scheduled date.

This gives:

```text
52 matches × 30 years × 5 date offsets = 7800 schedule-conditioned realisations
```

For each realisation, the match window is defined from one hour before to two
hours after kick-off. The following event indicators are computed:

* maximum WBGT during the match window
* accumulated rainfall during the match window
* 24-h antecedent rainfall before kick-off
* 5-day antecedent rainfall before kick-off
* maximum 10-m gust during the match window
* mean wind speed during the match window

## Match-weather typology

Each realisation is assigned to hazard families using the following baseline
thresholds:

```text
Thermal: WBGT_max >= 26 °C
Wet:     match rain >= 1 mm OR 24-h rain >= 5 mm OR 5-day rain >= 25 mm
Wind:    gust_max >= 12 m s-1
```

The final match-weather type is defined as:

```text
Compound: at least two hazard families active
Hot:      only the thermal family active
Wet:      only the wet family active
Windy:    only the wind family active
Fair:     no hazard family active
```

The combined `Wet + Compound` class is used as the main moisture-affected
response variable.

## Climate-mode classification

ENSO phase is classified from the OND Oceanic Niño Index:

```text
El Niño: ONI >= +0.5 °C
La Niña: ONI <= -0.5 °C
Neutral: otherwise
```

IOD phase is classified from SON Dipole Mode Index values:

```text
Positive IOD: DMI >= +0.4 °C
Negative IOD: DMI <= -0.4 °C
Neutral IOD: otherwise
```

## Uncertainty

Uncertainty is estimated using a year-level block bootstrap. Complete years are
resampled with replacement so that the dependence among matches, venues and
date offsets within a given year is preserved.

For ENSO and IOD contrasts, years are resampled independently within the two
contrasted phase groups. The Wet+Compound share is recomputed in each group,
and the contrast distribution is used to derive percentile 95% confidence
intervals.

## Main outputs

The workflow produces:

* event-level historical analogue samples
* ENSO and IOD season classifications
* tournament-wide type distributions
* phase-conditioned type distributions
* venue-level phase sensitivities
* threshold-sensitivity diagnostics
* Wet-family precipitation-trigger decomposition
* figures used in the manuscript

## Repository structure

```text
.
├── data/
│   ├── raw/                 # user-provided raw or downloaded inputs
│   ├── processed/           # processed event samples and classifications
│   └── external/            # ONI, DMI, schedule and metadata
├── scripts/
│   ├── build_hourly_master.R
│   ├── build_rugby_event_samples.R
│   ├── rwc2027_analysis.py
│   └── make_rwc2027_figures.py
└── README.md
```

## Reproducibility note

The repository is intended to document the analysis workflow used in the
manuscript. Because raw reanalysis files are large and subject to provider
terms of use, they are not included. Reproduction requires downloading the
corresponding ERA5 and ERA5-Land fields and adapting local file paths in the
processing scripts.

## Citation

If you use this code, please cite the associated manuscript or preprint once
available.

## Licence

Licence to be added before public release.
