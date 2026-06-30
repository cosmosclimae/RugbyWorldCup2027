import argparse
from pathlib import Path

import cdsapi

# [N, W, S, E]
AREAS = {
    "france": [51.5, -5.5, 41.0, 10.0],
    # Domaine Sénégal (tu avais mis ça dans AREA_FRANCE par erreur)
    "senegal": [17.0, -18.0, 12.0, -11.0],
}

VARS_BASELINE_TP_T2M = [
    "2m_temperature",
    "total_precipitation",
]

VARS_MONTHLY_FULL = [
    "2m_dewpoint_temperature",
    "2m_temperature",
    "surface_solar_radiation_downwards",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "total_precipitation",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--country",
        choices=sorted(AREAS.keys()),
        required=True,
        help="Domaine géographique (détermine le 'area' CDS).",
    )
    ap.add_argument("--year", required=True)
    ap.add_argument("--month", required=True)  # "01".."12"
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--profile",
        choices=["baseline", "monthly"],
        required=True,
        help="baseline = t2m+tp only, monthly = full variables",
    )
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    variables = VARS_BASELINE_TP_T2M if args.profile == "baseline" else VARS_MONTHLY_FULL
    area = AREAS[args.country]

    request = {
        "variable": variables,
        "year": args.year,
        "month": args.month,
        "day": [f"{d:02d}" for d in range(1, 32)],
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": area,
        "data_format": "netcdf",
        "download_format": "unarchived",
    }

    print(
        f"[CDS] ERA5-Land {args.profile} download {args.country} {args.year}-{args.month} -> {out}"
    )
    c = cdsapi.Client()
    c.retrieve("reanalysis-era5-land", request, str(out))


if __name__ == "__main__":
    main()
