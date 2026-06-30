import argparse
from pathlib import Path
import shutil
import zipfile

import cdsapi

VARS_ERA5 = [
    "10m_wind_gust_since_previous_post_processing",
]


def extract_zip_to_csv(zip_path: Path, final_csv: Path) -> None:
    with zipfile.ZipFile(zip_path, "r") as zf:
        members = zf.namelist()
        csv_members = [m for m in members if m.lower().endswith(".csv")]
        if not csv_members:
            raise RuntimeError(f"No CSV found inside ZIP archive: {zip_path}\nMembers: {members}")

        member = csv_members[0]
        with zf.open(member) as src, open(final_csv, "wb") as dst:
            shutil.copyfileobj(src, dst)


def main() -> None:
    ap = argparse.ArgumentParser(description="Download ERA5 single-levels gust timeseries for one stadium point.")
    ap.add_argument("--stadium", required=True)
    ap.add_argument("--lat", required=True, type=float)
    ap.add_argument("--lon", required=True, type=float)
    ap.add_argument("--start-date", required=True)
    ap.add_argument("--end-date", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    final_out = Path(args.output)
    final_out.parent.mkdir(parents=True, exist_ok=True)

    tmp_out = final_out.with_suffix(".zip")

    dataset = "reanalysis-era5-single-levels-timeseries"
    request = {
        "variable": VARS_ERA5,
        "location": {
            "longitude": args.lon,
            "latitude": args.lat,
        },
        "date": [f"{args.start_date}/{args.end_date}"],
        "data_format": "csv",
    }

    print(
        f"[CDS] ERA5 gust timeseries | stadium={args.stadium} "
        f"| lat={args.lat} lon={args.lon} "
        f"| {args.start_date} -> {args.end_date} "
        f"| output={final_out}"
    )

    client = cdsapi.Client()
    client.retrieve(dataset, request).download(str(tmp_out))

    if zipfile.is_zipfile(tmp_out):
        extract_zip_to_csv(tmp_out, final_out)
        tmp_out.unlink()
    else:
        tmp_out.rename(final_out)


if __name__ == "__main__":
    main()