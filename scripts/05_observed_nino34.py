"""
Observed Nino 3.4 anomalies for the lead-in segment of the Nino 3.4 plume.

Each recent GiOCEAN monthly ocean field is turned into an index the same
way a forecast is: the sea-surface temperature averaged over the Nino 3.4
box (5S-5N, 170W-120W), as an anomaly against the 1998-2024 mean for that
calendar month. Months not yet in the main GiOCEAN directory are read
from the near-real-time stream.

Runs on Discover (the GiOCEAN ocean data lives there) and writes a small
CSV that is committed, so the plume script on PFE can draw the observed
segment:

    outputs/observed_recent/nino34/observed_nino34.csv

To run (on Discover):
    module load python/GEOSpyD
    python scripts/05_observed_nino34.py
"""

import os
import csv
import glob

import numpy as np
import xarray as xr

# ===== settings ============================================================
# the months to compute, newest last (YYYYMM)
MONTHS      = ["202511", "202512", "202601", "202602", "202603", "202604",
               "202605", "202606", "202607"]

REANALYSIS_DIR = "/discover/nobackup/projects/gmao/geos-s2s-3/GiOCEAN_e1/ocn_tavg_1mo_glo_L720x361_slv"
NRT_DIR        = "/gpfsm/dnb07/projects/p236/GiOcean-NRT/ocn_tavg_1mo_glo_L720x361_slv"
CLIM_YEARS     = range(1998, 2025)

# candidate variable names for sea-surface temperature, tried in order
SST_CANDIDATES = ["TSFG", "BULK_OCEANTEMP", "TS", "SST"]

# the Nino 3.4 box
LAT_MIN, LAT_MAX = -5.0, 5.0
LON_MIN, LON_MAX = -170.0, -120.0
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))


def month_file(tag):
    """The monthly ocean file for YYYYMM, from the main directory where
    present and the near-real-time stream otherwise."""
    for root in (REANALYSIS_DIR, NRT_DIR):
        matches = glob.glob(os.path.join(root, f"*.{tag}01_0000z.nc4"))
        if matches:
            return matches[0]
    return None


def sst_name(dataset):
    for name in SST_CANDIDATES:
        if name in dataset:
            return name
    raise SystemExit(f"no SST variable found; data_vars: "
                     f"{list(dataset.data_vars)}")


def box_mean(path):
    with xr.open_dataset(path) as data:
        field = data[sst_name(data)]
        if "time" in field.dims:
            field = field.isel(time=0)
        field = field.sel(lat=slice(LAT_MIN, LAT_MAX),
                          lon=slice(LON_MIN, LON_MAX))
        weights = np.cos(np.deg2rad(field["lat"]))
        return float(field.weighted(weights).mean(("lat", "lon")))


def main():
    # the 1998-2024 climatology of the box mean, per calendar month
    climatology = {}
    for month in range(1, 13):
        values = []
        for year in CLIM_YEARS:
            path = month_file(f"{year}{month:02d}")
            if path is not None and REANALYSIS_DIR in path:
                values.append(box_mean(path))
        if values:
            climatology[month] = float(np.mean(values))
            print(f"  climatology {month:02d}: {len(values)} years, "
                  f"{climatology[month]:.2f} K")

    rows = []
    for tag in MONTHS:
        month = int(tag[4:])
        path = month_file(tag)
        if path is None or month not in climatology:
            print(f"  {tag}: file or climatology missing, skipped")
            continue
        anomaly = box_mean(path) - climatology[month]
        print(f"  {tag}: {anomaly:+.2f} K")
        rows.append([tag, f"{anomaly:.4f}"])

    out = os.path.normpath(os.path.join(
        script_dir, "..", "outputs", "observed_recent", "nino34"))
    os.makedirs(out, exist_ok=True)
    csv_path = os.path.join(out, "observed_nino34.csv")
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["verifying", "nino34_anom_K"])
        writer.writerows(rows)
    print("  saved:", csv_path)
    print("\nDone.")


if __name__ == "__main__":
    main()
