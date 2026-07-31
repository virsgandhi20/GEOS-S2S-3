"""
Nino 3.4 forecast plume from the GEOS-S2S-3 near-real-time forecasts.

The same chain as the teleconnection plumes, applied to the tropical
Pacific: for every ensemble member and verifying month, the sea-surface
temperature anomaly against the archived drift climatology (same
initialization month and verifying month), area-averaged over the
Nino 3.4 box (5S-5N, 170W-120W). No pattern projection is involved; the
index is the box-mean anomaly itself, in kelvin.

The monthly ocean fields come from the ocn_tavg_1mo collection, which has
a drift climatology in the same archive as the height fields. The SST
variable name is auto-detected from a small candidate list, so the script
adapts if the collection names it TS or SST.

Outputs:

    outputs/plumes/2026jun/nino34/forecast_nino34.csv
    outputs/plumes/2026jun/nino34/plume_nino34.png

To run (on PFE):
    module load python3/3.11.5
    python3 scripts/04_nino34_plume.py
"""

import os
import csv
import glob
import re

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ===== settings ============================================================
# the initialization month to process (all start dates inside it)
INIT_YEAR   = 2026
INIT_MONTH  = 6

NRT_ROOT    = "/nobackupp28/knakada/GEOSS2S3/GEOS_fcst4nrt"
DRIFT_ROOT  = "/nobackupp28/knakada/GEOSS2S3/GEOS_fcst/data/DRFT/DRFT_2001_2020/monthly/ocn_tavg_1mo_glo_L720x361_slv"
COLLECTION  = "ocn_tavg_1mo_glo_L720x361_slv"

# candidate variable names for sea-surface temperature, tried in order
SST_CANDIDATES = ["TS", "SST", "ts", "sst", "TSKINW"]

# the Nino 3.4 box
LAT_MIN, LAT_MAX = -5.0, 5.0
LON_MIN, LON_MAX = -170.0, -120.0
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))
_MONTH = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
_MONTH_DIR = ["", "jan", "feb", "mar", "apr", "may", "jun",
              "jul", "aug", "sep", "oct", "nov", "dec"]
_TAG = re.compile(r"\.monthly\.(\d{6})\.nc4$")


def find_members():
    """Every (init date, member) pair for the initialization month, with the
    monthly ocean file per verifying month."""
    records = []
    month_glob = os.path.join(NRT_ROOT, f"{INIT_YEAR}{INIT_MONTH:02d}??")
    for init_dir in sorted(glob.glob(month_glob)):
        init = os.path.basename(init_dir)
        for member_dir in sorted(glob.glob(os.path.join(init_dir, "ens*"))):
            files = {}
            for path in sorted(glob.glob(os.path.join(
                    member_dir, COLLECTION,
                    f"{init}.{COLLECTION}.monthly.*.nc4"))):
                match = _TAG.search(os.path.basename(path))
                if match:
                    files[match.group(1)] = path
            if files:
                records.append({"init": init,
                                "member": os.path.basename(member_dir),
                                "files": files})
    return records


def drift_file(v_month):
    """The drift climatology for this initialization month and verifying
    calendar month, or None where it does not exist."""
    mon = _MONTH_DIR[INIT_MONTH]
    name = f"{mon}.{COLLECTION}.monthly.drift.{v_month:02d}.nc4"
    for candidate in (os.path.join(DRIFT_ROOT, mon, name),
                      os.path.join(DRIFT_ROOT, name)):
        if os.path.exists(candidate):
            return candidate
    return None


def sst_name(dataset):
    for name in SST_CANDIDATES:
        if name in dataset:
            return name
    raise SystemExit(f"no SST variable found; data_vars: "
                     f"{list(dataset.data_vars)}")


def box_mean(dataset, name):
    """Area-weighted mean of the variable over the Nino 3.4 box."""
    field = dataset[name]
    if "time" in field.dims:
        field = field.isel(time=0)
    field = field.sel(lat=slice(LAT_MIN, LAT_MAX),
                      lon=slice(LON_MIN, LON_MAX))
    weights = np.cos(np.deg2rad(field["lat"]))
    return float(field.weighted(weights).mean(("lat", "lon")))


def main():
    members = find_members()
    if not members:
        raise SystemExit("no forecasts found for the initialization month")
    tags = sorted({t for r in members for t in r["files"]})
    print(f"{len(members)} members across "
          f"{len({r['init'] for r in members})} start dates; "
          f"target months {tags[0]}..{tags[-1]}")

    # drift box means per verifying calendar month (shared by all members)
    drift_mean = {}
    for tag in tags:
        v_month = int(tag[4:])
        path = drift_file(v_month)
        if path is None:
            print(f"  {tag}: no drift climatology for month {v_month:02d}, "
                  f"skipping that month")
            continue
        with xr.open_dataset(path) as drift:
            drift_mean[tag] = box_mean(drift, sst_name(drift))
    usable = [t for t in tags if t in drift_mean]

    rows = []
    for record in members:
        values = {}
        for tag, path in record["files"].items():
            if tag not in drift_mean:
                continue
            with xr.open_dataset(path) as data:
                values[tag] = box_mean(data, sst_name(data)) - drift_mean[tag]
        print(f"  {record['init']} {record['member']}: {len(values)} months")
        for tag, value in sorted(values.items()):
            rows.append([record["init"], record["member"], tag,
                         f"{value:.4f}"])

    out = os.path.normpath(os.path.join(
        script_dir, "..", "outputs", "plumes",
        f"{INIT_YEAR}{_MONTH_DIR[INIT_MONTH]}", "nino34"))
    os.makedirs(out, exist_ok=True)
    csv_path = os.path.join(out, "forecast_nino34.csv")
    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["init", "member", "verifying", "nino34_anom_K"])
        writer.writerows(rows)
    print("  saved:", csv_path)

    # ---- the plume, in the layout of the other plume figures --------------
    table = {}
    for init, member, tag, value in rows:
        table.setdefault((init, member), {})[tag] = float(value)
    months = usable
    labels = [f"{_MONTH[int(t[4:])]}\n{t[:4]}" for t in months]
    mean = [np.mean([v[t] for v in table.values() if t in v])
            for t in months]

    fig, ax = plt.subplots(figsize=(11, 6))
    colors = plt.get_cmap("tab10")
    inits = sorted({k[0] for k in table})
    for i, init in enumerate(inits):
        label = f"{_MONTH[int(init[4:6])].lower()}{init[6:]}"
        first = True
        for (r_init, member), values in sorted(table.items()):
            if r_init != init:
                continue
            curve = [values.get(t, np.nan) for t in months]
            ax.plot(labels, curve, "--", color=colors(i % 10),
                    linewidth=0.9, label=label if first else None)
            first = False
    ax.plot(labels, mean, "o-", color="red", linewidth=2.5, label="Ensmean")
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_ylabel("Nino 3.4 SST anomaly (K)")
    ax.set_xlabel("Forecast Month")
    r_month = INIT_MONTH % 12 + 1
    r_year = INIT_YEAR + (1 if INIT_MONTH == 12 else 0)
    ax.set_title(f"GEOS-S2S-3 Nino 3.4 (5S-5N, 170W-120W): {r_year} "
                 f"{_MONTH[r_month].lower()} released forecast")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    png_path = os.path.join(out, "plume_nino34.png")
    fig.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("  saved:", png_path)
    print("\nDone.")


if __name__ == "__main__":
    main()
