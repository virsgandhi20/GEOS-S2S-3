"""
Helpers for turning GEOS-S2S forecasts into teleconnection indices.

The idea: the GiOCEAN work saved its rotated patterns (patterns.nc, one per
level/season/region). Here each forecast height anomaly is fitted against those
patterns by least squares, giving an index per mode for every forecast. The
scripts set their options at the top and call the functions here.

Forecast quirks these functions absorb:
  - initializations sit in directories named like feb10 or nov07, one every
    five days, under runx/<year>/
  - ensemble member numbering varies by initialization (some have ens1, some
    ens2..ens5), so members are discovered by glob, never assumed
  - a forecast's climatology depends on the initialization date and the lead
    (models drift with lead), so anomalies are taken against the mean over
    years for the same init date and lead
"""

import os
import csv
import glob

import numpy as np
import xarray as xr

_MONTH_DIR = ["", "jan", "feb", "mar", "apr", "may", "jun",
              "jul", "aug", "sep", "oct", "nov", "dec"]


def verifying_month(init_year, init_month, lead):
    """The (year, month) a forecast verifies in: lead 1 is the month after the
    initialization month."""
    total = init_month + lead
    return init_year + (total - 1) // 12, (total - 1) % 12 + 1


def find_forecasts(root, years, init_months, lead,
                   collection="geosgcm_vis2d", members="*"):
    """Find the monthly forecast files for the given initialization months and
    lead. Returns a list of records with the init date, member and file path.

    members: glob for the member number, "*" for every member found, "1" for
    just ens1, "[2-5]" for ens2 to ens5, and so on."""
    records = []
    for year in years:
        for month in init_months:
            v_year, v_month = verifying_month(year, month, lead)
            tag = f"{v_year}{v_month:02d}"
            pattern = os.path.join(
                str(root), str(year), f"{_MONTH_DIR[month]}??",
                f"ens{members}", collection,
                f"*.{collection}.monthly.{tag}.nc4")
            for path in sorted(glob.glob(pattern)):
                parts = path.split(os.sep)
                records.append({"year": year,
                                "init": parts[-4],      # e.g. feb10
                                "member": parts[-3],    # e.g. ens1
                                "verifying": tag,
                                "path": path})
    return records


def read_field(path, variable, lat, lon, level=None):
    """Read one field and put it on the pattern grid. Give `level` (hPa) when
    the file carries the variable on pressure levels (the GEOS-S2S-3 style
    files); leave it out for ready-made 2D fields like H500 in vis2d.

    The GEOS-S2S atmospheric grids match the GiOCEAN grid exactly, so this is
    normally just a subset; if a grid ever differs (ERA5, for example), it
    falls back to interpolation."""
    with xr.open_dataset(path) as data:
        field = data[variable]
        if level is not None:
            field = field.sel(lev=level, method="nearest")
        field = field.isel(time=0).load()
    sub = field.sel(lat=slice(float(lat.min()), float(lat.max())),
                    lon=slice(float(lon.min()), float(lon.max())))
    if (sub.lat.size == lat.size and sub.lon.size == lon.size
            and np.allclose(sub.lat, lat) and np.allclose(sub.lon, lon)):
        return sub
    return field.interp(lat=lat, lon=lon)


def reanalysis_month_mean(data_dir, variable, level, month, years, lat, lon):
    """The long-term mean of one calendar month from the reanalysis the
    patterns were built on, on the pattern grid. Used as the baseline the
    real-time forecast anomaly is taken against."""
    fields = []
    for year in years:
        for path in sorted(glob.glob(os.path.join(
                data_dir, f"*.monthly.{year}{month:02d}.nc4"))):
            fields.append(read_field(path, variable, lat, lon, level=level))
    if not fields:
        raise SystemExit(f"no reanalysis files found for month {month:02d} "
                         f"under {data_dir}")
    return xr.concat(fields, "case").mean("case")


def area_weight(field):
    """The same sqrt(cos(latitude)) weighting the patterns were found with."""
    coslat = np.cos(np.deg2rad(field.lat.values)).clip(0.0, 1.0)
    return field * np.sqrt(coslat)[:, np.newaxis]


def fit_indices(weighted_anomaly, patterns):
    """Fit one weighted anomaly map against the historical patterns by least
    squares and put the result on the historical index scale.

    Solves y = E b for b, where E holds the pattern loadings as columns and y
    is the anomaly map, using only the points that are valid in both. b is then
    centred and scaled with the historical index mean and spread saved in the
    pattern file."""
    n_modes = patterns.sizes["mode"]
    E = patterns["loadings"].values.reshape(n_modes, -1).T   # (space, modes)
    y = np.asarray(weighted_anomaly.values, dtype=float).ravel()
    good = np.isfinite(y) & np.all(np.isfinite(E), axis=1)
    b, *_ = np.linalg.lstsq(E[good], y[good], rcond=None)
    return ((b - patterns["index_mean"].values)
            / patterns["index_std"].values)


def write_indices(rows, n_modes, out_csv):
    """Write one row per forecast: init date, member, verifying month, and the
    index for each mode."""
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["init", "member", "verifying"]
                        + [f"REOF{i+1}" for i in range(n_modes)])
        for row in rows:
            writer.writerow(row[:3] + [f"{v:.5f}" for v in row[3:]])
    print("  saved:", out_csv)
