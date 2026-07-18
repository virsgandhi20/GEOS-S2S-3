"""
Helpers for turning GEOS-S2S-3 forecasts into teleconnection indices.

The idea: the GiOCEAN work saved its rotated patterns (patterns.nc, one per
level/season/region). Here each forecast height anomaly is measured against
those patterns, giving an index per mode for every forecast. The anomaly
baseline is the drift climatology read from the archive, one file per
initialization month and verifying month, so the model drift is removed
before the index is computed. The scripts set their options at the top and
call the functions here.
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


def fit_indices(weighted_anomaly, patterns, method="least_squares"):
    """Turn one weighted anomaly map into an index per pattern, on the
    historical scale. Two conventions:

    "least_squares" solves y = E b for all patterns at once, where E holds the
    pattern loadings as columns and y is the anomaly map, using only the
    points valid in both. Each coefficient is the amount of that pattern
    after accounting for the others; applied to a historical map this
    reproduces the historical rotated index exactly. Scaled with the saved
    index_mean / index_std.

    "projection" projects y onto each unit-length pattern one at a time (the
    per-pattern convention, as at CPC). Overlap between patterns is not
    separated. Scaled with proj_mean / proj_std, the statistics of the same
    projection applied to the historical maps, so forecast and observed
    indices share a scale."""
    n_modes = patterns.sizes["mode"]
    E = patterns["loadings"].values.reshape(n_modes, -1).T   # (space, modes)
    y = np.asarray(weighted_anomaly.values, dtype=float).ravel()
    good = np.isfinite(y) & np.all(np.isfinite(E), axis=1)
    if method == "projection":
        if "proj_mean" not in patterns:
            raise SystemExit(
                "patterns.nc has no projection statistics; rerun the GiOCEAN "
                "analysis (02_rotated_modes.py) to regenerate it")
        unit = E[good] / np.linalg.norm(E[good], axis=0)
        p = y[good] @ unit
        return ((p - patterns["proj_mean"].values)
                / patterns["proj_std"].values)
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
