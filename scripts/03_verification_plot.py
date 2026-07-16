"""
Forecast verification figures, one per mode: the pattern map on top, and below
it the forecast index (black) against the historical GiOCEAN index (blue,
dashed) over the hindcast years, with their correlation printed on the panel.

The forecast curve is the average over everything verifying in a given month
(all start dates and members for the "members" source; the archive has already
averaged for "ensemble_mean"), so it plays the role of the ensemble-average curve in
the group's ENSO verification plots. The observed curve is the GiOCEAN rotated
index for the same verifying months. Only in-season months exist, so the
curves break between years.

The map panel needs the patterns file the indices were fitted against
(patterns.nc, which lives on Discover); where it is missing the time series
are drawn on their own.

Outputs, alongside the forecast CSVs:

    outputs/hindcasts/500hPa/DJF/20N_90N/lead1/members/projection/verification_REOF1.png ...

To run (on Discover):
    module load python/GEOSpyD
    python scripts/03_verification_plot.py
"""

import os
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

# ===== settings ============================================================
LEVELS      = [500, 250]
SEASONS     = ["DJF", "JJA"]
LEAD        = 1
SOURCES     = ["members", "ensemble_mean"]
DOMAIN      = "20N_90N"

# the index convention to verify: "projection" (the group's convention) or
# "least_squares". the observed indices are read in the same convention, so both
# sides of the comparison are computed the same way.
METHOD      = "projection"
_INDEX_FILE = {"projection": "projection_indices.csv",
               "least_squares": "teleconnection_indices.csv"}
N_PLOT      = 10         # how many leading modes to plot per case

# the GiOCEAN repository (for the observed indices, relative to the root of
# this repository) and its patterns file (for the map panel; skipped where
# the file is absent)
GIOCEAN_INDICES = os.path.join(
    "..", "GiOCEAN", "outputs", "{level}hPa", "{season}", DOMAIN,
    "regression", "{index_file}")
PATTERNS = "/discover/nobackup/vgandhi2/GiOCEAN/outputs/{level}hPa/{season}/{domain}/regression/patterns.nc"
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))


def forecast_csv(level, season, source):
    return os.path.normpath(os.path.join(
        script_dir, "..", "outputs", "hindcasts", f"{level}hPa", season,
        DOMAIN,
        f"lead{LEAD}", source, METHOD, "forecast_indices.csv"))


def read_forecast(path):
    """Average the forecast indices over everything verifying in the same
    month. Returns {verifying YYYYMM: array of indices per mode}."""
    by_month = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            by_month.setdefault(row[2], []).append(
                np.array([float(v) for v in row[3:]]))
    return {tag: np.mean(rows, axis=0) for tag, rows in by_month.items()}


def read_observed(path):
    """Read the GiOCEAN indices as {YYYYMM: array of indices per mode}."""
    table = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            table[row[0][:4] + row[0][5:7]] = np.array(
                [float(v) for v in row[1:]])
    return table


def monthly_axis(tags):
    """A continuous month axis spanning the tags, with NaN placeholders so the
    plotted line breaks across the months outside the season."""
    years = [int(t[:4]) for t in tags]
    start, end = min(years), max(years) + 1
    axis = [f"{y}{m:02d}" for y in range(start, end) for m in range(1, 13)]
    times = np.array([np.datetime64(f"{t[:4]}-{t[4:]}-01") for t in axis])
    return axis, times


def curve(axis, table, mode):
    return np.array([table[t][mode] if t in table else np.nan for t in axis])


def pattern_map(level, season, mode):
    """The mode's loading map from patterns.nc, or None where the file is not
    available (it stays on Discover)."""
    path = PATTERNS.format(level=level, season=season, domain=DOMAIN)
    if not os.path.exists(path):
        return None
    import xarray as xr
    with xr.open_dataset(path) as data:
        return data["loadings"].isel(mode=mode).load()


def _white_centre(levels):
    n_bins = len(levels) - 1
    base = plt.get_cmap("RdBu_r")
    shades = list(base(np.linspace(0.0, 1.0, n_bins)))
    shades[n_bins // 2 - 1] = shades[n_bins // 2] = (1.0, 1.0, 1.0, 1.0)
    return mcolors.ListedColormap(shades)


def plot_case(loading, times, forecast, observed, correlation, heading,
              out_png):
    """The map (where available) with the two half-period index panels
    underneath, forecast in black and observation in dashed blue."""
    n_rows = 3 if loading is not None else 2
    fig = plt.figure(figsize=(11, 3.2 * n_rows))
    ratios = ([2.0, 1.0, 1.0] if loading is not None else [1.0, 1.0])
    grid = fig.add_gridspec(n_rows, 1, height_ratios=ratios, hspace=0.45)
    row = 0

    if loading is not None:
        import cartopy.crs as ccrs
        import cartopy.util as cutil
        ax = fig.add_subplot(grid[row],
                             projection=ccrs.PlateCarree(central_longitude=180))
        filled, lon = cutil.add_cyclic_point(loading.values,
                                             coord=loading["lon"].values)
        top = np.nanpercentile(np.abs(filled), 99)
        contour_levels = np.linspace(-top, top, 21)
        shaded = ax.contourf(lon, loading["lat"].values, filled,
                             levels=contour_levels,
                             cmap=_white_centre(contour_levels),
                             transform=ccrs.PlateCarree(), extend="both")
        ax.coastlines(linewidth=0.5)
        fig.colorbar(shaded, ax=ax, orientation="horizontal", shrink=0.7,
                     pad=0.06, aspect=40)
        ax.set_title(heading)
        row += 1

    half = len(times) // 2
    for lo, hi in ((0, half), (half, len(times))):
        ax = fig.add_subplot(grid[row])
        ax.plot(times[lo:hi], forecast[lo:hi], "k-", linewidth=1.4)
        ax.plot(times[lo:hi], observed[lo:hi], "b--", linewidth=1.2)
        ax.axhline(0, color="limegreen", linewidth=0.8, linestyle=":")
        ax.set_ylim(-3, 3)
        ax.set_title(f"{LEAD}-month lead forecast (ensemble average, black) "
                     "& GiOCEAN (blue)", loc="left", fontsize=10)
        if row == (1 if loading is not None else 0):
            ax.set_title(f"Corr.={correlation:.2f}", loc="right",
                         fontweight="bold")
        ax.margins(x=0.01)
        row += 1

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  saved:", out_png)


def main():
    for level in LEVELS:
        for season in SEASONS:
            observed = read_observed(os.path.normpath(os.path.join(
                script_dir, "..",
                GIOCEAN_INDICES.format(level=level, season=season,
                                       index_file=_INDEX_FILE[METHOD]))))
            for source in SOURCES:
                path = forecast_csv(level, season, source)
                if not os.path.exists(path):
                    print(f"missing {path}, skipping")
                    continue
                print(f"\n=== {level}hPa {season} lead {LEAD} {source} ===")
                forecast = read_forecast(path)
                axis, times = monthly_axis(sorted(forecast))
                common = sorted(set(forecast) & set(observed))
                print(f"  {len(common)} common verifying months")
                for mode in range(N_PLOT):
                    f_curve = curve(axis, forecast, mode)
                    o_curve = curve(axis, observed, mode)
                    f = np.array([forecast[t][mode] for t in common])
                    o = np.array([observed[t][mode] for t in common])
                    correlation = float(np.corrcoef(f, o)[0, 1])
                    print(f"  REOF{mode+1}: corr {correlation:+.2f}")
                    plot_case(
                        pattern_map(level, season, mode), times,
                        f_curve, o_curve, correlation,
                        f"REOF{mode+1} H {level}hPa {season} ({source})",
                        os.path.join(os.path.dirname(path),
                                     f"verification_REOF{mode+1}.png"))
    print("\nDone.")


if __name__ == "__main__":
    main()
