"""
Teleconnection indices for one real forecast.

Takes a single forecast (ensemble-mean files, one per verifying month, in the
GiOCEAN-style pressure-level format) and turns each lead into teleconnection
indices: anomaly against a baseline, the usual latitude weighting, a
least-squares fit to the historical patterns, and scaling to the historical
index spread.

Two baselines are available (the BASELINE setting), each written to its own
folder:
  "drift"   - the GEOS-S2S-3 drift climatology (2001-2020), one file per
              initialization month and verifying month, which removes the
              model's drift with lead. This is the proper baseline.
  "giocean" - the GiOCEAN monthly climatology, the stopgap used before the
              drift files were identified; it ignores the drift.

Writes a CSV (one row per lead) and a small figure of the leading indices:

    outputs/2026feb/500hPa/20N_90N/drift/forecast_indices.csv
    outputs/2026feb/500hPa/20N_90N/drift/indices.png

To run (on Discover):
    module load python/GEOSpyD
    python scripts/02_realtime_indices.py
"""

import os
import glob
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analysis

# ===== settings ============================================================
# the forecast: where its ensemble-mean monthly files live, and when it starts
FORECAST_DIR = "/gpfsm/dnb33/pyadav1/Ylim_codes/PY/2026feb_fcst"
INIT_YEAR    = 2026
INIT_MONTH   = 2
LABEL        = "2026feb"       # output folder name
LEADS        = [1, 2, 3, 4, 5, 6, 7, 8]

VARIABLE     = "H"             # pressure-level height in these files
LEVELS       = [500, 250]
DOMAIN       = "20N_90N"

# historical patterns (written by the GiOCEAN repository) and which pattern season to
# use for each verifying month. MAM and SON need those seasons added to the
# GiOCEAN run first; leads whose pattern file is missing are skipped with a note.
PATTERNS     = "/discover/nobackup/vgandhi2/GiOCEAN/outputs/{level}hPa/{season}/{domain}/regression/patterns.nc"
SEASON_OF    = {12: "DJF", 1: "DJF", 2: "DJF",
                3: "MAM", 4: "MAM", 5: "MAM",
                6: "JJA", 7: "JJA", 8: "JJA",
                9: "SON", 10: "SON", 11: "SON"}

# the anomaly baseline: "drift" (the GEOS-S2S-3 drift climatology) or
# "giocean" (the reanalysis monthly mean, the earlier stopgap)
BASELINE     = "drift"

# the drift climatology: <init mon>/<init mon>.<collection>.monthly.drift.<verifying MM>.nc4
DRIFT_DIR    = "/gpfsm/dnb07/projects/p236/GEOSS2S3/GEOS_fcst/data/DRFT/DRFT_2001_2020/atm_inst_6hr_glo_L720x361_p49"
DRIFT_COLLECTION = "atm_inst_6hr_glo_L720x361_p49"

# the reanalysis the "giocean" baseline comes from (same data the patterns used)
REANALYSIS_DIR = "/discover/nobackup/projects/gmao/geos-s2s-3/GiOCEAN_e1/atm_inst_6hr_glo_L720x361_p49"
CLIM_YEARS     = range(1998, 2025)

N_PLOT       = 4               # modes shown in the figure
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))
_MONTH = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def out_dir(level):
    return os.path.normpath(os.path.join(
        script_dir, "..", "outputs", LABEL, f"{int(level)}hPa", DOMAIN,
        BASELINE))


def drift_file(v_month):
    mon = analysis._MONTH_DIR[INIT_MONTH]
    return os.path.join(DRIFT_DIR, mon,
                        f"{mon}.{DRIFT_COLLECTION}.monthly.drift."
                        f"{v_month:02d}.nc4")


def forecast_file(v_year, v_month):
    matches = glob.glob(os.path.join(
        FORECAST_DIR, f"*.monthly.{v_year}{v_month:02d}*.nc4"))
    return matches[0] if matches else None


def main():
    for level in LEVELS:
        rows, plotted = [], []
        for lead in LEADS:
            v_year, v_month = analysis.verifying_month(INIT_YEAR, INIT_MONTH, lead)
            season = SEASON_OF[v_month]
            pattern_file = PATTERNS.format(level=int(level), season=season,
                                           domain=DOMAIN)
            if not os.path.exists(pattern_file):
                print(f"  lead {lead} ({_MONTH[v_month]}): no {season} patterns "
                      f"yet, skipping (add {season} to the GiOCEAN run)")
                continue
            path = forecast_file(v_year, v_month)
            if path is None:
                print(f"  lead {lead} ({_MONTH[v_month]}): forecast file "
                      f"missing, skipping")
                continue

            patterns = xr.open_dataset(pattern_file)
            lat, lon = patterns["lat"], patterns["lon"]
            if BASELINE == "drift":
                baseline = analysis.read_field(
                    drift_file(v_month), VARIABLE, lat, lon, level=level)
            else:
                baseline = analysis.reanalysis_month_mean(
                    REANALYSIS_DIR, VARIABLE, level, v_month, CLIM_YEARS,
                    lat, lon)
            forecast = analysis.read_field(path, VARIABLE, lat, lon, level=level)
            anomaly = analysis.area_weight(forecast - baseline)
            indices = analysis.fit_indices(anomaly, patterns)
            n_modes = patterns.sizes["mode"]
            patterns.close()

            print(f"  lead {lead} ({_MONTH[v_month]} {v_year}, {season} "
                  f"patterns): " + "  ".join(
                      f"REOF{i+1} {v:+.2f}" for i, v in enumerate(indices[:N_PLOT])))
            rows.append([f"{v_year}-{v_month:02d}", str(lead), season]
                        + list(indices))
            plotted.append((lead, f"{_MONTH[v_month]}", indices))

        if not rows:
            print(f"nothing produced for {level} hPa")
            continue

        out = out_dir(level)
        os.makedirs(out, exist_ok=True)
        import csv
        with open(os.path.join(out, "forecast_indices.csv"), "w",
                  newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["verifying", "lead", "patterns"]
                            + [f"REOF{i+1}" for i in range(len(rows[0]) - 3)])
            for row in rows:
                writer.writerow(row[:3] + [f"{v:.5f}" for v in row[3:]])
        print("  saved:", os.path.join(out, "forecast_indices.csv"))

        fig, ax = plt.subplots(figsize=(9, 5))
        months = [name for _, name, _ in plotted]
        for i in range(min(N_PLOT, len(plotted[0][2]))):
            ax.plot(months, [vals[i] for _, _, vals in plotted],
                    marker="o", label=f"REOF{i+1}")
        ax.axhline(0, color="k", linewidth=0.5)
        ax.set_ylabel("index")
        ax.set_xlabel("verifying month")
        ax.set_title(f"{VARIABLE} ({int(level)}hPa) forecast indices, "
                     f"{_MONTH[INIT_MONTH]} {INIT_YEAR} start")
        ax.legend()
        fig.savefig(os.path.join(out, "indices.png"), dpi=150,
                    bbox_inches="tight")
        plt.close(fig)
        print("  saved:", os.path.join(out, "indices.png"))
    print("\nDone.")


if __name__ == "__main__":
    main()
