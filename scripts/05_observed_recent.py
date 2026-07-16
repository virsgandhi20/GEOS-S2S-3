"""
Observed teleconnection indices for recent months, for the lead-in segment of
the plume figures.

Each recent GiOCEAN monthly map is turned into an index the same way a
forecast is: anomaly against the 1998-2024 mean for that calendar month, the
usual latitude weighting, then the index by the chosen method against the
fixed historical patterns. The patterns are never refit; recent months are
only measured against them.

Runs on Discover (the GiOCEAN data and patterns live there) and writes a
small CSV that is committed, so the plume script on PFE can draw the
observed segment:

    outputs/observed_recent/500hPa/20N_90N/projection/observed_indices.csv

To run (on Discover):
    module load python/GEOSpyD
    python scripts/05_observed_recent.py
"""

import os
import csv

import xarray as xr
import analysis

# ===== settings ============================================================
# the months to compute, newest last (YYYYMM)
MONTHS      = ["202511", "202512", "202601"]

VARIABLE    = "H"
LEVELS      = [500, 250]
DOMAIN      = "20N_90N"
METHOD      = "projection"     # match the plume's METHOD

REANALYSIS_DIR = "/discover/nobackup/projects/gmao/geos-s2s-3/GiOCEAN_e1/atm_inst_6hr_glo_L720x361_p49"
CLIM_YEARS     = range(1998, 2025)

PATTERNS    = "/discover/nobackup/vgandhi2/GiOCEAN/outputs/{level}hPa/{season}/{domain}/regression/patterns.nc"
SEASON_OF   = {12: "DJF", 1: "DJF", 2: "DJF",
               3: "MAM", 4: "MAM", 5: "MAM",
               6: "JJA", 7: "JJA", 8: "JJA",
               9: "SON", 10: "SON", 11: "SON"}
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))


def month_file(tag):
    import glob
    matches = glob.glob(os.path.join(REANALYSIS_DIR,
                                     f"*.monthly.{tag}.nc4"))
    return matches[0] if matches else None


def main():
    for level in LEVELS:
        rows = []
        for tag in MONTHS:
            v_month = int(tag[4:])
            season = SEASON_OF[v_month]
            pattern_file = PATTERNS.format(level=int(level), season=season,
                                           domain=DOMAIN)
            path = month_file(tag)
            if path is None or not os.path.exists(pattern_file):
                print(f"  {tag}: file or {season} patterns missing, skipped")
                continue
            patterns = xr.open_dataset(pattern_file)
            lat, lon = patterns["lat"], patterns["lon"]
            baseline = analysis.reanalysis_month_mean(
                REANALYSIS_DIR, VARIABLE, level, v_month, CLIM_YEARS,
                lat, lon)
            field = analysis.read_field(path, VARIABLE, lat, lon,
                                        level=level)
            anomaly = analysis.area_weight(field - baseline)
            indices = analysis.fit_indices(anomaly, patterns, method=METHOD)
            n_modes = patterns.sizes["mode"]
            patterns.close()
            print(f"  {level}hPa {tag} ({season}): " + "  ".join(
                f"REOF{i+1} {v:+.2f}" for i, v in enumerate(indices[:4])))
            rows.append([tag, season] + list(indices))

        if not rows:
            print(f"nothing produced for {level} hPa")
            continue
        out_csv = os.path.normpath(os.path.join(
            script_dir, "..", "outputs", "observed_recent",
            f"{int(level)}hPa", DOMAIN, METHOD, "observed_indices.csv"))
        os.makedirs(os.path.dirname(out_csv), exist_ok=True)
        with open(out_csv, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["verifying", "patterns"]
                            + [f"REOF{i+1}" for i in range(n_modes)])
            for row in rows:
                writer.writerow(row[:2] + [f"{v:.5f}" for v in row[2:]])
        print("  saved:", out_csv)
    print("\nDone.")


if __name__ == "__main__":
    main()
