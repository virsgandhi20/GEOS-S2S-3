"""
Teleconnection indices from GEOS-S2S forecasts.

For each level, initialization season and lead, this reads every forecast of
geopotential height, takes its anomaly against the lead-dependent forecast
climatology, and fits it against the historical rotated patterns saved by the
GiOCEAN analysis. The result is one index per mode per forecast, written to a
CSV organised the same way as the other repos:

    outputs/500hPa/DJF/20N_90N/lead1/forecast_indices.csv

Settings are at the top; the machinery lives in analysis.py.

To run (on Discover):
    module load python/GEOSpyD
    python scripts/01_forecast_indices.py
"""

import os
import numpy as np
import xarray as xr
import analysis

# ===== settings ============================================================
LEVELS      = [500, 250]       # uses the ready-made H500 / H250 fields
YEAR_START  = 1991             # initialization years to include
YEAR_END    = 2023

# initialization seasons: (init months, label). lead 1 is the month after the
# init month, so DJF inits at lead 1 verify in JFM.
SEASONS     = [([12, 1, 2], "DJF"),
               ([6, 7, 8],  "JJA")]
LEADS       = [1]              # lead months to process, e.g. [1, 2, 3]

# which ensemble members: "*" for all found, "1" for just ens1
MEMBERS     = "1"

# where the forecasts live. GEOS-S2S-2 on Discover for now; point this at the
# GEOS-S2S-3 tree on NAS once access comes through.
DATA_ROOT   = "/discover/nobackup/projects/gmao/m2oasf/aogcm/g5fcst/forecast/production/geos-s2s/runx"
COLLECTION  = "geosgcm_vis2d"

# the historical patterns to fit against (written by the GiOCEAN repo's
# regression display) and the region they were found on
DOMAIN      = "20N_90N"
PATTERNS    = "/discover/nobackup/vgandhi2/GiOCEAN/outputs/{level}hPa/{season}/{domain}/regression/patterns.nc"
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))


def out_path(level, season, lead):
    return os.path.normpath(os.path.join(
        script_dir, "..", "outputs", f"{int(level)}hPa", season, DOMAIN,
        f"lead{lead}", "forecast_indices.csv"))


def main():
    years = range(YEAR_START, YEAR_END + 1)
    for level in LEVELS:
        variable = f"H{int(level)}"
        for months, season in SEASONS:
            pattern_file = PATTERNS.format(level=int(level), season=season,
                                           domain=DOMAIN)
            patterns = xr.open_dataset(pattern_file)
            lat, lon = patterns["lat"], patterns["lon"]
            n_modes = patterns.sizes["mode"]

            for lead in LEADS:
                print(f"\n=== {variable}  {season} inits  lead {lead}  "
                      f"({YEAR_START}-{YEAR_END}, ens {MEMBERS}) ===")
                records = analysis.find_forecasts(
                    DATA_ROOT, years, months, lead,
                    collection=COLLECTION, members=MEMBERS)
                if not records:
                    print("  no forecast files found, skipping")
                    continue
                print(f"  {len(records)} forecasts found; reading...")

                for count, record in enumerate(records, 1):
                    record["field"] = analysis.read_field(
                        record["path"], variable, lat, lon)
                    if count % 100 == 0:
                        print(f"  ...{count}/{len(records)}")

                # lead-dependent climatology: mean over years (and members)
                # for each initialization date
                by_init = {}
                for record in records:
                    by_init.setdefault(record["init"], []).append(record["field"])
                climatology = {init: xr.concat(fields, "case").mean("case")
                               for init, fields in by_init.items()}
                for init, fields in sorted(by_init.items()):
                    print(f"  climatology {init}: {len(fields)} forecasts")

                rows = []
                for record in records:
                    anomaly = analysis.area_weight(
                        record["field"] - climatology[record["init"]])
                    indices = analysis.fit_indices(anomaly, patterns)
                    rows.append([f"{record['year']}-{record['init']}",
                                 record["member"], record["verifying"]]
                                + list(indices))
                analysis.write_indices(rows, n_modes,
                                       out_path(level, season, lead))
            patterns.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
