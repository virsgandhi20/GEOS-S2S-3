"""
Teleconnection indices from GEOS-S2S forecasts.

For each level, initialization season and lead, this reads the forecasts of
geopotential height, takes each one's anomaly against the lead-dependent
forecast climatology, and fits it against the historical rotated patterns
saved by the GiOCEAN analysis. The result is one index per mode per forecast.

Two input sources are supported, each written to its own folder so the results
sit side by side:

    outputs/500hPa/DJF/20N_90N/lead1/members/forecast_indices.csv
    outputs/500hPa/DJF/20N_90N/lead1/ensmean/forecast_indices.csv

"members" processes the individual forecasts (every five-day start date, one
row per member); "ensmean" processes the pre-averaged product, one forecast per
initialization month, where the averaging over start dates and members has
already been done. Settings are at the top; the machinery lives in analysis.py.

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

# forecast inputs to process, each written to its own folder:
#   "members" - individual forecasts, every five-day start date
#   "ensmean" - the pre-averaged product, one forecast per initialization month
SOURCES     = ["members", "ensmean"]

# which ensemble members (the "members" source): "*" for all found, "1" for ens1
MEMBERS     = "1"

# where the forecasts live. GEOS-S2S-2 on Discover for now; point these at the
# GEOS-S2S-3 trees on NAS once access comes through.
DATA_ROOT    = "/discover/nobackup/projects/gmao/m2oasf/aogcm/g5fcst/forecast/production/geos-s2s/runx"
ENSMEAN_ROOT = "/gpfsm/dnb10/projects/p71/aogcm/g5fcst/forecast/production/geos-s2s/runx/ensmean"
COLLECTION   = "geosgcm_vis2d"

# the historical patterns to fit against (written by the GiOCEAN repo's
# regression display) and the region they were found on
DOMAIN      = "20N_90N"
PATTERNS    = "/discover/nobackup/vgandhi2/GiOCEAN/outputs/{level}hPa/{season}/{domain}/regression/patterns.nc"
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))


def out_path(level, season, lead, source):
    return os.path.normpath(os.path.join(
        script_dir, "..", "outputs", f"{int(level)}hPa", season, DOMAIN,
        f"lead{lead}", source, "forecast_indices.csv"))


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
                for source in SOURCES:
                    print(f"\n=== {variable}  {season} inits  lead {lead}  "
                          f"{source}  ({YEAR_START}-{YEAR_END}) ===")
                    if source == "ensmean":
                        records = analysis.find_ensemble_mean(
                            ENSMEAN_ROOT, years, months, lead,
                            collection=COLLECTION)
                    else:
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

                    # lead-dependent climatology: for each initialization (a
                    # start date for members, a month for ensmean), the mean
                    # over years
                    by_init = {}
                    for record in records:
                        by_init.setdefault(record["init"],
                                           []).append(record["field"])
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
                                           out_path(level, season, lead, source))
            patterns.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
