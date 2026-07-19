"""
Forecast plume from the GEOS-S2S-3 near-real-time forecasts, one figure per
mode: every ensemble member's index as a thin dashed line over the target
months, and their average as a thick solid line (the layout of the GMAO
Nino 3.4 plume plots).

One initialization month is processed at a time. The archive holds an
initialization every five days with five members each, except the last
initialization of the month, which has fifteen (so a January start month
carries 45 members in total). All members from all start dates go into one
plume, and every available forecast month is plotted, including the partial
initialization month itself (a January 1 start contributes January onward).

The chain per member and lead: anomaly against the drift climatology for the
same initialization month and verifying month, the usual latitude weighting,
then the index by projection onto each standardized observed pattern
(METHOD = "projection", the per-pattern convention; "least_squares" gives
the simultaneous fit instead). Each lead is measured against the pattern set of
the season its verifying month falls in, so plume curves should be read
within a season.

The data lives on PFE (NAS):

    <NRT_ROOT>/YYYYMMDD/ens*/<collection>/YYYYMMDD.<collection>.monthly.YYYYMM.nc4

Only the "monthly" files are used, never the "dailymean" ones (the exact-name
match takes care of that). Outputs:

    outputs/plumes/2026jan/500hPa/20N_90N/projection/forecast_indices.csv
    outputs/plumes/2026jan/500hPa/20N_90N/projection/plume_REOF1.png ...

To run (on PFE):
    python scripts/01_plume_forecasts.py
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
import analysis

# ===== settings ============================================================
# the initialization month to process (all start dates inside it)
INIT_YEAR   = 2026
INIT_MONTH  = 1

# where the near-real-time forecasts and the drift climatology live (PFE).
# the collection is a plain string so any collection can be used.
NRT_ROOT    = "/nobackupp28/knakada/GEOSS2S3/GEOS_fcst4nrt"
DRIFT_ROOT  = "/nobackupp28/knakada/GEOSS2S3/GEOS_fcst/data/DRFT/DRFT_2001_2020/monthly/atm_inst_6hr_glo_L720x361_p49"
COLLECTION  = "atm_inst_6hr_glo_L720x361_p49"

VARIABLE    = "H"              # pressure-level height in these files
LEVELS      = [500, 250]
DOMAIN      = "20N_90N"

# how each index is computed: "projection" (each map onto one standardized
# observed pattern at a time) or "least_squares" (all patterns fitted at once)
METHOD      = "projection"

# historical patterns (written by the GiOCEAN repository). patterns.nc lives
# on Discover and is not in git, so running on PFE means copying the pattern
# files over first and pointing this at the copy.
PATTERNS    = "/nobackupp28/vgandhi2/patterns/{level}hPa/{season}/{domain}/regression/patterns.nc"
SEASON_OF   = {12: "DJF", 1: "DJF", 2: "DJF",
               3: "MAM", 4: "MAM", 5: "MAM",
               6: "JJA", 7: "JJA", 8: "JJA",
               9: "SON", 10: "SON", 11: "SON"}

N_PLOT      = 10               # modes to draw a plume for

# recent observed indices for the lead-in segment of the plume (written by
# 02_observed_recent.py on Discover and committed, since the GiOCEAN data
# does not exist on PFE). drawn where the file is present.
OBSERVED    = os.path.join("..", "outputs", "observed_recent",
                           "{level}hPa", "{domain}", "{method}",
                           "observed_indices.csv")
# ===========================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))
_MONTH = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
_TAG = re.compile(r"\.monthly\.(\d{6})\.nc4$")


def find_members():
    """Every (init date, member) pair for the initialization month, with the
    monthly file per verifying month. Returns records like
    {init, member, files: {YYYYMM tag: path}}."""
    records = []
    month_glob = os.path.join(NRT_ROOT,
                              f"{INIT_YEAR}{INIT_MONTH:02d}??")
    for init_dir in sorted(glob.glob(month_glob)):
        init = os.path.basename(init_dir)
        for member_dir in sorted(glob.glob(os.path.join(init_dir, "ens*"))):
            files = {}
            for path in sorted(glob.glob(os.path.join(
                    member_dir, COLLECTION,
                    f"{init}.{COLLECTION}.monthly.*.nc4"))):
                match = _TAG.search(os.path.basename(path))
                if match:
                    files[match.group(1)] = path   # every month, including
                                                   # the partial init month
            if files:
                records.append({"init": init,
                                "member": os.path.basename(member_dir),
                                "files": files})
    return records


def drift_file(v_month):
    """The drift climatology for this initialization month and verifying
    calendar month, or None where it does not exist (the drift sets cover
    nine months from the initialization, so the last target month of the
    longest forecast can fall outside them). The drift tree may or may not
    nest a month directory."""
    mon = analysis._MONTH_DIR[INIT_MONTH]
    name = f"{mon}.{COLLECTION}.monthly.drift.{v_month:02d}.nc4"
    for candidate in (os.path.join(DRIFT_ROOT, mon, name),
                      os.path.join(DRIFT_ROOT, name)):
        if os.path.exists(candidate):
            return candidate
    return None


def read_observed(level, mode):
    """The recent observed indices for the lead-in segment, as
    {YYYYMM: value}, or an empty dict where the file is absent."""
    path = os.path.normpath(os.path.join(script_dir, OBSERVED.format(
        level=int(level), domain=DOMAIN, method=METHOD)))
    if not os.path.exists(path):
        return {}
    table = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            table[row[0]] = float(row[2 + mode])
    return table


def plume_plot(months, groups, mean, observed, heading, out_png):
    """The layout of the GMAO plume plots: the observed index as a solid
    black lead-in, each start date's members dashed in one colour, and the
    ensemble mean as a heavy red line with markers."""
    obs_months = sorted(observed)
    axis = ([f"{_MONTH[int(t[4:])]}\n{t[:4]}" for t in obs_months]
            + list(months))
    pad = [np.nan] * len(obs_months)

    fig, ax = plt.subplots(figsize=(11, 6))
    if obs_months:
        ax.plot(axis[:len(obs_months)],
                [observed[t] for t in obs_months],
                "-", color="black", linewidth=2.5, label="GiOCEAN")
    colors = plt.get_cmap("tab10")
    for i, (init, curves) in enumerate(sorted(groups.items())):
        label = f"{_MONTH[int(init[4:6])].lower()}{init[6:]}"
        for j, values in enumerate(curves):
            ax.plot(axis, pad + list(values), "--", color=colors(i % 10),
                    linewidth=0.9, label=label if j == 0 else None)
    ax.plot(axis, pad + list(mean), "o-", color="red", linewidth=2.5,
            label="Ensmean")
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_ylim(-3, 3)
    ax.set_ylabel("index")
    ax.set_xlabel("Forecast Month")
    ax.set_title(heading)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(linewidth=0.3, alpha=0.5)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  saved:", out_png)


def main():
    records = find_members()
    if not records:
        raise SystemExit(f"no forecasts found under {NRT_ROOT} for "
                         f"{INIT_YEAR}-{INIT_MONTH:02d}")
    tags = sorted({tag for r in records for tag in r["files"]})
    print(f"{len(records)} members across "
          f"{len({r['init'] for r in records})} start dates; "
          f"target months {tags[0]}..{tags[-1]}")

    label = f"{INIT_YEAR}{analysis._MONTH_DIR[INIT_MONTH]}"
    for level in LEVELS:
        print(f"\n=== {VARIABLE} {level}hPa ===")
        # patterns and drift per verifying month, loaded once
        month_setup = {}
        for tag in tags:
            v_month = int(tag[4:])
            season = SEASON_OF[v_month]
            pattern_file = PATTERNS.format(level=int(level), season=season,
                                           domain=DOMAIN)
            if not os.path.exists(pattern_file):
                print(f"  {tag}: no {season} patterns, skipping that month")
                continue
            drift_path = drift_file(v_month)
            if drift_path is None:
                print(f"  {tag}: no drift climatology for month "
                      f"{v_month:02d}, skipping that month")
                continue
            patterns = xr.open_dataset(pattern_file).load()
            drift = analysis.read_field(drift_path, VARIABLE,
                                        patterns["lat"], patterns["lon"],
                                        level=level)
            month_setup[tag] = (patterns, drift, season)

        rows, indexed = [], {}
        for record in records:
            for tag, path in sorted(record["files"].items()):
                if tag not in month_setup:
                    continue
                patterns, drift, season = month_setup[tag]
                forecast = analysis.read_field(path, VARIABLE,
                                               patterns["lat"],
                                               patterns["lon"], level=level)
                anomaly = analysis.area_weight(forecast - drift)
                indices = analysis.fit_indices(anomaly, patterns,
                                               method=METHOD)
                key = (record["init"], record["member"])
                indexed.setdefault(key, {})[tag] = indices
                rows.append([record["init"], record["member"], tag, season]
                            + list(indices))
            print(f"  {record['init']} {record['member']}: "
                  f"{len(record['files'])} months")

        out = os.path.normpath(os.path.join(
            script_dir, "..", "outputs", "plumes", label,
            f"{int(level)}hPa", DOMAIN, METHOD))
        os.makedirs(out, exist_ok=True)
        n_modes = len(rows[0]) - 4
        with open(os.path.join(out, "forecast_indices.csv"), "w",
                  newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["init", "member", "verifying", "patterns"]
                            + [f"REOF{i+1}" for i in range(n_modes)])
            for row in rows:
                writer.writerow(row[:4] + [f"{v:.5f}" for v in row[4:]])
        print("  saved:", os.path.join(out, "forecast_indices.csv"))

        months_axis = [tag for tag in tags if tag in month_setup]
        month_names = [f"{_MONTH[int(t[4:])]}\n{t[:4]}" for t in months_axis]
        r_year, r_month = analysis.verifying_month(INIT_YEAR, INIT_MONTH, 1)
        for mode in range(min(N_PLOT, n_modes)):
            groups, everything = {}, []
            for key in sorted(indexed):
                curve = [indexed[key][t][mode]
                         if t in indexed[key] else np.nan
                         for t in months_axis]
                groups.setdefault(key[0], []).append(curve)
                everything.append(curve)
            mean = np.nanmean(np.array(everything, dtype=float), axis=0)
            plume_plot(
                month_names, groups, mean, read_observed(level, mode),
                f"REOF{mode+1} {VARIABLE} ({int(level)}hPa): "
                f"{r_year} {_MONTH[r_month].lower()} released forecast",
                os.path.join(out, f"plume_REOF{mode+1}.png"))

        for setup in month_setup.values():
            setup[0].close()
    print("\nDone.")


if __name__ == "__main__":
    main()
