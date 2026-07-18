# Methodology

This document describes how the GEOS-S2S-3 forecast teleconnection indices
and the plume figures are computed, with the reasoning behind each step. Each
section names the function that implements it, so the document can be read
alongside the code. The settings that control a run are at the top of each
script.

## Overview

The GiOCEAN analysis identifies the recurring patterns of Northern Hemisphere
height variability (the North Atlantic Oscillation, the Pacific/North
American pattern, and so on) and saves them as data (`patterns.nc`). This
repository does not derive any new patterns. Each forecast height map is
measured against the saved patterns, giving one teleconnection index per
pattern per forecast, on the same scale as the historical record. A
complicated forecast map is thereby reduced to a short set of interpretable
quantities, such as a positive NAO phase and a negative PNA phase.

The method was developed and verified on the GEOS-S2S-2 hindcast archive
(the GEOS-S2S-2 repository), where the indices could be checked against the
observed record; the leading winter modes verify at about 0.6 correlation at
one month lead. The identical chain runs here on the GEOS-S2S-3
near-real-time forecasts.

## The dataset and its consequences for the code

The GEOS-S2S-3 near-real-time archive (PFE) holds an initialization every
five days, each carrying five ensemble members; a subset of members,
including several on the month's final start date (which carries fifteen),
extends to nine lead months while the rest stop at three. Every run includes
the partial initialization month itself, and all available forecast months
are processed. The height field `H` sits on pressure levels, on the same
half-degree grid as GiOCEAN, so no regridding is needed (`read_field` falls
back to interpolation should a dataset differ).

## The calculation, per member and forecast month

### 1. The drift climatology, read from the archive

Forecast models drift: the model's typical lead-1 January differs from its
typical lead-7 January, and both differ from the observed January climate.
Anomalies must therefore be taken against the model's own behaviour for the
same initialization month and lead. The archive provides this as a
precomputed drift climatology (2001-2020), one file per initialization
month and verifying month; the code reads these files (`drift_file`) and
does not compute the climatology itself. Forecast months without a drift
file are skipped with a note.

### 2. The anomaly

```
A = F - drift
```

The forecast's departure from the model drift. Subtracting the drift file
removes the model drift from every forecast before any index is computed.

### 3. Area weighting

```
A~ = A * sqrt(cos(latitude))     (analysis.area_weight)
```

Grid cells shrink towards the pole, and the square root of the cosine
restores equal-area contributions. Equally important is consistency: the
historical patterns were derived in this weighted space, so any field
measured against them must receive the same transformation first.

### 4. The index, against standardized patterns

Two conventions are implemented (`analysis.fit_indices`), and each writes to
its own folder:

- **Projection** (the default, the group's convention): the weighted anomaly
  is projected onto each pattern one at a time. The patterns are
  standardized: each is scaled to unit length before the projection, and
  the raw projection is then standardized with the mean and spread of the
  same projection applied to the historical maps (`proj_mean` and
  `proj_std` in the pattern file). The result measures how much the map
  resembles that single pattern.
- **Least squares**: all patterns are fitted at once (`y = E b + error`),
  so each coefficient is the amount of that pattern after accounting for
  the others. Scaled with `index_mean` and `index_std`.

In both conventions the final index is standardized on the historical
scale: a value of -1.5 means 1.5 historical standard deviations into the
negative phase, directly comparable with the observed record.

### 5. Season matching

Each forecast month is measured against the pattern set of the season it
falls in (the `SEASON_OF` mapping: MAM patterns for a March map, and so
on). REOF2 in May belongs to the spring pattern set and REOF2 in June to
the summer set, two different physical structures, so index curves should
be read within a season; a jump across a season boundary partly reflects
the change of pattern set rather than the forecast itself.

## Validity on data the patterns have not seen

- **No refitting.** The patterns, means and spreads are frozen from the
  reanalysis. Forecasts are only measured against them and never used to
  define anything, so there is no leakage, and the yardstick applies
  equally to any year.
- **Limitation.** An index represents only the part of a map that resembles
  its pattern; structure unlike all ten patterns is not summarised by the
  indices.

## The forecast plume

`01_plume_forecasts.py` processes one initialization month at a time: every
member from every start date, every available forecast month (45 members
for a January start month). Each member's maps become indices exactly as
above, and the figure follows the layout of the GMAO Nino 3.4 plume plots:

- the recent observed index as a solid black lead-in, computed by
  `02_observed_recent.py`, which measures recent GiOCEAN months against the
  fixed patterns with the same convention
- each member as a dashed line in its start date's colour
- the ensemble mean as a heavy red line, averaged at each forecast month
  over the members that reach it (all 45 for the first months, the
  extended subset beyond)

## The single-forecast case

`03_single_forecast.py` runs the same chain on one forecast delivered as a
set of ensemble-mean files (the February 2026 case). The anomaly baseline
is a setting: the drift climatology (the proper choice), or the GiOCEAN
monthly mean (an earlier alternative that ignores the model drift). Each
baseline and convention writes to its own folder.

## Outputs

```
outputs/plumes/<label>/<level>hPa/<region>/<method>/plume_REOF<n>.png
outputs/plumes/<label>/<level>hPa/<region>/<method>/forecast_indices.csv
outputs/observed_recent/<level>hPa/<region>/<method>/observed_indices.csv
outputs/realtime/<label>/<level>hPa/<region>/<baseline>/<method>/...
```

A plume CSV records the initialization date, the member, the verifying
month, the pattern season, and the ten indices, all in standard deviations
on the historical scale. Each initialization month writes to its own dated
folder, so successive months accumulate side by side.

## Checks

- The index computation recovers known indices from synthetic maps, handles
  missing grid points, and scales linearly.
- The lead arithmetic wraps the year correctly (a December start verifies
  January of the following year at lead 1).
- The verification of the same chain on the GEOS-S2S-2 hindcasts (that
  repository) is the evidence that the indices carry forecast skill.

## Reproducing

```bash
python scripts/01_plume_forecasts.py     # the forecast plume (PFE)
python scripts/02_observed_recent.py     # its observed lead-in (Discover)
python scripts/03_single_forecast.py     # one single forecast (Discover)
```

The GiOCEAN repository must have been run with the regression display first,
since that step writes `patterns.nc`; running on PFE requires a copy of the
pattern files, which stay on the machine that produced them.
