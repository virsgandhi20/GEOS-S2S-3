# Methodology

This document describes how the forecast teleconnection indices are computed
and the reasoning behind each step. Each section names the function that
implements it, so the document can be read alongside the code. The settings
that control a run are at the top of each script.

## Overview

The GiOCEAN analysis identifies the recurring patterns of Northern Hemisphere
height variability (the North Atlantic Oscillation, the Pacific/North American
pattern, and so on) and saves them as data (`patterns.nc`). This repository
does not derive any new patterns. Each forecast height map is instead measured
against the saved patterns: a least-squares fit determines how much of each
pattern is present in the map. The result is one number per pattern per
forecast, a teleconnection index, on the same scale as the historical record.
A complicated forecast map is thereby reduced to a short set of interpretable
quantities, such as a positive NAO phase and a negative PNA phase, which can be
compared directly with observations and with the historical indices.

## The dataset and its consequences for the code

GEOS-S2S is a forecast system, and its archive is organised accordingly:

- A new forecast is initialized every five days. Under `runx/<year>/` the
  initialization dates appear as directories such as `dec02`, `dec07`, `feb10`.
- Lead month 1 is the month after the initialization month: a February start
  verifies March at lead 1, April at lead 2, and so on out to roughly nine
  months, with one monthly file per lead.
- Each initialization carries ensemble members (`ens1`, `ens2`, ...). The
  numbering varies between initializations, so `find_forecasts` discovers the
  members by glob rather than assuming a fixed set.
- The height fields `H500` and `H250` are provided directly in the
  `geosgcm_vis2d` collection, in metres, on the same half-degree grid as
  GiOCEAN. The fit therefore requires no regridding; `read_field` falls back
  to interpolation should a future dataset differ.

One property of forecast models shapes the anomaly step: drift. A model drifts
away from reality as the lead grows, so its typical lead-1 January differs
from its typical lead-7 January, and both differ from the observed January
climate. Anomalies must therefore be taken against the model's own behaviour
at the same initialization date and lead; otherwise the drift enters the
indices as a spurious signal.

## The hindcast calculation

Implemented by `01_forecast_indices.py`, one case at a time (a level, a set of
initialization months, a lead, a region). Let `F(d, y)` denote the forecast
height map initialized on date `d` (for example `dec27`) of year `y`, at the
fixed lead.

### 1. Climatology per initialization date and lead

```
C(d) = mean over years y of F(d, y)
```

`C(d)` is the model's typical prediction at this lead for forecasts started on
date `d`. Because the average is taken at a fixed initialization date and
lead, the drift is contained in the baseline and cancels exactly when
subtracted. (In the code: the `by_init` grouping and the `climatology`
dictionary.)

### 2. Anomaly

```
A(d, y) = F(d, y) - C(d)
```

The forecast's departure from its own norm. By construction, the anomalies for
each initialization date sum to zero over the years, so the mean index over
all forecasts in a case must equal zero to rounding. The computed value is
0.000 over 448 winter forecasts, which follows from the algebra rather than
from chance and serves as a check that the implementation matches it.

### 3. Area weighting

```
A~ = A * sqrt(cos(latitude))     (analysis.area_weight)
```

Grid cells shrink towards the pole, and the square root of the cosine restores
equal-area contributions once quantities are squared downstream. Equally
important is consistency: the historical patterns were derived in this
weighted space, so any field fitted against them must receive the same
transformation first. The weighting is applied before the fit in all cases.

### 4. The least-squares fit

The rotated patterns are loaded from `patterns.nc`. Each pattern map is
flattened into a column, and the ten columns are stacked into a matrix `E`
(grid points by modes). The weighted anomaly is flattened into a vector `y`,
with grid points that are missing in either (heights below ground) removed
from both. The system

```
y = E b + error
```

is solved for `b` by least squares (`analysis.fit_indices`, via
`numpy.linalg.lstsq`), which minimises `|| y - E b ||^2`. The solution is the
combination of the ten known patterns that best reconstructs the forecast map:
`b1` measures how much of pattern 1 is present, `b2` how much of pattern 2,
and so on. All ten coefficients are fitted simultaneously. This matters
because the rotated patterns are not exactly orthogonal on the masked region;
the simultaneous fit is a multiple regression rather than ten independent
projections, and it accounts for the overlap between patterns.

### 4b. The projection index (the per-pattern convention)

An alternative to the simultaneous fit is used for the real-time products:
the weighted anomaly is projected onto each pattern one at a time, with the
pattern normalized to unit length. The result measures how much the map
resembles that single pattern, without separating the overlap between
patterns. The two conventions coincide for orthogonal patterns and differ
mildly for rotated ones; the projection is the convention used by CPC and by
the group's earlier forecast work, so the forecast indices become directly
comparable with those. The `method` argument of `fit_indices` selects
between them ("lstsq" and "projection"), and each writes to its own folder.

### 5. Scaling to the historical record

`patterns.nc` also carries the mean and standard deviation of each historical
index. The final index is

```
I_j = (b_j - mean_j) / std_j
```

The projection index is scaled the same way, with the mean and spread of the
same projection applied to the historical maps (`proj_mean` and `proj_std`
in the pattern file), so both conventions place forecast and observation on
one scale.

This places every forecast on the historical scale: an index of -1.7 in a 2009
forecast has the same meaning as -1.7 in the GiOCEAN record, namely 1.7
historical standard deviations into the negative phase. The historical mean is
close to zero and the spread close to one by construction, so the step is
nearly the identity, but applying it exactly removes any residual difference
in convention.

## Validity on data the patterns have not seen

Two properties support applying the fixed patterns to new data:

- **Consistency.** Applied to a historical GiOCEAN map, this procedure
  recovers the historical rotated indices exactly, because those indices are
  themselves the least-squares solution of the same equation; that is how
  rotated principal component analysis defines them. This was verified
  numerically on synthetic data (agreement to 1e-10). The forecast indices are
  therefore the out-of-sample extension of the historical ones: the same
  formula applied to new data.
- **No refitting.** The patterns, means and spreads are frozen from the
  reanalysis. Forecasts are only measured against them and never used to
  define anything, so there is no leakage, and the yardstick applies equally
  to 1981, to 2026, or to any other year.

One limitation should be stated: least squares represents only the part of a
forecast map that projects onto the ten patterns. Structure orthogonal to all
of them remains in the residual and is not summarised by the indices. This is
inherent to what an index is, but it bounds what the indices can be taken to
say about a forecast.

## The real-time case

`02_realtime_indices.py` runs the same chain on a single live forecast (the
February 2026 ensemble mean, which arrives in the GiOCEAN-style pressure-level
format, so the level is selected through the `lev` coordinate rather than read
as a 2D field). It differs from the hindcast processing in two respects:

- **Climatology.** The anomaly is taken against the GEOS-S2S-3 drift
  climatology (2001-2020), one file per initialization month and verifying
  month, which removes the model's drift with lead just as the per-init
  climatology does for the hindcasts. The GiOCEAN monthly mean remains
  available as an alternative baseline (the `BASELINE` setting); it was the
  stopgap used before the drift files were identified, and it ignores the
  drift. Each baseline writes to its own folder.
- **Season matching.** Each lead verifies in a different month, so each lead
  is fitted against the pattern set for the season containing its verifying
  month (the `SEASON_OF` mapping: MAM patterns for a March map, and so on).
  This has an interpretation consequence: REOF2 in May belongs to the spring
  pattern set and REOF2 in June to the summer set, two different physical
  structures. Index curves should be read within a season; a jump across a
  season boundary partly reflects the change of pattern set rather than the
  forecast itself. The `patterns` column in the CSV records the set used for
  each row.

## The forecast plume

`04_plume_s2s3.py` processes one initialization month of the GEOS-S2S-3
near-real-time archive on PFE. That archive holds an initialization every
five days, each carrying five ensemble members; a subset of members
(including several on the month's final start date, which carries fifteen)
extends to nine lead months, while the rest stop at three. A January start
month therefore holds 45 members in total, and all of them enter one plume.

Each member's monthly maps are turned into indices exactly as above: anomaly
against the drift climatology for the same initialization month and verifying
month, weighting, then the projection index. Target months without a drift
file are skipped. The figure follows the layout of the GMAO Nino 3.4 plume
plots: the recent observed index as a solid black lead-in (computed by
`05_observed_recent.py`, which measures recent GiOCEAN months against the
fixed patterns on Discover), each member dashed in its start date's colour,
and the ensemble mean as a heavy red line. The mean at each target month is
taken over the members that reach it, so it rests on all 45 members for the
first two target months and on the extended subset beyond.

Each lead is measured against the pattern set of the season its verifying
month falls in, so a plume crossing a season boundary changes yardstick
there; curves should be read within a season.

## Outputs

```
outputs/<level>hPa/<season>/<region>/lead<L>/<source>/<method>/forecast_indices.csv   hindcasts
outputs/<level>hPa/<season>/<region>/lead<L>/<source>/<method>/verification_REOF<n>.png  their verification figures
outputs/<label>/<level>hPa/<region>/<baseline>/<method>/forecast_indices.csv  one live forecast
outputs/<label>/<level>hPa/<region>/<baseline>/<method>/indices.png           its summary figure
outputs/plumes/<label>/<level>hPa/<region>/<method>/plume_REOF<n>.png       the forecast plumes
outputs/observed_recent/<level>hPa/<region>/<method>/observed_indices.csv   the plume's observed lead-in
```

Hindcasts are produced from two input sources, written to separate folders so
the results can be compared: `members` processes the individual forecasts (one
per five-day start date, one row per member), while `ensmean` processes the
archive's pre-averaged product, one forecast per initialization month.
Averaging over members retains the predictable part of a forecast and damps
the unpredictable part, so the ensemble-mean indices are smoother and fewer;
the per-member indices are noisier but show the spread between forecasts. The
climatology grouping adapts to the source: per start date for `members`, per
initialization month for `ensmean`.

The verification figures (`03_verification_plot.py`) place each mode's
pattern map above the forecast index and the historical GiOCEAN index over
the hindcast years, with their correlation printed. The forecast curve is the
average over everything verifying in a given month; the correlation runs over
the months both records cover (the GiOCEAN seasonal index has no March, so
lead-1 winter forecasts verify against January and February only).

A hindcast row records the initialization date, the member, the verifying
month, and the ten indices for that forecast. All values are in standard
deviations on the historical scale and are directly comparable with the
historical index files in the GiOCEAN repository. Each live forecast is
written to its own dated folder (the `LABEL` setting), so successive months
accumulate side by side rather than overwriting one another.

## Checks

- On synthetic data, the fit recovers known indices from maps built with them,
  handles missing grid points, and scales linearly (doubling the anomaly
  doubles the index).
- The lead arithmetic wraps the year correctly (a December start at lead 1
  verifies in January of the following year).
- The mean index over all 448 winter hindcasts is 0.000, as the anomaly
  construction requires.
- The winter of 2009/10, the most negative NAO winter in the record, appears
  clearly in the hindcast indices, including the transition from positive to
  strongly negative values as the initialization dates approached the event.

## Reproducing

```bash
module load python/GEOSpyD
python scripts/00_inspect_data.py        # survey the archive
python scripts/01_forecast_indices.py    # hindcast indices
python scripts/02_realtime_indices.py    # the live forecast
```

The GiOCEAN repository must have been run with the regression display first,
since that step writes `patterns.nc`.
