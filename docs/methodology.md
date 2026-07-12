# Methodology

How the forecast indices are computed and why each step is there. This is meant
to be readable by someone who has the code open next to it: every section names
the function that implements it. The settings that control a run sit at the top
of each script.

## The idea in one paragraph

The GiOCEAN work found the recurring patterns of Northern Hemisphere height
variability (the NAO, the PNA and so on) and saved them as data (`patterns.nc`).
This repo never finds new patterns. It takes forecast height maps and asks, for
each one: how much of each known pattern is in this map? The answer is one
number per pattern per forecast, a teleconnection index, on the same scale as
the historical record. Instead of "the model predicts this complicated map for
June", you get "the model predicts the NAO slightly positive and the PNA
negative", which can be reasoned about and verified.

## The data (and why the code is shaped the way it is)

GEOS-S2S is a forecast system, so the archive is organised very differently
from a reanalysis:

- A new forecast starts every five days. Under `runx/<year>/` the start dates
  are directories like `dec02`, `dec07`, `feb10`.
- Lead month 1 is the month after the start: a February start verifies March at
  lead 1, April at lead 2, out to about nine months. One monthly file per lead.
- Each start has ensemble members (`ens1`, `ens2`, ...). The numbering varies
  between starts, so `find_forecasts` discovers members with a glob and never
  assumes them.
- The height fields `H500` and `H250` come ready-made in the `geosgcm_vis2d`
  collection, in metres, on the same half-degree grid as GiOCEAN. That grid
  match means the fit needs no regridding; `read_field` still falls back to
  interpolation in case a future dataset differs.

One property of forecast models drives the whole anomaly step: **drift**. A
model slowly wanders from reality as the lead grows, so its typical lead-1
January differs from its typical lead-7 January, and both differ from the real
January. Any anomaly has to be taken against the model's own habits at that
same start date and lead, or the drift leaks into the indices as a fake signal.

## The hindcast calculation, step by step

Implemented by `01_forecast_indices.py`, one case at a time (a level, a set of
initialization months, a lead, a region). Write `F(d, y)` for the forecast
height map started on date `d` (say `dec27`) of year `y`, at the fixed lead.

### 1. Climatology per start date and lead

```
C(d) = mean over years y of F(d, y)
```

This is "what the model typically predicts, one month out, when started on
December 27". Averaging over years at a fixed `(d, lead)` makes the drift part
of the baseline, so subtracting it removes the bias exactly. (In the code: the
`by_init` grouping and the `climatology` dictionary.)

### 2. Anomaly

```
A(d, y) = F(d, y) - C(d)
```

The forecast's departure from its own norm; the signal. A useful built-in
check: for each start date the anomalies sum to zero over the years, so the
mean index over all forecasts in a case must come out 0 to rounding. It does
(0.000 over 448 winter forecasts), which is a construction guarantee rather
than luck, and a quick way to confirm the code matches the algebra.

### 3. Area weighting

```
A~ = A * sqrt(cos(latitude))     (analysis.area_weight)
```

Two reasons. Grid cells shrink towards the pole, and the square root of the
cosine restores equal-area contributions once things get squared downstream.
More importantly, the historical patterns were derived in this weighted space,
so anything fitted against them must be transformed the same way first.
Weight, then fit; never the other way round.

### 4. The least-squares fit (the heart of it)

Load the rotated patterns from `patterns.nc`, flatten each pattern map into a
column, and stack the ten columns into a matrix `E` (grid points by modes).
Flatten the weighted anomaly into a vector `y`, dropping the grid points that
are missing in either (heights below ground). Then solve

```
y = E b + error        =>       b = argmin || y - E b ||^2
```

by least squares (`analysis.fit_indices`, using `numpy.linalg.lstsq`). In
words: find the mix of the ten known patterns that best reconstructs this
forecast map. `b1` is how much of pattern 1 is present, `b2` how much of
pattern 2, and so on. All ten are fitted at once, which matters because the
rotated patterns are not exactly orthogonal on the masked region; fitting them
simultaneously is a multiple regression, not ten separate projections.

### 5. Scaling to the historical yardstick

`patterns.nc` also carries the mean and standard deviation of each historical
index. The final index is

```
I_j = (b_j - mean_j) / std_j
```

This puts every forecast on exactly the historical scale: an index of -1.7 in
a 2009 forecast means the same thing as -1.7 in the GiOCEAN record, 1.7
historical standard deviations into the negative phase. The mean is close to 0
and the spread close to 1 by construction, so this step is nearly the
identity, but applying it exactly removes any residual convention mismatch.

## Why this is legitimate on data the patterns never saw

Two facts carry the argument:

- **Consistency.** Feed this exact procedure a historical GiOCEAN map and the
  recovered indices equal the historical rotated indices exactly, because the
  historical indices are themselves the least-squares solution of the same
  equation; that is how rotated principal component analysis defines them.
  This was checked numerically (agreement to 1e-10 on synthetic data). The
  forecast indices are therefore the natural out-of-sample extension of the
  historical ones: same formula, new data.
- **Nothing is refitted.** The patterns, the mean and the spread are frozen
  from the reanalysis. Forecasts are only measured against them, never used to
  define anything, so there is no leakage, and applying the yardstick to 1981
  or to 2026 is as valid as to any other year.

One honest limitation: least squares represents only the part of a forecast
map that projects onto the ten patterns. Structure orthogonal to all of them
lands in the residual and is simply not summarised. That is what an index is,
but it is worth saying out loud.

## The real-time case

`02_realtime_indices.py` runs the same chain on a single live forecast (the
February 2026 ensemble mean, which arrives in the GiOCEAN-style pressure-level
format, so the level is selected with `lev` rather than read as a 2D field).
Two differences from the hindcasts:

- **Climatology.** There is no hindcast archive for this model version yet, so
  the anomaly is taken against the GiOCEAN monthly mean for the verifying
  month (`analysis.reanalysis_month_mean`). That ignores the model's drift; if
  a drift-corrected hindcast climatology becomes available it should replace
  this, and the script says so in its header.
- **Season matching.** Each lead verifies in a different month, so each lead is
  fitted against the pattern set for the season containing its verifying month
  (the `SEASON_OF` mapping: MAM patterns for a March map, and so on). This has
  an interpretation consequence: "REOF2" in May is a spring pattern and
  "REOF2" in June is a summer pattern, two different physical structures. Read
  the index curves within a season; a jump across a season boundary partly
  reflects the change of pattern set, not the forecast. The `patterns` column
  in the CSV records which set each row used.

## Outputs

```
outputs/<level>hPa/<season>/<region>/lead<L>/forecast_indices.csv   hindcasts
outputs/<label>/<level>hPa/<region>/forecast_indices.csv           one live forecast
outputs/<label>/<level>hPa/<region>/indices.png                    its summary figure
```

A hindcast row reads: the forecast started on this date, from this member,
verifying this month, predicted these ten indices. All values are in standard
deviations on the historical scale, directly comparable with the historical
index CSVs in the GiOCEAN repo. Each live forecast gets its own dated folder
(the `LABEL` setting), so successive months accumulate side by side instead of
overwriting.

## Checks

- The fit recovers known indices from synthetic maps built with them, handles
  missing grid points, and scales linearly (double the anomaly, double the
  index).
- The lead arithmetic wraps the year correctly (a December start at lead 1
  verifies in January of the next year).
- The mean index over all 448 winter hindcasts is 0.000, as the anomaly
  construction requires.
- The 2009/10 winter, the most negative NAO winter in the record, appears in
  the hindcast indices, including the model flipping from positive to strongly
  negative as the start dates approached the event.

## Reproducing

```bash
module load python/GEOSpyD
python scripts/00_inspect_data.py        # look at the archive
python scripts/01_forecast_indices.py    # hindcast indices
python scripts/02_realtime_indices.py    # the live forecast
```

The GiOCEAN repo must have been run with the regression display first, since
that is what writes `patterns.nc`.
