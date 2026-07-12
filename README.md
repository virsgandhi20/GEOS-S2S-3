# GEOS-S2S forecast teleconnection indices

Applying the rotated patterns found in the GiOCEAN work to GEOS-S2S forecast
fields. Instead of finding new modes here, each forecast height anomaly is
fitted against the historical rotated patterns by least squares, which gives a
teleconnection index for every forecast month and lead. That tells us which
phase of each pattern (NAO, PNA and so on) the model predicts, and lets us
compare forecast indices against the historical record on the same scale.

## The dataset (and why the path is a setting)

GEOS-S2S is a forecast system: a new forecast starts every five days, and each
run carries about nine lead months. The conventions used here:

- the initialization month is the month the forecast starts in
- lead month 1 is the month after the initialization month (a March start has
  April as lead 1, May as lead 2, and so on)
- picking a season at a fixed lead shifts the verifying months (initializations
  in D, J, F at lead 1 verify in J, F, M)

The version we ultimately want, GEOS-S2S-3, lives on NAS, which we don't have
access to yet. A close relative, GEOS-S2S-2, is on Discover. The code is being
built against the Discover copy, with the data root kept as a setting so
switching to the NAS dataset later is a one-line change.

## The method

1. Load the historical rotated patterns from the GiOCEAN analysis (saved by
   that repo as `patterns.nc`, together with each index's mean and spread).
2. Form the forecast anomaly: forecast minus the forecast climatology for the
   same initialization month and lead. Forecast models drift with lead time, so
   the climatology has to be lead-dependent.
3. Apply the same `sqrt(cos(latitude))` weighting the patterns were found with.
4. Fit the weighted anomaly against the patterns by least squares
   (`y = E b + error`, solving for `b`), one fit per forecast map.
5. Scale `b` by the historical index spread so the forecast indices sit on the
   same normalized scale as the historical ones.

If the forecast grid differs from the GiOCEAN grid, the forecast is
interpolated onto the pattern grid before the fit.

## Layout

```
scripts/     the code
data/        empty (the data lives on Discover / NAS)
outputs/     forecast indices and figures
docs/        notes
```

## What the archive looks like (GEOS-S2S-2 on Discover)

```
runx/<year>/<init, e.g. feb10>/ens<member>/geosgcm_vis2d/
    feb10.geosgcm_vis2d.monthly.<verifying YYYYMM>.nc4
```

- hindcast years from 1981; initializations every five days (dec02, dec07, ...)
- one monthly file per lead: the partial initialization month plus about nine
  full lead months
- `H500` and `H250` come ready-made in `geosgcm_vis2d`, in metres
- the atmospheric grid is the same half-degree grid as GiOCEAN, so the fit
  needs no regridding (the code falls back to interpolation if a future
  dataset differs)
- member numbering varies by initialization (most have `ens1`; some carry
  `ens2` to `ens5` instead), so members are found by glob, never assumed

## Running it

```bash
module load python/GEOSpyD
python scripts/00_inspect_data.py            # look at the data
python scripts/01_forecast_indices.py        # the projection
```

`01_forecast_indices.py` writes one CSV per level, initialization season and
lead (`outputs/500hPa/DJF/20N_90N/lead1/forecast_indices.csv`), with a row per
forecast: init date, member, verifying month, and the index for each mode. The
anomaly is taken against the mean over years for the same initialization date
and lead, so the model's drift is removed. Settings (levels, seasons, leads,
members, years, data root, pattern file) are at the top.

## Status

- [x] Repo scaffold
- [x] GiOCEAN saves its patterns (`patterns.nc`) for use here
- [x] GEOS-S2S-2 tree located and inspected (grid matches GiOCEAN)
- [x] Projection script (anomaly, weighting, least squares, scaling)
- [x] Real-time script for a single forecast (the February 2026 case, in the
      GiOCEAN-style pressure-level format; anomaly against the reanalysis
      monthly mean until a drift-corrected climatology is settled)
- [ ] First run, and a check of the indices against known winters
- [ ] Ensemble handling settled with Young-Kwon's approach (his runs used four members)
- [ ] Switch the data root to GEOS-S2S-3 on NAS once access comes through
