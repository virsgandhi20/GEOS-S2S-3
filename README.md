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

## Running it

```bash
module load python/GEOSpyD
python scripts/00_inspect_data.py /path/to/geos-s2s-2   # look at the data first
```

The projection script comes after the inspection settles the file layout,
naming, grid and how the initialization and lead information is stored.

## Status

- [x] Repo scaffold
- [x] GiOCEAN saves its patterns (`patterns.nc`) for use here
- [ ] Locate the GEOS-S2S-2 tree on Discover and inspect a file
- [ ] Forecast climatology found (per initialization month and lead) or computed
- [ ] Projection script (anomaly, weighting, least squares, scaling)
- [ ] Switch the data root to GEOS-S2S-3 on NAS once access comes through
