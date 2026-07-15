# GEOS-S2S forecast teleconnection indices

Applying the rotated patterns found in the GiOCEAN work to GEOS-S2S forecast
fields. No new modes are derived here; each forecast height anomaly is fitted
against the historical rotated patterns by least squares, giving a
teleconnection index for every forecast month and lead. The indices state which
phase of each pattern (NAO, PNA and so on) the model predicts, on the same
scale as the historical record, so forecast and observation are directly
comparable.

## The dataset (and why the path is a setting)

GEOS-S2S is a forecast system: a new forecast starts every five days, and each
run carries about nine lead months. The conventions used here:

- the initialization month is the month the forecast starts in
- lead month 1 is the month after the initialization month (a March start has
  April as lead 1, May as lead 2, and so on)
- picking a season at a fixed lead shifts the verifying months (initializations
  in D, J, F at lead 1 verify in J, F, M)

The target dataset, GEOS-S2S-3, resides on NAS, to which access is pending.
Its close relative, GEOS-S2S-2, is available on Discover. The code is built
against the Discover copy, with the data root kept as a setting so that
switching to the NAS dataset later is a one-line change.

## The method

1. Load the historical rotated patterns from the GiOCEAN analysis (saved by
   that repository as `patterns.nc`, together with each index's mean and spread).
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
python scripts/01_forecast_indices.py        # hindcast indices
python scripts/02_realtime_indices.py        # one live forecast
python scripts/03_verification_plot.py       # forecast vs observed figures
python scripts/04_plume_s2s3.py              # GEOS-S2S-3 forecast plume (PFE)
```

`01_forecast_indices.py` processes two input sources, each written to its own
folder (`outputs/500hPa/DJF/20N_90N/lead1/members/` and `.../lead1/ensmean/`),
with a row per forecast: init date, member, verifying month, and the index for
each mode. The "members" source reads the individual forecasts, one per
five-day start date; the "ensmean" source reads the pre-averaged product, one
forecast per initialization month, where the averaging over start dates and
members has already been done. In both cases the anomaly is taken against the
mean over years for the same initialization and lead, so the model's drift is
removed. Settings (sources, levels, seasons, leads, members, years, data roots,
pattern file) are at the top.

`03_verification_plot.py` draws one figure per mode alongside each hindcast
CSV: the pattern map on top, and below it the forecast index (the average
over everything verifying in a month, in black) against the historical
GiOCEAN index (blue), with their correlation printed. At lead 1 the leading
winter modes correlate at about 0.6.

`04_plume_s2s3.py` is the GEOS-S2S-3 near-real-time step and runs on PFE
(NAS), where that archive lives. For one initialization month it collects
every member from every five-day start date (five members each, fifteen on
the month's last start), computes each member's indices against the drift
climatology, and draws one plume per mode: members as thin dashed lines,
their average as a thick solid line, in the layout of the GMAO Nino 3.4
plume plots. The index is computed by projection onto each standardized
observed pattern (the `METHOD` setting; the simultaneous least-squares fit
remains available as `"lstsq"`). The GiOCEAN pattern files have to be copied
to PFE first, since they stay on the machine that produced them.

## More detail

`docs/methodology.md` walks through the whole calculation step by step, with the
reasoning behind each piece: the drift-corrected climatology, the weighting, the
least-squares fit, the scaling, why the method is valid on data the patterns
never saw, and the checks used to trust it.

## Status

- [x] Repository scaffold
- [x] GiOCEAN saves its patterns (`patterns.nc`) for use here
- [x] GEOS-S2S-2 tree located and inspected (grid matches GiOCEAN)
- [x] Projection script (anomaly, weighting, least squares, scaling)
- [x] Real-time script for a single forecast (the February 2026 case, in the
      GiOCEAN-style pressure-level format; anomaly against the GEOS-S2S-3
      drift climatology, with the reanalysis monthly mean kept as an
      alternative baseline)
- [x] First run: 448 winter and 396 summer hindcasts indexed; the mean index
      is zero as required, and the 2009/10 negative NAO winter appears
      clearly, including the transition from positive to negative indices as
      the initialization dates approached the event
- [x] Verification figures for all ten modes, both sources: forecast index
      against the GiOCEAN index with the correlation printed
- [ ] Ensemble handling settled with Young-Kwon's approach (his runs used
      four members, ens2 to ens5)
- [ ] Index computation switched from the simultaneous least-squares fit to
      projection onto each standardized observed pattern, per the group's
      convention (details to be settled)
- [ ] Switch the data root to GEOS-S2S-3 on NAS once access comes through
      (the February 2026 real-time chain already runs on GEOS-S2S-3: the
      forecast file and the drift climatology both come from it)
