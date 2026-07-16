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

The target dataset is GEOS-S2S-3. Its near-real-time forecasts and drift
climatology live on PFE (NAS) and feed the real-time products (the plume and
the single-forecast script). Its close relative, GEOS-S2S-2, is available on
Discover and was used to build and verify the hindcast side; the data roots
are settings, so the hindcast scripts point at a GEOS-S2S-3 hindcast archive
whenever one becomes reachable.

## The method

1. Load the historical rotated patterns from the GiOCEAN analysis (saved by
   that repository as `patterns.nc`, together with each index's mean and spread).
2. Form the forecast anomaly: forecast minus the forecast climatology for the
   same initialization month and lead. Forecast models drift with lead time, so
   the climatology has to be lead-dependent.
3. Apply the same `sqrt(cos(latitude))` weighting the patterns were found with.
4. Compute the index against the observed patterns. The default projects
   the weighted anomaly onto each standardized pattern one at a time (the
   group's convention); a simultaneous least-squares fit against all
   patterns (`y = E b + error`) remains available as an option, and each
   convention writes to its own folder.
5. Scale by the statistics of the same computation applied to the historical
   record, so the forecast indices sit on the same normalized scale as the
   historical ones.

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
python scripts/05_observed_recent.py         # observed lead-in for the plume
```

The hindcast scripts (`01`, `03`) and their outputs use GEOS-S2S-2, the
version available on Discover; they are kept as the working template for the
GEOS-S2S-3 hindcasts once that archive is reachable. Current work targets
GEOS-S2S-3 (`02`, `04`, `05`).

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
every member from every five-day start date, computes each member's indices
against the drift climatology, and draws one plume per mode in the layout of
the GMAO Nino 3.4 plume plots: the recent observed index as a solid black
lead-in, each member as a dashed line in its start date's colour, and the
ensemble mean as a heavy red line. The index is computed by projection onto
each standardized observed pattern (the `METHOD` setting; the simultaneous
least-squares fit remains available as `"lstsq"`). The GiOCEAN pattern files
have to be copied to PFE first, since they stay on the machine that produced
them.

`05_observed_recent.py` runs on Discover and measures recent GiOCEAN months
against the fixed patterns, writing the small CSV the plume draws its
observed lead-in from (committed, so it travels to PFE through git).

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
- [x] Projection index implemented (each map onto one standardized observed
      pattern at a time, the group's convention), alongside the simultaneous
      least-squares fit
- [x] GEOS-S2S-3 forecast plume on PFE: January 2026, all 45 members across
      the seven start dates, indices against the drift climatology, one
      figure per mode with the observed lead-in
- [ ] Standardization convention of the projection index confirmed (scaled
      here by the spread of the same projection over the historical maps)
- [ ] GEOS-S2S-3 hindcast archive located and processed (the S2S-2 hindcast
      scripts are the template)
