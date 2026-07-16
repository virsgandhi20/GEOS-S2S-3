# GEOS-S2S-3 forecast teleconnection indices

Teleconnection forecast products from the GEOS-S2S-3 forecast system, built on
the rotated patterns found in the GiOCEAN work. No new modes are derived here;
each forecast height anomaly is measured against the historical patterns,
giving a teleconnection index for every forecast month, member and lead. The
indices state which phase of each pattern (NAO, PNA and so on) the model
predicts, on the same scale as the historical record, so forecast and
observation are directly comparable. The main product is a forecast plume per
mode in the style of the GMAO Nino 3.4 plume plots.

## The datasets

GEOS-S2S is a forecast system: a new forecast starts every five days, and each
run carries up to about nine lead months. The conventions used here:

- the initialization month is the month the forecast starts in
- lead month 1 is the month after the initialization month (a March start has
  April as lead 1, May as lead 2, and so on)
- picking a season at a fixed lead shifts the verifying months (initializations
  in D, J, F at lead 1 verify in J, F, M)

**GEOS-S2S-3** is the target system. Its near-real-time forecasts and drift
climatology live on PFE (NAS) and feed the forecast products:

```
GEOS_fcst4nrt/YYYYMMDD/ens<member>/<collection>/
    YYYYMMDD.<collection>.monthly.<verifying YYYYMM>.nc4
DRFT/DRFT_2001_2020/monthly/<collection>/<init mon>/
    <init mon>.<collection>.monthly.drift.<verifying MM>.nc4
```

- one directory per five-day start date, five ensemble members each; a subset
  of members (fifteen on the month's final start date) runs to nine lead
  months, the rest to three
- the height field `H` sits on pressure levels in the
  `atm_inst_6hr_glo_L720x361_p49` collection, on the same half-degree grid as
  GiOCEAN; the collection name is a plain setting, so any collection applies
- the drift climatology (2001-2020) provides the anomaly baseline per
  initialization month and verifying month
- only the `monthly` files are used, never the `dailymean` ones

**GEOS-S2S-2**, the previous version of the system, has a long hindcast
archive on Discover (`runx/<year>/<init>/ens<member>/geosgcm_vis2d/`,
initializations from 1981). It served to develop and verify the method: with
decades of past forecasts, the forecast indices can be checked against the
observed record, which a near-real-time archive cannot provide. The hindcast
scripts remain the working template for GEOS-S2S-3 hindcasts.

## The method

1. Load the historical rotated patterns from the GiOCEAN analysis (saved by
   that repository as `patterns.nc`, together with the index statistics).
2. Form the forecast anomaly: forecast minus the forecast climatology for the
   same initialization month and lead (the drift climatology for GEOS-S2S-3,
   the mean over hindcast years for GEOS-S2S-2). Forecast models drift with
   lead time, so the climatology has to be lead-dependent.
3. Apply the same `sqrt(cos(latitude))` weighting the patterns were found with.
4. Compute the index against the observed patterns. The default projects
   the weighted anomaly onto each standardized pattern one at a time (the
   group's convention); a simultaneous least-squares fit against all
   patterns (`y = E b + error`) remains available as an option, and each
   convention writes to its own folder.
5. Scale by the statistics of the same computation applied to the historical
   record, so the forecast indices sit on the same normalized scale as the
   historical ones.

If a forecast grid differs from the GiOCEAN grid, the forecast is
interpolated onto the pattern grid before the fit.

## Layout

```
scripts/     the code
data/        empty (the data lives on PFE / Discover)
outputs/
  plumes/            the forecast plumes, one folder per initialization month
  observed_recent/   the plumes' observed lead-in
  realtime/          single forecasts, one folder per case
  hindcasts/         the GEOS-S2S-2 indices and verification figures
docs/        notes
```

## Running it

The GEOS-S2S-3 products (PFE, plus one helper on Discover):

```bash
python scripts/04_plume_s2s3.py              # the forecast plume (PFE)
python scripts/05_observed_recent.py         # observed lead-in (Discover)
python scripts/02_realtime_indices.py        # one single forecast (Discover)
```

The GEOS-S2S-2 hindcast side (Discover):

```bash
module load python/GEOSpyD
python scripts/00_inspect_data.py            # look at the data
python scripts/01_forecast_indices.py        # hindcast indices
python scripts/03_verification_plot.py       # forecast vs observed figures
```

`04_plume_s2s3.py` is the main product. For one initialization month it
collects every member from every five-day start date (45 for a January start
month), computes each member's indices against the drift climatology, and
draws one plume per mode in the layout of the GMAO Nino 3.4 plume plots: the
recent observed index as a solid black lead-in, each member as a dashed line
in its start date's colour, and the ensemble mean as a heavy red line. The
GiOCEAN pattern files have to be copied to PFE first, since they stay on the
machine that produced them.

`05_observed_recent.py` runs on Discover and measures recent GiOCEAN months
against the fixed patterns, writing the small CSV the plume draws its
observed lead-in from (committed, so it travels to PFE through git).

`02_realtime_indices.py` handles a single forecast delivered as one set of
ensemble-mean files (the February 2026 case), with a choice of anomaly
baseline (the drift climatology, or the reanalysis monthly mean).

`01_forecast_indices.py` processes the GEOS-S2S-2 hindcasts into
`outputs/hindcasts/`, from two input sources, each written to its own folder: "members" reads the individual
forecasts, one per five-day start date, and "ensmean" reads the archive's
pre-averaged product. The anomaly is taken against the mean over years for
the same initialization and lead. `03_verification_plot.py` then draws one
figure per mode: the pattern map on top, and the forecast index against the
historical GiOCEAN index with their correlation printed. At lead 1 the
leading winter modes correlate at about 0.6, which is the skill evidence
behind the forecast products.

## More detail

`docs/methodology.md` walks through the whole calculation step by step, with the
reasoning behind each piece: the drift-corrected climatology, the weighting, the
two index conventions, the scaling, why the method is valid on data the patterns
never saw, and the checks used to trust it.

## Status

- [x] Repository scaffold; GiOCEAN saves its patterns (`patterns.nc`)
- [x] Method developed and verified on the GEOS-S2S-2 hindcasts: 448 winter
      and 396 summer forecasts indexed, mean index zero as required, the
      2009/10 negative NAO winter recovered, and lead-1 correlations of
      about 0.6 for the leading winter modes
- [x] Projection index implemented (each map onto one standardized observed
      pattern at a time, the group's convention), alongside the simultaneous
      least-squares fit; both conventions available in every script
- [x] GEOS-S2S-3 single-forecast indices (February 2026) against the drift
      climatology
- [x] GEOS-S2S-3 forecast plume on PFE: January 2026, all 45 members across
      the seven start dates, one figure per mode with the observed lead-in
- [ ] Standardization convention of the projection index confirmed (scaled
      here by the spread of the same projection over the historical maps)
- [ ] GEOS-S2S-3 hindcast archive located and processed (the GEOS-S2S-2
      hindcast scripts are the template)
