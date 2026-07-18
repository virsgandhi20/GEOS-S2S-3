# GEOS-S2S-3 forecast teleconnection indices

Teleconnection forecast products from the GEOS-S2S-3 forecast system, built
on the rotated patterns found in the GiOCEAN work. No new modes are derived
here; each forecast height anomaly is measured against the historical
patterns, giving a teleconnection index for every forecast month, member and
lead. The indices state which phase of each pattern (NAO, PNA and so on) the
model predicts, on the same scale as the historical record, so forecast and
observation are directly comparable. The main product is a forecast plume
per mode in the style of the GMAO Nino 3.4 plume plots.

The method was developed and verified on the GEOS-S2S-2 hindcast archive,
where decades of past forecasts allow the indices to be checked against the
observed record; that work lives in the GEOS-S2S-2 repository, and the
leading winter modes verify at about 0.6 correlation one month ahead.

## The dataset

The GEOS-S2S-3 near-real-time forecasts and drift climatology live on PFE
(NAS):

```
GEOS_fcst4nrt/YYYYMMDD/ens<member>/<collection>/
    YYYYMMDD.<collection>.monthly.<verifying YYYYMM>.nc4
DRFT/DRFT_2001_2020/monthly/<collection>/<init mon>/
    <init mon>.<collection>.monthly.drift.<verifying MM>.nc4
```

- one directory per five-day start date, five ensemble members each; a
  subset of members (fifteen on the month's final start date) runs to nine
  lead months, the rest to three, and every run includes the partial
  initialization month itself
- the height field `H` sits on pressure levels in the
  `atm_inst_6hr_glo_L720x361_p49` collection, on the same half-degree grid
  as GiOCEAN; the collection name is a plain setting, so any collection
  applies
- the drift climatology (2001-2020) is read from the archive, one file per
  initialization month and verifying month; the code does not compute it
- only the `monthly` files are used, never the `dailymean` ones

## The method

1. Load the historical rotated patterns from the GiOCEAN analysis
   (`patterns.nc`, with the index statistics). The patterns are
   standardized: each is scaled to unit length before the projection, and
   the resulting indices are standardized with the statistics of the
   historical record.
2. Form the forecast anomaly: forecast minus the drift climatology for the
   same initialization month and verifying month. Forecast models drift
   with lead, and subtracting the drift climatology removes that model
   drift from every forecast.
3. Apply the same `sqrt(cos(latitude))` weighting the patterns were found
   with.
4. Compute the index per pattern: by default the weighted anomaly is
   projected onto each standardized pattern one at a time (the group's
   convention); a simultaneous least-squares fit against all patterns
   remains available as an option, and each convention writes to its own
   folder.

If a forecast grid differs from the GiOCEAN grid, the forecast is
interpolated onto the pattern grid first.

## Layout

```
scripts/     the code
data/        empty (the data lives on PFE / Discover)
outputs/
  plumes/            the forecast plumes, one folder per initialization month
  observed_recent/   the plumes' observed lead-in
  realtime/          single forecasts, one folder per case
docs/        notes
```

## Running it

```bash
python scripts/01_plume_forecasts.py         # the forecast plume (PFE)
python scripts/02_observed_recent.py         # observed lead-in (Discover)
python scripts/03_single_forecast.py         # one single forecast (Discover)
```

`01_plume_forecasts.py` is the main product. For one initialization month it
collects every member from every five-day start date (45 for a January start
month) and every available forecast month, including the partial
initialization month, computes each member's indices against the drift
climatology, and draws one plume per mode in the layout of the GMAO Nino 3.4
plume plots: the recent observed index as a solid black lead-in, each member
as a dashed line in its start date's colour, and the ensemble mean as a
heavy red line. The GiOCEAN pattern files have to be copied to PFE first,
since they stay on the machine that produced them.

`02_observed_recent.py` runs on Discover and measures recent GiOCEAN months
against the fixed patterns, writing the small CSV the plume draws its
observed lead-in from (committed, so it travels to PFE through git).

`03_single_forecast.py` handles a single forecast delivered as one set of
ensemble-mean files (the February 2026 case), with a choice of anomaly
baseline (the drift climatology, or the reanalysis monthly mean).

## More detail

`docs/methodology.md` walks through the whole calculation step by step, with
the reasoning behind each piece: the drift climatology, the weighting, the
two index conventions, the standardization, why the method is valid on data
the patterns never saw, and how the plume is assembled.

## Status

- [x] Repository scaffold; GiOCEAN saves its patterns (`patterns.nc`);
      method verified on the GEOS-S2S-2 hindcasts (see that repository)
- [x] Projection index (each map onto one standardized observed pattern at
      a time, the group's convention), alongside the simultaneous
      least-squares fit
- [x] Single-forecast indices (February 2026) against the drift climatology
- [x] Forecast plume on PFE: January 2026, all 45 members across the seven
      start dates, every available forecast month including the partial
      initialization month, one figure per mode with the observed lead-in
- [ ] Standardization convention of the projection index confirmed (scaled
      here by the spread of the same projection over the historical maps)
- [ ] GEOS-S2S-3 hindcast archive located and processed (the GEOS-S2S-2
      repository holds the template)
