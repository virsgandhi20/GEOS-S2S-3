"""
Have a look at the forecast data before writing any analysis code.

GEOS-S2S is a forecast dataset, so the files are organised very differently
from a reanalysis: there is a new initialization every five days, each with its
own set of lead months (and possibly ensemble members). Before the projection
code can be written, the directory layout, the file naming, the variables
inside, the grid, and the way the time/lead information is stored all have
to be seen.

This script does two things:
  - survey a directory: list its contents and find a few forecast files
  - inspect one file: print its variables, dimensions and coordinates

Usage (on Discover):
    module load python/GEOSpyD
    python scripts/00_inspect_data.py /path/to/the/geos-s2s-2/tree
    python scripts/00_inspect_data.py /full/path/to/one_file.nc4

To locate the GEOS-S2S-2 tree, try:
    ls /discover/nobackup/projects/gmao | grep -i s2s
or check the paths used inside Priyanka's Tele-V2-Retro scripts.
"""

import os
import sys
import fnmatch

import numpy as np
import xarray as xr

# ===== settings ============================================================
# roots to try when no path is given on the command line. the GEOS-S2S-2
# forecast archive sits under runx/, organised as
#   runx/<year>/<MMDD init>/ens<member>/<collection>/<MMDD>.<collection>.monthly.<YYYYMM>.nc4
CANDIDATE_ROOTS = [
    "/discover/nobackup/projects/gmao/m2oasf/aogcm/g5fcst/forecast/production/geos-s2s/runx",
]
FILE_PATTERN = "*.monthly.*.nc4"   # the monthly files, one per lead
MAX_DEPTH    = 8           # how deep to search below the root
MAX_MATCHES  = 20          # stop after finding this many files
# ===========================================================================


def survey(root):
    """Show what's under a directory and find a few data files."""
    print("=" * 78)
    print("SURVEYING:", root)
    print("=" * 78)

    entries = sorted(os.listdir(root))
    dirs = [e for e in entries if os.path.isdir(os.path.join(root, e))]
    files = [e for e in entries if not os.path.isdir(os.path.join(root, e))]
    print(f"\n{len(dirs)} directories, {len(files)} files at the top level")
    for d in dirs[:15]:
        print("  dir :", d)
    if len(dirs) > 15:
        print(f"  ... and {len(dirs) - 15} more directories")
    for f in files[:10]:
        print("  file:", f)

    print(f"\nsearching for {FILE_PATTERN} (depth <= {MAX_DEPTH})...")
    matches = []
    base_depth = root.rstrip("/").count("/")
    for dirpath, subdirs, filenames in os.walk(root, followlinks=True):
        if dirpath.count("/") - base_depth >= MAX_DEPTH:
            subdirs[:] = []
            continue
        for fn in fnmatch.filter(filenames, FILE_PATTERN):
            matches.append(os.path.join(dirpath, fn))
            if len(matches) >= MAX_MATCHES:
                break
        if len(matches) >= MAX_MATCHES:
            break

    if not matches:
        print("no matching files found; the data may sit deeper, or under a "
              "different extension. adjust FILE_PATTERN / MAX_DEPTH.")
        return
    print(f"first {len(matches)} matches:")
    for m in matches:
        print("  ", m)
    print("\ninspecting the first one:\n")
    inspect(matches[0])


def inspect(path):
    """Print what's inside one file."""
    print("=" * 78)
    print("SAMPLE FILE:", path)
    print("=" * 78)
    data = xr.open_dataset(path, chunks={})

    print("\n--- xarray repr ---------------------------------------------------")
    print(data)

    print("\n--- COORDINATES ----------------------------------------------------")
    for name in data.coords:
        c = data.coords[name]
        vals = np.asarray(c.values).ravel()
        if vals.size <= 12:
            preview = ", ".join(str(v) for v in vals)
        else:
            preview = (f"{vals[0]} ... {vals[-1]}  (n={vals.size})")
        print(f"  {name:16s} dims={c.dims} units={c.attrs.get('units', '')!r}")
        print(f"       values: {preview}")

    print("\n--- DATA VARIABLES -------------------------------------------------")
    for name in data.data_vars:
        v = data[name]
        print(f"  {name:16s} dims={v.dims} units={v.attrs.get('units', '')!r}  "
              f"{v.attrs.get('long_name', '')}")

    print("\n--- GLOBAL ATTRIBUTES ----------------------------------------------")
    for key, val in data.attrs.items():
        print(f"  {key}: {val}")
    data.close()


def main():
    if len(sys.argv) > 1:
        target = sys.argv[1]
        if os.path.isdir(target):
            survey(target)
        else:
            inspect(target)
        return

    for root in CANDIDATE_ROOTS:
        if os.path.isdir(root):
            survey(root)
            return
    sys.exit("None of the candidate roots exist. Pass the data directory as an "
             "argument:  python scripts/00_inspect_data.py /path/to/tree")


if __name__ == "__main__":
    main()
