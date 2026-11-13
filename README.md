# G08 GNSS Positioning Lab

AAE4203 Guidance and Navigation, Week 9, The Hong Kong Polytechnic University.
This repository reproduces the G08 own-data and UrbanNav benchmark submission
CSVs. All trajectory coordinates are WGS-84 ECEF in metres.

## Reproduce the CSVs on Windows

Use 64-bit Python 3.12 and PowerShell. The package includes the RTKLIB Explorer
2.5.1 command-line programs for Windows.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe reproduce.py
```

The script regenerates the two CSVs, runs the course format and epoch-membership
checks, and compares the generated files with `submission.sha256`. Run the
script again to replace generated results from an earlier run. Raw inputs are
read-only during processing.

| Submission file | Solved epochs | Observed epochs | Coverage |
| --- | ---: | ---: | ---: |
| `G08_owndata.csv` | 487 | 549 | 88.71% |
| `benchmark/G08_benchmark.csv` | 1,437 | 1,534 | 93.68% |

Each CSV uses the columns `gps_week,gps_tow_s,x_ecef_m,y_ecef_m,z_ecef_m`.
The course leaderboard determines UrbanNav accuracy and rank.

## Inputs and code

- `raw/G08_owndata_20261030.ubx` is the original G08 receiver log. CONVBIN extracts rover RINEX observations and navigation. RTKLIB produces the GPS SPP comparison.
- `code/process_g08.py` reads checksummed UBX messages and GPS broadcast ephemerides. It solves GPS L1 pseudoranges by iterative elevation-weighted least squares, with satellite clock, atmospheric, group-delay and Earth-rotation corrections. It writes the own-data CSV and validation plots. Failed epochs are recorded in `data/omitted_epochs.csv`.
- `course/benchmark/` contains the UrbanNav observation and navigation files and the GPS L1 matrices supplied by the course. The matrices, generated with `course/rinex2csv.m`, let the benchmark be reproduced in Python without running MATLAB. `course/WLSE_spp.py` solves those matrices and writes the benchmark CSV.
- `course/example/` contains the 2025 Block X worked example used to verify the solver.
- `raw/hksc/` contains the HKSC observations and GPS navigation used for the differential comparison in the report. To rerun it after regenerating rover RINEX, use the command below.

```powershell
$root = (Get-Location).Path
.\tools\rnx2rtkp.exe -k "$root\code\dgps.conf" -o "$root\data\G08_HKSC_DGNSS.pos" "$root\data\G08.obs" "$root\raw\hksc\hksc303i.26o" "$root\data\G08.nav" "$root\raw\hksc\hksc303i.26n"
```

## Sources

- [AAE4203 course code and teaching data](https://github.com/weisongwen/AAE4203-S12627/tree/main/Lab1). See `course/LICENSE`.
- [RTKLIB Explorer 2.5.1](https://github.com/rtklibexplorer/RTKLIB/releases/tag/v2.5.1). See `tools/RTKLIB_readme.txt`.
- Hong Kong Satellite Positioning Reference Station Network, HKSC, 30 October 2026, 08:00–08:59:59 GPST.

## Course verification

Run `python course/tools/verify_package.py` from the repository root. It checks the 186-epoch historical regression, synthetic GPS position and clock recovery, missing-epoch handling, malformed CSVs, and benchmark input hashes and coverage. The suite passed on this package.
