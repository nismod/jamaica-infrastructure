"""
The Snakemake workflow eventually needs to generate the files that are referenced
in CSVs from the irv-jamaica repository's etl folder.

Those CSVs are:

- adaptation_files.csv
- damage_exp_files.csv
- damage_rp_files.csv
- damage_rp_files_with_buildings.csv
- hazard_layers.csv
- network_layers.csv
- network_tilelayers.csv
- network_tilelayers_with_buildings.csv
- other_tilelayers.csv
- hotspot_layers.csv
- storm_layers.csv

"""

import os
import subprocess
import pandas

# Check the files that we need to generate
csv_path_refs = {
    # "adaptation_files": ["avoided_risk"],
    "damage_exp_files": ["expected"],
    "damage_rp_files": ["damage", "exposure", "loss"],
    "damage_rp_files_with_buildings": ["damage", "exposure", "loss"],
    # "hazard_layers": ["path"],
    # "hotspot_layers": ["path"],
    # "network_layers": ["path", "single_failure_scenarios"],
    # "network_tilelayers": [],
    # "network_tilelayers_with_buildings": [],
    # "other_tilelayers": [],
    # "storm_layers": ["path"],
}

required_files = set()
for f, cols in csv_path_refs.items():
    csv_file = pandas.read_csv(f"https://github.com/nismod/irv-jamaica/raw/refs/heads/main/etl/{f}.csv")
    for col in cols:
        required_files.update(csv_file[col].dropna().values)

# alphabetize the files for ease of comparison
required_files = sorted(required_files)

processed_data = "../processed_data"

missing_files = []
ambiguous_files = []

for f in required_files:
    matches = subprocess.run(
        ["bash", "-c", f"find {processed_data} -iwholename *{f}"],
        capture_output=True,
        text=True,
    )
    count = subprocess.run(
        ["bash", "-c", "wc -l"],
        input=matches.stdout,
        capture_output=True,
        text=True,
    )
    try:
        n = int(count.stdout.strip())
    except ValueError:
        print(f"Error finding {f}: {count.stdout}")
        continue

    if n == 0:
        missing_files.append(f)
    elif n > 1:
        ambiguous_files.append({"search": f, "matches": matches.stdout})

if missing_files:
    print("The following files are missing:")
    print("\n".join(missing_files))

if ambiguous_files:
    print("The following files are ambiguous:")
    for f in ambiguous_files:
        print(f['search'] + ":")
        print(f['matches'])

print(f"{len(required_files) - len(missing_files)} files present in {processed_data}.")
print(f"{len(ambiguous_files)} files matched multiple paths in {processed_data}.")
print(f"{len(missing_files)} files need to be generated by the workflow.")
