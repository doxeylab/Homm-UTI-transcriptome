# Quickly get quality information from fastp files for Nader

# Load packages
import json
import glob
import os
import pandas as pd


# Define paths to all batches (including wildcard for files)
batch_1 = glob.glob("/data/mchomm/pathogen_surveillance/data/summary_batch1_fastp/*.json")
batch_2 = glob.glob("/data/mchomm/pathogen_surveillance/data/summary_batch2_fastp/*.json")
batch_3 = glob.glob("/data/mchomm/pathogen_surveillance/data/summary_batch3_fastp/*.json")
batch_4 = glob.glob("/data/mchomm/pathogen_surveillance/data/summary_batch4_fastp/*.json")

# Combine all batches into a single flat list rather than a list of lists
batches = batch_1 + batch_2 + batch_3 + batch_4
print(len(batches))

# Load the JSON file
proportion_good_dict = {}
for filepath in batches:
    # Skip empty files
    if os.path.getsize(filepath) == 0:
        print(f"Skipping empty file: {filepath}")
        continue

    # Load JSON with error handling
    try:
        with open(filepath, "r") as f:
            fastp_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Skipping invalid JSON: {filepath}")
        continue

    # Extract ptid and metrics
    ptid = os.path.basename(filepath).split("_")[0]
    reads_before = fastp_data["summary"]["before_filtering"]["total_reads"]
    reads_after  = fastp_data["summary"]["after_filtering"]["total_reads"]
    proportion_good_dict[ptid] = reads_after / reads_before

proportion_good_dict = {
    (key[:-5] if key.endswith(".json") else key): value
    for key, value in proportion_good_dict.items()
}

# Create a DataFrame from the dict
df = pd.DataFrame.from_dict(
    proportion_good_dict,
    orient='index',
    columns=['proportion_good']
)

# Turn the index into a column
df.index.name = 'sample_id'
df.reset_index(inplace=True)

# Write to CSV
df.to_csv('proportions.csv', index=False)

