# Running fastp on all patient files

# Dependencies

import subprocess
import os
import glob
import pandas as pd
# fastp v0.23.4 must first be installed, see https://github.com/OpenGene/fastp

# Specify file paths

data = "" # Replace with path to fastq files

# Run fastp

# 1) -> Create a dictionary where the key is the ID, and the value is a 2-item list with location of both paired fastq files

paired_dict = {}
for file in sorted(glob.glob(data + "*.fastq.gz")):
    patient_full_id = file.split("/")[8].split(".")[0] # Must tailor to your file
    patient_partial_id = patient_full_id.split("_")[0]
    if patient_partial_id in paired_dict:
        paired_dict[patient_partial_id].append(file)
    else:
        paired_dict[patient_partial_id] = [file]

# 2) -> Iterate through this dict and run fastp

# To see progress throughout
iteration_count = 0
total_iterations = len(overall_dict)

# Loop through each patient
for k, v in paired_dict.items():

    # Update which iteration you're on
    iteration_count += 1
    print(f"fastp progress -> currently processing sample {k}: {iteration_count}/{total_iterations}...")

    # Variables specifying which read: may need to tailor to your own files
    # *_run specifies whether R1 or R2, and *_id specifies its ID
    v1_run = v[0].split("/")[-1].split(".")[0].split("_")[3]
    v2_run = v[1].split("/")[-1].split(".")[0].split("_")[3]
    v1_id = v[0].split("/")[-1].split(".")[0].split("_")[2]
    v2_id = v[1].split("/")[-1].split(".")[0].split("_")[2]

    # Orders them -> R1, R2
    if v1_run != "R1":
        v[0], v[1] = v[1], v[0]

    # Creates the output directory for new fastq files
    fastq_output_dir = "" # Replace with desired output directory location
    os.makedirs(fastq_output_dir, exist_ok=True)

    # Creates the output directory for qc information
    summary_output_dir = "" # Replace with location of where to put HTML and JSON summary data
    os.makedirs(summary_output_dir, exist_ok=True)

    # Runs fastp
    subprocess.run(f"fastp --in1 {v[0]} --in2 {v[1]} --out1 fastq_output_dir/{k}_{v1_id}_{v1_run}.fastq.gz "
                   f"--out2 fastq_output_dir/{k}_{v2_id}_{v2_run}.fastq.gz --html summary_output_dir//{k}.html "
                   f"--json summary_output_dir/{k}.json --report_title {k}_{v1_id}_report --trim_front1 8 --trim_front2 8", shell=True)
