# Quickly get RPM data (kraken and bracken) for validation batches 5 and 6

# Load packages

import json
import glob
import os
import pandas as pd

# Define filepaths

krakenV = glob.glob("/data/mchomm/Homm-UTI-transcriptome/data/validation/kraken_out/*")
krakenV.extend(glob.glob("/data/mchomm/Homm-UTI-transcriptome/data/validation/kraken_out_batch_5/*"))
kraken_total = krakenV
bracken_final = [f for f in kraken_total if f.endswith("bracken.out")]

print(len(bracken_final))

# Define function to derive read and percentage numbers

def extract_data(filepaths, tax, lvl):
    '''
    Returns a pd df with human-mapped reads and read%
    '''

    # Want a dictionary of dictionaries
    human_dict = {}

    # Loop through, keep track of human information
    for file in filepaths: 
        ptid = file.split("/")[-1][0:5]
        #print(f"ptid = {ptid}")
        human_dict[ptid] = {"percentage" : "0.00", "reads" : "0"}
        with open(file) as current_file: 
            for line in current_file: 
                #print(f"line = {line}")
                line_list = line.split("\t")
                #print(f"line list = {line_list}")
                taxon = line_list[0].strip()
                taxon_level = line_list[2].strip()
                percentage = line_list[-1].strip()
                reads = line_list[-2].strip()
                if taxon == tax and taxon_level == lvl: 
                    human_dict[ptid] = {"percentage" : percentage, "reads" : reads}
    
    # Add to pandas df
    df = pd.DataFrame.from_dict(human_dict, orient='index')
    return(df)

# Human Testing
bracken_human = extract_data(bracken_final, "Homo sapiens", "S")
print(bracken_human)

# Bacteria Testing
bracken_bac = extract_data(bracken_final, "Bacteria", "D")
print(bracken_bac)

dataframes = {
    "bracken_human": bracken_human,
    "bracken_bac": bracken_bac
}

# Write CSVs for all of these
output_dir = "/data/mchomm/Homm-UTI-transcriptome/misc/good_rpm/"
for name, df in dataframes.items():
    df.to_csv(f"{output_dir}{name}.csv")
    
# Define function to derive actual RPM dfs

def extract_rpm(filepaths, rpm_column_index=2):
    pass


# Define function to find percentages and raw read counts for a list of defined pathogens

def extract_data_multi(filepaths, tax_list, lvl):
    '''
    Returns a pandas DataFrame with the total percentage and read count
    for all taxa in tax_list at the specified taxonomic level.
    '''
    combined_dict = {}

    for file in filepaths:
        ptid = file.split("/")[-1][0:5]
        total_reads = 0
        total_percentage = 0.0

        with open(file) as current_file:
            for line in current_file:
                line_list = line.strip().split("\t")
                if len(line_list) < 4:
                    continue
                taxon = line_list[0].strip()
                taxon_level = line_list[2].strip()
                percentage = line_list[-1].strip()
                reads = line_list[-2].strip()

                if taxon_level == lvl and taxon in tax_list:
                    try:
                        total_reads += int(reads)
                        total_percentage += float(percentage)
                    except ValueError:
                        continue

        combined_dict[ptid] = {
            "percentage": round(total_percentage, 4),
            "reads": total_reads
        }

    df = pd.DataFrame.from_dict(combined_dict, orient='index')
    return df

# Testing

tax_list = ["Escherichia coli", "Klebsiella pneumoniae", "Enterococcus faecalis", "Proteus mirabilis", 
                     "Pseudomonas aeruginosa", "Staphylococcus saprophyticus", "Streptococcus agalactiae"]

# Pathogen Testing
bracken_pathogen = extract_data_multi(bracken_final, tax_list, lvl="S")
print(bracken_pathogen)

dataframes = {
    "bracken_pathogen": bracken_pathogen,
}

# Write CSVs for all of these
output_dir = "/data/mchomm/Homm-UTI-transcriptome/misc/good_rpm/"
for name, df in dataframes.items():
    df.to_csv(f"{output_dir}{name}.csv")


# Define function to actually compute the RPM values

import pandas as pd
import numpy as np

def extract_nonhuman_rpm(filepaths, tax_list, lvl="S", pseudocount=1e-6):
    """
    Returns a DataFrame with:
      - total, human, non-human, and total pathogen reads
      - per-pathogen read counts
      - per-pathogen non-human RPMs
      - summed total non-human RPM
      - log10(non-human RPM) values (clamped at 0)
    """

    results = {}

    for file in filepaths:
        ptid = file.split("/")[-1][0:5]

        total_reads = 0
        human_reads = 0
        nonhuman_reads = 0

        # initialize per-pathogen dict for reads
        pathogen_reads_dict = {tax: 0 for tax in tax_list}

        with open(file) as current_file:
            for line in current_file:
                line_list = line.strip().split("\t")
                if len(line_list) < 4:
                    continue

                taxon = line_list[0].strip()
                taxon_level = line_list[2].strip()
                reads_str = line_list[-2].strip()

                try:
                    reads = int(reads_str)
                except ValueError:
                    continue

                total_reads += reads

                # record human reads
                if taxon == "Homo sapiens" and taxon_level == "S":
                    human_reads += reads

                # record pathogen reads (species level match)
                elif taxon_level == lvl and taxon in tax_list:
                    pathogen_reads_dict[taxon] += reads

        nonhuman_reads = total_reads - human_reads

        # compute per-pathogen RPMs
        pathogen_rpms = {}
        if nonhuman_reads > 0:
            for taxon in tax_list:
                pathogen_rpms[taxon + "_RPM"] = (pathogen_reads_dict[taxon] / nonhuman_reads) * 1e6
        else:
            for taxon in tax_list:
                pathogen_rpms[taxon + "_RPM"] = 0.0

        total_pathogen_reads = sum(pathogen_reads_dict.values())
        total_pathogen_rpm = sum(pathogen_rpms.values())

        # log10 transformation (safe + clamped)
        log10_rpm_raw = np.log10(total_pathogen_rpm + pseudocount)
        log10_rpm_clamped = max(log10_rpm_raw, 0)

        # assemble row
        row = {
            "total_reads": total_reads,
            "human_reads": human_reads,
            "nonhuman_reads": nonhuman_reads,
            "pathogen_reads_total": total_pathogen_reads,
            "nonhuman_RPM_total": round(total_pathogen_rpm, 6),
            "log10_RPM": round(log10_rpm_clamped, 6)
        }

        # merge in pathogen-specific read counts and RPMs
        for taxon in tax_list:
            row[taxon + "_reads"] = pathogen_reads_dict[taxon]
            row[taxon + "_RPM"] = round(pathogen_rpms[taxon + "_RPM"], 6)

        results[ptid] = row

    return pd.DataFrame.from_dict(results, orient="index")




bracken_nonhuman_rpm = extract_nonhuman_rpm(bracken_final, tax_list, lvl="S")
print(bracken_nonhuman_rpm)

output_dir = "/data/mchomm/Homm-UTI-transcriptome/data/validation/"
bracken_nonhuman_rpm.to_csv(f"{output_dir}validation_bracken_pathogen_nonhumanRPM.csv")