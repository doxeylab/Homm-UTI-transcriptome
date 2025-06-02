# Quickly get RPM data (kraken and bracken) for Nader

# Load packages

import json
import glob
import os
import pandas as pd

# Define filepaths

kraken1 = glob.glob("/data/mchomm/pathogen_surveillance/data/kraken_output/*.kraken")
kraken3 = glob.glob("/data/mchomm/pathogen_surveillance/data/kraken_output_batch3/*.kraken")
kraken4 = glob.glob("/data/mchomm/pathogen_surveillance/data/kraken_output_batch4/*.kraken")
kraken_total = kraken1 + kraken3 + kraken4
kraken_final = [f for f in kraken_total if not f.endswith("species.kraken")]
bracken_final = [f for f in kraken_total if f.endswith("species.kraken")]

print(len(bracken_final))
print(len(kraken_final))

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
        human_dict[ptid] = {"percentage" : "0.00", "reads" : "0"}
        with open(file) as current_file: 
            for line in current_file: 
                line_list = line.split("\t")
                taxon = line_list[-1].strip()
                taxon_level = line_list[-3].strip()
                percentage = line_list[0].strip()
                reads = line_list[1].strip()
                if taxon == tax and taxon_level == lvl: 
                    human_dict[ptid] = {"percentage" : percentage, "reads" : reads}
    
    # Add to pandas df
    df = pd.DataFrame.from_dict(human_dict, orient='index')
    return(df)

# Human Testing
kraken_human = extract_data(kraken_final, "Homo sapiens", "S")
bracken_human = extract_data(bracken_final, "Homo sapiens", "S")
print(kraken_human)
print(bracken_human)

# Bacteria Testing
kraken_bac = extract_data(kraken_final, "Bacteria", "D")
bracken_bac = extract_data(bracken_final, "Bacteria", "D")
print(kraken_bac)
print(bracken_bac)

dataframes = {
    "kraken_human": kraken_human,
    "bracken_human": bracken_human,
    "kraken_bac": kraken_bac,
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
                percentage = line_list[0].strip()
                reads = line_list[1].strip()
                taxon_level = line_list[-3].strip()
                taxon = line_list[-1].strip()

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
kraken_pathogen = extract_data_multi(kraken_final, tax_list, lvl="S")
bracken_pathogen = extract_data_multi(bracken_final, tax_list, lvl="S")
print(kraken_pathogen)
print(bracken_pathogen)

dataframes = {
    "kraken_pathogen": kraken_pathogen,
    "bracken_pathogen": bracken_pathogen,
}

# Write CSVs for all of these
output_dir = "/data/mchomm/Homm-UTI-transcriptome/misc/good_rpm/"
for name, df in dataframes.items():
    df.to_csv(f"{output_dir}{name}.csv")


# Define function to actually compute the RPM values

def compute_rpm(filepaths):
    '''
    Returns a pandas DataFrame with the total percentage and read count
    for all taxa in tax_list at the specified taxonomic level.
    '''