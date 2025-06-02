######################################################################################
# Author: Max Homm
# Date: May 30, 2025
# Purpose: Actually run salmon on all the species
######################################################################################

from pathogen import *
import glob

# Good patients to include: 
good_patients = [
    "10005", "10006", "10010", "10016", "10021", "10022", "10043", "10048", "10050", "10065", "10108", "10117", "10136", "10145", "10160", "10166", "10168", "10199", "10201", "10206",
    "10220", "10224", "10254", "10257", "10263", "10286", "10292", "10308", "10315", "10318", "10346", "10354", "10365", "10369", "10377", "10380", "10383", "10385", "10389", "10393",
    "10395", "10397", "10399", "10422", "10423", "10428", "10430", "10436", "10443", "10445", "10448", "10454", "10458", "10459", "10466", "10468", "10470", "10487", "10490", "10491",
    "10494", "10506", "10511", "10514", "10519", "10523", "10524", "10530", "10532", "10533", "10539", "10542", "10550", "10554", "10555", "10581", "10584", "10596", "10608", "10622",
    "10625", "10627", "10644", "10650", "10656", "10659", "10660", "10667", "10674", "10675", "10676", "10678", "10685", "10687", "10703", "10709", "10717", "10720", "10759", "10763",
    "10803", "10811", "10815", "10817", "10821", "10831", "10833", "10836", "10852", "10867", "10869", "10871", "10877", "10883", "10891", "10892", "10896", "10897", "10900", "10901",
    "10913", "10927", "10940", "10967", "10969", "10973", "10992", "11003", "11007", "11009", "11010", "11016", "11023", "11045", "11047", "11053", "11077", "11084", "11104", "11113",
    "11114", "11116", "11118", "11121", "11122", "11149", "11162", "11170", "11173", "11196", "11200", "11206", "11222", "11225", "11227", "11230", "11233", "11248", "11256", "11272",
    "11275", "11281", "11285", "11297", "11311", "11314", "11317", "11323", "11329", "11331", "11351", "11363", "11370", "11394", "11397", "15860", "15887", "15900", "15906", "15927",
    "15963", "15985", "15991", "15994", "16003", "16018", "16050", "16089", "16133", "16148", "16151", "16190", "20002", "20004", "20005", "20014", "20018", "20023", "20036", "20042",
    "20048", "20054", "20062", "20068", "20070", "20077", "20079", "20082", "20090", "20097", "20103", "20115", "20122", "20123", "20134", "20140", "20142", "20147", "20149", "20150",
    "20159", "20161", "20163", "32194", "32253", "32270", "32271", "32279", "32300", "32388", "32471", "32521", "32536", "32546", "32558", "32566", "32579", "32607", "32644"
]

# Create the dictionary of 2-element lists for each patient
def generate_paths(paths, patients): 
    dict_of_patients = {}
    for batch in paths: 
        for file in glob.glob(f"{batch}*.fastq.gz"): 
            ptid = file.split("/")[-1].split("_")[0]
            if ptid not in patients: 
                continue
            if ptid in dict_of_patients: 
                dict_of_patients[ptid].append(file)
            else: 
                dict_of_patients[ptid] = [file]
    
    # Sort so it's {id: [R1, R2]}
    sorted_dict = {k: sorted(v) for k, v in dict_of_patients.items()}
    return(sorted_dict)

# Get a dict of all of the paths
paths = ["/data/mchomm/pathogen_surveillance/data/batch4_fastp/", 
         "/data/mchomm/pathogen_surveillance/data/batch3_fastp/", 
         "/data/mchomm/pathogen_surveillance/data/batch2_fastp/", 
         "/data/mchomm/pathogen_surveillance/data/batch1_fastp/"]
fastp_files = generate_paths(paths=paths, patients=good_patients)


# Create your pathogen objects

# K. pneumoniae
kpn_output = "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/kn/salmon_output"
kpn = Pathogen("Klebsiella pneumoniae", "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/kn/salmon_index")
kpn.build_index("/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/kn/cds/ncbi_dataset/data/GCF_000240185.1/cds_from_genomic.fna")
kpn.run_salmon(fastp_files, kpn_output)

# E. faecalis
efs_output = "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/ef/salmon_output"
efs = Pathogen("Enterococcus faecalis", "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/ef/salmon_index")
efs.build_index("/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/ef/cds/ncbi_dataset/data/GCF_000393015.1/cds_from_genomic.fna")
efs.run_salmon(fastp_files, efs_output)

# P. aeruginosa
paa_output = "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/pa/salmon_output"
paa = Pathogen("Pseudomonas aeruginosa", "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/pa/salmon_index")
paa.build_index("/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/pa/cds/ncbi_dataset/data/GCF_000006765.1/cds_from_genomic.fna")
paa.run_salmon(fastp_files, paa_output)

# P. mirabilis
pms_output = "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/pm/salmon_output"
pms = Pathogen("Proteus mirabilis", "/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/pm/salmon_index")
pms.build_index("/data/mchomm/Homm-UTI-transcriptome/data/bacterial_transcriptomics/pm/cds/ncbi_dataset/data/GCF_000069965.1/cds_from_genomic.fna")
pms.run_salmon(fastp_files, pms_output)