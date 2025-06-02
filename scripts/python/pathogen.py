######################################################################################
# Author: Max Homm
# Date: May 29, 2025
# Purpose: Create Pathogen class to run salmon on species of choice
######################################################################################

# Packages

import os
import subprocess

# Pathogen class

class Pathogen:

    def __init__(self, name, index_path):
        """
        Represents a pathogen and its reference index.

        Parameters:
        - name (str): Species name (e.g., 'Klebsiella pneumoniae')
        - index_path (str): Path to Salmon index for this species
        """
        self.name = name
        self.index_path = index_path

    def __str__(self):
        return f"Pathogen(name='{self.name}', index='{self.index_path}')"
    
    def build_index(self, transcriptome_fasta):
        """
        Build a Salmon index for this pathogen from a transcriptome FASTA file.

        Parameters:
        - transcriptome_fasta (str): Path to the FASTA file of transcript sequences
        - kmer_size (int): K-mer size to use for the index (default = 31)
        """
        if os.path.exists(self.index_path):
            print(f"Index already exists at {self.index_path}, skipping.")
            return
        
        os.makedirs(self.index_path, exist_ok=True)
        cmd = [
            "salmon", "index",
            "-t", transcriptome_fasta,
            "-i", self.index_path
        ]
        print(f"Building Salmon index for {self.name}...")
        subprocess.run(cmd, check=True)
        print(f"Index built at {self.index_path}")
    
    def run_salmon(self, patients, output, threads = 4): 
        """
        Run Salmon quantification for all patients against this pathogen.

        Parameters:
        - patients (dict): patient_id -> [R1, R2] FASTQ file list
        - output_root (str): Root output dir (will create one per patient)
        - threads (int): Number of threads to use per job
        """
        ptid_count = 0
        total_ptid_count = len(patients)

        for patient_id, fastq_files in patients.items():
            ptid_count += 1
            print(f"Processing patient {ptid_count}/{total_ptid_count}")

            output_dir = os.path.join(output, self.name.replace(" ", "_"), patient_id)
            quant_file = os.path.join(output_dir, "quant.sf")

            # Skip if quantification already exists
            if os.path.exists(quant_file):
                print(f"Skipping {patient_id} - quant.sf already exists.")
                continue

            os.makedirs(output_dir, exist_ok=True)

            cmd = [
                "salmon", "quant",
                "-i", self.index_path,
                "-l", "A",
                "-1", fastq_files[0],
                "-2", fastq_files[1],
                "-p", str(threads),
                "-o", output_dir, 
                "--validateMappings", 
                "--seqBias", 
                "--gcBias", 
                "--minAssignedFrags", "0"
            ]
            print(f"Running Salmon for {patient_id} vs {self.name}")

            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Salmon failed for {patient_id}: {e}")
                continue