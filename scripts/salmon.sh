# Running Salmon on all patient files

# Load all required packages

module load StdEnv/2020
module load gcc/9.3.0
module load openmpi/4.1.5
module load salmon/1.10.2

# Initialize file locations

salmon_out="" # Replace with desired salmon output location
salmon_index="" # Replace with location of salmon index (must create beforehand with desired reference transcriptome)
data="" # Replace with location of fastq files (fastp filtered)
reference_file="" # Path to a TSV file containing <sample ID>\t<file name for fastq r1>\t<file name for fastq r2>

# Iterate through each pair of fastq files, run Salmon

while IFS=$'\t' read -r sample r1 r2; do
    salmon quant -i "$salmon_index" -l A -1 "$dir/$r1" -2 "$dir/$r2" --validateMappings --seqBias --gcBias -p 12 -o "$salmon_out/$sample"

done < "reference_file"
