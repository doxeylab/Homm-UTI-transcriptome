# Running Kraken2 + Bracken on all patient files

# Load all required packages

module load StdEnv/2020
module load gcc/9.3.0
module load kraken2/2.1.3
module load bracken/2.7
module load python

# Initialize file locations

kraken_out="" # Replace with desired kraken output location
kraken_db="" # Replace with location of kraken db

data="" # Replace with location of fastq files (fastp filtered)
reference_file="" # Path to a TSV file containing <sample ID>\t<file name for fastq r1>\t<file name for fastq r2>

# Iterate through each pair of fastq files, perform kraken2 + bracken, and place output in kraken_out

while IFS=$'\t' read -r sample r1 r2

do
  kraken2 --report $kraken_out/${sample}.kraken --confidence 0.5 --threads 20 --paired --gzip-compressed -db $kraken_db $data/$r1 $data/$r2 > $kraken_out/$sample.krakenRaw
  bracken -d $kraken_db -i $kraken_out/${sample}.kraken -o $kraken_out/${sample}.braken.out -t 0

done < $reference_file
