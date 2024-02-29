#!/bin/bash
mamba init
mamba activate fastp
mkdir -p ./data/02.Clean_data
mkdir -p ./data/02.Clean_data/report
run_fastp() {
    file="$1"
    fastp -i "$file" -I "${file%_1.fastq.gz}_2.fastq.gz" -o "./data/02.Clean_data/$(basename $file)" -O "./data/02.Clean_data/$(basename ${file%_1.fastq.gz}_2.fastq.gz)" -h "./data/02.Clean_data/report/$(basename ${file%_1.fastq.gz}).html" -j "./data/02.Clean_data/report/$(basename ${file%_1.fastq.gz}).json"
}

export -f run_fastp

find ./data/01.Raw_data -name "*_1.fastq.gz" | parallel --progress --keep-order --line-buffer run_fastp
# run fastqc on the cleaned data, output to ./data/02.Clean_data/fastqc
mkdir -p ./data/02.Clean_data/fastqc
fastqc -t 20 ./data/02.Clean_data/*.fastq.gz -o ./data/02.Clean_data/fastqc
/data0/apps/anaconda3/bin/multiqc ./data/02.Clean_data/fastqc -o ./data/02.Clean_data/fastqc/report
