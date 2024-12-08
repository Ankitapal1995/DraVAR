#!/bin/bash

# Print usage
usage() {
    echo "Usage: $0 --ref <reference.gb> --R1 <reads_1.fastq.gz> --R2 <reads_2.fastq.gz> --prefix <prefix> --output <output_dir>"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --ref) REF="$2"; shift ;;
        --R1) READS1="$2"; shift ;;
        --R2) READS2="$2"; shift ;;
        --prefix) PREFIX="$2"; shift ;;
        --output) OUTPUT="$2"; shift ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Check if mandatory arguments are provided
if [[ -z "$REF" || -z "$READS1" || -z "$READS2" || -z "$PREFIX" || -z "$OUTPUT" ]]; then
    echo "Error: Missing arguments."
    usage
fi
START_TIME=$(date +%s)

echo "Checking all the files and folders for variant calling annotation....."

sh snpeff_file_checking.sh

# Create output directory
mkdir -p "$OUTPUT"

# Run the pipeline
set -e  # Exit if any command fails

echo "Running genebank_to_fasta.py..."
python genebank_to_fasta.py "$REF"

if [[ "$REF" == *.gbk ]]; then
    BASE_NAME=$(basename "$REF" .gbk)
elif [[ "$REF" == *.gb ]]; then
    BASE_NAME=$(basename "$REF" .gb)
else
    echo "Error: File must have .gb or .gbk extension."
    exit 1
fi


# Define the output file
REF_fasta="${BASE_NAME}.fasta"

echo "Running FastQC..."
fastqc "$READS1" "$READS2" -o "$OUTPUT"

echo "Running Trimmomatic..."
trimmomatic PE \
    "$READS1" "$READS2" \
    "$OUTPUT/${PREFIX}_trimmed_R1.fastq.gz" "$OUTPUT/${PREFIX}_unpaired_R1.fastq.gz" \
    "$OUTPUT/${PREFIX}_trimmed_R2.fastq.gz" "$OUTPUT/${PREFIX}_unpaired_R2.fastq.gz" \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

echo "Running BWA..."
bwa index "$REF_fasta"
bwa mem -M -R "@RG\tID:${PREFIX}\tLB:${PREFIX}\tPL:ILLUMINA\tPM:HISEQ\tSM:${PREFIX}" \
    "$REF_fasta" "$OUTPUT/${PREFIX}_trimmed_R1.fastq.gz" "$OUTPUT/${PREFIX}_trimmed_R2.fastq.gz" > "$OUTPUT/${PREFIX}_aligned_reads.sam"

echo "Converting SAM to BAM and sorting..."
samtools view -S -b "$OUTPUT/${PREFIX}_aligned_reads.sam" > "$OUTPUT/${PREFIX}_aligned.bam"
samtools sort "$OUTPUT/${PREFIX}_aligned.bam" -o "$OUTPUT/${PREFIX}_sorted.bam"
samtools sort -n "$OUTPUT/${PREFIX}_sorted.bam" -o "$OUTPUT/${PREFIX}_query_sorted.bam"

echo "Running Samtools Fixmate..."
samtools fixmate -m "$OUTPUT/${PREFIX}_query_sorted.bam" "$OUTPUT/${PREFIX}_fixed.bam"


samtools sort "$OUTPUT/${PREFIX}_fixed.bam" -o "$OUTPUT/${PREFIX}_sorted_fixed.bam"

echo "Marking duplicates..."
samtools markdup -r "$OUTPUT/${PREFIX}_sorted_fixed.bam" "$OUTPUT/${PREFIX}_deduplicated.bam"

echo "Calling variants with FreeBayes..."
#freebayes -p 2 -P 0 -C 2 -F 0.05 --min-coverage 10 --min-repeat-entropy 1.0 -q 13 -m 60 --strict-vcf -f "$REF_fasta" "$OUTPUT/${PREFIX}_deduplicated.bam" > "$OUTPUT/${PREFIX}_raw_variants.vcf"

samtools index "$OUTPUT/${PREFIX}_deduplicated.bam"


freebayes-parallel <(fasta_generate_regions.py data/AL123456.fasta.fai 100000) 36 -f "$REF_fasta" "$OUTPUT/${PREFIX}_deduplicated.bam" > "$OUTPUT/${PREFIX}_raw_variants.vcf"
echo "Filtering variants with BCFtools..."
bcftools filter --include 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0' "$OUTPUT/${PREFIX}_raw_variants.vcf" -o "$OUTPUT/${PREFIX}_filtered_variants.vcf"

echo "Annotating variants with SnpEff..."
snpEff -noStats -no-downstream -no-upstream -no-utr Mycobacterium_tuberculosis_h37rv "$OUTPUT/${PREFIX}_filtered_variants.vcf" > "$OUTPUT/${PREFIX}_annotated_variants.vcf"

echo "VCF to table...."

python vcf_to_table_format.py "$OUTPUT/${PREFIX}_annotated_variants.vcf"

echo "Compressing and indexing VCF..."
bgzip -c "$OUTPUT/${PREFIX}_annotated_variants.vcf" > "$OUTPUT/${PREFIX}_annotated_variants.vcf.gz"
bcftools index "$OUTPUT/${PREFIX}_annotated_variants.vcf.gz"

echo "Generating consensus sequence..."
bcftools consensus -f "$REF_fasta" -o "$OUTPUT/${PREFIX}_consensus.fasta" "$OUTPUT/${PREFIX}_annotated_variants.vcf.gz"

echo "Pipeline completed successfully. Results are in the $OUTPUT directory."
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Total time taken: $DURATION seconds"
