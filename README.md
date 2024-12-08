# DraVAR
DraVAR (Drug Resistance Variant Analysis and Reporting) is a pipeline for analyzing genomic data to identify variants associated with drug resistance. The pipeline integrates tools for quality control, alignment, variant calling, and annotation.

Features

    Automated variant calling and annotation.
    Integrates widely used bioinformatics tools.
    Flexible and customizable for different reference genomes and datasets.

Installation
Install Dependencies

The pipeline requires the following software. Install them via conda: <br/>
`conda create -n draVAR_env -c bioconda -c conda-forge \
    fastqc trimmomatic bwa samtools bcftools freebayes snpeff vcfpy pandas
conda activate draVAR_env`

`fastqc --version
trimmomatic -version
bwa
samtools --version
bcftools --version
freebayes --version
python -m vcfpy --help`

Usage
Input Files

The pipeline requires three input files: <br/>

    Reference Genome (`ref.gb`): A GenBank file for the reference genome.
    Read Files (`reads_1.fastq.gz` and `reads_2.fastq.gz`): Paired-end sequencing data.

Run the Pipeline

Run the pipeline with the following command:<br/>
`sh DraVAR.sh --ref 'ref.gb' --R1 reads_1.fastq.gz --R2 reads_2.fastq.gz --prefix prefix --output out`
Parameters <br/>

    --ref: Path to the reference genome in .gb or .gbk format.
    --R1: Path to the first paired-end read file (reads_1.fastq.gz).
    --R2: Path to the second paired-end read file (reads_2.fastq.gz).
    --prefix: Prefix for naming intermediate and result files.
    --output: Directory to store the results.
Contribution

Feel free to fork the repository and submit pull requests to improve the pipeline.

For questions or issues, contact Ankita Pal at ap578817@gmail.com.

Example <br/>
`DraVAR.sh --ref 'AL123456.gb' --R1 'reads_1.fastq.gz' --R2 'reads_2.fastq.gz' --prefix 'example' --output 'results'`

    

