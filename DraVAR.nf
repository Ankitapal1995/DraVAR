#!/usr/bin/env nextflow

// Define parameters
params.ref = null
params.reads1 = null
params.reads2 = null
params.prefix = 'output'
params.output_dir = './results'
params.truseq_adapters = '/data1/ankita/variant_calling/TruSeq3-PE-2.fa'

// Workflow definition
workflow {
    // Create output directory
    Channel.fromPath(params.output_dir).ifEmpty { file(params.output_dir).mkdirs() }

    // Generate FASTA from GenBank
    processGenBankToFasta(params.ref) | generateFastaIndex

    // Perform quality checks and trimming
    processFastQC(params.reads1, params.reads2) | processTrimmomatic(params.reads1, params.reads2)

    // Align reads and process alignments
    processBwaMem(params.reads1, params.reads2, params.ref) | processSamtools

    // Variant calling and filtering
    processFreeBayes | processBcftoolsFilter

    // Variant annotation
    processSnpEff | processAnnotateCSV

    // Final outputs: consensus sequence and compressed VCF
    processBgzipAndIndex | processConsensus
}

// Process definitions

process processGenBankToFasta {
    input:
    path ref_file from Channel.fromPath(params.ref)

    output:
    path "*.fasta" into fasta_channel

    script:
    """
    python genebank_to_fasta.py ${ref_file}
    """
}

process generateFastaIndex {
    input:
    path fasta_file from fasta_channel

    output:
    path "*.fai"

    script:
    """
    samtools faidx ${fasta_file}
    """
}

process processFastQC {
    input:
    path reads1 from Channel.fromPath(params.reads1)
    path reads2 from Channel.fromPath(params.reads2)

    output:
    path "*.html", emit: fastqc_output

    script:
    """
    fastqc ${reads1} ${reads2} -o ${params.output_dir}
    """
}

process processTrimmomatic {
    input:
    path reads1 from Channel.fromPath(params.reads1)
    path reads2 from Channel.fromPath(params.reads2)

    output:
    path "${params.prefix}_trimmed_R1.fastq.gz"
    path "${params.prefix}_trimmed_R2.fastq.gz"

    script:
    """
    trimmomatic PE \\
        ${reads1} ${reads2} \\
        ${params.output_dir}/${params.prefix}_trimmed_R1.fastq.gz \\
        ${params.output_dir}/${params.prefix}_unpaired_R1.fastq.gz \\
        ${params.output_dir}/${params.prefix}_trimmed_R2.fastq.gz \\
        ${params.output_dir}/${params.prefix}_unpaired_R2.fastq.gz \\
        ILLUMINACLIP:${params.truseq_adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process processBwaMem {
    input:
    path reads1 from Channel.fromPath(params.reads1)
    path reads2 from Channel.fromPath(params.reads2)
    path ref_fasta from fasta_channel

    output:
    path "*.sam"

    script:
    """
    bwa index ${ref_fasta}
    bwa mem -M -R "@RG\\tID:${params.prefix}\\tLB:${params.prefix}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${params.prefix}" \\
        ${ref_fasta} ${reads1} ${reads2} > ${params.output_dir}/${params.prefix}_aligned_reads.sam
    """
}

process processSamtools {
    input:
    path sam_file from processBwaMem.out.collect()

    output:
    path "${params.prefix}_deduplicated.bam"

    script:
    """
    samtools view -S -b ${sam_file} > ${params.output_dir}/${params.prefix}_aligned.bam
    samtools sort ${params.output_dir}/${params.prefix}_aligned.bam -o ${params.output_dir}/${params.prefix}_sorted.bam
    samtools sort -n ${params.output_dir}/${params.prefix}_sorted.bam -o ${params.output_dir}/${params.prefix}_query_sorted.bam
    samtools fixmate -m ${params.output_dir}/${params.prefix}_query_sorted.bam ${params.output_dir}/${params.prefix}_fixed.bam
    samtools sort ${params.output_dir}/${params.prefix}_fixed.bam -o ${params.output_dir}/${params.prefix}_sorted_fixed.bam
    samtools markdup -r ${params.output_dir}/${params.prefix}_sorted_fixed.bam ${params.output_dir}/${params.prefix}_deduplicated.bam
    """
}

process processFreeBayes {
    input:
    path bam_file from processSamtools.out.collect()

    output:
    path "${params.prefix}_raw_variants.vcf"

    script:
    """
    freebayes-parallel <(fasta_generate_regions.py ${params.ref}.fai 100000) 36 -f ${params.ref} ${bam_file} > ${params.output_dir}/${params.prefix}_raw_variants.vcf
    """
}

process processBcftoolsFilter {
    input:
    path vcf_file from processFreeBayes.out.collect()

    output:
    path "${params.prefix}_filtered_variants.vcf"

    script:
    """
    bcftools filter --include 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0' ${vcf_file} > ${params.output_dir}/${params.prefix}_filtered_variants.vcf
    """
}

process processSnpEff {
    input:
    path vcf_file from processBcftoolsFilter.out.collect()

    output:
    path "${params.prefix}_annotated_variants.vcf"

    script:
    """
    snpEff -noStats -no-downstream -no-upstream -no-utr Mycobacterium_tuberculosis_h37rv ${vcf_file} > ${params.output_dir}/${params.prefix}_annotated_variants.vcf
    """
}

process processAnnotateCSV {
    input:
    path vcf_file from processSnpEff.out.collect()

    output:
    path "${params.prefix}_annotated_AMR_graded_variants.csv"

    script:
    """
    python vcf_to_table_format.py ${vcf_file}
    python vcf_to_amr.py ${params.output_dir}/${params.prefix}_annotated_variants.csv > ${params.output_dir}/${params.prefix}_annotated_AMR_graded_variants.csv
    """
}

process processBgzipAndIndex {
    input:
    path vcf_file from processSnpEff.out.collect()

    output:
    path "${params.prefix}_annotated_variants.vcf.gz"

    script:
    """
    bgzip -c ${vcf_file} > ${params.output_dir}/${params.prefix}_annotated_variants.vcf.gz
    bcftools index ${params.output_dir}/${params.prefix}_annotated_variants.vcf.gz
    """
}

process processConsensus {
    input:
    path vcf_file from processBgzipAndIndex.out.collect()

    output:
    path "${params.prefix}_consensus.fasta"

    script:
    """
    bcftools consensus -f ${params.ref} -o ${params.output_dir}/${params.prefix}_consensus.fasta ${vcf_file}
    """
}

