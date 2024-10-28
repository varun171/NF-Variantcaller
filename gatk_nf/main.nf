#!/usr/bin/env nextflow

params.reads = '/path/to/raw/reads/*_R{1,2}.fastq.gz'
params.reference = '/path/to/reference/genome.fasta'
params.dbSNP = '/path/to/dbsnp.vcf'
params.output = './results'

process fastqc {
    input:
    file fastq from channel.fromFilePairs(params.reads)

    output:
    file("*.zip") into qc

    script:
    """
    fastqc ${fastq}
    """
}

process bwa_mem {
    input:
    file read1 from reads.first
    file read2 from reads.second

    output:
    file("*.bam") into bams

    script:
    """
    bwa mem -t 16 ${params.reference} ${read1} ${read2} | samtools view -Sb - > ${sample}.bam
    samtools sort -@ 8 -o ${sample}_sorted.bam ${sample}.bam
    """
}

process mark_duplicates {
    input:
    file bam from bams

    output:
    file("*.marked.bam") into dedup_bams

    script:
    """
    gatk MarkDuplicates -I ${bam} -O ${sample}_marked.bam -M ${sample}_marked.metrics
    """
}

process base_recalibration {
    input:
    file bam from dedup_bams

    output:
    file("*.recal.bam") into recal_bams

    script:
    """
    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${params.reference} \
        --known-sites ${params.dbSNP} \
        -O ${sample}_recal_data.table
    gatk ApplyBQSR \
        -I ${bam} \
        -R ${params.reference} \
        -O ${sample}_recal.bam \
        -bqsr ${sample}_recal_data.table
    """
}

process haplotype_caller {
    input:
    file recal_bam from recal_bams

    output:
    file("*.g.vcf") into gvcfs

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.reference} \
        -I ${recal_bam} \
        -O ${sample}.g.vcf.gz \
        -ERC GVCF
    """
}

process joint_genotyping {
    input:
    file gvcf from gvcfs.collect()

    output:
    file("*.vcf") into raw_vcfs

    script:
    """
    gatk GenotypeGVCFs \
        -R ${params.reference} \
        -V ${gvcf} \
        -O ${sample}.vcf.gz
    """
}

process variant_filtering {
    input:
    file vcf from raw_vcfs

    output:
    file("*.filtered.vcf.gz") into filtered_vcfs

    script:
    """
    gatk VariantFiltration \
        -R ${params.reference} \
        -V ${vcf} \
        -O ${sample}_filtered.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
        --filter-name "basic_snp_filter"
    """
}

process variant_annotation {
    input:
    file vcf from filtered_vcfs

    output:
    file("*.annotated.vcf") into annotated_vcfs

    script:
    """
    gatk Funcotator \
        -R ${params.reference} \
        -V ${vcf} \
        --output ${sample}_annotated.vcf.gz \
        --output-file-format VCF
    """
}

workflow {
    reads = Channel.fromFilePairs(params.reads)
    fastqc(reads)
    bwa_mem(reads)
    mark_duplicates(bams)
    base_recalibration(dedup_bams)
    haplotype_caller(recal_bams)
    joint_genotyping(gvcfs)
    variant_filtering(raw_vcfs)
    variant_annotation(filtered_vcfs)
}
