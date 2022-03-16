# OpenPBTA-workflows
[![DOI](https://zenodo.org/badge/225902612.svg)](https://zenodo.org/badge/latestdoi/225902612)

This is a central repository to organize publication-related workflows for maximum reproducibility of results.

## Quick summary of usages

### cwl/kfdrc_combined_somatic_wgs_cnv_wf.cwl
 - Assay types: DNA, WGS
 - Analysis type: CNV
 - Major tools: ControlFreeC v11.6, CNVKit, v0.9.3
 - Major inputs: fasta reference, aligned cram/bam, b allele frequency file

### cwl/kfdrc_strelka2_mutect2_manta_workflow.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations, structural variants
 - Major tools: strelka2 v2.9.3, mutect2 (GATK) v4.1.1.0, manta v1.4.0
 - Major inputs: fasta reference, aligned cram/bam

### cwl/kfdrc_RNAseq_workflow.cwl
 - Assay types: DNA
 - Analysis types: RNA expression, RNA fusion
 - Major tools: STAR aligner v2.6.1d, rsem v1.3.1, STAR-Fusion v1.5.0, Arriba v1.1.0
 - Major inputs: fasta reference, GENCODE gtf, Star Index, Fusion Reference, input reads as fastq or bam

### cwl/kfdrc-alignment-cram-only-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *tumor* reads as *BAM*

### cwl/kfdrc-alignment-fq-input.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *normal* reads as *fastq*

### cwl/kfdrc-alignment-fq-input-cram-only-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *tumor* reads as *fastq*

### cwl/kfdrc-alignment-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *normal* reads as *BAM*

### cwl/kfdrc-lancet-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations
 - Major tools: lancet v1.0.7
 - Major inputs: fasta reference, aligned cram/bam

### cwl/kfdrc-mendqc-wf.cwl
 - Assay types: RNA
 - Analysis types: QC
 - Major tools: sambamba v0.6.7, UCSC-Treehouse MendQC
 - Major inputs: UCSC bed file

### cwl/kfdrc-mutect2_strelka2-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations
 - Major tools: strelka2 v2.9.3, mutect2 (GATK) v4.1.1.0
 - Major inputs: fasta reference, aligned cram/bam

### cwl/kfdrc-single-genotype-basic.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Germline mutations
 - Major tools: GATK v4.0.12.0
 - Major inputs: fasta reference, aligned cram/bam

### cwl/kfdrc-vardict-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations
 - Major tools: VarDict Java v1.5.8
 - Major inputs: fasta reference, aligned cram/bam
 
 ### cwl/kfdrc_annot_vcf_sub_wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations
 - Major tools: VEP, bcftools
 - Major inputs: pass VCF

### bash/run-gistic.sh
 - Assay types: DNA, WGS
 - Analysis types: gene-level copy number
 - Major tools: GISTIC v.2.0.23
 - Major inputs: CNVkit seg file, hg38.UCSC.add_miR.160920.refgene.mat

### bash/run_gistic_consensus.sh
 - Assay types: DNA, WGS
 - Analysis types: gene-level copy number
 - Major tools: GISTIC v.2.0.23
 - Major inputs: OpenPBTA consensus seg file, hg38.UCSC.add_miR.160920.refgene.mat
