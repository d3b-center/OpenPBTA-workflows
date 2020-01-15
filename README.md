# publication_workflows
This is a central repository to organize publication-related workflows for maximum reproducibility of results.

## Quick summary of usages

### openPBTA/kfdrc_combined_somatic_wgs_cnv_wf.cwl
 - Assay types: DNA, WGS
 - Analysis type: CNV
 - Major tools: ControlFreeC v11.6, CNVKit, v0.9.3
 - Major inputs: fasta reference, aligned cram/bam, b allele frequency file

### openPBTA/kfdrc_strelka2_mutect2_manta_workflow.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations, structural variants
 - Major tools: strelka2 v2.9.3, mutect2 (GATK) v4.1.1.0, manta v1.4.0
 - Major inputs: fasta reference, aligned cram/bam

### openPBTA/kfdrc_RNAseq_workflow.cwl
 - Assay types: DNA
 - Analysis types: RNA expression, RNA fusion
 - Major tools: STAR aligner v2.6.1d, rsem v1.3.1, STAR-Fusion v1.5.0, Arriba v1.1.0
 - Major inputs: fasta reference, GENCODE gtf, Star Index, Fusion Reference, input reads as fastq or bam

### openPBTA/kfdrc-alignment-cram-only-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *tumor* reads as *BAM*

### openPBTA/kfdrc-alignment-fq-input.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *normal* reads as *fastq*

### openPBTA/kfdrc-alignment-fq-input-cram-only-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *tumor* reads as *fastq*

### openPBTA/kfdrc-alignment-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: NGS read alignment
 - Major tools: bwa v0.7.17, sambamba v0.6.3, GATK v4.0.3.0
 - Major inputs: fasta reference, bwa indices, input *normal* reads as *BAM*

### openPBTA/kfdrc-lancet-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations
 - Major tools: lancet v1.0.7
 - Major inputs: fasta reference, aligned cram/bam

### openPBTA/kfdrc-mendqc-wf.cwl
 - Assay types: RNA
 - Analysis types: QC
 - Major tools: sambamba v0.6.7, UCSC-Treehouse MendQC
 - Major inputs: UCSC bed file

### openPBTA/kfdrc-mutect2_strelka2-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations
 - Major tools: strelka2 v2.9.3, mutect2 (GATK) v4.1.1.0
 - Major inputs: fasta reference, aligned cram/bam

### openPBTA/kfdrc-single-genotype-basic.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Germline mutations
 - Major tools: GATK v4.0.12.0
 - Major inputs: fasta reference, aligned cram/bam

### openPBTA/kfdrc-vardict-wf.cwl
 - Assay types: DNA, WGS, WXS
 - Analysis types: Somatic mutations
 - Major tools: VarDict Java v1.5.8
 - Major inputs: fasta reference, aligned cram/bam
