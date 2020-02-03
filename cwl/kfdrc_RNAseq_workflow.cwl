class: Workflow
cwlVersion: v1.0
id: kfdrc_rnaseq_wf
doc: >-
  ![data service
  logo](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS9BnbvIsTkK3QlSGMDvlgu0tZQJ1q4crMvA-S3fcWfIq6y2d2Y) 


  This is the Kids First Data Resource Center RNA-Seq Workflow, which includes
  fusion and expression detection.  This workflow takes either BAM or FASTQ
  input and uses STAR for alignment to the reference.  STAR gives aligned
  genomic, chimeric, transcriptomic, and junction read outputs.  Arriba and STAR
  fusion mode are run for fusion estimations on STAR alignment chimeric output. 
  RSEM is run for gene expression estimations on STAR transcriptomic output.
  Additional gene abundance estimations are gathered from Kallisto which is run
  on trimmed reads of the raw read input. RNAseQC is run to generate a number of
  metrics including mapping rates, transcript counts, gene counts, sense and
  antisense mapping along with others. 


  ### Gene Expression Abundance Estimation:

  [STAR v2.6.1d](https://doi.org/f4h523) is used to align paired-end RNA-seq
  reads. This output is used for all subsequent RNA analysis. The reference we
  used, and is recommended for use, was that of ENSEMBL's [GENCODE
  27](https://www.gencodegenes.org/human/release_27.html), "Comprehensive gene
  annotation." [RSEM v1.3.1](https://doi:10/cwg8n5) is used for transcript- and
  gene-level quantification. A second method of quantification was added using
  [Kallisto v0.43.1](https://doi:10.1038/nbt.3519). This method differs in that
  it uses pseudoaligments using fastq reads directly to the aforementioned
  GENCODE 27 reference.

  ### RNA Fusion Calling:

  [Arriba v1.1.0](https://github.com/suhrig/arriba/) and [STAR-Fusion
  1.5.0](https://doi:10.1101/120295) fusion detection tools are set up for
  fusion calling. For both of these tools, aligned BAM and chimeric SAM files
  are used from STAR as inputs and `GRCh38_gencode_v27` GTF for gene annotation
  is recommended. STAR-Fusion is set with default parameters and we recommend to
  annotate all fusion calls with
  `GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz` provided in the STAR-fusion
  release. For Arriba, we used a blacklist file
  `blacklist_hg38_GRCh38_2018-11-04.tsv.gz` from the Arriba release tarballs to
  remove recurrent fusion artifacts and transcripts present in healthy tissue.

  ### Tips To Run:

  1) For fastq input, please enter the reads 1 file in reads1 and the reads 2
  file in reads2. For bam input, please enter the reads file in reads1 and leave
  reads2 empty as it is optional.


  2) r1_adapter and r2_adapter are OPTIONAL.  If the input reads have already
  been trimmed, leave these as null and cutadapt step will simple pass on the
  fastq files to STAR.  If they do need trimming, supply the adapters and the
  cutadapt step will trim, and pass trimmed fastqs along


  3) `wf_strand_param` is a workflow convenience param so that, if you input the
  following, the equivalent will propagate to the four tools that use that
  parameter:

    - `default`: 'rsem_std': null, 'kallisto_std': null, 'rnaseqc_std': null, 'arriba_std': null. This means unstranded or auto in the case of arriba.
    - `rf_stranded`: 'rsem_std': 0, 'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std': 'reverse'.  This means if read1 in the input fastq/bam is reverse complement to the transcript that it maps to.
    - `fr-stranded`: 'rsem_std': 1, 'kallisto_std': 'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'. This means if read1 in the input fastq/bam is the same sense (maps 5' to 3') to the transcript that it maps to.

  4) Suggested `STAR_outSAMattrRGline`, with **TABS SEPARATING THE TAGS**, 
  format is:

    - `ID:sample_name LB:aliquot_id   PL:platform SM:BSID`, for example:
    - `ID:7316-242   LB:750189 PL:ILLUMINA SM:BS_W72364MN`

  5) Suggested inputs are:


    * FusionGenome: [GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/)
    * gtf_anno: gencode.v27.primary_assembly.annotation.gtf
    * RNAseQC_GTF: gencode.v27.primary_assembly.RNAseQC.gtf
    * RSEMgenome: RSEM_GENCODE27.tar.gz
    * STARgenome: STAR_GENCODE27.tar.gz
    * reference_fasta: GRCh38.primary_assembly.genome.fa
    * kallisto_idx: gencode.v27.kallisto.index

    
  ### Links/Resources:

  The related Github branch for this app is located
  [here](https://github.com/kids-first/kf-rnaseq-workflow/tree/be-publish).
    
label: Kids First DRC RNA-Seq Workflow
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: FusionGenome
    type: File
    doc: GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
    'sbg:suggestedValue':
      class: File
      name: GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz
      path: 5d8bb21fe4b0950c4028f854
    'sbg:x': 844.7615356445312
    'sbg:y': 723.5
  - id: RNAseQC_GTF
    type: File
    doc: gencode.v27.primary_assembly.RNAseQC.gtf
    'sbg:suggestedValue':
      class: File
      name: gencode.v27.primary_assembly.RNAseQC.gtf
      path: 5d8bb21fe4b0950c4028f852
    'sbg:x': 0
    'sbg:y': 428
  - id: RSEMgenome
    type: File
    doc: RSEM_GENCODE27.tar.gz
    'sbg:suggestedValue':
      class: File
      name: RSEM_GENCODE27.tar.gz
      path: 5d8bb21fe4b0950c4028f851
    'sbg:x': 0
    'sbg:y': 321
  - id: STAR_outSAMattrRGline
    type: string
    doc: >-
      Suggested setting, with TABS SEPARATING THE TAGS, format is:
      ID:sample_name LB:aliquot_id PL:platform SM:BSID for example ID:7316-242
      LB:750189 PL:ILLUMINA SM:BS_W72364MN
    'sbg:x': 0
    'sbg:y': 214
  - id: STARgenome
    type: File
    doc: STAR_GENCODE27.tar.gz
    'sbg:suggestedValue':
      class: File
      name: STAR_GENCODE27.tar.gz
      path: 5d8bb21fe4b0950c4028f853
    'sbg:x': 0
    'sbg:y': 107
  - id: gtf_anno
    type: File
    doc: gencode.v27.primary_assembly.annotation.gtf
    'sbg:suggestedValue':
      class: File
      name: gencode.v27.primary_assembly.annotation.gtf
      path: 5d8bb21fe4b0950c4028f84f
    'sbg:x': 0
    'sbg:y': 1284
  - id: input_type
    type:
      type: enum
      symbols:
        - BAM
        - FASTQ
      name: input_type
    doc: 'Please select one option for input file type, BAM or FASTQ.'
    'sbg:x': 0
    'sbg:y': 1177
  - id: kallisto_idx
    type: File
    doc: gencode.v27.kallisto.index
    'sbg:suggestedValue':
      class: File
      name: gencode.v27.kallisto.index
      path: 5d8bb21fe4b0950c4028f850
    'sbg:x': 0
    'sbg:y': 1070
  - id: r1_adapter
    type: string?
    doc: >-
      Optional input. If the input reads have already been trimmed, leave these
      as null. If they do need trimming, supply the adapters.
    'sbg:x': 0
    'sbg:y': 963
  - id: r2_adapter
    type: string?
    doc: >-
      Optional input. If the input reads have already been trimmed, leave these
      as null. If they do need trimming, supply the adapters.
    'sbg:x': 0
    'sbg:y': 856
  - id: reads1
    type: File
    doc: >-
      For FASTQ input, please enter reads 1 here. For BAM input, please enter
      reads here.
    'sbg:x': 0
    'sbg:y': 749
  - id: reads2
    type: File?
    doc: 'For FASTQ input, please enter reads 2 here. For BAM input, leave empty.'
    'sbg:x': 0
    'sbg:y': 642
  - id: reference_fasta
    type: File
    doc: GRCh38.primary_assembly.genome.fa
    'sbg:suggestedValue':
      class: File
      name: GRCh38.primary_assembly.genome.fa
      path: 5d8bb21fe4b0950c4028f855
    'sbg:x': 0
    'sbg:y': 535
  - id: runThread
    type: int
    doc: Amount of threads for analysis.
    'sbg:x': 539.3709106445312
    'sbg:y': 614
  - id: sample_name
    type: string
    doc: Sample name used for file base name of all outputs
    'sbg:x': 539.3709106445312
    'sbg:y': 507
  - id: wf_strand_param
    type:
      type: enum
      symbols:
        - default
        - rf-stranded
        - fr-stranded
      name: wf_strand_param
    doc: >-
      use 'default' for unstranded/auto, 'rf-stranded' if read1 in the fastq
      read pairs is reverse complement to the transcript, 'fr-stranded' if read1
      same sense as transcript
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: RNASeQC_Metrics
    outputSource:
      - rna_seqc/Metrics
    type: File
    doc: 'Metrics on mapping, intronic, exonic rates, count information, etc'
    'sbg:x': 1894.18310546875
    'sbg:y': 695.5
  - id: RNASeQC_counts
    outputSource:
      - supplemental/RNASeQC_counts
    type: File
    doc: 'Contains gene tpm, gene read, and exon counts'
    'sbg:x': 2209.035400390625
    'sbg:y': 642
  - id: RSEM_gene
    outputSource:
      - rsem/gene_out
    type: File
    doc: RSEM gene expression estimates
    'sbg:x': 1599.3170166015625
    'sbg:y': 781.5
  - id: RSEM_isoform
    outputSource:
      - rsem/isoform_out
    type: File
    doc: RSEM isoform expression estimates
    'sbg:x': 1599.3170166015625
    'sbg:y': 674.5
  - id: STAR-Fusion_results
    outputSource:
      - star_fusion/abridged_coding
    type: File
    doc: STAR fusion detection from chimeric reads
    'sbg:x': 1599.3170166015625
    'sbg:y': 139.5
  - id: STAR_chimeric_bam_out
    outputSource:
      - samtools_sort/chimeric_bam_out
    type: File
    doc: STAR bam output of chimeric reads
    'sbg:x': 1599.3170166015625
    'sbg:y': 567.5
  - id: STAR_chimeric_junctions
    outputSource:
      - star_fusion/chimeric_junction_compressed
    type: File
    doc: STAR chimeric junctions
    'sbg:x': 1599.3170166015625
    'sbg:y': 460.5
  - id: STAR_final_log
    outputSource:
      - star/log_final_out
    type: File
    doc: >-
      STAR metrics log file of unique, multi-mapping, unmapped, and chimeric
      reads
    'sbg:x': 1237.117431640625
    'sbg:y': 579
  - id: STAR_gene_count
    outputSource:
      - star/gene_counts
    type: File
    doc: STAR gene counts
    'sbg:x': 1237.117431640625
    'sbg:y': 337
  - id: STAR_junctions_out
    outputSource:
      - star/junctions_out
    type: File
    doc: STAR junction reads
    'sbg:x': 1237.117431640625
    'sbg:y': 230
  - id: STAR_sorted_genomic_bai
    outputSource:
      - samtools_sort/sorted_bai
    type: File
    doc: STAR index for sorted aligned bam
    'sbg:x': 1599.3170166015625
    'sbg:y': 353.5
  - id: STAR_sorted_genomic_bam
    outputSource:
      - samtools_sort/sorted_bam
    type: File
    doc: STAR sorted alignment bam
    'sbg:x': 1599.3170166015625
    'sbg:y': 246.5
  - id: STAR_transcriptome_bam
    outputSource:
      - star/transcriptome_bam_out
    type: File
    doc: STAR bam of transcriptome reads
    'sbg:x': 1237.117431640625
    'sbg:y': 123
  - id: arriba_fusion_results
    outputSource:
      - arriba_fusion/arriba_fusions
    type: File
    doc: Fusion output from Arriba
    'sbg:x': 1599.3170166015625
    'sbg:y': 1144.5
  - id: arriba_fusion_viz
    outputSource:
      - arriba_fusion/arriba_pdf
    type: File
    doc: pdf output from Arriba
    'sbg:x': 1599.3170166015625
    'sbg:y': 1037.5
  - id: cutadapt_stats
    outputSource:
      - cutadapt/cutadapt_stats
    type: File
    doc: 'Cutadapt stats output, only if adapter is supplied.'
    'sbg:x': 844.7615356445312
    'sbg:y': 830.5
  - id: kallisto_Abundance
    outputSource:
      - kallisto/abundance_out
    type: File
    doc: Gene abundance output from STAR genomic bam file
    'sbg:x': 1237.117431640625
    'sbg:y': 970
steps:
  - id: arriba_fusion
    in:
      - id: arriba_strand_flag
        source: strand_parse/arriba_std
      - id: chimeric_sam_out
        source: star/chimeric_sam_out
      - id: genome_aligned_bai
        source: samtools_sort/sorted_bai
      - id: genome_aligned_bam
        source: samtools_sort/sorted_bam
      - id: gtf_anno
        source: gtf_anno
      - id: outFileNamePrefix
        source: sample_name
      - id: reference_fasta
        source: reference_fasta
    out:
      - id: arriba_fusions
      - id: arriba_pdf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: arriba_fusion
      baseCommand:
        - /arriba_v1.1.0/arriba
      inputs:
        - id: arriba_strand_flag
          type: string?
        - id: chimeric_sam_out
          type: File
        - id: genome_aligned_bai
          type: File
        - id: genome_aligned_bam
          type: File
        - id: gtf_anno
          type: File
        - id: outFileNamePrefix
          type: string
        - id: reference_fasta
          type: File
      outputs:
        - id: arriba_fusions
          type: File
          outputBinding:
            glob: $(inputs.outFileNamePrefix).arriba.fusions.tsv
        - id: arriba_pdf
          type: File
          outputBinding:
            glob: $(inputs.outFileNamePrefix).arriba.fusions.pdf
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            -c $(inputs.chimeric_sam_out.path) -x
            $(inputs.genome_aligned_bam.path) -a $(inputs.reference_fasta.path)
            -g $(inputs.gtf_anno.path) -o
            $(inputs.outFileNamePrefix).arriba.fusions.tsv -O
            $(inputs.outFileNamePrefix).arriba.discarded_fusions.tsv -b
            /arriba_v1.1.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz -T
            -P ${
              if(inputs.arriba_strand_flag == null){
                return "-s auto";
              }
              else{
                return "-s " + inputs.arriba_strand_flag;
              }
            } && /arriba_v1.1.0/draw_fusions.R
            --annotation=$(inputs.gtf_anno.path)
            --fusions=$(inputs.outFileNamePrefix).arriba.fusions.tsv
            --alignments=$(inputs.genome_aligned_bam.path)
            --cytobands=/arriba_v1.1.0/database/cytobands_hg38_GRCh38_2018-02-23.tsv
            --proteinDomains=/arriba_v1.1.0/database/protein_domains_hg38_GRCh38_2018-03-06.gff3
            --output=$(inputs.outFileNamePrefix).arriba.fusions.pdf
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 64000
          coresMin: 8
        - class: DockerRequirement
          dockerPull: 'kfdrc/arriba:1.1.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 1237.117431640625
    'sbg:y': 1119
  - id: bam2fastq
    in:
      - id: SampleID
        source: sample_name
      - id: input_reads_1
        source: reads1
      - id: input_reads_2
        source: reads2
      - id: input_type
        source: input_type
      - id: runThreadN
        source: runThread
    out:
      - id: fq1
      - id: fq2
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: bam2fastq
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: SampleID
          type: string
        - id: input_reads_1
          type: File
          doc: >-
            For FASTQ input, please enter reads 1 here. For BAM input, please
            enter reads here.
        - id: input_reads_2
          type: File?
          doc: >-
            For FASTQ input, please enter reads 2 here. For BAM input, leave
            empty.
        - id: input_type
          type:
            type: enum
            symbols:
              - BAM
              - FASTQ
            name: input_type
          doc: 'Please select one option for input file type, BAM or FASTQ.'
        - id: runThreadN
          type: int
      outputs:
        - id: fq1
          type: File
          outputBinding:
            glob: '*.converted_1.fastq.gz'
        - id: fq2
          type: File
          outputBinding:
            glob: '*.converted_2.fastq.gz'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: |-
            set -eo pipefail
            ${
                if(inputs.input_type == "BAM"){
                    var command = "samtools sort -m 1G -n -O SAM -@ " + inputs.runThreadN + " " + inputs.input_reads_1.path + " | samtools fastq -c 2 -1 " + inputs.SampleID + ".converted_1.fastq.gz -2 " + inputs.SampleID + ".converted_2.fastq.gz -@ " + inputs.runThreadN+ " -"
                    return command
                }
                
                if(inputs.input_type == "FASTQ"){
                    var command =  "cp " + inputs.input_reads_1.path + " " + inputs.input_reads_1.nameroot + ".converted_1.fastq.gz && cp " + inputs.input_reads_2.path + " " + inputs.input_reads_2.nameroot + ".converted_2.fastq.gz"
                    return command
                }
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 30000
          coresMin: 36
        - class: DockerRequirement
          dockerPull: 'kfdrc/samtools:1.9'
        - class: InlineJavascriptRequirement
    'sbg:x': 250.96875
    'sbg:y': 695.5
  - id: cutadapt
    in:
      - id: r1_adapter
        source: r1_adapter
      - id: r2_adapter
        source: r2_adapter
      - id: readFilesIn1
        source: bam2fastq/fq1
      - id: readFilesIn2
        source: bam2fastq/fq2
      - id: sample_name
        source: sample_name
    out:
      - id: cutadapt_stats
      - id: trimmedReadsR1
      - id: trimmedReadsR2
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: cutadapter
      baseCommand: []
      inputs:
        - id: r1_adapter
          type: string?
        - id: r2_adapter
          type: string?
        - id: readFilesIn1
          type: File
        - id: readFilesIn2
          type: File
        - id: sample_name
          type: string
      outputs:
        - id: cutadapt_stats
          type: File?
          outputBinding:
            glob: $(inputs.sample_name).cutadapt_results.txt
        - id: trimmedReadsR1
          type: File
          outputBinding:
            glob: $("*TRIMMED." + inputs.readFilesIn1.basename)
        - id: trimmedReadsR2
          type: File
          outputBinding:
            glob: $("*TRIMMED." + inputs.readFilesIn2.basename)
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: |-
            ${ if (inputs.r1_adapter == null){
              var cmd = "cp " + inputs.readFilesIn1.path + " ./UNTRIMMED." + inputs.readFilesIn1.basename
              + ";cp " + inputs.readFilesIn2.path + " ./UNTRIMMED." + inputs.readFilesIn2.basename;
              return cmd;
            } else{
              var cmd = "cutadapt -j 16 -m 20 --quality-base=33 -q 20 -a " + inputs.r1_adapter
              + " -A " + inputs.r2_adapter + " -o TRIMMED." + inputs.readFilesIn1.basename
              + " -p TRIMMED." + inputs.readFilesIn2.basename + " " + inputs.readFilesIn1.path + " " +
              inputs.readFilesIn2.path + " > " + inputs.sample_name + ".cutadapt_results.txt";
              return cmd;
            } }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 24000
          coresMin: 16
        - class: DockerRequirement
          dockerPull: 'kfdrc/cutadapt:latest'
        - class: InlineJavascriptRequirement
    'sbg:x': 539.3709106445312
    'sbg:y': 749
  - id: kallisto
    in:
      - id: SampleID
        source: sample_name
      - id: reads1
        source: cutadapt/trimmedReadsR1
      - id: reads2
        source: cutadapt/trimmedReadsR2
      - id: strand
        source: strand_parse/kallisto_std
      - id: transcript_idx
        source: kallisto_idx
    out:
      - id: abundance_out
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: kallisto
      baseCommand:
        - kallisto
      inputs:
        - id: SampleID
          type: string
        - id: reads1
          type: File
        - id: reads2
          type: File
        - id: strand
          type: string?
          doc: 'input none if unstranded, otherwise rf-stranded or fr-stranded'
        - id: transcript_idx
          type: File
      outputs:
        - id: abundance_out
          type: File
          outputBinding:
            glob: '*.abundance.tsv.gz'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            ${
              var cmd = "quant -i " + inputs.transcript_idx.path + " -o output -b 10 -t 8";
              if (inputs.strand != null && inputs.strand != "default"){
                cmd += " --" + inputs.strand;
              }
              cmd += " " + inputs.reads1.path + " " + inputs.reads2.path;
              return cmd;
            }

            mv output/abundance.tsv $(inputs.SampleID).kallisto.abundance.tsv &&
            gzip $(inputs.SampleID).kallisto.abundance.tsv
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 10000
          coresMin: 8
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/uros_sipetic/kallisto:0.43.1'
        - class: InlineJavascriptRequirement
    'sbg:x': 844.7615356445312
    'sbg:y': 588.5
  - id: rna_seqc
    in:
      - id: Aligned_sorted_bam
        source: samtools_sort/sorted_bam
      - id: collapsed_gtf
        source: RNAseQC_GTF
      - id: strand
        source: strand_parse/rnaseqc_std
    out:
      - id: Exon_count
      - id: Gene_TPM
      - id: Gene_count
      - id: Metrics
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: rnaseqc
      baseCommand:
        - rnaseqc
      inputs:
        - id: Aligned_sorted_bam
          type: File
        - id: collapsed_gtf
          type: File
        - id: strand
          type: string?
      outputs:
        - id: Exon_count
          type: File
          outputBinding:
            glob: output/*.exon_reads.gct
        - id: Gene_TPM
          type: File
          outputBinding:
            glob: output/*.gene_tpm.gct
        - id: Gene_count
          type: File
          outputBinding:
            glob: output/*.gene_reads.gct
        - id: Metrics
          type: File
          outputBinding:
            glob: output/*.metrics.tsv
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.collapsed_gtf.path) $(inputs.Aligned_sorted_bam.path)
            output/ ${
              var cmd = "--legacy";
              if (inputs.strand != null && inputs.strand != "default"){
                cmd += " --stranded=" + inputs.strand;
              }
              return cmd;
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 10000
          coresMin: 8
        - class: DockerRequirement
          dockerPull: 'gcr.io/broad-cga-aarong-gtex/rnaseqc:latest'
        - class: InlineJavascriptRequirement
    'sbg:x': 1599.3170166015625
    'sbg:y': 909.5
  - id: rsem
    in:
      - id: bam
        source: star/transcriptome_bam_out
      - id: genomeDir
        source: RSEMgenome
      - id: outFileNamePrefix
        source: sample_name
      - id: strandedness
        source: strand_parse/rsem_std
    out:
      - id: gene_out
      - id: isoform_out
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: rsem_calculate_expression
      baseCommand:
        - tar
      inputs:
        - id: bam
          type: File
        - id: genomeDir
          type: File
        - id: outFileNamePrefix
          type: string
        - id: strandedness
          type: string?
          doc: 'Options relative to upstream reads - none, forward, reverse'
      outputs:
        - id: gene_out
          type: File
          outputBinding:
            glob: '*genes.results.gz'
        - id: isoform_out
          type: File
          outputBinding:
            glob: '*isoforms.results.gz'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: |-
            -zxf $(inputs.genomeDir.path) && ${
              var cmd = "rsem-calculate-expression --paired-end --alignments --append-names --no-bam-output -p 16";
              if (inputs.strandedness != null && inputs.strandedness != "default"){
                cmd += " --strandedness " + inputs.strandedness;
              }
              cmd += " " + inputs.bam.path + " ./" + inputs.genomeDir.nameroot.split('.')[0] + "/" + inputs.genomeDir.nameroot.split('.')[0] + " " +  inputs.outFileNamePrefix + ".rsem";
              return cmd
            } && gzip *results
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 24000
          coresMin: 16
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/uros_sipetic/rsem:1.3.1'
        - class: InlineJavascriptRequirement
    'sbg:x': 1237.117431640625
    'sbg:y': 842
  - id: samtools_sort
    in:
      - id: chimeric_sam_out
        source: star/chimeric_sam_out
      - id: unsorted_bam
        source: star/genomic_bam_out
    out:
      - id: chimeric_bam_out
      - id: sorted_bai
      - id: sorted_bam
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: samtools_sort
      baseCommand:
        - samtools
      inputs:
        - id: chimeric_sam_out
          type: File
        - id: unsorted_bam
          type: File
      outputs:
        - id: chimeric_bam_out
          type: File
          outputBinding:
            glob: $(inputs.chimeric_sam_out.nameroot).bam
        - id: sorted_bai
          type: File
          outputBinding:
            glob: '*.sorted.bai'
        - id: sorted_bam
          type: File
          outputBinding:
            glob: '*.sorted.bam'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            sort $(inputs.unsorted_bam.path) -@ 16 -m 1G -O bam >
            $(inputs.unsorted_bam.nameroot).sorted.bam && samtools index -@ 16
            $(inputs.unsorted_bam.nameroot).sorted.bam
            $(inputs.unsorted_bam.nameroot).sorted.bai && samtools view -bh -@
            16 $(inputs.chimeric_sam_out.path) -o
            $(inputs.chimeric_sam_out.nameroot).bam
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 24000
          coresMin: 16
        - class: DockerRequirement
          dockerPull: 'kfdrc/samtools:1.9'
        - class: InlineJavascriptRequirement
    'sbg:x': 1237.117431640625
    'sbg:y': 700
  - id: star
    in:
      - id: genomeDir
        source: STARgenome
      - id: outFileNamePrefix
        source: sample_name
      - id: outSAMattrRGline
        source: STAR_outSAMattrRGline
      - id: readFilesIn1
        source: cutadapt/trimmedReadsR1
      - id: readFilesIn2
        source: cutadapt/trimmedReadsR2
      - id: runThreadN
        source: runThread
    out:
      - id: chimeric_junctions
      - id: chimeric_sam_out
      - id: gene_counts
      - id: genomic_bam_out
      - id: junctions_out
      - id: log_final_out
      - id: log_out
      - id: log_progress_out
      - id: transcriptome_bam_out
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: star_align_reads
      baseCommand:
        - tar
        - '-xzf'
      inputs:
        - id: genomeDir
          type: File
        - id: outFileNamePrefix
          type: string
        - id: outSAMattrRGline
          type: string
        - id: readFilesIn1
          type: File
        - id: readFilesIn2
          type: File
        - id: runThreadN
          type: int
      outputs:
        - id: chimeric_junctions
          type: File
          outputBinding:
            glob: '*Chimeric.out.junction'
        - id: chimeric_sam_out
          type: File
          outputBinding:
            glob: '*Chimeric.out.sam'
        - id: gene_counts
          type: File
          outputBinding:
            glob: '*ReadsPerGene.out.tab.gz'
        - id: genomic_bam_out
          type: File
          outputBinding:
            glob: '*Aligned.out.bam'
        - id: junctions_out
          type: File
          outputBinding:
            glob: '*SJ.out.tab.gz'
        - id: log_final_out
          type: File
          outputBinding:
            glob: '*Log.final.out'
        - id: log_out
          type: File
          outputBinding:
            glob: '*Log.out'
        - id: log_progress_out
          type: File
          outputBinding:
            glob: '*Log.progress.out'
        - id: transcriptome_bam_out
          type: File
          outputBinding:
            glob: '*Aligned.toTranscriptome.out.bam'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.genomeDir.path) && STAR --outSAMattrRGline
            $(inputs.outSAMattrRGline) --genomeDir
            ./$(inputs.genomeDir.nameroot.split('.')[0])/ --readFilesIn
            $(inputs.readFilesIn1.path) $(inputs.readFilesIn2.path)
            --readFilesCommand zcat --runThreadN $(inputs.runThreadN)
            --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin
            8 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5
            --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1
            --alignIntronMax 100000 --chimSegmentReadGapMax 3
            --chimOutJunctionFormat 1 --alignMatesGapMax 100000 --outFilterType
            BySJout --outFilterScoreMinOverLread 0.33
            --outFilterMatchNminOverLread 0.33 --outReadsUnmapped None
            --limitSjdbInsertNsj 1200000 --outFileNamePrefix
            $(inputs.outFileNamePrefix). --outSAMstrandField intronMotif
            --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes
            --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted
            --outSAMunmapped Within --genomeLoad NoSharedMemory --chimSegmentMin
            12 --chimJunctionOverhangMin 12 --chimOutType Junctions
            SeparateSAMold WithinBAM SoftClip --chimMainSegmentMultNmax 1
            --outSAMattributes NH HI AS nM NM ch && gzip *ReadsPerGene.out.tab 
            *SJ.out.tab
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 60000
          coresMin: 16
        - class: DockerRequirement
          dockerPull: 'kfdrc/star:latest'
        - class: InlineJavascriptRequirement
    'sbg:x': 844.7615356445312
    'sbg:y': 397.5
  - id: star_fusion
    in:
      - id: Chimeric_junction
        source: star/chimeric_junctions
      - id: SampleID
        source: sample_name
      - id: genomeDir
        source: FusionGenome
    out:
      - id: abridged_coding
      - id: chimeric_junction_compressed
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: star_fusion
      baseCommand:
        - tar
      inputs:
        - id: Chimeric_junction
          type: File
        - id: SampleID
          type: string
        - id: genomeDir
          type: File
      outputs:
        - id: abridged_coding
          type: File
          outputBinding:
            glob: '*.fusion_predictions.abridged.coding_effect.tsv'
        - id: chimeric_junction_compressed
          type: File
          outputBinding:
            glob: $(inputs.Chimeric_junction.basename).gz
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            -zxf $(inputs.genomeDir.path) && /usr/local/STAR-Fusion/STAR-Fusion
            --genome_lib_dir
            ./GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir -J
            $(inputs.Chimeric_junction.path) --output_dir STAR-Fusion_outdir
            --examine_coding_effect --CPU 16 && mv
            STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv
            $(inputs.SampleID).STAR.fusion_predictions.abridged.coding_effect.tsv
            && gzip -c $(inputs.Chimeric_junction.path) >
            $(inputs.Chimeric_junction.basename).gz
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 64000
          coresMin: 16
        - class: DockerRequirement
          dockerPull: 'kfdrc/star-fusion:1.5.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 1237.117431640625
    'sbg:y': 458
  - id: strand_parse
    in:
      - id: wf_strand_param
        source: wf_strand_param
    out:
      - id: arriba_std
      - id: kallisto_std
      - id: rnaseqc_std
      - id: rsem_std
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: >-
        ${ var strand = 'default'; if (inputs.wf_strand_param != null){ strand =
        inputs.wf_strand_param; } var parse_dict = { 'default': {'rsem_std':
        'none', 'kallisto_std': 'default', 'rnaseqc_std': 'default',
        'arriba_std': 'auto'}, 'rf-stranded': {'rsem_std': 'reverse',
        'kallisto_std': 'rf-stranded', 'rnaseqc_std': 'rf', 'arriba_std':
        'reverse'}, 'fr-stranded': {'rsem_std': 'forward', 'kallisto_std':
        'fr-stranded', 'rnaseqc_std': 'fr', 'arriba_std': 'yes'} }; if (strand
        in parse_dict){ return parse_dict[strand];

        } else{ throw new Error(strand + ' is a not a valid strand param. Use
        one of default, rf-stranded, fr-stranded'); } }
      id: expression_strand_params
      inputs:
        wf_strand_param:
          doc: >-
            use 'default' for unstranded/auto, rf_stranded if read1 in the fastq
            read pairs is reverse complement to the transcript, fr-stranded if
            read1 same sense as transcript
          type:
            - name: wf_strand_param
              symbols:
                - default
                - rf-stranded
                - fr-stranded
              type: enum
      outputs:
        arriba_std:
          type: string
        kallisto_std:
          type: string
        rnaseqc_std:
          type: string
        rsem_std:
          type: string
      requirements:
        - class: InlineJavascriptRequirement
    'sbg:x': 250.96875
    'sbg:y': 539.5
  - id: supplemental
    in:
      - id: Exon_count
        source: rna_seqc/Exon_count
      - id: Gene_TPM
        source: rna_seqc/Gene_TPM
      - id: Gene_count
        source: rna_seqc/Gene_count
      - id: outFileNamePrefix
        source: sample_name
    out:
      - id: RNASeQC_counts
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: supplemental_tar_gz
      baseCommand:
        - mkdir
      inputs:
        - id: Exon_count
          type: File
        - id: Gene_TPM
          type: File
        - id: Gene_count
          type: File
        - id: outFileNamePrefix
          type: string
      outputs:
        - id: RNASeQC_counts
          type: File
          outputBinding:
            glob: $(inputs.outFileNamePrefix).RNASeQC.counts.tar.gz
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.outFileNamePrefix)_RNASeQC_counts

            cp $(inputs.Gene_TPM.path) $(inputs.Gene_count.path)
            $(inputs.Exon_count.path) $(inputs.outFileNamePrefix)_RNASeQC_counts

            tar -czf $(inputs.outFileNamePrefix).RNASeQC.counts.tar.gz
            $(inputs.outFileNamePrefix)_RNASeQC_counts
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 1600
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'ubuntu:18.04'
        - class: InlineJavascriptRequirement
    'sbg:x': 1894.18310546875
    'sbg:y': 567.5
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 3
requirements: []
