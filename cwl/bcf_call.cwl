cwlVersion: v1.2
class: Workflow
label: KFDRC NGS Checkmate Preprocess
doc: |-
  # BCF Filter Tool
  Preprocessing workflow to use bcftools to subset bams and create a bcftools-called vcf

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  ## bcf_call.cwl

  Creates input vcfs for ngs checkmate. Especially useful to run when inputs are large WGS bam files.

  ### inputs
  ```yaml
  inputs:
    input_align: File[]
    chr_list: File
    reference_fasta: File
    snp_bed: File
  ```
  Suggested inputs:
  ```text
  chr_list: chr_list.txt
  snp_bed: SNP_hg38_liftover_wChr.bed
  reference_fasta: Homo_sapiens_assembly38.fasta
  ```
  ### outputs
  ```yaml
  bcf_called_vcf: {type: File[], outputSource: bcf_filter/bcf_call}
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement

inputs:
- id: input_align
  type:
    type: array
    items: File
  secondaryFiles:
  - pattern: .crai
    required: false
  - pattern: .bai
    required: false
  - pattern: ^.bai
    required: false
- id: chr_list
  type: File
  sbg:suggestedValue:
    name: chr_list.txt
    class: File
    path: 5f50018fe4b054958bc8d2e2
- id: reference_fasta
  type: File
  secondaryFiles:
  - pattern: .fai
    required: true
  sbg:suggestedValue:
    name: Homo_sapiens_assembly38.fasta
    class: File
    secondaryFiles:
    - name: Homo_sapiens_assembly38.fasta.fai
      class: File
      path: 60639016357c3a53540ca7af
    path: 60639014357c3a53540ca7a3
- id: snp_bed
  type: File
  sbg:suggestedValue:
    name: SNP_hg38_liftover_wChr.bed
    class: File
    path: 5f50018fe4b054958bc8d2e4

outputs:
- id: bcf_called_vcf
  type:
    type: array
    items: File
  outputSource: bcf_filter/bcf_call

steps:
- id: bcf_filter
  in:
  - id: input_align
    source: input_align
  - id: chr_list
    source: chr_list
  - id: reference_fasta
    source: reference_fasta
  - id: snp_bed
    source: snp_bed
  scatter:
  - input_align
  run:
    cwlVersion: v1.2
    class: CommandLineTool

    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/ngscheckmate:1.3
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 4
      ramMin: 4000

    inputs:
    - id: input_align
      type: File
      secondaryFiles:
      - pattern: .crai
        required: false
      - pattern: .bai
        required: false
      - pattern: ^.bai
        required: false
    - id: chr_list
      type: File
    - id: reference_fasta
      type: File
      secondaryFiles:
      - .fai
    - id: snp_bed
      type: File

    outputs:
    - id: bcf_call
      type: File
      outputBinding:
        glob: $(inputs.input_align.nameroot).bcf.called.vcf

    baseCommand:
    - cat
    arguments:
    - position: 1
      valueFrom: |-
        $(inputs.chr_list.path) | xargs -ICH -P 4 sh -c 'bcftools mpileup --output /dev/stdout --fasta-ref $(inputs.reference_fasta.path) -r CH -T $(inputs.snp_bed.path) $(inputs.input_align.path) > $(inputs.input_align.nameroot).CH.pileup.vcf' && find . -not -empty -name '*.pileup.vcf' > pileup_list.txt && vcf-concat -f pileup_list.txt > $(inputs.input_align.nameroot).merged.vcf && bcftools call -c $(inputs.input_align.nameroot).merged.vcf > $(inputs.input_align.nameroot).bcf.called.vcf
      shellQuote: false
    id: bcf_filter
  out:
  - bcf_call

hints:
- class: sbg:AWSInstanceType
  value: c5.9xlarge;ebs-gp2;850
- class: sbg:maxNumberOfParallelInstances
  value: 4
id: |-
  https://cavatica-api.sbgenomics.com/v2/apps/brownm28/ngs-checkmate-dev/bcf-call/8/raw/
sbg:appVersion:
- v1.2
sbg:categories:
- NGSCHECKMATE
- PREPROCESS
sbg:content_hash: a5ee1b8f0cae5032bcc8348e6795bc2ab02735cda1099b55c6cc74cddec7c30c8
sbg:contributors:
- brownm28
sbg:createdBy: brownm28
sbg:createdOn: 1552682992
sbg:id: brownm28/ngs-checkmate-dev/bcf-call/8
sbg:image_url: |-
  https://cavatica.sbgenomics.com/ns/brood/images/brownm28/ngs-checkmate-dev/bcf-call/8.png
sbg:latestRevision: 8
sbg:license: Apache License 2.0
sbg:links:
- id: https://github.com/kids-first/ngs_checkmate_wf/releases/tag/v1.0.0
  label: github-release
sbg:modifiedBy: brownm28
sbg:modifiedOn: 1650466217
sbg:project: brownm28/ngs-checkmate-dev
sbg:projectName: ngs_checkmate_dev
sbg:publisher: sbg
sbg:revision: 8
sbg:revisionNotes: |-
  Uploaded using sbpack v2021.10.07. 
  Source: 
  repo: https://github.com/kids-first/ngs_checkmate_wf
  file: workflows/bcf_call.cwl
  commit: (uncommitted file)
sbg:revisionsInfo:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1552682992
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650399002
  sbg:revision: 1
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: b18a59e
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650399779
  sbg:revision: 2
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: 3a45a50
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650400341
  sbg:revision: 3
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: (uncommitted file)
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650400542
  sbg:revision: 4
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: (uncommitted file)
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650403480
  sbg:revision: 5
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: 38c4aaa
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650404050
  sbg:revision: 6
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: (uncommitted file)
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650465779
  sbg:revision: 7
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: (uncommitted file)
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650466217
  sbg:revision: 8
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/bcf_call.cwl
    commit: (uncommitted file)
sbg:sbgMaintained: false
sbg:validationErrors: []
sbg:workflowLanguage: CWL
