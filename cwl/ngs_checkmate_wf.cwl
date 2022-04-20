cwlVersion: v1.0
class: Workflow
label: KFDRC NGS Checkmate Sample QC
doc: |
  # ngs checkmate workflow

  ## Introduction
  Based on the tool from https://github.com/parklab/NGSCheckMate, "NGSCheckMate uses depth-dependent correlation models of allele fractions of known single-nucleotide polymorphisms (SNPs) to identify samples from the same individual." 

  ![data service logo](https://github.com/d3b-center/d3b-research-workflows/raw/master/doc/kfdrc-logo-sm.png)

  ### ngs_checkmate_wf.cwl
  Runs ngscheckmate in vcf mode - requires output from bcf_call step.

  #### inputs
  ```yaml
  inputs:
    input_vcf:
      type:
          type: array
          items:
              type: array
              items: File
    
    snp_bed: File
    output_basename: string[]
    ram: 
      type: ['null', int]
      default: 4000
  ```
  1) Ram input param optional - use if you plan on batching ~20+ vcfs.
  2) `input_vcf` is an array of arrays - basically an array of groups of vcfs that you'd like ot see evaluated together.
  3) `output_basename` is an array of file output prefixes - should line up with the first level of array elements from `input_vcf`

  #### outputs
  ```yaml
  outputs:
    match_results: {type: 'File[]', outputSource: ngs_checkmate/match_results}
    correlation_matrix: {type: 'File[]', outputSource: ngs_checkmate/correlation_matrix}
  ```
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement

inputs:
- id: input_vcf
  type:
    type: array
    items:
      type: array
      items: File
- id: snp_bed
  type: File
  sbg:suggestedValue:
    name: SNP_hg38_liftover_wChr.bed
    class: File
    path: 5f50018fe4b054958bc8d2e4
- id: output_basename
  type:
    type: array
    items: string
- id: ram
  type:
  - 'null'
  - int
  default: 4000

outputs:
- id: match_results
  type:
    type: array
    items: File
  outputSource: ngs_checkmate/match_results
- id: correlation_matrix
  type:
    type: array
    items: File
  outputSource: ngs_checkmate/correlation_matrix

steps:
- id: ngs_checkmate
  in:
  - id: input_vcf
    source: input_vcf
  - id: snp_bed
    source: snp_bed
  - id: output_basename
    source: output_basename
  - id: ram
    source: ram
  scatter:
  - input_vcf
  - output_basename
  scatterMethod: dotproduct
  run:
    cwlVersion: v1.0
    class: CommandLineTool

    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/ngscheckmate:1.3
    - class: InitialWorkDirRequirement
      listing: $(inputs.input_vcf)
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 2
      ramMin: $(inputs.ram)

    inputs:
    - id: input_vcf
      type:
        type: array
        items: File
    - id: snp_bed
      type: File
    - id: output_basename
      type: string
    - id: ram
      type:
      - 'null'
      - int
      default: 4000

    outputs:
    - id: match_results
      type: File
      outputBinding:
        glob: $(inputs.output_basename)_all.txt
    - id: correlation_matrix
      type: File
      outputBinding:
        glob: '$(inputs.output_basename)_corr_matrix.txt '

    baseCommand:
    - python
    arguments:
    - position: 1
      valueFrom: |-
        /NGSCheckMate/ncm.py -V -d ./ -bed $(inputs.snp_bed.path) -O ./ -N $(inputs.output_basename) && mv output_corr_matrix.txt $(inputs.output_basename)_corr_matrix.txt
      shellQuote: false
    id: ngs_checkmate
  out:
  - match_results
  - correlation_matrix

hints:
- class: sbg:AWSInstanceType
  value: c5.9xlarge;ebs-gp2;400
- class: sbg:maxNumberOfParallelInstances
  value: 4
id: |-
  https://cavatica-api.sbgenomics.com/v2/apps/brownm28/ngs-checkmate-dev/ngs-checkmate-wf/10/raw/
sbg:appVersion:
- v1.0
sbg:categories:
- NGSCHECKMATE
- QC
sbg:content_hash: a9f93af553e5928d6ad0b8d359a66a315d501c2322e639b221c706cb73d317c21
sbg:contributors:
- brownm28
sbg:createdBy: brownm28
sbg:createdOn: 1552920917
sbg:id: brownm28/ngs-checkmate-dev/ngs-checkmate-wf/10
sbg:image_url: |-
  https://cavatica.sbgenomics.com/ns/brood/images/brownm28/ngs-checkmate-dev/ngs-checkmate-wf/10.png
sbg:latestRevision: 10
sbg:license: Apache License 2.0
sbg:links:
- id: https://github.com/kids-first/ngs_checkmate_wf/releases/tag/v1.0.0
  label: github-release
sbg:modifiedBy: brownm28
sbg:modifiedOn: 1650402050
sbg:project: brownm28/ngs-checkmate-dev
sbg:projectName: ngs_checkmate_dev
sbg:publisher: sbg
sbg:revision: 10
sbg:revisionNotes: |-
  Uploaded using sbpack v2021.10.07. 
  Source: 
  repo: https://github.com/kids-first/ngs_checkmate_wf
  file: workflows/ngs_checkmate_wf.cwl
  commit: 38c4aaa
sbg:revisionsInfo:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1552920917
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1552921257
  sbg:revision: 1
  sbg:revisionNotes: test output basename string
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1552922064
  sbg:revision: 2
  sbg:revisionNotes: testing nested dot product
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1552923118
  sbg:revision: 3
  sbg:revisionNotes:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1552923494
  sbg:revision: 4
  sbg:revisionNotes:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1552950740
  sbg:revision: 5
  sbg:revisionNotes:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1553007277
  sbg:revision: 6
  sbg:revisionNotes: troubleshoot ram defaults
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1553007825
  sbg:revision: 7
  sbg:revisionNotes: swtich to dotproduct for scatter
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650394406
  sbg:revision: 8
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/ngs_checkmate_wf.cwl
    commit: 017a66c
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650398680
  sbg:revision: 9
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/ngs_checkmate_wf.cwl
    commit: b18a59e
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650402050
  sbg:revision: 10
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/ngs_checkmate_wf
    file: workflows/ngs_checkmate_wf.cwl
    commit: 38c4aaa
sbg:sbgMaintained: false
sbg:validationErrors: []
sbg:workflowLanguage: CWL
