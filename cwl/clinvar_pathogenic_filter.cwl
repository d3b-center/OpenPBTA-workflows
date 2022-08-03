cwlVersion: v1.2
class: Workflow
label: clinvar-pathogenic-filter
doc: |-
  Small workflow to run snpEff to annotate a file with a vcf input, then use bcftools to filter output based on supplied include and exclude criterion
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

inputs:
- id: input_vcf
  doc: VCF file (with TBI) to be annotated
  type: File
  sbg:x: 0
  sbg:y: 221.0333251953125
- id: db_file
  doc: Reference database for annotating the input
  type: File
  sbg:x: 0
  sbg:y: 773.6166381835938
- id: mode
  doc: Mode of SnpSift to run
  type:
    name: mode
    type: enum
    symbols:
    - annotate
    - dbnsfp
    - gwasCat
  sbg:x: 0
  sbg:y: 110.51666259765625
- id: db_name
  doc: Name of the database being used
  type: string
  sbg:x: 0
  sbg:y: 663.0999755859375
- id: fields
  doc: |-
    Comma-separated list of fields from the database that will be used as annotations
  type: string?
  sbg:x: 0
  sbg:y: 442.066650390625
- id: include_expression
  type: string?
  sbg:x: 0
  sbg:y: 331.54998779296875
- id: exclude_expression
  type: string?
  sbg:x: 0
  sbg:y: 552.5833129882812
- id: output_basename
  doc: String that will be used in the output filenames
  type: string
  sbg:x: 0
  sbg:y: 0

outputs:
- id: filtered_vcf
  type: File
  outputSource:
  - bcftools_filter_vcf/filtered_vcf
  sbg:x: 789.61669921875
  sbg:y: 386.8083190917969

steps:
- id: snpsift_annotate
  label: SnpSift
  in:
  - id: mode
    source: mode
  - id: db_file
    source: db_file
  - id: db_name
    source: db_name
  - id: fields
    source: fields
  - id: input_vcf
    source: input_vcf
  - id: output_basename
    source: output_basename
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: SnpSift
    doc: |
      Simplified descrition of what this tool does:
        1. Run SnpEff on input VCF
        2. BGZIP output VCF
        3. TABIX ouptut VCF

      SnpSift Parameters:
        1. v: Verbose logging
        2. a: In cases where variants do not have annotations in the database the tool will report an empty value (FIELD=.) rather than no annotation
        3. info: Comma separated list of INFO fields from the reference database the tool will use to annotate matching listings
        4. f: Same as "info" but specficially used when the tool is used running in dbnsfp mode
        5. db: The reference database file used for annotating the input
        6. tabix: VCF database is tabix-indexed

      An example run of this tool will use a command like this:
        /bin/bash -c
        set -eo pipefail 
        java -jar /snpEff/SnpSift.jar 
          annotate 
          -v 
          -a 
          -info fields-string-value 
          -db /path/to/db_file.ext 
          -tabix
          /path/to/input_vcf.ext | 
        bgzip -c > output_basename-string-value.SnpSift.db_name-string-value.snpEff.vcf.gz && 
        tabix output_basename-string-value.SnpSift.db_name-string-value.snpEff.vcf.gz
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: kfdrc/snpeff:4_3t
    - class: InlineJavascriptRequirement

    inputs:
    - id: mode
      doc: Mode of SnpSift to run
      type:
        name: mode
        type: enum
        symbols:
        - annotate
        - dbnsfp
        - gwasCat
    - id: db_file
      doc: Reference database for annotating the input
      type: File
      secondaryFiles:
      - .tbi
    - id: db_name
      doc: Name of the database being used
      type: string
    - id: fields
      doc: |-
        Comma-separated list of fields from the database that will be used as annotations
      type: string?
    - id: input_vcf
      doc: VCF file (with TBI) to be annotated
      type: File
      secondaryFiles:
      - .tbi
    - id: output_basename
      doc: String that will be used in the output filenames
      type: string

    outputs:
    - id: output_vcf
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.vcf.gz'

    baseCommand:
    - /bin/bash
    - -c
    arguments:
    - position: 1
      valueFrom: |-
        set -eo pipefail
        java -jar /snpEff/SnpSift.jar $(inputs.mode) -v ${ if (inputs.mode == 'gwasCat') {return ''} else {return '-a -tabix'}} ${ if (inputs.fields && inputs.mode == 'annotate') {
          return '-info ' + inputs.fields
        } else if (inputs.fields && inputs.mode == 'dbnsfp') {
          return '-f ' + inputs.fields
        } else {
          return ''
        }} -db $(inputs.db_file.path) $(inputs.input_vcf.path) | bgzip -c > $(inputs.output_basename).SnpSift.$(inputs.db_name).snpEff.vcf.gz && tabix $(inputs.output_basename).SnpSift.$(inputs.db_name).snpEff.vcf.gz
      shellQuote: false
    id: d3b-bixu/dev-wgsa/snpsift-annotate/5
    sbg:appVersion:
    - v1.0
    sbg:content_hash: a70d600e752ae35519dc63cab16ed055cd1ad347228a8679d4fec6d0ecce950a0
    sbg:contributors:
    - danmiller
    - brownm28
    sbg:createdBy: danmiller
    sbg:createdOn: 1585074760
    sbg:id: d3b-bixu/dev-wgsa/snpsift-annotate/5
    sbg:image_url:
    sbg:latestRevision: 5
    sbg:modifiedBy: brownm28
    sbg:modifiedOn: 1659463516
    sbg:project: d3b-bixu/dev-wgsa
    sbg:projectName: KF Annotation
    sbg:publisher: sbg
    sbg:revision: 5
    sbg:revisionNotes: |-
      Uploaded using sbpack v2022.03.16. 
      Source: 
      repo: https://github.com/kids-first/kf-annotation
      file: tools/snpsift_annotate.cwl
      commit: 6d2b7a2
    sbg:revisionsInfo:
    - sbg:modifiedBy: danmiller
      sbg:modifiedOn: 1585074760
      sbg:revision: 0
      sbg:revisionNotes:
    - sbg:modifiedBy: danmiller
      sbg:modifiedOn: 1585074987
      sbg:revision: 1
      sbg:revisionNotes: fix stuff
    - sbg:modifiedBy: danmiller
      sbg:modifiedOn: 1585075169
      sbg:revision: 2
      sbg:revisionNotes: question mark
    - sbg:modifiedBy: danmiller
      sbg:modifiedOn: 1585335952
      sbg:revision: 3
      sbg:revisionNotes: documented
    - sbg:modifiedBy: danmiller
      sbg:modifiedOn: 1585592565
      sbg:revision: 4
      sbg:revisionNotes: tabix flag
    - sbg:modifiedBy: brownm28
      sbg:modifiedOn: 1659463516
      sbg:revision: 5
      sbg:revisionNotes: |-
        Uploaded using sbpack v2022.03.16. 
        Source: 
        repo: https://github.com/kids-first/kf-annotation
        file: tools/snpsift_annotate.cwl
        commit: 6d2b7a2
    sbg:sbgMaintained: false
    sbg:validationErrors: []
    sbg:workflowLanguage: CWL
  out:
  - id: output_vcf
  sbg:x: 213.5166778564453
  sbg:y: 351.67498779296875
- id: bcftools_filter_vcf
  label: bcftools-filter-vcf
  in:
  - id: input_vcf
    source: snpsift_annotate/output_vcf
  - id: include_expression
    source: include_expression
  - id: exclude_expression
    source: exclude_expression
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: bcftools-filter-vcf
    doc: |-
      More generic tool to take in an include expression and optionally an exclude expresssion to filter a vcf
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 1000
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest
    - class: InlineJavascriptRequirement

    inputs:
    - id: input_vcf
      type: File
    - id: include_expression
      type: string?
    - id: exclude_expression
      type: string?
    - id: output_basename
      type: string?

    outputs:
    - id: filtered_vcf
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.vcf.gz'

    baseCommand: []
    arguments:
    - position: 1
      valueFrom: |-
        ${
          var out_base = inputs.output_basename;
          if (out_base == null){
            out_base = inputs.input_vcf.nameroot + ".bcf_filtered"
          }
          var cmd = "bcftools view ";
          if (inputs.include_expression != null){
              cmd += "--include '" + inputs.include_expression + "' " + inputs.input_vcf.path;
              if (inputs.exclude_expression != null){
                  cmd += " | bcftools view --exclude '" + inputs.exclude_expression + "' -O z > " + out_base + ".vcf.gz;";
              } else {
                  cmd += " -O z > " + out_base + ".vcf.gz;";
              }
          } else if (inputs.include_expression == null && inputs.exclude_expression != null){
              cmd += "--exclude '" + inputs.exclude_expression + "' " + inputs.input_vcf.path + " -O z > " + out_base + ".vcf.gz;";
          } else if (inputs.include_expression == null && inputs.exclude_expression == null){
              cmd = "cp " + inputs.input_vcf.path + " ./" + out_base + ".vcf.gz;";
          }
          cmd += "tabix " + out_base + ".vcf.gz;"
          return cmd;
        }
      shellQuote: false
    id: d3b-bixu/dev-wgsa/bcftools-filter-vcf/1
    sbg:appVersion:
    - v1.0
    sbg:content_hash: a33e9a472508d42a1215eba6e90a7fc7f072421e7d06e1107b34af42ee1271c5a
    sbg:contributors:
    - brownm28
    sbg:createdBy: brownm28
    sbg:createdOn: 1583350187
    sbg:id: d3b-bixu/dev-wgsa/bcftools-filter-vcf/1
    sbg:image_url:
    sbg:latestRevision: 1
    sbg:modifiedBy: brownm28
    sbg:modifiedOn: 1659550506
    sbg:project: d3b-bixu/dev-wgsa
    sbg:projectName: KF Annotation
    sbg:publisher: sbg
    sbg:revision: 1
    sbg:revisionNotes: |-
      Uploaded using sbpack v2022.03.16. 
      Source: 
      repo: https://github.com/kids-first/kf-somatic-workflow
      file: tools/bcftools_filter_vcf.cwl
      commit: 1.1.1-306-gf6fd590
    sbg:revisionsInfo:
    - sbg:modifiedBy: brownm28
      sbg:modifiedOn: 1583350187
      sbg:revision: 0
      sbg:revisionNotes:
    - sbg:modifiedBy: brownm28
      sbg:modifiedOn: 1659550506
      sbg:revision: 1
      sbg:revisionNotes: |-
        Uploaded using sbpack v2022.03.16. 
        Source: 
        repo: https://github.com/kids-first/kf-somatic-workflow
        file: tools/bcftools_filter_vcf.cwl
        commit: 1.1.1-306-gf6fd590
    sbg:sbgMaintained: false
    sbg:validationErrors: []
    sbg:workflowLanguage: CWL
  out:
  - id: filtered_vcf
  sbg:x: 505.16668701171875
  sbg:y: 372.8083190917969
id: |-
  https://cavatica-api.sbgenomics.com/v2/apps/d3b-bixu/dev-wgsa/clinvar-pathogenic-filter/3/raw/
sbg:appVersion:
- v1.2
- v1.0
sbg:content_hash: a5819db524533c30d8f95839ac7369dbc929a19df333cf572e10fbe8dfd14df7c
sbg:contributors:
- brownm28
sbg:createdBy: brownm28
sbg:createdOn: 1659551966
sbg:id: d3b-bixu/dev-wgsa/clinvar-pathogenic-filter/3
sbg:image_url: |-
  https://cavatica.sbgenomics.com/ns/brood/images/d3b-bixu/dev-wgsa/clinvar-pathogenic-filter/3.png
sbg:latestRevision: 3
sbg:modifiedBy: brownm28
sbg:modifiedOn: 1659552623
sbg:project: d3b-bixu/dev-wgsa
sbg:projectName: KF Annotation
sbg:publisher: sbg
sbg:revision: 3
sbg:revisionNotes: added output basename, some descriptors
sbg:revisionsInfo:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1659551966
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1659552319
  sbg:revision: 1
  sbg:revisionNotes: initial creation
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1659552446
  sbg:revision: 2
  sbg:revisionNotes: open needed ports
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1659552623
  sbg:revision: 3
  sbg:revisionNotes: added output basename, some descriptors
sbg:sbgMaintained: false
sbg:toolAuthor: Miguel Brown
sbg:validationErrors: []
sbg:workflowLanguage: CWL
