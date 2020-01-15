class: Workflow
cwlVersion: v1.0
id: kfdrc-harmonization/sd-bhjxbdqk/kfdrc-mendqc-wf/0
doc: >-
  RNAseq QC UCSC Treehouse. Adapted from
  https://github.com/UCSC-Treehouse/mend_qc. Calculates the number of Mapped
  Exonic Non-Duplicate (MEND) reads in a bam file
label: kfdrc-mendqc-wf
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: gencode_bed
    type: File
    doc: >-
      Basic gencode bed file matching ref used for aligned bam. Obtain from
      UCSC.
    'sbg:x': 192.703125
    'sbg:y': 160.5
  - id: input_bam
    type: File
    doc: RNA aligned to genome bam
    'sbg:x': 0
    'sbg:y': 160.5
  - id: output_basename
    type: string
    'sbg:x': 0
    'sbg:y': 53.5
outputs:
  - id: mend_qc_json
    outputSource:
      - mend_qc/output_json
    type: File
    'sbg:x': 748.10009765625
    'sbg:y': 214
  - id: mend_qc_readDist
    outputSource:
      - mend_qc/output_readDist
    type: File
    'sbg:x': 748.10009765625
    'sbg:y': 107
  - id: mend_qc_tsv
    outputSource:
      - mend_qc/output_tsv
    type: File
    'sbg:x': 748.10009765625
    'sbg:y': 0
  - id: out_md_bam
    outputSource:
      - preprocess_bam/out_md_bam
    type: File
    doc: 'Sorted by coord, dup marked bam'
    'sbg:x': 459.0577697753906
    'sbg:y': 39.5
steps:
  - id: mend_qc
    in:
      - id: gencode_bed
        source: gencode_bed
      - id: input_md_bam
        source: preprocess_bam/out_md_bam
      - id: output_basename
        source: output_basename
    out:
      - id: output_json
      - id: output_readDist
      - id: output_tsv
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: mend_qc
      baseCommand:
        - /usr/local/bin/read_distribution.py
      inputs:
        - id: gencode_bed
          type: File
          doc: >-
            Basic gencode bed file matching ref used for aligned bam. Obtain
            from UCSC.
        - id: input_md_bam
          type: File
          secondaryFiles:
            - .bai
        - id: output_basename
          type: string
      outputs:
        - id: output_json
          type: File
          outputBinding:
            glob: '*_qc.json'
        - id: output_readDist
          type: File
          outputBinding:
            glob: '*readDist.txt'
        - id: output_tsv
          type: File
          outputBinding:
            glob: '*qc.tsv'
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            -i $(inputs.input_md_bam.path) -r $(inputs.gencode_bed.path) >
            $(inputs.output_basename).readDist.txt && /usr/bin/Rscript --vanilla
            /mend_qc/parseReadDist.R $(inputs.output_basename).readDist.txt

            mv bam_umend_qc.json $(inputs.output_basename).bam_umend_qc.json

            mv bam_umend_qc.tsv $(inputs.output_basename).bam_umend_qc.tsv
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 16000
          coresMin: 8
        - class: DockerRequirement
          dockerPull: kfdrc/mend_qc
        - class: InlineJavascriptRequirement
    'sbg:x': 459.0577697753906
    'sbg:y': 160.5
  - id: preprocess_bam
    in:
      - id: input_bam
        source: input_bam
      - id: output_basename
        source: output_basename
    out:
      - id: out_md_bam
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: mend_qc
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: input_bam
          type: File
          doc: RNA aligned to genome bam
        - id: output_basename
          type: string
      outputs:
        - id: out_md_bam
          type: File
          outputBinding:
            glob: '*.sortedByCoord.md.bam'
          secondaryFiles:
            - .bai
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            set -eo pipefail

            mkdir TMP

            sambamba sort -t 4 -m 4GB --tmpdir TMP --sort-by-name --out
            $(inputs.output_basename).sortedByName.bam $(inputs.input_bam.path)

            sambamba view -h $(inputs.output_basename).sortedByName.bam |
            samblaster | sambamba view --sam-input --format bam /dev/stdin >
            $(inputs.output_basename).sortedByName.md.bam

            sambamba sort --tmpdir TMP --show-progress -t 4 -m 4GB
            --out=$(inputs.output_basename).sortedByCoord.md.bam
            $(inputs.output_basename).sortedByName.md.bam
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: kfdrc/mend_qc
        - class: InlineJavascriptRequirement
    'sbg:x': 192.703125
    'sbg:y': 46.5
requirements: []
'sbg:appVersion':
  - v1.0
'sbg:id': kfdrc-harmonization/sd-bhjxbdqk/kfdrc-mendqc-wf/0
'sbg:revision': 0
'sbg:revisionNotes': Copy of kfdrc-harmonization/sd-bhjxbdqk-06/kfdrc-mendqc-wf/0
'sbg:modifiedOn': 1577550399
'sbg:modifiedBy': brownm28
'sbg:createdOn': 1577550399
'sbg:createdBy': brownm28
'sbg:project': kfdrc-harmonization/sd-bhjxbdqk
'sbg:projectName': CBTTC BGI RNAseq
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - brownm28
'sbg:latestRevision': 0
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': brownm28
    'sbg:modifiedOn': 1577550399
    'sbg:revisionNotes': Copy of kfdrc-harmonization/sd-bhjxbdqk-06/kfdrc-mendqc-wf/0
'sbg:image_url': >-
  https://cavatica.sbgenomics.com/ns/brood/images/kfdrc-harmonization/sd-bhjxbdqk/kfdrc-mendqc-wf/0.png
'sbg:publisher': sbg
'sbg:content_hash': a0c70016b5e78f1c99ce5005b8196b5eaa09595c37f50d503a80f3d1a9b88bb73
'sbg:copyOf': kfdrc-harmonization/sd-bhjxbdqk-06/kfdrc-mendqc-wf/0
