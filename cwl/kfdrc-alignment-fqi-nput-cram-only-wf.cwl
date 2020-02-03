class: Workflow
cwlVersion: v1.0
id: kfdrc-harmonization/sd-bhjxbdqk-05/kfdrc-alignment-fqinput-cram-only-wf/5
label: kfdrc-alignment-fqinput-cram-only-wf
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: files_R1
    type: 'File[]'
    'sbg:x': 0
    'sbg:y': 535
  - id: files_R2
    type: 'File[]'
    'sbg:x': 0
    'sbg:y': 428
  - id: indexed_reference_fasta
    type: File
    'sbg:x': 1892.3577880859375
    'sbg:y': 321
  - id: knownsites
    type: 'File[]'
    'sbg:x': 0
    'sbg:y': 321
  - id: output_basename
    type: string
    'sbg:x': 265.5
    'sbg:y': 246.5
  - id: reference_dict
    type: File
    'sbg:x': 0
    'sbg:y': 214
  - id: rgs
    type: 'string[]'
    'sbg:x': 0
    'sbg:y': 107
  - id: wgs_coverage_interval_list
    type: File
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: aggregation_metrics
    outputSource:
      - picard_collectaggregationmetrics/output
    type: 'File[]'
    'sbg:x': 2352.254638671875
    'sbg:y': 374.5
  - id: bqsr_report
    outputSource:
      - gatk_gatherbqsrreports/output
    type: File
    'sbg:x': 1585.91162109375
    'sbg:y': 321
  - id: cram
    outputSource:
      - samtools_coverttocram/output
    type: File
    'sbg:x': 2352.254638671875
    'sbg:y': 267.5
  - id: wgs_metrics
    outputSource:
      - picard_collectwgsmetrics/output
    type: File
    'sbg:x': 2352.254638671875
    'sbg:y': 160.5
steps:
  - id: bwa_mem
    in:
      - id: file_R1
        source: files_R1
      - id: file_R2
        source: files_R2
      - id: ref
        source: indexed_reference_fasta
      - id: rg
        source: rgs
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: bwa_mem_fq
      baseCommand: []
      inputs:
        - id: file_R1
          type: File
        - id: file_R2
          type: File
        - id: ref
          type: File
          secondaryFiles:
            - .64.amb
            - .64.ann
            - .64.bwt
            - .64.pac
            - .64.sa
            - .64.alt
            - ^.dict
            - .amb
            - .ann
            - .bwt
            - .pac
            - .sa
        - id: rg
          type: string
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: '*.bam'
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            bwa mem -K 100000000 -v 3 -t 36 -Y $(inputs.ref.path) -R
            '$(inputs.rg)' $(inputs.file_R1.path) $(inputs.file_R2.path) |
            /opt/samblaster/samblaster -i /dev/stdin -o /dev/stdout |
            /opt/sambamba_0.6.3/sambamba_v0.6.3 view -t 17 -f bam -l 0 -S
            /dev/stdin | /opt/sambamba_0.6.3/sambamba_v0.6.3 sort -t 17
            --natural-sort -m 15GiB --tmpdir ./ -o
            $(inputs.file_R1.nameroot).unsorted.bam -l 5 /dev/stdin
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 50000
          coresMin: 36
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/bogdang/bwa-kf-bundle:0.1.17'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: c5.9xlarge;ebs-gp2;750
    scatter:
      - file_R1
      - file_R2
      - rg
    scatterMethod: dotproduct
    'sbg:x': 265.5
    'sbg:y': 374.5
  - id: gatk_applybqsr
    in:
      - id: bqsr_report
        source: gatk_gatherbqsrreports/output
      - id: input_bam
        source: sambamba_sort/sorted_bam
      - id: reference
        source: indexed_reference_fasta
      - id: sequence_interval
        source: python_createsequencegroups/sequence_intervals_with_unmapped
    out:
      - id: recalibrated_bam
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_applybqsr
      baseCommand:
        - /gatk
        - ApplyBQSR
      inputs:
        - id: bqsr_report
          type: File
        - id: input_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
        - id: sequence_interval
          type: File
      outputs:
        - id: recalibrated_bam
          type: File
          outputBinding:
            glob: '*bam'
          secondaryFiles:
            - ^.bai
            - .md5
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --java-options "-Xms3000m -XX:+PrintFlagsFinal
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails
            -Xloggc:gc_log.log -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
            --create-output-bam-md5 --add-output-sam-program-record -R
            $(inputs.reference.path) -I $(inputs.input_bam.path)
            --use-original-qualities -O
            $(inputs.input_bam.nameroot).aligned.duplicates_marked.recalibrated.bam
            -bqsr $(inputs.bqsr_report.path) --static-quantized-quals 10
            --static-quantized-quals 20 --static-quantized-quals 30 -L
            $(inputs.sequence_interval.path)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 4500
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.3.0'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: c5.4xlarge;ebs-gp2;500
    scatter:
      - sequence_interval
    'sbg:x': 1585.91162109375
    'sbg:y': 193
  - id: gatk_baserecalibrator
    in:
      - id: input_bam
        source: sambamba_sort/sorted_bam
      - id: knownsites
        source:
          - knownsites
      - id: reference
        source: indexed_reference_fasta
      - id: sequence_interval
        source: python_createsequencegroups/sequence_intervals
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatkv4_baserecalibrator
      baseCommand:
        - /gatk
        - BaseRecalibrator
      inputs:
        - id: input_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: knownsites
          type:
            type: array
            items: File
            inputBinding:
              prefix: '--known-sites'
          inputBinding:
            position: 1
          secondaryFiles:
            - .tbi
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
        - id: sequence_interval
          type: File
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: '*.recal_data.csv'
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            --java-options "-Xms4000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10
            -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps
            -XX:+PrintGCDetails -Xloggc:gc_log.log" -R $(inputs.reference.path)
            -I $(inputs.input_bam.path) --use-original-qualities -O
            $(inputs.input_bam.nameroot).recal_data.csv -L
            $(inputs.sequence_interval.path)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.3.0'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: c5.4xlarge;ebs-gp2;500
    scatter:
      - sequence_interval
    'sbg:x': 1102.87646484375
    'sbg:y': 246.5
  - id: gatk_gatherbqsrreports
    in:
      - id: input_brsq_reports
        source:
          - gatk_baserecalibrator/output
      - id: output_basename
        source: output_basename
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_gatherbqsrreports
      baseCommand:
        - /gatk
        - GatherBQSRReports
      inputs:
        - id: input_brsq_reports
          type:
            type: array
            items: File
            inputBinding:
              prefix: '-I'
              separate: true
        - id: output_basename
          type: string
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: '*.csv'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --java-options "-Xms3000m" -O
            $(inputs.output_basename).GatherBqsrReports.recal_data.csv
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.3.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 1351.447509765625
    'sbg:y': 260.5
  - id: picard_collectaggregationmetrics
    in:
      - id: input_bam
        source: picard_gatherbamfiles/output
      - id: reference
        source: indexed_reference_fasta
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: picard_collectaggregationmetrics
      baseCommand:
        - java
        - '-Xms5000m'
        - '-jar'
        - /picard.jar
        - CollectMultipleMetrics
      inputs:
        - id: input_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
      outputs:
        - id: output
          type: 'File[]'
          outputBinding:
            glob: $(inputs.input_bam.nameroot).*
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            INPUT=$(inputs.input_bam.path)
            REFERENCE_SEQUENCE=$(inputs.reference.path)
            OUTPUT=$(inputs.input_bam.nameroot) ASSUME_SORTED=true
            PROGRAM="null" PROGRAM="CollectAlignmentSummaryMetrics" 
            PROGRAM="CollectInsertSizeMetrics" 
            PROGRAM="CollectSequencingArtifactMetrics" 
            PROGRAM="CollectGcBiasMetrics"  PROGRAM="QualityScoreDistribution" 
            METRIC_ACCUMULATION_LEVEL="null" 
            METRIC_ACCUMULATION_LEVEL="SAMPLE" 
            METRIC_ACCUMULATION_LEVEL="LIBRARY"
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 12000
        - class: DockerRequirement
          dockerPull: 'kfdrc/picard-r:latest-dev'
        - class: InlineJavascriptRequirement
    'sbg:x': 2152.21240234375
    'sbg:y': 388.5
  - id: picard_collectwgsmetrics
    in:
      - id: input_bam
        source: picard_gatherbamfiles/output
      - id: intervals
        source: wgs_coverage_interval_list
      - id: reference
        source: indexed_reference_fasta
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: picard_collectwgsmetrics
      baseCommand:
        - java
        - '-Xms2000m'
        - '-jar'
        - /picard.jar
        - CollectWgsMetrics
      inputs:
        - id: input_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: intervals
          type: File
        - id: reference
          type: File
          secondaryFiles:
            - .fai
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: $(inputs.input_bam.nameroot).wgs_metrics
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            INPUT=$(inputs.input_bam.path) VALIDATION_STRINGENCY=SILENT
            REFERENCE_SEQUENCE=$(inputs.reference.path)
            INCLUDE_BQ_HISTOGRAM=true INTERVALS=$(inputs.intervals.path)
            OUTPUT=$(inputs.input_bam.nameroot).wgs_metrics
            USE_FAST_ALGORITHM=true READ_LENGTH=250
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
        - class: DockerRequirement
          dockerPull: 'kfdrc/picard:2.18.2-dev'
        - class: InlineJavascriptRequirement
    'sbg:x': 2152.21240234375
    'sbg:y': 260.5
  - id: picard_gatherbamfiles
    in:
      - id: input_bam
        source:
          - gatk_applybqsr/recalibrated_bam
      - id: output_bam_basename
        source: output_basename
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: picard_gatherbamfiles
      baseCommand: []
      inputs:
        - id: input_bam
          type: 'File[]'
          secondaryFiles:
            - ^.bai
        - id: output_bam_basename
          type: string
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: '*.bam'
          secondaryFiles:
            - ^.bai
            - .md5
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            rm_bams="${
              var arr = [];
              for (var i=0; i<inputs.input_bam.length; i++)
                  arr = arr.concat(inputs.input_bam[i].path)
              return (arr.join(' '))
            }"

            input_bams="${
              var arr = [];
              for (var i=0; i<inputs.input_bam.length; i++)
                  arr = arr.concat(inputs.input_bam[i].path)
              return (arr.join(' INPUT='))
            }"

            java -Xms2000m -jar /picard.jar GatherBamFiles
            OUTPUT=$(inputs.output_bam_basename).bam INPUT=$input_bams
            CREATE_INDEX=true CREATE_MD5_FILE=true && rm $rm_bams
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
        - class: DockerRequirement
          dockerPull: 'kfdrc/picard:2.18.2-dev'
        - class: InlineJavascriptRequirement
    'sbg:x': 1892.3577880859375
    'sbg:y': 207
  - id: python_createsequencegroups
    in:
      - id: ref_dict
        source: reference_dict
    out:
      - id: sequence_intervals
      - id: sequence_intervals_with_unmapped
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: python_createsequencegroups
      baseCommand:
        - python
        - '-c'
      inputs:
        - id: ref_dict
          type: File
      outputs:
        - id: sequence_intervals
          type: 'File[]'
          outputBinding:
            glob: sequence_grouping_*.intervals
        - id: sequence_intervals_with_unmapped
          type: 'File[]'
          outputBinding:
            glob: '*.intervals'
      arguments:
        - position: 0
          valueFrom: |-
            def main():
                with open("$(inputs.ref_dict.path)", "r") as ref_dict_file:
                    sequence_tuple_list = []
                    longest_sequence = 0
                    for line in ref_dict_file:
                        if line.startswith("@SQ"):
                            line_split = line.split("\t")
                            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
                    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
                hg38_protection_tag = ":1+"
                tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
                temp_size = sequence_tuple_list[0][1]
                i = 0
                for sequence_tuple in sequence_tuple_list[1:]:
                    if temp_size + sequence_tuple[1] <= longest_sequence:
                        temp_size += sequence_tuple[1]
                        tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
                    else:
                        i += 1
                        pad = "{:0>2d}".format(i)
                        tsv_file_name = "sequence_grouping_" + pad + ".intervals"
                        with open(tsv_file_name, "w") as tsv_file:
                            tsv_file.write(tsv_string)
                            tsv_file.close()
                        tsv_string = sequence_tuple[0] + hg38_protection_tag
                        temp_size = sequence_tuple[1]
                i += 1
                pad = "{:0>2d}".format(i)
                tsv_file_name = "sequence_grouping_" + pad + ".intervals"
                with open(tsv_file_name, "w") as tsv_file:
                    tsv_file.write(tsv_string)
                    tsv_file.close()

                with open("unmapped.intervals", "w") as tsv_file:
                    tsv_file.write("unmapped")
                    tsv_file.close()

            if __name__ == "__main__":
                main()
      requirements:
        - class: DockerRequirement
          dockerPull: 'kfdrc/python:2.7.13'
        - class: InlineJavascriptRequirement
    'sbg:x': 265.5
    'sbg:y': 132.5
  - id: sambamba_merge
    in:
      - id: bams
        source:
          - bwa_mem/output
      - id: base_file_name
        source: output_basename
    out:
      - id: merged_bam
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: sambamba_merge
      baseCommand: []
      inputs:
        - id: bams
          type: 'File[]'
          inputBinding:
            position: 1
            shellQuote: false
        - id: base_file_name
          type: string
      outputs:
        - id: merged_bam
          type: File
          outputBinding:
            glob: '*.bam'
            outputEval: |-
              ${
                  if(inputs.bams.length > 1) return self
                  else return inputs.bams
              }
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: |-
            ${
                if (inputs.bams.length != 1)
                 return "/opt/sambamba_0.6.3/sambamba_v0.6.3 merge -t 36 " + inputs.base_file_name + ".aligned.duplicates_marked.unsorted.bam"
                else return "echo "
             }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 2024
          coresMin: 36
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/bogdang/sambamba:0.6.3'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: c5.9xlarge;ebs-gp2;1500
    'sbg:x': 598.1827392578125
    'sbg:y': 260.5
  - id: sambamba_sort
    in:
      - id: bam
        source: sambamba_merge/merged_bam
      - id: base_file_name
        source: output_basename
    out:
      - id: sorted_bam
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: sambamba_sort
      baseCommand: []
      inputs:
        - id: bam
          type: File
        - id: base_file_name
          type: string
        - default: aligned.duplicates_marked.sorted
          id: suffix
          type: string
      outputs:
        - id: sorted_bam
          type: File
          outputBinding:
            glob: '*.bam'
          secondaryFiles:
            - ^.bai
          format: BAM
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /opt/sambamba_0.6.3/sambamba_v0.6.3 sort -t 36 -m 10G -o
            $(inputs.base_file_name).$(inputs.suffix).bam $(inputs.bam.path)

            mv $(inputs.base_file_name).$(inputs.suffix).bam.bai
            $(inputs.base_file_name).$(inputs.suffix).bai

            rm $(inputs.bam.path)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 15000
          coresMin: 36
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/bogdang/sambamba:0.6.3'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: c5.9xlarge;ebs-gp2;1500
    'sbg:x': 854.099853515625
    'sbg:y': 260.5
  - id: samtools_coverttocram
    in:
      - id: input_bam
        source: picard_gatherbamfiles/output
      - id: reference
        source: indexed_reference_fasta
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: samtools_convert_to_cram
      baseCommand:
        - samtools
        - view
      inputs:
        - id: input_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: reference
          type: File
          secondaryFiles:
            - .fai
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: '*.cram'
          secondaryFiles:
            - .crai
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            -C -T $(inputs.reference.path) -o $(inputs.input_bam.nameroot).cram
            $(inputs.input_bam.path) && samtools index
            $(inputs.input_bam.nameroot).cram
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 4000
        - class: DockerRequirement
          dockerPull: 'kfdrc/samtools:1.8-dev'
        - class: InlineJavascriptRequirement
    'sbg:x': 2152.21240234375
    'sbg:y': 132.5
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
requirements:
  - class: ScatterFeatureRequirement
'sbg:projectName': GH-CBTTC-FQ-INPUT-Tumor
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': jinhaoxuan
    'sbg:modifiedOn': 1529692474
    'sbg:revisionNotes': >-
      Copy of
      kfdrc-harmonization/sd-bhjxbdqk-04/kfdrc-alignment-fqinput-cram-only-wf/1
  - 'sbg:revision': 1
    'sbg:modifiedBy': jinhaoxuan
    'sbg:modifiedOn': 1530131796
    'sbg:revisionNotes': c48 1.2T
  - 'sbg:revision': 2
    'sbg:modifiedBy': jinhaoxuan
    'sbg:modifiedOn': 1530131827
    'sbg:revisionNotes': r4.8 1.2T
  - 'sbg:revision': 3
    'sbg:modifiedBy': jinhaoxuan
    'sbg:modifiedOn': 1530195946
    'sbg:revisionNotes': r4.8 1.2T cores32
  - 'sbg:revision': 4
    'sbg:modifiedBy': zhangb1
    'sbg:modifiedOn': 1550160711
    'sbg:revisionNotes': add unmapped and change instance
  - 'sbg:revision': 5
    'sbg:modifiedBy': zhangb1
    'sbg:modifiedOn': 1550504911
    'sbg:revisionNotes': add 1500G to sambamba
'sbg:image_url': >-
  https://cavatica.sbgenomics.com/ns/brood/images/kfdrc-harmonization/sd-bhjxbdqk-05/kfdrc-alignment-fqinput-cram-only-wf/5.png
'sbg:appVersion':
  - v1.0
'sbg:id': kfdrc-harmonization/sd-bhjxbdqk-05/kfdrc-alignment-fqinput-cram-only-wf/5
'sbg:revision': 5
'sbg:revisionNotes': add 1500G to sambamba
'sbg:modifiedOn': 1550504911
'sbg:modifiedBy': zhangb1
'sbg:createdOn': 1529692474
'sbg:createdBy': jinhaoxuan
'sbg:project': kfdrc-harmonization/sd-bhjxbdqk-05
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - jinhaoxuan
  - zhangb1
'sbg:latestRevision': 5
'sbg:publisher': sbg
'sbg:content_hash': a1fe60b20a1fbb5a002d5fb188abab830e39d02ff574de8e428b070737a54a63f
