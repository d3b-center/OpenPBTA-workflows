class: Workflow
cwlVersion: v1.0
"$namespaces":
  sbg: https://sevenbridges.com
inputs:
- id: biospecimen_name
  type: string
  sbg:x: 265.5
  sbg:y: 328
- id: indexed_reference_fasta
  type: File
  sbg:x: 0
  sbg:y: 428
- id: input_reads
  type: File
  sbg:x: 0
  sbg:y: 321
- id: knownsites
  type: File[]
  sbg:x: 0
  sbg:y: 214
- id: output_basename
  type: string
  sbg:x: 598.1827392578125
  sbg:y: 146.5
- id: reference_dict
  type: File
  sbg:x: 0
  sbg:y: 107
- id: wgs_coverage_interval_list
  type: File
  sbg:x: 0
  sbg:y: 0
outputs:
- id: aggregation_metrics
  outputSource:
  - picard_collectaggregationmetrics/output
  type: File[]
  sbg:x: 2666.65625
  sbg:y: 321
- id: bqsr_report
  outputSource:
  - gatk_gatherbqsrreports/output
  type: File
  sbg:x: 1900.3133544921875
  sbg:y: 267.5
- id: cram
  outputSource:
  - samtools_coverttocram/output
  type: File
  sbg:x: 2666.65625
  sbg:y: 214
- id: wgs_metrics
  outputSource:
  - picard_collectwgsmetrics/output
  type: File
  sbg:x: 2666.65625
  sbg:y: 107
steps:
- id: bwa_mem
  in:
  - id: indexed_reference_fasta
    source: indexed_reference_fasta
  - id: input_reads
    source: samtools_split/bam_files
  - id: sample_name
    source: biospecimen_name
  out:
  - id: aligned_bams
  run:
    class: Workflow
    cwlVersion: v1.0
    id: bwa_mem_wf
    inputs:
    - id: indexed_reference_fasta
      type: File
    - id: input_reads
      type: File
    - id: sample_name
      type: string
    outputs:
    - id: aligned_bams
      outputSource:
      - bwa_mem_split/output
      type: File[]
    steps:
    - id: bwa_input_prepare
      in:
      - id: input_bam
        source: input_reads
      out:
      - id: output
      - id: rg
      run:
        class: CommandLineTool
        cwlVersion: v1.0
        id: bwa_input_prepare
        baseCommand: []
        inputs:
        - id: input_bam
          type: File
        - default: 20000000000
          id: max_siz
          type: int
        outputs:
        - id: output
          type: File[]
          outputBinding:
            glob: "*.fq"
            outputEval: |-
              ${
                if( inputs.input_bam.size < inputs.max_siz ) return [inputs.input_bam]
                else return self
              }
        - id: rg
          type: File
          outputBinding:
            glob: rg.txt
        arguments:
        - position: 0
          shellQuote: false
          valueFrom: |-
            samtools view -H $(inputs.input_bam.path) | grep ^@RG > rg.txt
            if [ $(inputs.input_bam.size) -gt $(inputs.max_siz) ]; then
              bamtofastq tryoq=1 filename=$(inputs.input_bam.path) | split -dl 680000000 - reads-
              ls reads-* | xargs -i mv {} {}.fq
              rm $(inputs.input_bam.path)
            fi
        requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 1000
        - class: DockerRequirement
          dockerPull: kfdrc/bwa-bundle:dev
        - class: InlineJavascriptRequirement
    - id: bwa_mem_split
      in:
      - id: reads
        source: bwa_input_prepare/output
      - id: ref
        source: indexed_reference_fasta
      - id: rg
        source: expression_updatergsample/rg_str
      out:
      - id: output
      run:
        class: CommandLineTool
        cwlVersion: v1.0
        id: bwa_mem_split
        baseCommand:
        - "/bin/bash"
        - "-c"
        inputs:
        - id: reads
          type: File
        - id: ref
          type: File
          secondaryFiles:
          - ".64.amb"
          - ".64.ann"
          - ".64.bwt"
          - ".64.pac"
          - ".64.sa"
          - ".64.alt"
          - "^.dict"
          - ".amb"
          - ".ann"
          - ".bwt"
          - ".pac"
          - ".sa"
        - id: rg
          type: string
        outputs:
        - id: output
          type: File
          outputBinding:
            glob: "*.bam"
        arguments:
        - position: 0
          valueFrom: |-
            set -eo pipefail
            if [ $(inputs.reads.nameext) = ".bam" ]; then
              CMD="/opt/biobambam2/2.0.87-release-20180301132713/x86_64-etch-linux-gnu/bin/bamtofastq tryoq=1 filename=$(inputs.reads.path)"
            else
              CMD="cat $(inputs.reads.path)"
            fi
            $CMD | bwa mem -K 100000000 -p -v 3 -t 36 -Y $(inputs.ref.path) -R "$(inputs.rg)" - | /opt/samblaster/samblaster -i /dev/stdin -o /dev/stdout | /opt/sambamba_0.6.3/sambamba_v0.6.3 view -t 36 -f bam -l 0 -S /dev/stdin | /opt/sambamba_0.6.3/sambamba_v0.6.3 sort -t 36 --natural-sort -m 15GiB --tmpdir ./ -o $(inputs.reads.nameroot).unsorted.bam -l 5 /dev/stdin
            rm $(inputs.reads.path)
        requirements:
        - class: ResourceRequirement
          ramMin: 50000
          coresMin: 36
        - class: DockerRequirement
          dockerPull: images.sbgenomics.com/bogdang/bwa-kf-bundle:0.1.17
        - class: InlineJavascriptRequirement
        hints:
        - class: sbg:AWSInstanceType
          value: c5.9xlarge;ebs-gp2;750
      scatter:
      - reads
    - id: expression_updatergsample
      in:
      - id: rg
        source: bwa_input_prepare/rg
      - id: sample
        source: sample_name
      out:
      - id: rg_str
      run:
        class: ExpressionTool
        cwlVersion: v1.0
        expression: "${ var arr = inputs.rg.contents.split('\\n')[0].split('\\t');
          for (var i=1; i<arr.length; i++){ if (arr[i].startsWith('SM')){ arr[i] =
          'SM:' + inputs.sample; } } return {rg_str: arr.join('\\\\t')}; }"
        id: expression_preparerg
        inputs:
          rg:
            inputBinding:
              loadContents: true
            type: File
          sample: string
        outputs:
          rg_str: string
        requirements:
        - class: InlineJavascriptRequirement
    requirements:
    - class: ScatterFeatureRequirement
  scatter:
  - input_reads
  sbg:x: 598.1827392578125
  sbg:y: 267.5
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
    - "/gatk"
    - ApplyBQSR
    inputs:
    - id: bqsr_report
      type: File
    - id: input_bam
      type: File
      secondaryFiles:
      - "^.bai"
    - id: reference
      type: File
      secondaryFiles:
      - "^.dict"
      - ".fai"
    - id: sequence_interval
      type: File
    outputs:
    - id: recalibrated_bam
      type: File
      outputBinding:
        glob: "*bam"
      secondaryFiles:
      - "^.bai"
      - ".md5"
    arguments:
    - position: 1
      shellQuote: false
      valueFrom: --java-options "-Xms3000m -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps
        -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -XX:GCTimeLimit=50
        -XX:GCHeapFreeLimit=10" --create-output-bam-md5 --add-output-sam-program-record
        -R $(inputs.reference.path) -I $(inputs.input_bam.path) --use-original-qualities
        -O $(inputs.input_bam.nameroot).aligned.duplicates_marked.recalibrated.bam
        -bqsr $(inputs.bqsr_report.path) --static-quantized-quals 10 --static-quantized-quals
        20 --static-quantized-quals 30 -L $(inputs.sequence_interval.path)
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 4500
    - class: DockerRequirement
      dockerPull: kfdrc/gatk:4.0.3.0
    - class: InlineJavascriptRequirement
    hints:
    - class: sbg:AWSInstanceType
      value: c5.4xlarge;ebs-gp2;500
  scatter:
  - sequence_interval
  sbg:x: 1900.3133544921875
  sbg:y: 139.5
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
    - "/gatk"
    - BaseRecalibrator
    inputs:
    - id: input_bam
      type: File
      secondaryFiles:
      - "^.bai"
    - id: knownsites
      type:
        type: array
        items: File
        inputBinding:
          prefix: "--known-sites"
      inputBinding:
        position: 1
      secondaryFiles:
      - ".tbi"
    - id: reference
      type: File
      secondaryFiles:
      - "^.dict"
      - ".fai"
    - id: sequence_interval
      type: File
    outputs:
    - id: output
      type: File
      outputBinding:
        glob: "*.recal_data.csv"
    arguments:
    - position: 0
      shellQuote: false
      valueFrom: --java-options "-Xms4000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10
        -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails
        -Xloggc:gc_log.log" -R $(inputs.reference.path) -I $(inputs.input_bam.path)
        --use-original-qualities -O $(inputs.input_bam.nameroot).recal_data.csv -L
        $(inputs.sequence_interval.path)
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: kfdrc/gatk:4.0.3.0
    - class: InlineJavascriptRequirement
    hints:
    - class: sbg:AWSInstanceType
      value: c5.4xlarge;ebs-gp2;500
  scatter:
  - sequence_interval
  sbg:x: 1417.278076171875
  sbg:y: 193
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
    - "/gatk"
    - GatherBQSRReports
    inputs:
    - id: input_brsq_reports
      type:
        type: array
        items: File
        inputBinding:
          prefix: "-I"
          separate: true
    - id: output_basename
      type: string
    outputs:
    - id: output
      type: File
      outputBinding:
        glob: "*.csv"
    arguments:
    - position: 1
      shellQuote: false
      valueFrom: --java-options "-Xms3000m" -O $(inputs.output_basename).GatherBqsrReports.recal_data.csv
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: kfdrc/gatk:4.0.3.0
    - class: InlineJavascriptRequirement
  sbg:x: 1665.8492431640625
  sbg:y: 207
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
    - "-Xms5000m"
    - "-jar"
    - "/picard.jar"
    - CollectMultipleMetrics
    inputs:
    - id: input_bam
      type: File
      secondaryFiles:
      - "^.bai"
    - id: reference
      type: File
      secondaryFiles:
      - "^.dict"
      - ".fai"
    outputs:
    - id: output
      type: File[]
      outputBinding:
        glob: "$(inputs.input_bam.nameroot).*"
    arguments:
    - position: 1
      shellQuote: false
      valueFrom: INPUT=$(inputs.input_bam.path) REFERENCE_SEQUENCE=$(inputs.reference.path)
        OUTPUT=$(inputs.input_bam.nameroot) ASSUME_SORTED=true PROGRAM="null" PROGRAM="CollectAlignmentSummaryMetrics"  PROGRAM="CollectInsertSizeMetrics"  PROGRAM="CollectSequencingArtifactMetrics"  PROGRAM="CollectGcBiasMetrics"  PROGRAM="QualityScoreDistribution"  METRIC_ACCUMULATION_LEVEL="null"  METRIC_ACCUMULATION_LEVEL="SAMPLE"  METRIC_ACCUMULATION_LEVEL="LIBRARY"
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 12000
    - class: DockerRequirement
      dockerPull: kfdrc/picard-r:latest-dev
    - class: InlineJavascriptRequirement
  sbg:x: 2466.614013671875
  sbg:y: 335
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
    - "-Xms2000m"
    - "-jar"
    - "/picard.jar"
    - CollectWgsMetrics
    inputs:
    - id: input_bam
      type: File
      secondaryFiles:
      - "^.bai"
    - id: intervals
      type: File
    - id: reference
      type: File
      secondaryFiles:
      - ".fai"
    outputs:
    - id: output
      type: File
      outputBinding:
        glob: "$(inputs.input_bam.nameroot).wgs_metrics"
    arguments:
    - position: 1
      shellQuote: false
      valueFrom: INPUT=$(inputs.input_bam.path) VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$(inputs.reference.path)
        INCLUDE_BQ_HISTOGRAM=true INTERVALS=$(inputs.intervals.path) OUTPUT=$(inputs.input_bam.nameroot).wgs_metrics
        USE_FAST_ALGORITHM=true READ_LENGTH=250
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: kfdrc/picard:2.18.2-dev
    - class: InlineJavascriptRequirement
  sbg:x: 2466.614013671875
  sbg:y: 207
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
      type: File[]
      secondaryFiles:
      - "^.bai"
    - id: output_bam_basename
      type: string
    outputs:
    - id: output
      type: File
      outputBinding:
        glob: "*.bam"
      secondaryFiles:
      - "^.bai"
      - ".md5"
    arguments:
    - position: 1
      shellQuote: false
      valueFrom: |-
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
        java -Xms2000m -jar /picard.jar GatherBamFiles OUTPUT=$(inputs.output_bam_basename).bam INPUT=$input_bams CREATE_INDEX=true CREATE_MD5_FILE=true && rm $rm_bams
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: kfdrc/picard:2.18.2-dev
    - class: InlineJavascriptRequirement
  sbg:x: 2206.759521484375
  sbg:y: 207
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
    - "-c"
    inputs:
    - id: ref_dict
      type: File
    outputs:
    - id: sequence_intervals
      type: File[]
      outputBinding:
        glob: sequence_grouping_*.intervals
    - id: sequence_intervals_with_unmapped
      type: File[]
      outputBinding:
        glob: "*.intervals"
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
      dockerPull: kfdrc/python:2.7.13
    - class: InlineJavascriptRequirement
  sbg:x: 265.5
  sbg:y: 214
- id: sambamba_merge
  in:
  - id: bams
    source:
    - bwa_mem/aligned_bams
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
      type:
        type: array
        items:
          items: File
          type: array
    - id: base_file_name
      type: string
    outputs:
    - id: merged_bam
      type: File
      outputBinding:
        glob: "*.bam"
      secondaryFiles:
      - ".bai"
      - "^.bai"
      format: BAM
    arguments:
    - position: 0
      shellQuote: false
      valueFrom: |-
        bams="${
          var arr = [];
          for (var i=0; i<inputs.bams.length; i++)
            for (var j=0; j<inputs.bams[i].length; j++)
              arr = arr.concat(inputs.bams[i][j].path)
          return (arr.join(' '))
        }"
        /opt/sambamba_0.6.3/sambamba_v0.6.3 merge -t 36 $(inputs.base_file_name).aligned.duplicates_marked.unsorted.bam $bams && rm $bams
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 2024
      coresMin: 36
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/bogdang/sambamba:0.6.3
    - class: InlineJavascriptRequirement
  sbg:x: 912.58447265625
  sbg:y: 207
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
        glob: "*.bam"
      secondaryFiles:
      - "^.bai"
      format: BAM
    arguments:
    - position: 0
      shellQuote: false
      valueFrom: |-
        /opt/sambamba_0.6.3/sambamba_v0.6.3 sort -t 36 -m 10G -o $(inputs.base_file_name).$(inputs.suffix).bam $(inputs.bam.path)
        mv $(inputs.base_file_name).$(inputs.suffix).bam.bai $(inputs.base_file_name).$(inputs.suffix).bai
        rm $(inputs.bam.path)
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 15000
      coresMin: 36
    - class: DockerRequirement
      dockerPull: images.sbgenomics.com/bogdang/sambamba:0.6.3
    - class: InlineJavascriptRequirement
  sbg:x: 1168.5015869140625
  sbg:y: 207
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
      - "^.bai"
    - id: reference
      type: File
      secondaryFiles:
      - ".fai"
    outputs:
    - id: output
      type: File
      outputBinding:
        glob: "*.cram"
      secondaryFiles:
      - ".crai"
    arguments:
    - position: 1
      shellQuote: false
      valueFrom: "-C -T $(inputs.reference.path) -o $(inputs.input_bam.nameroot).cram
        $(inputs.input_bam.path) && samtools index $(inputs.input_bam.nameroot).cram"
    requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: kfdrc/samtools:1.8-dev
    - class: InlineJavascriptRequirement
  sbg:x: 2466.614013671875
  sbg:y: 79
- id: samtools_split
  in:
  - id: input_bam
    source: input_reads
  - id: reference
    source: indexed_reference_fasta
  out:
  - id: bam_files
  run:
    class: CommandLineTool
    cwlVersion: v1.0
    id: samtools_split
    baseCommand: []
    inputs:
    - id: input_bam
      type: File
    - id: reference
      type: File
    outputs:
    - id: bam_files
      type: File[]
      outputBinding:
        glob: "*.bam"
        outputEval: |-
          ${
            if (self.length == 0) return [inputs.input_bam]
            else return self
          }
    arguments:
    - position: 0
      shellQuote: false
      valueFrom: |-
        RG_NUM=`samtools view -H $(inputs.input_bam.path) | grep -c ^@RG`
        if [ $RG_NUM != 1 ]; then
          samtools split -f '%!.bam' -@ 36 --reference $(inputs.reference.path) $(inputs.input_bam.path)
          rm $(inputs.input_bam.path)
        fi
    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: kfdrc/samtools:1.8-dev
    - class: InlineJavascriptRequirement
  sbg:x: 265.5
  sbg:y: 93
hints:
- class: sbg:maxNumberOfParallelInstances
  value: 4
requirements:
- class: SubworkflowFeatureRequirement
- class: ScatterFeatureRequirement
label: kf-alignment-cram-only-wf
sbg:appVersion:
- v1.0
id: https://cavatica-api.sbgenomics.com/v2/apps/kfdrc-harmonization/sd-bhjxbdqk-02/kf-alignment-cram-only-wf/6/raw/
sbg:id: kfdrc-harmonization/sd-bhjxbdqk-02/kf-alignment-cram-only-wf/6
sbg:revision: 6
sbg:revisionNotes: add instance
sbg:modifiedOn: 1550074882
sbg:modifiedBy: zhangb1
sbg:createdOn: 1528749498
sbg:createdBy: zhangb1
sbg:project: kfdrc-harmonization/sd-bhjxbdqk-02
sbg:projectName: GH-CBTTC-Tumor
sbg:sbgMaintained: false
sbg:validationErrors: []
sbg:contributors:
- zhangb1
sbg:latestRevision: 6
sbg:revisionsInfo:
- sbg:revision: 0
  sbg:modifiedBy: zhangb1
  sbg:modifiedOn: 1528749498
  sbg:revisionNotes: 
- sbg:revision: 1
  sbg:modifiedBy: zhangb1
  sbg:modifiedOn: 1528749813
  sbg:revisionNotes: 
- sbg:revision: 2
  sbg:modifiedBy: zhangb1
  sbg:modifiedOn: 1550073531
  sbg:revisionNotes: add_unmapped_and_instance change
- sbg:revision: 3
  sbg:modifiedBy: zhangb1
  sbg:modifiedOn: 1550074375
  sbg:revisionNotes: bwa mem add instacne
- sbg:revision: 4
  sbg:modifiedBy: zhangb1
  sbg:modifiedOn: 1550074604
  sbg:revisionNotes: add bwa mem split instance
- sbg:revision: 5
  sbg:modifiedBy: zhangb1
  sbg:modifiedOn: 1550074705
  sbg:revisionNotes: bwa mem split add instance
- sbg:revision: 6
  sbg:modifiedBy: zhangb1
  sbg:modifiedOn: 1550074882
  sbg:revisionNotes: add instance
sbg:image_url: https://cavatica.sbgenomics.com/ns/brood/images/kfdrc-harmonization/sd-bhjxbdqk-02/kf-alignment-cram-only-wf/6.png
sbg:publisher: sbg
sbg:content_hash: a1c6e785416e640d6fc5febc5a0eafd8a042b800fe5ef7f9f2f23564b4b4bb22c
