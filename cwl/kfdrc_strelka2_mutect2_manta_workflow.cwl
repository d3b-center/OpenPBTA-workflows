class: Workflow
cwlVersion: v1.0
id: kfdrc_somatic_wf
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: af_only_gnomad_vcf
    type: File
    secondaryFiles:
      - .tbi
    'sbg:x': 247.3125
    'sbg:y': 736.6875
  - id: exac_common_vcf
    type: File
    secondaryFiles:
      - .tbi
    'sbg:x': 567.62158203125
    'sbg:y': 650.8125
  - id: exome_flag
    type: string?
    doc: insert 'Y' if exome mode
    'sbg:x': 0
    'sbg:y': 855
  - id: hg38_strelka_bed
    type: File
    secondaryFiles:
      - .tbi
    'sbg:x': 0
    'sbg:y': 748.125
  - id: indexed_reference_fasta
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    'sbg:x': 0
    'sbg:y': 641.25
  - id: input_normal_aligned
    type: File
    doc: normal BAM or CRAM
    secondaryFiles:
      - |
        ${
          var dpath = self.location.replace(self.basename, "")
          if(self.nameext == '.bam'){
            return {"location": dpath+self.nameroot+".bai", "class": "File"}
          }
          else{
            return {"location": dpath+self.basename+".crai", "class": "File"}
          }
        }
    'sbg:x': 0
    'sbg:y': 534.375
  - id: input_normal_name
    type: string
    'sbg:x': 247.3125
    'sbg:y': 494.9375
  - id: input_tumor_aligned
    type: File
    doc: tumor BAM or CRAM
    secondaryFiles:
      - |
        ${
          var dpath = self.location.replace(self.basename, "")
          if(self.nameext == '.bam'){
            return {"location": dpath+self.nameroot+".bai", "class": "File"}
          }
          else{
            return {"location": dpath+self.basename+".crai", "class": "File"}
          }
        }
    'sbg:x': 0
    'sbg:y': 427.5
  - id: input_tumor_name
    type: string
    'sbg:x': 247.3125
    'sbg:y': 388.0625
  - id: output_basename
    type: string
    'sbg:x': 1334.072265625
    'sbg:y': 264.625
  - id: reference_dict
    type: File
    'sbg:x': 0
    'sbg:y': 320.625
  - id: select_vars_mode
    type:
      - 'null'
      - type: enum
        symbols:
          - gatk
          - grep
        name: select_vars_mode
    doc: 'Choose ''gatk'' for SelectVariants tool, or ''grep'' for grep expression'
    default: gatk
    'sbg:x': 0
    'sbg:y': 213.75
  - id: vep_cache
    type: File
    doc: tar gzipped cache from ensembl/local converted cache
    'sbg:x': 0
    'sbg:y': 106.875
  - id: wgs_calling_interval_list
    type: File
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: manta_pass_vcf
    outputSource:
      - gatk_selectvariants_manta/pass_vcf
    type: File
    'sbg:x': 1334.072265625
    'sbg:y': 371.5
  - id: manta_prepass_vcf
    outputSource:
      - rename_manta_samples/reheadered_vcf
    type: File
    'sbg:x': 931.1585693359375
    'sbg:y': 643.8125
  - id: mutect2_prepass_vcf
    outputSource:
      - filter_mutect2_vcf/filtered_vcf
    type: File
    'sbg:x': 1652.5955810546875
    'sbg:y': 413.5
  - id: mutect2_vep_maf
    outputSource:
      - vep_annot_mutect2/output_maf
    type: File
    'sbg:x': 2599.4423828125
    'sbg:y': 534.375
  - id: mutect2_vep_tbi
    outputSource:
      - vep_annot_mutect2/output_tbi
    type: File
    'sbg:x': 2599.4423828125
    'sbg:y': 427.5
  - id: mutect2_vep_vcf
    outputSource:
      - vep_annot_mutect2/output_vcf
    type: File
    'sbg:x': 2599.4423828125
    'sbg:y': 320.625
  - id: strelka2_prepass_vcf
    outputSource:
      - rename_strelka_samples/reheadered_vcf
    type: File
    'sbg:x': 1334.072265625
    'sbg:y': 157.75
  - id: strelka2_vep_maf
    outputSource:
      - vep_annot_strelka2/output_maf
    type: File
    'sbg:x': 2126.01904296875
    'sbg:y': 587.8125
  - id: strelka2_vep_tbi
    outputSource:
      - vep_annot_strelka2/output_tbi
    type: File
    'sbg:x': 2126.01904296875
    'sbg:y': 480.9375
  - id: strelka2_vep_vcf
    outputSource:
      - vep_annot_strelka2/output_vcf
    type: File
    'sbg:x': 2126.01904296875
    'sbg:y': 374.0625
steps:
  - id: filter_mutect2_vcf
    in:
      - id: contamination_table
        source: mutect2_filter_support/contamination_table
      - id: mutect_stats
        source: merge_mutect2_stats/merged_stats
      - id: mutect_vcf
        source: merge_mutect2_vcf/merged_vcf
      - id: ob_priors
        source: mutect2_filter_support/f1r2_bias
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
      - id: segmentation_table
        source: mutect2_filter_support/segmentation_table
    out:
      - id: filtered_vcf
      - id: stats_table
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_filtermutect2calls
      baseCommand:
        - /gatk
        - FilterMutectCalls
      inputs:
        - id: contamination_table
          type: File
        - id: mutect_stats
          type: File
        - id: mutect_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: ob_priors
          type: File
        - id: output_basename
          type: string
        - id: reference
          type: File
        - id: segmentation_table
          type: File
      outputs:
        - id: filtered_vcf
          type: File
          outputBinding:
            glob: '*.mutect2_filtered.vcf.gz'
          secondaryFiles:
            - .tbi
        - id: stats_table
          type: File
          outputBinding:
            glob: '*.mutect2_filtered.txt'
      label: GATK Filter Mutect2
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            --java-options "-Xmx4000m" -V $(inputs.mutect_vcf.path) -O
            $(inputs.output_basename).mutect2_filtered.vcf.gz -R
            $(inputs.reference.path) --contamination-table
            $(inputs.contamination_table.path) --tumor-segmentation
            $(inputs.segmentation_table.path) --ob-priors
            $(inputs.ob_priors.path) --filtering-stats
            $(inputs.output_basename).mutect2_filtered.txt --stats
            $(inputs.mutect_stats.path)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 4000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: GATK Filter Mutect2
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 1334.072265625
    'sbg:y': 655.25
  - id: gatk_intervallisttools
    in:
      - id: bands
        valueFrom: '${return 80000000}'
      - id: exome_flag
        source: exome_flag
      - id: interval_list
        source: wgs_calling_interval_list
      - id: reference_dict
        source: reference_dict
      - id: scatter_ct
        valueFrom: '${return 50}'
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_intervallist2bed
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: bands
          type: int
        - default: 'N'
          id: break_by_chr
          type: string?
          doc: >-
            If Y, break up files by chr.  If creating smaler intervals,
            recommend scatter_ct=1
        - default: 'N'
          id: exome_flag
          type: string?
          doc: 'If ''Y'', will set bands to 0 to prevent breaking up of intervals'
        - id: interval_list
          type: File
        - 'sbg:suggestedValue':
            name: >-
              Provide only if input is bed file instead of gatk style
              .interval_list
          id: reference_dict
          type: File?
          doc: >-
            Provide only if input is bed file instead of gatk style
            .interval_list
        - id: scatter_ct
          type: int
      outputs:
        - id: output
          type: 'File[]'
          outputBinding:
            glob: '*.bed'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: |-
            set -eo pipefail
            ${
              var cmd = "";
              if (inputs.interval_list.nameext == '.interval_list'){
                cmd = "LIST=" + inputs.interval_list.path + ";";
              }
              else{
                cmd = "/gatk BedToIntervalList -I " + inputs.interval_list.path + " -O " + inputs.interval_list.nameroot 
                + ".interval_list -SD " + inputs.reference_dict.path + "; LIST=" + inputs.interval_list.nameroot 
                + ".interval_list;";

              }
              if (inputs.exome_flag == "Y"){
                  cmd += "BANDS=0;";
                }
                
              else{
                cmd += "BANDS=" + inputs.bands + ";";
              }
              return cmd;
            }
            ${
              if (inputs.break_by_chr == "N"){
                var cmd = "/gatk IntervalListTools --java-options \"-Xmx2000m\" --SCATTER_COUNT=" + inputs.scatter_ct + " --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --UNIQUE=true --SORT=true --BREAK_BANDS_AT_MULTIPLES_OF=$BANDS --INPUT=$LIST --OUTPUT=.;"
                cmd += "CT=`find . -name 'temp_0*' | wc -l`;";
                cmd += "seq -f \"%04g\" $CT | xargs -I N -P 4 /gatk IntervalListToBed --java-options -Xmx100m -I temp_N_of_$CT/scattered.interval_list -O temp_N_of_$CT/scattered.interval_list.N.bed;";
                cmd += "mv temp_0*/*.bed .;";
              }
              else{
                cmd = "mkdir intvl_by_chr;"
                cmd += "/gatk IntervalListTools --java-options \"-Xmx2000m\" --SCATTER_COUNT=" + inputs.scatter_ct + " --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW --UNIQUE=true --SORT=true --BREAK_BANDS_AT_MULTIPLES_OF=$BANDS --INPUT=$LIST --OUTPUT=intvl_by_chr/scattered.interval_list;";
                cmd += "/gatk IntervalListToBed --java-options -Xmx100m -I intvl_by_chr/scattered.interval_list -O intvl_by_chr/scattered.interval_list.bed;"
                cmd += "cut -f 1 intvl_by_chr/scattered.interval_list.bed | uniq | xargs -ICM sh -c 'grep -P \"CM\\t\" intvl_by_chr/scattered.interval_list.bed > CM_intervals.bed';";
              }
              return cmd;
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 2000
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 247.3125
    'sbg:y': 615.8125
  - id: gatk_selectvariants_manta
    in:
      - id: input_vcf
        source: rename_manta_samples/reheadered_vcf
      - id: mode
        source: select_vars_mode
      - id: output_basename
        source: output_basename
      - id: tool_name
        valueFrom: '${return "manta"}'
    out:
      - id: pass_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_selectvariants
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - default: gatk
          id: mode
          type:
            - 'null'
            - type: enum
              symbols:
                - gatk
                - grep
              name: select_vars_mode
          doc: 'Choose ''gatk'' for SelectVariants tool, or ''grep'' for grep expression'
        - id: output_basename
          type: string
        - id: tool_name
          type: string
      outputs:
        - id: pass_vcf
          type: File
          outputBinding:
            glob: '*.PASS.vcf.gz'
          secondaryFiles:
            - .tbi
      label: GATK Select PASS
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: |-
            set -eo pipefail
            ${
              var run_mode = inputs.mode;
              if (run_mode == 'grep' || run_mode == 'gatk'){
                var in_vcf = inputs.input_vcf.path;
                var out_vcf = inputs.output_basename + '.' + inputs.tool_name + '.PASS.vcf.gz';
                var cmd = '/gatk SelectVariants --java-options "-Xmx8000m" -V ' + in_vcf +  ' -O ' + out_vcf + ' --exclude-filtered TRUE';
                if (run_mode == 'grep'){
                  cmd = 'zcat ' + in_vcf + ' | grep -E "^#|PASS" | bgzip > ' + out_vcf + '; tabix ' + out_vcf;
                }
                return cmd;
              }
              else{
                throw new Error(run_mode + ' is a not a valid mode.  Choices are gatk or grep.');
              }
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: GATK Select Manta PASS
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 931.1585693359375
    'sbg:y': 764.6875
  - id: gatk_selectvariants_mutect2
    in:
      - id: input_vcf
        source: filter_mutect2_vcf/filtered_vcf
      - id: mode
        source: select_vars_mode
      - id: output_basename
        source: output_basename
      - id: tool_name
        valueFrom: '${return "mutect2"}'
    out:
      - id: pass_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_selectvariants
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - default: gatk
          id: mode
          type:
            - 'null'
            - type: enum
              symbols:
                - gatk
                - grep
              name: select_vars_mode
          doc: 'Choose ''gatk'' for SelectVariants tool, or ''grep'' for grep expression'
        - id: output_basename
          type: string
        - id: tool_name
          type: string
      outputs:
        - id: pass_vcf
          type: File
          outputBinding:
            glob: '*.PASS.vcf.gz'
          secondaryFiles:
            - .tbi
      label: GATK Select PASS
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: |-
            set -eo pipefail
            ${
              var run_mode = inputs.mode;
              if (run_mode == 'grep' || run_mode == 'gatk'){
                var in_vcf = inputs.input_vcf.path;
                var out_vcf = inputs.output_basename + '.' + inputs.tool_name + '.PASS.vcf.gz';
                var cmd = '/gatk SelectVariants --java-options "-Xmx8000m" -V ' + in_vcf +  ' -O ' + out_vcf + ' --exclude-filtered TRUE';
                if (run_mode == 'grep'){
                  cmd = 'zcat ' + in_vcf + ' | grep -E "^#|PASS" | bgzip > ' + out_vcf + '; tabix ' + out_vcf;
                }
                return cmd;
              }
              else{
                throw new Error(run_mode + ' is a not a valid mode.  Choices are gatk or grep.');
              }
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: GATK Select Mutect2 PASS
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 1652.5955810546875
    'sbg:y': 534.375
  - id: gatk_selectvariants_strelka2
    in:
      - id: input_vcf
        source: rename_strelka_samples/reheadered_vcf
      - id: mode
        source: select_vars_mode
      - id: output_basename
        source: output_basename
      - id: tool_name
        valueFrom: '${return "strelka2"}'
    out:
      - id: pass_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_selectvariants
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - default: gatk
          id: mode
          type:
            - 'null'
            - type: enum
              symbols:
                - gatk
                - grep
              name: select_vars_mode
          doc: 'Choose ''gatk'' for SelectVariants tool, or ''grep'' for grep expression'
        - id: output_basename
          type: string
        - id: tool_name
          type: string
      outputs:
        - id: pass_vcf
          type: File
          outputBinding:
            glob: '*.PASS.vcf.gz'
          secondaryFiles:
            - .tbi
      label: GATK Select PASS
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: |-
            set -eo pipefail
            ${
              var run_mode = inputs.mode;
              if (run_mode == 'grep' || run_mode == 'gatk'){
                var in_vcf = inputs.input_vcf.path;
                var out_vcf = inputs.output_basename + '.' + inputs.tool_name + '.PASS.vcf.gz';
                var cmd = '/gatk SelectVariants --java-options "-Xmx8000m" -V ' + in_vcf +  ' -O ' + out_vcf + ' --exclude-filtered TRUE';
                if (run_mode == 'grep'){
                  cmd = 'zcat ' + in_vcf + ' | grep -E "^#|PASS" | bgzip > ' + out_vcf + '; tabix ' + out_vcf;
                }
                return cmd;
              }
              else{
                throw new Error(run_mode + ' is a not a valid mode.  Choices are gatk or grep.');
              }
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: GATK Select Strelka2 PASS
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 1334.072265625
    'sbg:y': 492.375
  - id: manta
    in:
      - id: hg38_strelka_bed
        source: hg38_strelka_bed
      - id: input_normal_cram
        source: input_normal_aligned
      - id: input_tumor_cram
        source: input_tumor_aligned
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
    out:
      - id: output_sv
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: kfdrc_manta_sv
      baseCommand:
        - /manta-1.4.0.centos6_x86_64/bin/configManta.py
      inputs:
        - default: 18
          id: cores
          type: int?
        - id: hg38_strelka_bed
          type: File
          secondaryFiles:
            - .tbi
        - id: input_normal_cram
          type: File?
          secondaryFiles:
            - .crai
        - id: input_tumor_cram
          type: File?
          secondaryFiles:
            - .crai
        - id: output_basename
          type: string
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
      outputs:
        - id: output_sv
          type: File
          outputBinding:
            glob: '*SV.vcf.gz'
          secondaryFiles:
            - .tbi
      doc: >-
        Calls structural variants.  Tool designed to pick correct run mode based
        on if tumor, normal, or both crams are given
      label: Manta sv caller
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: |-
            ${
              var std = " --ref " + inputs.reference.path + " --callRegions " + inputs.hg38_strelka_bed.path + " --runDir=./ && ./runWorkflow.py -m local -j " + inputs.cores + " --quiet ";
              var mv = " && mv results/variants/";
              if (typeof inputs.input_tumor_cram === 'undefined' || inputs.input_tumor_cram === null){
                var mv_cmd = mv + "diploidSV.vcf.gz " +  inputs.output_basename + ".manta.diploidSV.vcf.gz" + mv + "diploidSV.vcf.gz.tbi " + inputs.output_basename + ".manta.diploidSV.vcf.gz.tbi";
                return "--bam ".concat(inputs.input_normal_cram.path, std, mv_cmd);
              }
              else if (typeof inputs.input_normal_cram === 'undefined' || inputs.input_normal_cram === null){
                var mv_cmd = mv + "tumorSV.vcf.gz " + inputs.output_basename + ".manta.tumorSV.vcf.gz" + mv + "tumorSV.vcf.gz.tbi " + inputs.output_basename + ".manta.tumorSV.vcf.gz.tbi";
                return "--tumorBam " + inputs.input_tumor_cram.path + std + mv_cmd;
              }
              else{
                var mv_cmd = mv + "somaticSV.vcf.gz " + inputs.output_basename + ".manta.somaticSV.vcf.gz" + mv + "somaticSV.vcf.gz.tbi " + inputs.output_basename + ".manta.somaticSV.vcf.gz.tbi";
                return "--tumorBam " + inputs.input_tumor_cram.path + " --normalBam " + inputs.input_normal_cram.path + std + mv_cmd;
              }
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 10000
          coresMin: $(inputs.cores)
        - class: DockerRequirement
          dockerPull: 'kfdrc/manta:1.4.0'
        - class: InlineJavascriptRequirement
    label: Manta sv caller
    'sbg:x': 247.3125
    'sbg:y': 253.1875
  - id: merge_mutect2_stats
    in:
      - id: input_stats
        source:
          - mutect2/mutect_stats
      - id: output_basename
        source: output_basename
    out:
      - id: merged_stats
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_mergepileup
      baseCommand:
        - /gatk
        - MergeMutectStats
      inputs:
        - id: input_stats
          type:
            type: array
            items: File
            inputBinding:
              prefix: '--stats'
          inputBinding:
            position: 1
        - id: output_basename
          type: string
      outputs:
        - id: merged_stats
          type: File
          outputBinding:
            glob: '*.merged.stats'
      label: GATK Merge Stats
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            --java-options "-Xmx3000m" -O
            $(inputs.output_basename).Mutect2.merged.stats 
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 4000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: Merge mutect2 stats
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 931.1585693359375
    'sbg:y': 529.9375
  - id: merge_mutect2_vcf
    in:
      - id: input_vcfs
        source:
          - mutect2/mutect2_vcf
      - id: output_basename
        source: output_basename
      - id: reference_dict
        source: reference_dict
      - id: tool_name
        valueFrom: '${return "mutect2"}'
    out:
      - id: merged_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_mergevcfs
      baseCommand:
        - /gatk
        - MergeVcfs
      inputs:
        - id: input_vcfs
          type:
            type: array
            items: File
            inputBinding:
              prefix: '-I'
          inputBinding:
            position: 1
          secondaryFiles:
            - .tbi
        - id: output_basename
          type: string
        - id: reference_dict
          type: File
        - id: silent_flag
          type: int?
        - id: tool_name
          type: string
      outputs:
        - id: merged_vcf
          type: File
          outputBinding:
            glob: '*.merged.vcf.gz'
          secondaryFiles:
            - .tbi
      doc: Merge input vcfs
      label: GATK Merge VCF
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            --java-options "-Xmx2000m" --TMP_DIR=./TMP --CREATE_INDEX=true
            --SEQUENCE_DICTIONARY=$(inputs.reference_dict.path) ${
              var cmd = "--OUTPUT=" + inputs.output_basename + "." + inputs.tool_name + ".merged.vcf.gz "
              if (typeof inputs.silent_flag !== 'undefined' && inputs.silent_flag == 1){
                cmd += "--VALIDATION_STRINGENCY SILENT"
              }
              return cmd
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 4000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: Merge & pass filter mutect2
    'sbg:x': 931.1585693359375
    'sbg:y': 402.0625
  - id: merge_strelka2_vcf
    in:
      - id: input_vcfs
        source:
          - strelka2/output_snv
          - strelka2/output_indel
      - id: output_basename
        source: output_basename
      - id: reference_dict
        source: reference_dict
      - id: tool_name
        valueFrom: '${ return "strelka2"}'
    out:
      - id: merged_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_mergevcfs
      baseCommand:
        - /gatk
        - MergeVcfs
      inputs:
        - id: input_vcfs
          type:
            type: array
            items: File
            inputBinding:
              prefix: '-I'
          inputBinding:
            position: 1
          secondaryFiles:
            - .tbi
        - id: output_basename
          type: string
        - id: reference_dict
          type: File
        - id: silent_flag
          type: int?
        - id: tool_name
          type: string
      outputs:
        - id: merged_vcf
          type: File
          outputBinding:
            glob: '*.merged.vcf.gz'
          secondaryFiles:
            - .tbi
      doc: Merge input vcfs
      label: GATK Merge VCF
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            --java-options "-Xmx2000m" --TMP_DIR=./TMP --CREATE_INDEX=true
            --SEQUENCE_DICTIONARY=$(inputs.reference_dict.path) ${
              var cmd = "--OUTPUT=" + inputs.output_basename + "." + inputs.tool_name + ".merged.vcf.gz "
              if (typeof inputs.silent_flag !== 'undefined' && inputs.silent_flag == 1){
                cmd += "--VALIDATION_STRINGENCY SILENT"
              }
              return cmd
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 4000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: Merge & pass filter strekla2
    'sbg:x': 567.62158203125
    'sbg:y': 529.9375
  - id: mutect2
    in:
      - id: af_only_gnomad_vcf
        source: af_only_gnomad_vcf
      - id: exome_flag
        source: exome_flag
      - id: input_normal_aligned
        source: input_normal_aligned
      - id: input_normal_name
        source: input_normal_name
      - id: input_tumor_aligned
        source: input_tumor_aligned
      - id: input_tumor_name
        source: input_tumor_name
      - id: interval_list
        source: gatk_intervallisttools/output
      - id: reference
        source: indexed_reference_fasta
    out:
      - id: f1r2_counts
      - id: mutect2_vcf
      - id: mutect_stats
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4__mutect2
      baseCommand:
        - /gatk
        - Mutect2
      inputs:
        - id: af_only_gnomad_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: exome_flag
          type: string?
          doc: 'Y if exome/capture, defaults to WGS'
        - id: input_normal_aligned
          type: File
          doc: normal BAM or CRAM
          secondaryFiles:
            - |
              ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                  return {"location": dpath+self.nameroot+".bai", "class": "File"}
                }
                else{
                  return {"location": dpath+self.basename+".crai", "class": "File"}
                }
              }
        - id: input_normal_name
          type: string
        - id: input_tumor_aligned
          type: File
          doc: tumor BAM or CRAM
          secondaryFiles:
            - |
              ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                  return {"location": dpath+self.nameroot+".bai", "class": "File"}
                }
                else{
                  return {"location": dpath+self.basename+".crai", "class": "File"}
                }
              }
        - id: input_tumor_name
          type: string
        - id: interval_list
          type: File
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
      outputs:
        - id: f1r2_counts
          type: File
          outputBinding:
            glob: '*.f1r2_counts.tar.gz'
        - id: mutect2_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
        - id: mutect_stats
          type: File
          outputBinding:
            glob: '*.stats'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --java-options "-Xmx6000m" -R $(inputs.reference.path) -I
            $(inputs.input_tumor_aligned.path) -I
            $(inputs.input_normal_aligned.path) -tumor
            $(inputs.input_tumor_name) -normal $(inputs.input_normal_name)
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -L
            $(inputs.interval_list.path) --germline-resource
            $(inputs.af_only_gnomad_vcf.path) --f1r2-tar-gz
            $(inputs.input_tumor_aligned.nameroot).$(inputs.interval_list.nameroot).f1r2_counts.tar.gz
            ${
              var arg = "-O " + inputs.input_tumor_aligned.nameroot + "." + inputs.interval_list.nameroot + ".Mutect2.vcf.gz"
              if (inputs.exome_flag == 'Y'){
                arg += " --disable-adaptive-pruning"
              }
              return arg
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 6000
          coresMin: 3
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    scatter:
      - interval_list
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;500
    'sbg:x': 567.62158203125
    'sbg:y': 360.0625
  - id: mutect2_filter_support
    in:
      - id: exac_common_vcf
        source: exac_common_vcf
      - id: f1r2_counts
        source:
          - mutect2/f1r2_counts
      - id: indexed_reference_fasta
        source: indexed_reference_fasta
      - id: input_normal_aligned
        source: input_normal_aligned
      - id: input_tumor_aligned
        source: input_tumor_aligned
      - id: output_basename
        source: output_basename
      - id: reference_dict
        source: reference_dict
      - id: wgs_calling_interval_list
        source:
          - gatk_intervallisttools/output
    out:
      - id: contamination_table
      - id: f1r2_bias
      - id: segmentation_table
    run:
      class: Workflow
      cwlVersion: v1.0
      id: mutect2_support
      $namespaces:
        sbg: 'https://sevenbridges.com'
      inputs:
        - id: exac_common_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: f1r2_counts
          type: 'File[]'
          doc: orientation counts from mutect2 outputs
        - id: indexed_reference_fasta
          type: File
          secondaryFiles:
            - .fai
            - ^.dict
        - id: input_normal_aligned
          type: File
          doc: normal BAM or CRAM
          secondaryFiles:
            - |
              ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                  return {"location": dpath+self.nameroot+".bai", "class": "File"}
                }
                else{
                  return {"location": dpath+self.basename+".crai", "class": "File"}
                }
              }
        - id: input_tumor_aligned
          type: File
          doc: tumor BAM or CRAM
          secondaryFiles:
            - |
              ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                  return {"location": dpath+self.nameroot+".bai", "class": "File"}
                }
                else{
                  return {"location": dpath+self.basename+".crai", "class": "File"}
                }
              }
        - id: output_basename
          type: string
        - id: reference_dict
          type: File
        - id: wgs_calling_interval_list
          type: 'File[]'
      outputs:
        - id: contamination_table
          outputSource:
            - gatk_calculate_contamination/contamination_table
          type: File
        - id: f1r2_bias
          outputSource:
            - gatk_learn_orientation_bias/f1r2_bias
          type: File
        - id: segmentation_table
          outputSource:
            - gatk_calculate_contamination/segmentation_table
          type: File
      steps:
        - id: gatk_calculate_contamination
          in:
            - id: normal_pileup
              source: gatk_gather_normal_pileup_summaries/merged_table
            - id: output_basename
              source: output_basename
            - id: tumor_pileup
              source: gatk_gather_tumor_pileup_summaries/merged_table
          out:
            - id: contamination_table
            - id: segmentation_table
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            id: gatk4_calulcate_contamination
            baseCommand:
              - /gatk
              - CalculateContamination
            inputs:
              - id: normal_pileup
                type: File
              - id: output_basename
                type: string
              - id: tumor_pileup
                type: File
            outputs:
              - id: contamination_table
                type: File
                outputBinding:
                  glob: '*.contamination.table'
              - id: segmentation_table
                type: File
                outputBinding:
                  glob: '*.segmentation.table'
            label: GATK Calculate Contamination
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: >-
                  --java-options "-Xmx4000m" -I $(inputs.tumor_pileup.path)
                  --matched-normal $(inputs.normal_pileup.path) -O
                  $(inputs.output_basename).contamination.table
                  --tumor-segmentation
                  $(inputs.output_basename).segmentation.table
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 4000
                coresMin: 2
              - class: DockerRequirement
                dockerPull: 'kfdrc/gatk:4.1.1.0'
              - class: InlineJavascriptRequirement
          label: GATK Calculate Contamination
        - id: gatk_gather_normal_pileup_summaries
          in:
            - id: input_tables
              source:
                - gatk_get_normal_pileup_summaries/pileup_table
            - id: output_basename
              source: output_basename
            - id: reference_dict
              source: reference_dict
            - id: tool_name
              valueFrom: '${return "mutect2"}'
          out:
            - id: merged_table
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            id: gatk4_mergepileup
            baseCommand:
              - /gatk
              - GatherPileupSummaries
            inputs:
              - id: input_tables
                type:
                  type: array
                  items: File
                  inputBinding:
                    prefix: '-I'
                inputBinding:
                  position: 1
              - id: output_basename
                type: string
              - id: reference_dict
                type: File
              - id: tool_name
                type: string
            outputs:
              - id: merged_table
                type: File
                outputBinding:
                  glob: '*.merged.pileup.table'
            label: GATK Merge Pileups
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: >-
                  --java-options "-Xmx3000m" --sequence-dictionary
                  $(inputs.reference_dict.path) -O
                  $(inputs.output_basename).$(inputs.tool_name).merged.pileup.table 
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 4000
                coresMin: 2
              - class: DockerRequirement
                dockerPull: 'kfdrc/gatk:4.1.1.0'
              - class: InlineJavascriptRequirement
          label: GATK merge normal pileup tables
        - id: gatk_gather_tumor_pileup_summaries
          in:
            - id: input_tables
              source:
                - gatk_get_tumor_pileup_summaries/pileup_table
            - id: output_basename
              source: output_basename
            - id: reference_dict
              source: reference_dict
            - id: tool_name
              valueFrom: '${return "mutect2"}'
          out:
            - id: merged_table
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            id: gatk4_mergepileup
            baseCommand:
              - /gatk
              - GatherPileupSummaries
            inputs:
              - id: input_tables
                type:
                  type: array
                  items: File
                  inputBinding:
                    prefix: '-I'
                inputBinding:
                  position: 1
              - id: output_basename
                type: string
              - id: reference_dict
                type: File
              - id: tool_name
                type: string
            outputs:
              - id: merged_table
                type: File
                outputBinding:
                  glob: '*.merged.pileup.table'
            label: GATK Merge Pileups
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: >-
                  --java-options "-Xmx3000m" --sequence-dictionary
                  $(inputs.reference_dict.path) -O
                  $(inputs.output_basename).$(inputs.tool_name).merged.pileup.table 
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 4000
                coresMin: 2
              - class: DockerRequirement
                dockerPull: 'kfdrc/gatk:4.1.1.0'
              - class: InlineJavascriptRequirement
          label: GATK merge tumor pileup tables
        - id: gatk_get_normal_pileup_summaries
          in:
            - id: aligned_reads
              source: input_normal_aligned
            - id: exac_common_vcf
              source: exac_common_vcf
            - id: interval_list
              source: wgs_calling_interval_list
            - id: reference
              source: indexed_reference_fasta
          out:
            - id: pileup_table
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            id: gatk4_getpileupsummary
            baseCommand:
              - /gatk
              - GetPileupSummaries
            inputs:
              - id: aligned_reads
                type: File
                secondaryFiles:
                  - .crai
              - id: exac_common_vcf
                type: File
                secondaryFiles:
                  - .tbi
              - id: interval_list
                type: File
              - id: reference
                type: File
            outputs:
              - id: pileup_table
                type: File
                outputBinding:
                  glob: '*.pileupsummary.table'
            label: GATK Pileup
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: >-
                  --java-options "-Xmx2000m" -I $(inputs.aligned_reads.path) -V
                  $(inputs.exac_common_vcf.path) -L $(inputs.interval_list.path)
                  -R $(inputs.reference.path) -O
                  $(inputs.aligned_reads.nameroot).pileupsummary.table
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 2000
                coresMin: 2
              - class: DockerRequirement
                dockerPull: 'kfdrc/gatk:4.1.1.0'
              - class: InlineJavascriptRequirement
          label: GATK normal pileup scatter
          scatter:
            - interval_list
        - id: gatk_get_tumor_pileup_summaries
          in:
            - id: aligned_reads
              source: input_tumor_aligned
            - id: exac_common_vcf
              source: exac_common_vcf
            - id: interval_list
              source: wgs_calling_interval_list
            - id: reference
              source: indexed_reference_fasta
          out:
            - id: pileup_table
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            id: gatk4_getpileupsummary
            baseCommand:
              - /gatk
              - GetPileupSummaries
            inputs:
              - id: aligned_reads
                type: File
                secondaryFiles:
                  - .crai
              - id: exac_common_vcf
                type: File
                secondaryFiles:
                  - .tbi
              - id: interval_list
                type: File
              - id: reference
                type: File
            outputs:
              - id: pileup_table
                type: File
                outputBinding:
                  glob: '*.pileupsummary.table'
            label: GATK Pileup
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: >-
                  --java-options "-Xmx2000m" -I $(inputs.aligned_reads.path) -V
                  $(inputs.exac_common_vcf.path) -L $(inputs.interval_list.path)
                  -R $(inputs.reference.path) -O
                  $(inputs.aligned_reads.nameroot).pileupsummary.table
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 2000
                coresMin: 2
              - class: DockerRequirement
                dockerPull: 'kfdrc/gatk:4.1.1.0'
              - class: InlineJavascriptRequirement
          label: GATK tumor pileup scatter
          scatter:
            - interval_list
        - id: gatk_learn_orientation_bias
          in:
            - id: input_tgz
              source:
                - f1r2_counts
            - id: output_basename
              source: output_basename
            - id: tool_name
              valueFrom: '${return "mutect2"}'
          out:
            - id: f1r2_bias
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            id: gatk4_learn_oritentation_bias
            baseCommand:
              - /gatk
              - LearnReadOrientationModel
            inputs:
              - id: input_tgz
                type:
                  type: array
                  items: File
                  inputBinding:
                    prefix: '-I'
                inputBinding:
                  position: 1
              - id: output_basename
                type: string
              - id: tool_name
                type: string
            outputs:
              - id: f1r2_bias
                type: File
                outputBinding:
                  glob: '*.f1r2_bias.tar.gz'
            label: GATK Learn Bias
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: >-
                  --java-options "-Xmx4000m" -O
                  $(inputs.output_basename).$(inputs.tool_name).f1r2_bias.tar.gz 
            requirements:
              - class: ShellCommandRequirement
              - class: ResourceRequirement
                ramMin: 4000
                coresMin: 2
              - class: DockerRequirement
                dockerPull: 'kfdrc/gatk:4.1.1.0'
              - class: InlineJavascriptRequirement
          label: Gatk learn bias
      requirements:
        - class: ScatterFeatureRequirement
    'sbg:x': 931.1585693359375
    'sbg:y': 232.1875
  - id: rename_manta_samples
    in:
      - id: input_normal_name
        source: input_normal_name
      - id: input_tumor_name
        source: input_tumor_name
      - id: input_vcf
        source: manta/output_sv
    out:
      - id: reheadered_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: bcftools_reheader_vcf
      baseCommand:
        - echo
      inputs:
        - id: input_normal_name
          type: string
        - id: input_tumor_name
          type: string
        - id: input_vcf
          type: File
      outputs:
        - id: reheadered_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.input_normal_name) > sample_list.txt && echo
            $(inputs.input_tumor_name) >> sample_list.txt && bcftools reheader
            -s sample_list.txt $(inputs.input_vcf.path) >
            $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz"))
            && tabix $(inputs.input_vcf.nameroot.replace(".vcf",
            ".reheadered.vcf.gz"))
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 1000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'kfdrc/bvcftools:latest'
        - class: InlineJavascriptRequirement
    'sbg:x': 567.62158203125
    'sbg:y': 190.1875
  - id: rename_strelka_samples
    in:
      - id: input_normal_name
        source: input_normal_name
      - id: input_tumor_name
        source: input_tumor_name
      - id: input_vcf
        source: merge_strelka2_vcf/merged_vcf
    out:
      - id: reheadered_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: bcftools_reheader_vcf
      baseCommand:
        - echo
      inputs:
        - id: input_normal_name
          type: string
        - id: input_tumor_name
          type: string
        - id: input_vcf
          type: File
      outputs:
        - id: reheadered_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.input_normal_name) > sample_list.txt && echo
            $(inputs.input_tumor_name) >> sample_list.txt && bcftools reheader
            -s sample_list.txt $(inputs.input_vcf.path) >
            $(inputs.input_vcf.nameroot.replace(".vcf", ".reheadered.vcf.gz"))
            && tabix $(inputs.input_vcf.nameroot.replace(".vcf",
            ".reheadered.vcf.gz"))
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 1000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'kfdrc/bvcftools:latest'
        - class: InlineJavascriptRequirement
    'sbg:x': 931.1585693359375
    'sbg:y': 62.3125
  - id: strelka2
    in:
      - id: exome_flag
        source: exome_flag
      - id: hg38_strelka_bed
        source: hg38_strelka_bed
      - id: input_normal_aligned
        source: input_normal_aligned
      - id: input_tumor_aligned
        source: input_tumor_aligned
      - id: reference
        source: indexed_reference_fasta
    out:
      - id: output_indel
      - id: output_snv
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: strelka2
      baseCommand:
        - /strelka-2.9.3.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py
      inputs:
        - default: 18
          id: cores
          type: int?
        - id: exome_flag
          type: string?
          doc: 'Y if exome/capture, defaults to WGS'
        - id: hg38_strelka_bed
          type: File
          label: gzipped bed file
          secondaryFiles:
            - .tbi
        - id: input_normal_aligned
          type: File
          doc: normal BAM or CRAM
          secondaryFiles:
            - |
              ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                  return {"location": dpath+self.nameroot+".bai", "class": "File"}
                }
                else{
                  return {"location": dpath+self.basename+".crai", "class": "File"}
                }
              }
        - id: input_tumor_aligned
          type: File
          doc: tumor BAM or CRAM
          secondaryFiles:
            - |
              ${
                var dpath = self.location.replace(self.basename, "")
                if(self.nameext == '.bam'){
                  return {"location": dpath+self.nameroot+".bai", "class": "File"}
                }
                else{
                  return {"location": dpath+self.basename+".crai", "class": "File"}
                }
              }
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
      outputs:
        - id: output_indel
          type: File
          outputBinding:
            glob: results/variants/*.indels.vcf.gz
          secondaryFiles:
            - .tbi
        - id: output_snv
          type: File
          outputBinding:
            glob: results/variants/*.snvs.vcf.gz
          secondaryFiles:
            - .tbi
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --normalBam $(inputs.input_normal_aligned.path) --tumorBam
            $(inputs.input_tumor_aligned.path) --ref $(inputs.reference.path)
            --callRegions $(inputs.hg38_strelka_bed.path) ${
              var arg = "--runDir=./";
              if (inputs.exome_flag == 'Y'){
                arg += " --exome"
              }
              return arg
            } && ./runWorkflow.py -m local -j $(inputs.cores)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 10000
          coresMin: $(inputs.cores)
        - class: DockerRequirement
          dockerPull: obenauflab/strelka
        - class: InlineJavascriptRequirement
    'sbg:x': 247.3125
    'sbg:y': 90.3125
  - id: vep_annot_mutect2
    in:
      - id: cache
        source: vep_cache
      - id: input_vcf
        source: gatk_selectvariants_mutect2/pass_vcf
      - id: normal_id
        source: input_normal_name
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
      - id: tool_name
        valueFrom: '${return "mutect2_somatic"}'
      - id: tumor_id
        source: input_tumor_name
    out:
      - id: output_maf
      - id: output_tbi
      - id: output_vcf
      - id: warn_txt
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: kfdrc_vep_somatic_annotate_maf
      baseCommand:
        - tar
        - '-xzf'
      inputs:
        - id: cache
          type: File
          label: tar gzipped cache from ensembl/local converted cache
        - default: 93
          id: cache_version
          type: int?
          doc: 'Version being used, should match build version'
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: normal_id
          type: string
        - id: output_basename
          type: string
        - default: GRCh38
          id: ref_build
          type: string?
          doc: 'Genome ref build used, should line up with cache.'
        - id: reference
          type: File
          label: Fasta genome assembly with index
          secondaryFiles:
            - .fai
        - id: tool_name
          type: string
        - id: tumor_id
          type: string
      outputs:
        - id: output_maf
          type: File
          outputBinding:
            glob: '*.maf'
        - id: output_tbi
          type: File
          outputBinding:
            glob: '*.vcf.gz.tbi'
        - id: output_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
        - id: warn_txt
          type: File?
          outputBinding:
            glob: '*.txt'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.cache.path) && gunzip -c $(inputs.input_vcf.path) >
            input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf
            input_file.vcf --output-maf
            $(inputs.output_basename).$(inputs.tool_name).vep.maf --filter-vcf 0
            --vep-path /ensembl-vep/ --vep-data $PWD --vep-forks 16 --ncbi-build
            $(inputs.ref_build) --cache-version $(inputs.cache_version)
            --ref-fasta $(inputs.reference.path) --tumor-id $(inputs.tumor_id)
            --normal-id $(inputs.normal_id) && mv input_file.vep.vcf
            $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf &&
            /ensembl-vep/htslib/bgzip
            $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf &&
            /ensembl-vep/htslib/tabix
            $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf.gz
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 24000
          coresMin: 16
        - class: DockerRequirement
          dockerPull: 'kfdrc/vep:r93_v2'
        - class: InlineJavascriptRequirement
    'sbg:x': 2126.01904296875
    'sbg:y': 232.1875
  - id: vep_annot_strelka2
    in:
      - id: cache
        source: vep_cache
      - id: input_vcf
        source: gatk_selectvariants_strelka2/pass_vcf
      - id: normal_id
        source: input_normal_name
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
      - id: tool_name
        valueFrom: '${return "strelka2_somatic"}'
      - id: tumor_id
        source: input_tumor_name
    out:
      - id: output_maf
      - id: output_tbi
      - id: output_vcf
      - id: warn_txt
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: kfdrc_vep_somatic_annotate_maf
      baseCommand:
        - tar
        - '-xzf'
      inputs:
        - id: cache
          type: File
          label: tar gzipped cache from ensembl/local converted cache
        - default: 93
          id: cache_version
          type: int?
          doc: 'Version being used, should match build version'
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: normal_id
          type: string
        - id: output_basename
          type: string
        - default: GRCh38
          id: ref_build
          type: string?
          doc: 'Genome ref build used, should line up with cache.'
        - id: reference
          type: File
          label: Fasta genome assembly with index
          secondaryFiles:
            - .fai
        - id: tool_name
          type: string
        - id: tumor_id
          type: string
      outputs:
        - id: output_maf
          type: File
          outputBinding:
            glob: '*.maf'
        - id: output_tbi
          type: File
          outputBinding:
            glob: '*.vcf.gz.tbi'
        - id: output_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
        - id: warn_txt
          type: File?
          outputBinding:
            glob: '*.txt'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.cache.path) && gunzip -c $(inputs.input_vcf.path) >
            input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf
            input_file.vcf --output-maf
            $(inputs.output_basename).$(inputs.tool_name).vep.maf --filter-vcf 0
            --vep-path /ensembl-vep/ --vep-data $PWD --vep-forks 16 --ncbi-build
            $(inputs.ref_build) --cache-version $(inputs.cache_version)
            --ref-fasta $(inputs.reference.path) --tumor-id $(inputs.tumor_id)
            --normal-id $(inputs.normal_id) && mv input_file.vep.vcf
            $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf &&
            /ensembl-vep/htslib/bgzip
            $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf &&
            /ensembl-vep/htslib/tabix
            $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf.gz
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 24000
          coresMin: 16
        - class: DockerRequirement
          dockerPull: 'kfdrc/vep:r93_v2'
        - class: InlineJavascriptRequirement
    'sbg:x': 1652.5955810546875
    'sbg:y': 271.625
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
