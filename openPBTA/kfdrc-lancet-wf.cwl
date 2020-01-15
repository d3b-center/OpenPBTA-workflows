class: Workflow
cwlVersion: v1.0
id: cavatica/openpbta-tcga/kfdrc-lancet-wf-baminput/0
label: kfdrc-lancet-wf
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: calling_interval_list
    type: File
    'sbg:x': 0
    'sbg:y': 963
  - id: exome_flag
    type: string?
    doc: 'set to ''Y'' for exome mode, most likely given run time'
    'sbg:x': 314.328125
    'sbg:y': 521
  - id: indexed_reference_fasta
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    'sbg:x': 0
    'sbg:y': 856
  - id: input_normal_aligned
    type: File
    'sbg:x': 0
    'sbg:y': 749
  - id: input_normal_name
    type: string
    'sbg:x': 1406.600830078125
    'sbg:y': 521
  - id: input_tumor_aligned
    type: File
    'sbg:x': 0
    'sbg:y': 642
  - id: input_tumor_name
    type: string
    'sbg:x': 1406.600830078125
    'sbg:y': 414
  - id: mutect2_vcf
    type: File?
    doc: >-
      PASS vcf from mutect2 run for the sample to be analyzed, Optional and
      recommneded to augment an exome interval list
    'sbg:x': 0
    'sbg:y': 535
  - id: output_basename
    type: string
    'sbg:x': 804.6476440429688
    'sbg:y': 439.5
  - id: padding
    type: int
    doc: 'If WGS (less likely), default 25, if exome+, recommend hald window size'
    'sbg:x': 0
    'sbg:y': 428
  - id: reference_dict
    type: File
    'sbg:x': 804.6476440429688
    'sbg:y': 332.5
  - id: select_vars_mode
    type: string
    doc: 'Choose ''gatk'' for SelectVariants tool, or ''grep'' for grep expression'
    'sbg:x': 0
    'sbg:y': 321
  - id: strelka2_vcf
    type: File?
    doc: >-
      PASS vcf from strelka2 run for the sample to be analyzed. Optional and
      recommneded to augment an exome interval list
    'sbg:x': 0
    'sbg:y': 214
  - id: vep_cache
    type: File
    label: tar gzipped cache from ensembl/local converted cache
    'sbg:x': 0
    'sbg:y': 107
  - id: window
    type: int
    doc: >-
      window size for lancet.  default is 600, recommend 500 for WGS, 600 for
      exome+
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: lancet_prepass_vcf
    outputSource:
      - sort_merge_lancet_vcf/merged_vcf
    type: File
    'sbg:x': 1406.600830078125
    'sbg:y': 307
  - id: lancet_vep_maf
    outputSource:
      - vep_annot_lancet/output_maf
    type: File
    'sbg:x': 2133.03515625
    'sbg:y': 588.5
  - id: lancet_vep_tbi
    outputSource:
      - vep_annot_lancet/output_tbi
    type: File
    'sbg:x': 2133.03515625
    'sbg:y': 481.5
  - id: lancet_vep_vcf
    outputSource:
      - vep_annot_lancet/output_vcf
    type: File
    'sbg:x': 2133.03515625
    'sbg:y': 374.5
steps:
  - id: bedops_gen_lancet_intervals
    in:
      - id: mutect2_vcf
        source: mutect2_vcf
      - id: output_basename
        source: output_basename
      - id: ref_bed
        source: calling_interval_list
      - id: strelka2_vcf
        source: strelka2_vcf
    out:
      - id: run_bed
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: bcftools_reheader_vcf
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: mutect2_vcf
          type: File?
          doc: PASS vcf from mutect2 run for the sample to be analyzed
        - id: output_basename
          type: string
        - id: ref_bed
          type: File
          doc: >-
            Exonic bed file recommended.  Will be augmented with strelka2 and/or
            mutect2 sites if supplied
        - id: strelka2_vcf
          type: File?
          doc: PASS vcf from strelka2 run for the sample to be analyzed
      outputs:
        - id: run_bed
          type: File
          outputBinding:
            glob: '*.lancet_intvervals.bed'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: |-
            set -eo pipefail
            ${
                var flag = 0
                var cmd = "";
                var bed = []
                if(inputs.strelka2_vcf != null){
                    flag = 1;
                    cmd += "gunzip -c " + inputs.strelka2_vcf.path + " > " + inputs.strelka2_vcf.basename + ";";
                    cmd += "vcf2bed --insertions < " + inputs.strelka2_vcf.basename + " | cut -f 1-3 > strelka2.insertions.bed;";
                    bed.push("strelka2.insertions.bed");
                    cmd += "vcf2bed --deletions < " + inputs.strelka2_vcf.basename + " | cut -f 1-3 > strelka2.deletions.bed;";
                    bed.push("strelka2.deletions.bed");
                    cmd += "vcf2bed --snvs < " + inputs.strelka2_vcf.basename + " | cut -f 1-3 > strelka2.snvs.bed;";
                    bed.push("strelka2.snvs.bed");
                }
                if(inputs.mutect2_vcf != null){
                    flag = 1;
                    cmd += "gunzip -c " + inputs.mutect2_vcf.path + " > " + inputs.mutect2_vcf.basename + ";";
                    cmd += "vcf2bed --insertions < " + inputs.mutect2_vcf.basename +  " | cut -f 1-3 > mutect2.insertions.bed;";
                    bed.push("mutect2.insertions.bed");
                    cmd += "vcf2bed --deletions < " + inputs.mutect2_vcf.basename + " | cut -f 1-3 > mutect2.deletions.bed;";
                    bed.push("mutect2.deletions.bed");
                    cmd += "vcf2bed --snvs < " + inputs.mutect2_vcf.basename +  " | cut -f 1-3 > mutect2.snvs.bed;";
                    bed.push("mutect2.snvs.bed");
                }
              if(flag == 0){
                  cmd += "echo \"No input vcfs found to convert.  Returning ref bed\"; cp " + inputs.ref_bed.path + " " + inputs.output_basename + ".lancet_intvervals.bed;";
              }
              else{
                  cmd += "cat " + bed.join(" ") + " " + inputs.ref_bed.path + " | bedtools sort | bedtools merge > " + inputs.output_basename + ".lancet_intvervals.bed;";
              }
              return cmd;
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/bedops:2.4.36'
        - class: InlineJavascriptRequirement
    'sbg:x': 314.328125
    'sbg:y': 649
  - id: gatk_intervallisttools
    in:
      - id: bands
        valueFrom: '${return 80000000}'
      - id: exome_flag
        source: exome_flag
      - id: interval_list
        source: bedops_gen_lancet_intervals/run_bed
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
        - id: exome_flag
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
            glob: temp*/*.bed
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
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
            } /gatk IntervalListTools --java-options "-Xmx2000m"
            --SCATTER_COUNT=$(inputs.scatter_ct)
            --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
            --UNIQUE=true --SORT=true --BREAK_BANDS_AT_MULTIPLES_OF=$BANDS
            --INPUT=$LIST --OUTPUT=. && CT=`find . -name 'temp_0*' | wc -l` &&
            seq -f "%04g" $CT | xargs -I N -P 4 /gatk IntervalListToBed
            --java-options -Xmx100m -I temp_N_of_$CT/scattered.interval_list -O
            temp_N_of_$CT/scattered.interval_list.N.bed
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 2000
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 585.2146606445312
    'sbg:y': 467.5
  - id: gatk_selectvariants_lancet
    in:
      - id: input_vcf
        source: sort_merge_lancet_vcf/merged_vcf
      - id: mode
        source: select_vars_mode
      - id: output_basename
        source: output_basename
      - id: tool_name
        valueFrom: '${return "lancet"}'
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
        - id: mode
          type: string
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
    label: GATK Select Lancet PASS
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 1406.600830078125
    'sbg:y': 642
  - id: lancet
    in:
      - id: bed
        source: gatk_intervallisttools/output
      - id: input_normal_bam
        source: samtools_cram2bam_plus_calmd_normal/bam_file
      - id: input_tumor_bam
        source: samtools_cram2bam_plus_calmd_tumor/bam_file
      - id: output_basename
        source: output_basename
      - id: padding
        source: padding
      - id: reference
        source: indexed_reference_fasta
      - id: window
        source: window
    out:
      - id: lancet_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: lancet
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: bed
          type: File
        - id: input_normal_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: input_tumor_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: output_basename
          type: string
        - id: padding
          type: int
          doc: >-
            If WGS (less likely), recommend 25, if exome+, recommend half window
            size
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
        - id: window
          type: int
          doc: >-
            window size for lancet.  default is 600, recommend 500 for WGS, 600
            for exome+
      outputs:
        - id: lancet_vcf
          type: File
          outputBinding:
            glob: '*.vcf'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            set -eo pipefail

            /lancet-1.0.7/lancet --tumor $(inputs.input_tumor_bam.path) --normal
            $(inputs.input_normal_bam.path) --ref $(inputs.reference.path) --bed
            $(inputs.bed.path) --num-threads 6 --window-size $(inputs.window)
            --padding $(inputs.padding) --max-indel-len 50 >
            $(inputs.input_tumor_bam.nameroot).$(inputs.bed.nameroot).vcf ||
            (echo 'active region filter failed, trying without' &&
            /lancet-1.0.7/lancet --tumor $(inputs.input_tumor_bam.path) --normal
            $(inputs.input_normal_bam.path) --ref $(inputs.reference.path) --bed
            $(inputs.bed.path) --num-threads 6 --window-size $(inputs.window)
            --active-region-off --padding $(inputs.padding) --max-indel-len 50 >
            $(inputs.input_tumor_bam.nameroot).$(inputs.bed.nameroot).vcf)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 12000
          coresMin: 6
        - class: DockerRequirement
          dockerPull: 'kfdrc/lancet:1.0.7'
        - class: InlineJavascriptRequirement
    scatter:
      - bed
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.9xlarge;ebs-gp2;500
    'sbg:x': 804.6476440429688
    'sbg:y': 588.5
  - id: samtools_cram2bam_plus_calmd_normal
    in:
      - id: input_reads
        source: input_normal_aligned
      - id: reference
        source: indexed_reference_fasta
      - id: threads
        default: 16
    out:
      - id: bam_file
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: samtools_cram2bam_plus_calmd
      baseCommand:
        - samtools
        - view
      inputs:
        - id: input_reads
          type: File
        - id: reference
          type: File
          secondaryFiles:
            - .fai
        - default: 16
          id: threads
          type: int?
      outputs:
        - id: bam_file
          type: File
          outputBinding:
            glob: '*.bam'
          secondaryFiles:
            - ^.bai
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            -@ $(inputs.threads) -h -T $(inputs.reference.path)
            $(inputs.input_reads.path) | samtools calmd -@ 16 -b --reference
            $(inputs.reference.path) - > $(inputs.input_reads.nameroot).bam &&
            samtools index $(inputs.input_reads.nameroot).bam
            $(inputs.input_reads.nameroot).bai
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 12000
          coresMin: $(inputs.threads)
        - class: DockerRequirement
          dockerPull: 'kfdrc/samtools:1.9'
        - class: InlineJavascriptRequirement
    'sbg:x': 314.328125
    'sbg:y': 407
  - id: samtools_cram2bam_plus_calmd_tumor
    in:
      - id: input_reads
        source: input_tumor_aligned
      - id: reference
        source: indexed_reference_fasta
      - id: threads
        default: 16
    out:
      - id: bam_file
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: samtools_cram2bam_plus_calmd
      baseCommand:
        - samtools
        - view
      inputs:
        - id: input_reads
          type: File
        - id: reference
          type: File
          secondaryFiles:
            - .fai
        - default: 16
          id: threads
          type: int?
      outputs:
        - id: bam_file
          type: File
          outputBinding:
            glob: '*.bam'
          secondaryFiles:
            - ^.bai
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            -@ $(inputs.threads) -h -T $(inputs.reference.path)
            $(inputs.input_reads.path) | samtools calmd -@ 16 -b --reference
            $(inputs.reference.path) - > $(inputs.input_reads.nameroot).bam &&
            samtools index $(inputs.input_reads.nameroot).bam
            $(inputs.input_reads.nameroot).bai
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 12000
          coresMin: $(inputs.threads)
        - class: DockerRequirement
          dockerPull: 'kfdrc/samtools:1.9'
        - class: InlineJavascriptRequirement
    'sbg:x': 314.328125
    'sbg:y': 286
  - id: sort_merge_lancet_vcf
    in:
      - id: input_vcfs
        source:
          - lancet/lancet_vcf
      - id: output_basename
        source: output_basename
      - id: reference_dict
        source: reference_dict
      - id: tool_name
        valueFrom: '${return "lancet"}'
    out:
      - id: merged_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk4_mergevcfs
      baseCommand:
        - /gatk
        - SortVcf
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
            --java-options "-Xmx6g" -O
            $(inputs.output_basename).$(inputs.tool_name).merged.vcf.gz
            --SEQUENCE_DICTIONARY $(inputs.reference_dict.path)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 6000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    label: GATK Sort & Merge lancet
    'sbg:x': 1129.397705078125
    'sbg:y': 467.5
  - id: vep_annot_lancet
    in:
      - id: cache
        source: vep_cache
      - id: input_vcf
        source: gatk_selectvariants_lancet/pass_vcf
      - id: normal_id
        source: input_normal_name
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
      - id: tool_name
        valueFrom: '${return "lancet_somatic"}'
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
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: normal_id
          type: string
        - id: output_basename
          type: string
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
            GRCh38 --ref-fasta $(inputs.reference.path) --tumor-id
            $(inputs.tumor_id) --normal-id $(inputs.normal_id) && mv
            input_file.vep.vcf
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
    'sbg:x': 1659.6119384765625
    'sbg:y': 446.5
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
requirements:
  - class: ScatterFeatureRequirement
'sbg:image_url': >-
  https://cavatica.sbgenomics.com/ns/brood/images/cavatica/openpbta-tcga/kfdrc-lancet-wf-baminput/0.png
'sbg:projectName': openPBTA-TCGA
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': kogantit
    'sbg:modifiedOn': 1578935248
    'sbg:revisionNotes': Copy of kogantit/teja-workspace/kfdrc-lancet-wf/7
'sbg:appVersion':
  - v1.0
'sbg:id': cavatica/openpbta-tcga/kfdrc-lancet-wf-baminput/0
'sbg:revision': 0
'sbg:revisionNotes': Copy of kogantit/teja-workspace/kfdrc-lancet-wf/7
'sbg:modifiedOn': 1578935248
'sbg:modifiedBy': kogantit
'sbg:createdOn': 1578935248
'sbg:createdBy': kogantit
'sbg:project': cavatica/openpbta-tcga
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - kogantit
'sbg:latestRevision': 0
'sbg:publisher': sbg
'sbg:content_hash': a50fcf86276a75bbd1edac6114f9727786da6848aafd32377e137d56b0d798bdf
