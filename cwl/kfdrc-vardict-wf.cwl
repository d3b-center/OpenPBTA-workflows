class: Workflow
cwlVersion: v1.0
id: kfdrc-harmonization/pbta-lancet-vardict-analysis/kfdrc-vardict-wf/0
label: kfdrc-vardict-wf
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: exome_flag
    type: string?
    doc: set to 'Y' for exome mode
    'sbg:x': 0
    'sbg:y': 670
  - id: indexed_reference_fasta
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    'sbg:x': 0
    'sbg:y': 563
  - id: input_normal_aligned
    type: File
    secondaryFiles:
      - .crai
    'sbg:x': 0
    'sbg:y': 456
  - id: input_normal_name
    type: string
    'sbg:x': 314.328125
    'sbg:y': 556
  - id: input_tumor_aligned
    type: File
    secondaryFiles:
      - .crai
    'sbg:x': 0
    'sbg:y': 349
  - id: input_tumor_name
    type: string
    'sbg:x': 314.328125
    'sbg:y': 449
  - id: min_vaf
    type: float?
    doc: Min variant allele frequency for vardict to consider.  Recommend 0.05
    default: 0.05
    'sbg:x': 314.328125
    'sbg:y': 342
  - id: output_basename
    type: string
    'sbg:x': 314.328125
    'sbg:y': 235
  - id: reference_dict
    type: File
    'sbg:x': 533.7611083984375
    'sbg:y': 402.5
  - id: select_vars_mode
    type: string
    doc: 'Choose ''gatk'' for SelectVariants tool, or ''grep'' for grep expression'
    'sbg:x': 0
    'sbg:y': 242
  - id: vep_cache
    type: File
    label: tar gzipped cache from ensembl/local converted cache
    'sbg:x': 0
    'sbg:y': 135
  - id: wgs_calling_interval_list
    type: File
    'sbg:x': 0
    'sbg:y': 28
outputs:
  - id: vardict_prepass_vcf
    outputSource:
      - sort_merge_vardict_vcf/merged_vcf
    type: File
    'sbg:x': 1156.032470703125
    'sbg:y': 288.5
  - id: vardict_vep_somatic_only_maf
    outputSource:
      - vep_annot_vardict/output_maf
    type: File
    'sbg:x': 2137.509033203125
    'sbg:y': 456
  - id: vardict_vep_somatic_only_tbi
    outputSource:
      - vep_annot_vardict/output_tbi
    type: File
    'sbg:x': 2137.509033203125
    'sbg:y': 349
  - id: vardict_vep_somatic_only_vcf
    outputSource:
      - vep_annot_vardict/output_vcf
    type: File
    'sbg:x': 2137.509033203125
    'sbg:y': 242
steps:
  - id: bcbio_filter_fp_somatic
    in:
      - id: input_vcf
        source: sort_merge_vardict_vcf/merged_vcf
      - id: output_basename
        source: output_basename
    out:
      - id: filtered_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: bcbio_vardict_fp_somatic_filter
      baseCommand:
        - python
      inputs:
        - id: input_vcf
          type: File
        - id: output_basename
          type: string
      outputs:
        - id: filtered_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /bcbio_vardict_filter.py $(inputs.input_vcf.path) | grep -E
            "^#|STATUS=StrongSomatic" >
            $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf

            bgzip $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf
            && tabix
            $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf.gz
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: kfdrc/bcbio_vardict_filter
        - class: InlineJavascriptRequirement
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 1156.032470703125
    'sbg:y': 402.5
  - id: gatk_intervallisttools
    in:
      - id: bands
        valueFrom: '${return 1000000}'
      - id: exome_flag
        source: exome_flag
      - id: interval_list
        source: wgs_calling_interval_list
      - id: reference_dict
        source: reference_dict
      - id: scatter_ct
        valueFrom: '${return 200}'
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
    'sbg:x': 314.328125
    'sbg:y': 677
  - id: gatk_selectvariants_vardict
    in:
      - id: input_vcf
        source: bcbio_filter_fp_somatic/filtered_vcf
      - id: mode
        source: select_vars_mode
      - id: output_basename
        source: output_basename
      - id: tool_name
        valueFrom: '${return "vardict"}'
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
    label: GATK Select Vardict PASS
    hints:
      - class: 'sbg:AWSInstanceType'
        value: c5.xlarge;ebs-gp2;250
    'sbg:x': 1411.0745849609375
    'sbg:y': 335
  - id: samtools_normal_cram2bam
    in:
      - id: input_reads
        source: input_normal_aligned
      - id: reference
        source: indexed_reference_fasta
      - id: threads
        valueFrom: '${return 16}'
    out:
      - id: bam_file
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: samtools_cram2bam
      baseCommand:
        - samtools
        - view
      inputs:
        - id: input_reads
          type: File
        - id: reference
          type: File
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
            -@ $(inputs.threads) -bh $(inputs.input_reads.path) -T
            $(inputs.reference.path) > $(inputs.input_reads.nameroot).bam &&
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
    'sbg:y': 121
  - id: samtools_tumor_cram2bam
    in:
      - id: input_reads
        source: input_tumor_aligned
      - id: reference
        source: indexed_reference_fasta
      - id: threads
        valueFrom: '${return 16}'
    out:
      - id: bam_file
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: samtools_cram2bam
      baseCommand:
        - samtools
        - view
      inputs:
        - id: input_reads
          type: File
        - id: reference
          type: File
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
            -@ $(inputs.threads) -bh $(inputs.input_reads.path) -T
            $(inputs.reference.path) > $(inputs.input_reads.nameroot).bam &&
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
    'sbg:y': 0
  - id: sort_merge_vardict_vcf
    in:
      - id: input_vcfs
        source:
          - vardict/vardict_vcf
      - id: output_basename
        source: output_basename
      - id: reference_dict
        source: reference_dict
      - id: tool_name
        valueFrom: '${return "vardict"}'
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
    label: GATK Sort & merge vardict
    'sbg:x': 878.829345703125
    'sbg:y': 335
  - id: vardict
    in:
      - id: bed
        source: gatk_intervallisttools/output
      - id: input_normal_bam
        source: samtools_normal_cram2bam/bam_file
      - id: input_normal_name
        source: input_normal_name
      - id: input_tumor_bam
        source: samtools_tumor_cram2bam/bam_file
      - id: input_tumor_name
        source: input_tumor_name
      - id: min_vaf
        source: min_vaf
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
    out:
      - id: vardict_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: vardictjava
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
        - id: input_normal_name
          type: string
        - id: input_tumor_bam
          type: File
          secondaryFiles:
            - ^.bai
        - id: input_tumor_name
          type: string
        - default: 0.05
          id: min_vaf
          type: float?
          doc: Recommend 0.05
        - id: output_basename
          type: string
        - id: reference
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
      outputs:
        - id: vardict_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            set -eo pipefail; export VAR_DICT_OPTS='"-Xms768m" "-Xmx16g"';
            /VarDict-1.5.8/bin/VarDict -G $(inputs.reference.path) -f
            $(inputs.min_vaf) -th 4 --nosv -N $(inputs.output_basename) -b
            '$(inputs.input_tumor_bam.path)|$(inputs.input_normal_bam.path)' -z
            -c 1 -S 2 -E 3 -g 4 -y -F 0x700 -Q 10 -V 0.01 -Y 100
            $(inputs.bed.path) | /VarDict-1.5.8/bin/testsomatic.R |
            /VarDict-1.5.8/bin/var2vcf_paired.pl -N
            '$(inputs.input_tumor_name)|$(inputs.input_normal_name)' -f
            $(inputs.min_vaf) -M -m 4.25 > $(inputs.output_basename).result.vcf
            && cat $(inputs.output_basename).result.vcf | perl -e 'while(<>){if
            ($_ =~ /^#/){print $_;} else{@a = split /\t/,$_; if($a[3] =~
            /[KMRYSWBVHDXkmryswbvhdx]/){$a[3] = "N";} if($a[4] =~
            /[KMRYSWBVHDXkmryswbvhdx]/){$a[4] = "N";} if($a[3] ne $a[4]){print
            join("\t", @a);}}}' >
            $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf
            && bgzip 
            $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf
            && tabix 
            $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf.gz
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 18
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/vardict:1.5.8'
        - class: InlineJavascriptRequirement
    scatter:
      - bed
    hints:
      - class: 'sbg:AWSInstanceType'
        value: r4.8xlarge;ebs-gp2;500
    'sbg:x': 533.7611083984375
    'sbg:y': 246.5
  - id: vep_annot_vardict
    in:
      - id: cache
        source: vep_cache
      - id: input_vcf
        source: gatk_selectvariants_vardict/pass_vcf
      - id: normal_id
        source: input_normal_name
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
      - id: tool_name
        valueFrom: '${return "vardict_somatic"}'
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
    'sbg:x': 1664.085693359375
    'sbg:y': 314
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
requirements:
  - class: ScatterFeatureRequirement
'sbg:appVersion':
  - v1.0
'sbg:id': kfdrc-harmonization/pbta-lancet-vardict-analysis/kfdrc-vardict-wf/0
'sbg:revision': 0
'sbg:revisionNotes': Copy of zhangb1/kf-somatic-tools-test/kfdrc-vardict-wf/4
'sbg:modifiedOn': 1567782742
'sbg:modifiedBy': ennisb
'sbg:createdOn': 1567782742
'sbg:createdBy': ennisb
'sbg:project': kfdrc-harmonization/pbta-lancet-vardict-analysis
'sbg:projectName': PBTA-Lancet-Vardict-analysis
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - ennisb
'sbg:latestRevision': 0
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': ennisb
    'sbg:modifiedOn': 1567782742
    'sbg:revisionNotes': Copy of zhangb1/kf-somatic-tools-test/kfdrc-vardict-wf/4
'sbg:image_url': >-
  https://cavatica.sbgenomics.com/ns/brood/images/kfdrc-harmonization/pbta-lancet-vardict-analysis/kfdrc-vardict-wf/0.png
'sbg:publisher': sbg
'sbg:content_hash': a871f19cd5f043d50feca6928c08d53fe4a8cc0ac045d4996d4eedf114d0d80b0
'sbg:copyOf': zhangb1/kf-somatic-tools-test/kfdrc-vardict-wf/4
'sbg:update': zhangb1/kf-somatic-tools-test/kfdrc-vardict-wf/32
'sbg:updateRevisionNotes': test new splitting method
'sbg:updateModifiedBy': brownm28
