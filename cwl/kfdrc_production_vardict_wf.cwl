cwlVersion: v1.0
class: Workflow
label: kfdrc_production_vardict_wf
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement

inputs:
- id: reference_fasta
  type: File
- id: reference_fai
  type:
  - 'null'
  - File
- id: reference_dict
  type:
  - 'null'
  - File
- id: input_tumor_aligned
  doc: tumor BAM or CRAM
  type: File
  secondaryFiles: |
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
- id: input_normal_aligned
  doc: normal BAM or CRAM
  type: File
  secondaryFiles: |
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
- id: vep_cache
  doc: tar gzipped cache from ensembl/local converted cache
  type: File
- id: output_basename
  doc: String value to use as basename for outputs
  type: string
- id: wgs_or_wxs
  doc: Select if this run is WGS or WXS
  type:
    name: sex
    type: enum
    symbols:
    - WGS
    - WXS
- id: select_vars_mode
  doc: Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression
  type:
  - 'null'
  - name: select_vars_mode
    type: enum
    symbols:
    - gatk
    - grep
  default: gatk
- id: vardict_cpus
  doc: Number of CPUs for Vardict to use
  type:
  - 'null'
  - int
  default: 9
- id: vardict_min_vaf
  doc: Min variant allele frequency for vardict to consider. Recommend 0.05
  type:
  - 'null'
  - float
  default: 0.05
- id: vardict_ram
  doc: GB of RAM to allocate to Vardict
  type:
  - 'null'
  - int
  default: 18
- id: vep_ref_build
  doc: Genome ref build used, should line up with cache
  type:
  - 'null'
  - string
  default: GRCh38
- id: exome_flag
  doc: Whether to run in exome mode for callers. Y for WXS, N for WGS
  type:
  - 'null'
  - string
- id: vardict_padding
  doc: |-
    Padding to add to input intervals, recommend 0 if intervals already padded such as in WXS, 150 if not such as in WGS
  type:
  - 'null'
  - int
- id: wgs_calling_interval_list
  doc: |-
    GATK intervals list-style, or bed file.  Recommend canocical chromosomes with N regions removed
  type:
  - 'null'
  - File
- id: padded_capture_regions
  doc: Recommend 100bp pad, for somatic variant
  type:
  - 'null'
  - File

outputs:
- id: vardict_vep_somatic_only_vcf
  type: File
  outputSource: run_vardict/vardict_vep_somatic_only_vcf
- id: vardict_vep_somatic_only_tbi
  type: File
  outputSource: run_vardict/vardict_vep_somatic_only_tbi
- id: vardict_vep_somatic_only_maf
  type: File
  outputSource: run_vardict/vardict_vep_somatic_only_maf
- id: vardict_prepass_vcf
  type: File
  outputSource: run_vardict/vardict_prepass_vcf

steps:
- id: choose_defaults
  in:
  - id: input_mode
    source: wgs_or_wxs
  - id: exome_flag
    source: exome_flag
  - id: vardict_padding
    source: vardict_padding
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: Selects the appropriate defaults based on given the mode

    requirements:
    - class: InlineJavascriptRequirement

    inputs:
    - id: input_mode
      doc: Select if this run is WGS or WXS
      type:
        name: input_mode
        type: enum
        symbols:
        - WGS
        - WXS
    - id: exome_flag
      doc: Whether to run in exome mode for callers.
      type:
      - 'null'
      - string
    - id: cnvkit_wgs_mode
      doc: |-
        Entering Y will run cnvkit in WGS mode, otherwise it will run in hybrid mode. Defaults to Y in wgs mode.
      type:
      - 'null'
      - string
    - id: i_flag
      doc: |-
        Flag to intersect germline calls on padded regions.  Use N if you want to skip this. Defaults to N in WGS mode.
      type:
      - 'null'
      - string
    - id: lancet_padding
      doc: |-
        Recommend 0 if interval file padded already, half window size if not. Recommended: 0 for WXS; 300 for WGS
      type:
      - 'null'
      - int
    - id: lancet_window
      doc: window size for lancet.  default is 600, recommend 500 for WGS, 600 for
        exome+
      type:
      - 'null'
      - int
    - id: vardict_padding
      doc: |-
        Padding to add to input intervals, recommend 0 if intervals already padded such as in WXS, 150 if not such as in WGS
      type:
      - 'null'
      - int

    outputs:
    - id: out_exome_flag
      type:
      - 'null'
      - string
      outputBinding:
        outputEval: |-
          ${
            if (inputs.exome_flag) { return inputs.exome_flag }
            else if (inputs.input_mode == 'WGS') { return 'N' }
            else if (inputs.input_mode == 'WXS') { return 'Y' }
          }
    - id: out_cnvkit_wgs_mode
      type:
      - 'null'
      - string
      outputBinding:
        outputEval: |-
          ${
            if (inputs.cnvkit_wgs_mode) { return inputs.cnvkit_wgs_mode }
            else if (inputs.input_mode == 'WGS') { return 'Y' }
            else if (inputs.input_mode == 'WXS') { return 'N' }
          }
    - id: out_i_flag
      type:
      - 'null'
      - string
      outputBinding:
        outputEval: |-
          ${
            if (inputs.i_flag) { return inputs.i_flag }
            else if (inputs.input_mode == 'WGS') { return 'N' }
            else if (inputs.input_mode == 'WXS') { return null }
          }
    - id: out_lancet_padding
      type:
      - 'null'
      - int
      outputBinding:
        outputEval: |-
          ${
            if (inputs.lancet_padding) { return inputs.lancet_padding }
            else if (inputs.input_mode == 'WGS') { return 300 }
            else if (inputs.input_mode == 'WXS') { return 0 }
          }
    - id: out_lancet_window
      type:
      - 'null'
      - int
      outputBinding:
        outputEval: |-
          ${
            if (inputs.lancet_window) { return inputs.lancet_window }
            else if (inputs.input_mode == 'WGS') { return 600 }
            else if (inputs.input_mode == 'WXS') { return 600 }
          }
    - id: out_vardict_padding
      type:
      - 'null'
      - int
      outputBinding:
        outputEval: |-
          ${
            if (inputs.vardict_padding) { return inputs.vardict_padding }
            else if (inputs.input_mode == 'WGS') { return 150 }
            else if (inputs.input_mode == 'WXS') { return 0 }
          }

    baseCommand: echo
    arguments:
    - position: 0
      valueFrom: Selecting $(inputs.input_mode) default values
    id: mode_defaults
  out:
  - out_exome_flag
  - out_cnvkit_wgs_mode
  - out_i_flag
  - out_lancet_padding
  - out_lancet_window
  - out_vardict_padding
- id: prepare_reference
  in:
  - id: input_fasta
    source: reference_fasta
  - id: input_fai
    source: reference_fai
  - id: input_dict
    source: reference_dict
  run:
    cwlVersion: v1.0
    class: Workflow
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement

    inputs:
    - id: input_fasta
      type: File
    - id: input_fai
      type:
      - 'null'
      - File
    - id: input_dict
      type:
      - 'null'
      - File
    - id: input_alt
      type:
      - 'null'
      - File
    - id: input_amb
      type:
      - 'null'
      - File
    - id: input_ann
      type:
      - 'null'
      - File
    - id: input_bwt
      type:
      - 'null'
      - File
    - id: input_pac
      type:
      - 'null'
      - File
    - id: input_sa
      type:
      - 'null'
      - File
    - id: generate_bwa_indexes
      type:
      - 'null'
      - boolean

    outputs:
    - id: indexed_fasta
      doc: Reference fasta with all available indexes as secondaryFiles
      type: File
      outputSource: bundle_secondaries/output
    - id: reference_dict
      doc: Standalone reference dict
      type: File
      outputSource: picard_create_sequence_dictionary/dict

    steps:
    - id: samtools_faidx
      in:
      - id: input_fasta
        source: input_fasta
      - id: input_index
        source: input_fai
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        doc: |-
          This tool takes an input fasta and optionally a input index for the input fasta.
          If the index is not provided this tool will generate one.
          Finally the tool will return the input reference file with the index (generated or provided) as a secondaryFile.

        requirements:
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
        - class: InitialWorkDirRequirement
          listing:
          - $(inputs.input_fasta)
          - $(inputs.input_index)
        - class: InlineJavascriptRequirement
        - class: ResourceRequirement
        - class: ShellCommandRequirement

        inputs:
        - id: input_fasta
          doc: Input fasta file
          type: File
          inputBinding:
            position: 1
        - id: input_index
          doc: Input fasta index
          type:
          - 'null'
          - File

        outputs:
        - id: fai
          type: File
          outputBinding:
            glob: '*.fai'

        baseCommand: []
        arguments:
        - position: 0
          valueFrom: "$(inputs.input_index ? 'echo samtools faidx' : 'samtools faidx'\
            \ )"
          shellQuote: false
        id: samtools_faidx
      out:
      - fai
    - id: picard_create_sequence_dictionary
      in:
      - id: input_fasta
        source: input_fasta
      - id: input_dict
        source: input_dict
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        doc: |-
          This tool conditionally creats a sequence dictionary from an input fasta using Picard CreateSequenceDictionary.
          The tool will only generate the index if an the input_dict is not passed.
          The tool returnts the dict as its only output.

        requirements:
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.7.0R
        - class: InitialWorkDirRequirement
          listing:
          - $(inputs.input_fasta)
          - $(inputs.input_dict)
        - class: InlineJavascriptRequirement
        - class: ResourceRequirement
        - class: ShellCommandRequirement

        inputs:
        - id: input_fasta
          type: File
          inputBinding:
            prefix: -R
            position: 2
        - id: input_dict
          type:
          - 'null'
          - File

        outputs:
        - id: dict
          type: File
          outputBinding:
            glob: '*.dict'

        baseCommand: []
        arguments:
        - position: 0
          valueFrom: |-
            $(inputs.input_dict ? 'echo java -jar /gatk-package-4.1.7.0-local.jar' : 'java -jar /gatk-package-4.1.7.0-local.jar' )
          shellQuote: false
        - position: 1
          valueFrom: CreateSequenceDictionary
          shellQuote: false
        id: picard_createsequencedictionary
      out:
      - dict
    - id: bwa_index
      in:
      - id: generate_bwa_indexes
        source: generate_bwa_indexes
      - id: input_fasta
        source: input_fasta
      - id: input_alt
        source: input_alt
      - id: input_amb
        source: input_amb
      - id: input_ann
        source: input_ann
      - id: input_bwt
        source: input_bwt
      - id: input_pac
        source: input_pac
      - id: input_sa
        source: input_sa
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        doc: |-
          This tool conditionally generates the bwa 64 indexes for an input fasta file using bwa index.
          The tool will generate the indexes only of generate_bwa_indexes is set to true AND any of the alt,
          amb, ann, bwt, pac, or sa files is missing.
          The tool returns the six indexes as its output.

        requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
        - class: InitialWorkDirRequirement
          listing:
          - $(inputs.input_fasta)
          - $(inputs.input_alt)
          - $(inputs.input_amb)
          - $(inputs.input_ann)
          - $(inputs.input_bwt)
          - $(inputs.input_pac)
          - $(inputs.input_sa)
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/bwa:0.7.17-dev
        - class: InlineJavascriptRequirement

        inputs:
        - id: generate_bwa_indexes
          type:
          - 'null'
          - boolean
        - id: input_fasta
          type: File
          inputBinding:
            position: 2
            valueFrom: $(self.basename)
        - id: input_alt
          type:
          - 'null'
          - File
        - id: input_amb
          type:
          - 'null'
          - File
        - id: input_ann
          type:
          - 'null'
          - File
        - id: input_bwt
          type:
          - 'null'
          - File
        - id: input_pac
          type:
          - 'null'
          - File
        - id: input_sa
          type:
          - 'null'
          - File

        outputs:
        - id: alt
          type:
          - 'null'
          - File
          outputBinding:
            glob: '*.64.alt'
        - id: amb
          type:
          - 'null'
          - File
          outputBinding:
            glob: '*.64.amb'
        - id: ann
          type:
          - 'null'
          - File
          outputBinding:
            glob: '*.64.ann'
        - id: bwt
          type:
          - 'null'
          - File
          outputBinding:
            glob: '*.64.bwt'
        - id: pac
          type:
          - 'null'
          - File
          outputBinding:
            glob: '*.64.pac'
        - id: sa
          type:
          - 'null'
          - File
          outputBinding:
            glob: '*.64.sa'

        baseCommand: []
        arguments:
        - position: 0
          valueFrom: |-
            $(inputs.input_alt && inputs.input_amb && inputs.input_ann && inputs.input_bwt && inputs.input_pac && inputs.input_sa ? 'echo bwa' : inputs.generate_bwa_indexes ? 'bwa' : 'echo bwa')
          shellQuote: false
        - position: 1
          valueFrom: 'index -6 -a bwtsw '
          shellQuote: false
        id: bwa_index
      out:
      - alt
      - amb
      - ann
      - bwt
      - pac
      - sa
    - id: bundle_secondaries
      in:
      - id: primary_file
        source: input_fasta
      - id: secondary_files
        source:
        - samtools_faidx/fai
        - picard_create_sequence_dictionary/dict
        - bwa_index/alt
        - bwa_index/amb
        - bwa_index/ann
        - bwa_index/bwt
        - bwa_index/pac
        - bwa_index/sa
        linkMerge: merge_flattened
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        doc: |-
          This tool takes a primary file and list of secondary files as input and passes the primary_file as
          the output with the secondary_files as secondaryFiles.

        requirements:
        - class: InlineJavascriptRequirement
        - class: ResourceRequirement
        - class: ShellCommandRequirement
        - class: InitialWorkDirRequirement
          listing:
          - $(inputs.primary_file)
          - $(inputs.secondary_files)

        inputs:
        - id: primary_file
          doc: Primary File
          type: File
        - id: secondary_files
          doc: List of secondary files
          type:
            type: array
            items: File

        outputs:
        - id: output
          type: File
          secondaryFiles: |-
            ${var arr = []; for (i = 0; i < inputs.secondary_files.length; i++) { if (inputs.secondary_files[i]) { arr.push(inputs.secondary_files[i].basename) } }; return arr}
          outputBinding:
            glob: $(inputs.primary_file.basename)

        baseCommand:
        - echo
        id: bundle_secondaryfiles
      out:
      - output

    hints:
    - class: sbg:maxNumberOfParallelInstances
      value: 2
    id: kfdrc_prepare_reference
  out:
  - indexed_fasta
  - reference_dict
- id: select_interval_list
  in:
  - id: input_mode
    source: wgs_or_wxs
  - id: wgs_input
    source: wgs_calling_interval_list
  - id: wxs_input
    source: padded_capture_regions
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: Selects the appropriate input to serve as the output given the mode

    requirements:
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/ubuntu:18.04
    - class: InlineJavascriptRequirement

    inputs:
    - id: input_mode
      doc: Select if this run is WGS or WXS
      type:
        name: input_mode
        type: enum
        symbols:
        - WGS
        - WXS
    - id: wgs_input
      doc: Input that should be passed when mode is WGS
      type:
      - 'null'
      - Any
    - id: wxs_input
      doc: Input that should be passed when mode is WXS
      type:
      - 'null'
      - Any

    outputs:
    - id: output
      type: Any
      outputBinding:
        outputEval: |-
          ${
            if (inputs.input_mode == 'WGS') { return inputs.wgs_input }
            else if (inputs.input_mode == 'WXS') { return inputs.wxs_input }
          }

    baseCommand:
    - /bin/bash
    - -c
    arguments:
    - position: 0
      valueFrom: |-
        set -eo pipefail
        ${
            var cmd = " >&2 echo Choosing inputs based on mode;";
            if (inputs.input_mode == 'WGS' && inputs.wgs_input == null){
              return "echo WGS run requires wgs_input >&2 && exit 1;"
            }
            else if (inputs.input_mode == 'WXS' && inputs.wxs_input == null){
              return "echo WXS run requires wxs_input >&2 && exit 1;"
            }
            return cmd;
        }
         
    id: mode_selector
  out:
  - output
- id: python_vardict_interval_split
  doc: |-
    Custom interval list generation for vardict input. Briefly, ~60M bp per interval list, 20K bp intervals, lists break on chr and N reginos only
  in:
  - id: wgs_bed_file
    source: select_interval_list/output
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    label: Create intervals for VarDict
    doc: |-
      This tool takes in an interval list with the WGS coords split by N regions. It splits the bed files into bed files with a max total number of base bairs, unless the regions is already larger, it stays in it's own file. Then within the split lists, intervals are split the specified chunks for easier processing by vardict.  This method prevents FP calls caused by regions with valid ACGT bases from being split between interval lists.  For example, for hg38 canonical chromosomes, using bp_target=60000000 and intvl_target_size=20000 will yield about 55 bed files, each with about 60M bp worth of coverage (unless the interval was already larger, it will be in it's own list), split into 20kb chunks.

    requirements:
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/python:2.7.13
    - class: InlineJavascriptRequirement

    inputs:
    - id: wgs_bed_file
      doc: |-
        Should be a bed file of WGS regions with N's removed.  GATK calling regions is a good source.
      type: File
    - id: bp_target
      doc: |-
        Intended max number of base pairs per file.  Existing intervals large than this will NOT be split into another file.
      type:
      - 'null'
      - int
      default: 60000000
    - id: intvl_target_size
      doc: For each file, split each interval into chuck of this size
      type:
      - 'null'
      - int
      default: 20000

    outputs:
    - id: split_intervals_bed
      type:
        type: array
        items: File
      outputBinding:
        glob: '*.bed'

    baseCommand:
    - python
    - -c
    arguments:
    - position: 0
      valueFrom: |-
        def main():
            import sys
            bp_target = $(inputs.bp_target)
            intvl_target_size = $(inputs.intvl_target_size)
            bed_file = open("$(inputs.wgs_bed_file.path)")

            i=0
            intvl_set = {}
            cur_size = 0
            for cur_intvl in bed_file:
                f = 0
                if i not in intvl_set:
                    intvl_set[i] = []
                # chr(9) is ASCII code for tab; chr(10) is ASCII code for newline
                data = cur_intvl.rstrip(chr(10)).split(chr(9))
                (chrom, start, end) = (data[0], data[1], data[2])
                intvl_size = int(end) - int(start)
                if intvl_size >= bp_target:
                    if len(intvl_set[i]) != 0:
                        i += 1
                        intvl_set[i] = []
                        f = 1
                elif cur_size + intvl_size > bp_target:
                    if len(intvl_set[i]) != 0:
                        i += 1
                        intvl_set[i] = []
                        cur_size = intvl_size
                else:
                    cur_size += intvl_size
                intvl_set[i].append([chrom, start, end])
                if f == 1:
                    i += 1
                    cur_size = 0
            bed_file.close()

            for set_i, invtl_list in sorted(intvl_set.items()):
                set_size = 0
                out = open("set_" + str(set_i) + ".bed", "w")
                for intervals in invtl_list:
                    (chrom, start, end) = (intervals[0], intervals[1], intervals[2])
                    intvl_size = int(end) - int(start)
                    set_size += intvl_size
                    for j in range(int(start), int(end), intvl_target_size):
                        new_end = j + intvl_target_size
                        if new_end > int(end):
                            new_end = end
                        out.write(chrom + chr(9) + str(j) + chr(9) + str(new_end) + chr(10))
                sys.stderr.write("Set " + str(set_i) + " size:" + chr(9) + str(set_size) + chr(10))
                out.close()

        if __name__ == "__main__":
            main()
      shellQuote: true
    id: python_vardict_interval_split
  out:
  - split_intervals_bed
- id: gatk_intervallisttools
  in:
  - id: interval_list
    source: select_interval_list/output
  - id: reference_dict
    source: prepare_reference/reference_dict
  - id: exome_flag
    source: choose_defaults/out_exome_flag
  - id: scatter_ct
    valueFrom: ${return 50}
  - id: bands
    valueFrom: ${return 80000000}
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
    - class: ResourceRequirement
      ramMin: 2000

    inputs:
    - id: interval_list
      type:
      - 'null'
      - File
    - id: bands
      type: int
    - id: scatter_ct
      type: int
    - id: reference_dict
      doc: Provide only if input is bed file instead of gatk style .interval_list
      type:
      - 'null'
      - File
      sbg:suggestedValue:
        name: Provide only if input is bed file instead of gatk style .interval_list
    - id: exome_flag
      doc: If 'Y', will set bands to 0 to prevent breaking up of intervals
      type:
      - 'null'
      - string
      default: N
    - id: break_by_chr
      doc: |-
        If Y, break up files by chr.  If creating smaler intervals, recommend scatter_ct=1
      type:
      - 'null'
      - string
      default: N

    outputs:
    - id: output
      type:
      - 'null'
      - type: array
        items: File
      outputBinding:
        glob: '*.bed'

    baseCommand:
    - /bin/bash
    - -c
    arguments:
    - position: 1
      valueFrom: |-
        set -eo pipefail
        ${
          if (inputs.interval_list == null) {
            return "echo No interval list exiting without input >&2 && exit 0;";
          }
          else {
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
      shellQuote: false
    id: gatk4_intervallist2bed
  out:
  - output
- id: samtools_cram2bam_plus_calmd_tumor
  in:
  - id: input_reads
    source: input_tumor_aligned
  - id: threads
    valueFrom: ${return 16;}
  - id: reference
    source: prepare_reference/indexed_fasta
  run:
    cwlVersion: v1.0
    class: CommandLineTool

    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: $(inputs.threads)
      ramMin: 12000

    inputs:
    - id: input_reads
      type: File
    - id: threads
      type:
      - 'null'
      - int
      default: 16
    - id: reference
      type: File
      secondaryFiles:
      - .fai

    outputs:
    - id: bam_file
      type: File
      secondaryFiles:
      - ^.bai
      outputBinding:
        glob: '*.bam'

    baseCommand:
    - /bin/bash -c
    arguments:
    - position: 1
      valueFrom: |-
        set -eo pipefail
        ${
          var bam_name = inputs.input_reads.nameroot + ".bam";
          var cmd = "samtools view -@ " + inputs.threads + " -h -T " + inputs.reference.path + " " + inputs.input_reads.path
          + " | samtools calmd -@ " + inputs.threads + " -b --reference " + inputs.reference.path + " - > " + bam_name + ";";
          if(inputs.input_reads.basename == bam_name){
            cmd = ">&2 echo input reads already have bam extension, indexing and passing through; cp " + inputs.input_reads.path
            + " " + bam_name + ";"
          }
          cmd += "samtools index -@ " + inputs.threads + " " + bam_name + " " + inputs.input_reads.nameroot + ".bai;"
          return cmd;
        }
      shellQuote: false
    id: samtools_cram2bam_plus_calmd
  out:
  - bam_file
- id: samtools_cram2bam_plus_calmd_normal
  in:
  - id: input_reads
    source: input_normal_aligned
  - id: threads
    valueFrom: ${return 16;}
  - id: reference
    source: prepare_reference/indexed_fasta
  run:
    cwlVersion: v1.0
    class: CommandLineTool

    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/samtools:1.9
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: $(inputs.threads)
      ramMin: 12000

    inputs:
    - id: input_reads
      type: File
    - id: threads
      type:
      - 'null'
      - int
      default: 16
    - id: reference
      type: File
      secondaryFiles:
      - .fai

    outputs:
    - id: bam_file
      type: File
      secondaryFiles:
      - ^.bai
      outputBinding:
        glob: '*.bam'

    baseCommand:
    - /bin/bash -c
    arguments:
    - position: 1
      valueFrom: |-
        set -eo pipefail
        ${
          var bam_name = inputs.input_reads.nameroot + ".bam";
          var cmd = "samtools view -@ " + inputs.threads + " -h -T " + inputs.reference.path + " " + inputs.input_reads.path
          + " | samtools calmd -@ " + inputs.threads + " -b --reference " + inputs.reference.path + " - > " + bam_name + ";";
          if(inputs.input_reads.basename == bam_name){
            cmd = ">&2 echo input reads already have bam extension, indexing and passing through; cp " + inputs.input_reads.path
            + " " + bam_name + ";"
          }
          cmd += "samtools index -@ " + inputs.threads + " " + bam_name + " " + inputs.input_reads.nameroot + ".bai;"
          return cmd;
        }
      shellQuote: false
    id: samtools_cram2bam_plus_calmd
  out:
  - bam_file
- id: select_vardict_bed_interval
  in:
  - id: input_mode
    source: wgs_or_wxs
  - id: wgs_input
    source: python_vardict_interval_split/split_intervals_bed
  - id: wxs_input
    source: gatk_intervallisttools/output
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: Selects the appropriate input to serve as the output given the mode

    requirements:
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/ubuntu:18.04
    - class: InlineJavascriptRequirement

    inputs:
    - id: input_mode
      doc: Select if this run is WGS or WXS
      type:
        name: input_mode
        type: enum
        symbols:
        - WGS
        - WXS
    - id: wgs_input
      doc: Input that should be passed when mode is WGS
      type:
      - 'null'
      - Any
    - id: wxs_input
      doc: Input that should be passed when mode is WXS
      type:
      - 'null'
      - Any

    outputs:
    - id: output
      type: Any
      outputBinding:
        outputEval: |-
          ${
            if (inputs.input_mode == 'WGS') { return inputs.wgs_input }
            else if (inputs.input_mode == 'WXS') { return inputs.wxs_input }
          }

    baseCommand:
    - /bin/bash
    - -c
    arguments:
    - position: 0
      valueFrom: |-
        set -eo pipefail
        ${
            var cmd = " >&2 echo Choosing inputs based on mode;";
            if (inputs.input_mode == 'WGS' && inputs.wgs_input == null){
              return "echo WGS run requires wgs_input >&2 && exit 1;"
            }
            else if (inputs.input_mode == 'WXS' && inputs.wxs_input == null){
              return "echo WXS run requires wxs_input >&2 && exit 1;"
            }
            return cmd;
        }
         
    id: mode_selector
  out:
  - output
- id: run_vardict
  in:
  - id: indexed_reference_fasta
    source: prepare_reference/indexed_fasta
  - id: input_tumor_aligned
    source: samtools_cram2bam_plus_calmd_tumor/bam_file
  - id: input_tumor_name
    source: input_tumor_name
  - id: input_normal_aligned
    source: samtools_cram2bam_plus_calmd_normal/bam_file
  - id: input_normal_name
    source: input_normal_name
  - id: output_basename
    source: output_basename
  - id: reference_dict
    source: prepare_reference/reference_dict
  - id: bed_invtl_split
    source: select_vardict_bed_interval/output
  - id: padding
    source: choose_defaults/out_vardict_padding
  - id: min_vaf
    source: vardict_min_vaf
  - id: select_vars_mode
    source: select_vars_mode
  - id: cpus
    source: vardict_cpus
  - id: ram
    source: vardict_ram
  - id: vep_cache
    source: vep_cache
  - id: vep_ref_build
    source: vep_ref_build
  run:
    cwlVersion: v1.0
    class: Workflow
    $namespaces:
      sbg: https://sevenbridges.com

    requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement

    inputs:
    - id: indexed_reference_fasta
      type: File
      secondaryFiles:
      - .fai
      - ^.dict
    - id: input_tumor_aligned
      type: File
      secondaryFiles:
      - ^.bai
    - id: input_tumor_name
      type: string
    - id: input_normal_aligned
      type: File
      secondaryFiles:
      - ^.bai
    - id: input_normal_name
      type: string
    - id: output_basename
      type: string
    - id: reference_dict
      type: File
    - id: padding
      doc: |-
        Padding to add to input intervals, recommened 0 if intervals already padded, 150 if not
      type:
      - 'null'
      - int
      default: 150
    - id: bed_invtl_split
      doc: Bed file intervals passed on from and outside pre-processing step
      type:
        type: array
        items: File
    - id: min_vaf
      doc: Min variant allele frequency for vardict to consider.  Recommend 0.05
      type:
      - 'null'
      - float
      default: 0.05
    - id: select_vars_mode
      doc: Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression
      type:
      - 'null'
      - name: select_vars_mode
        type: enum
        symbols:
        - gatk
        - grep
      default: gatk
    - id: cpus
      type:
      - 'null'
      - int
      default: 9
    - id: ram
      doc: In GB
      type:
      - 'null'
      - int
      default: 18
    - id: vep_cache
      label: tar gzipped cache from ensembl/local converted cache
      type: File
    - id: vep_ref_build
      doc: Genome ref build used, should line up with cache.
      type:
      - 'null'
      - string
      default: GRCh38

    outputs:
    - id: vardict_vep_somatic_only_vcf
      type: File
      outputSource: vep_annot_vardict/output_vcf
    - id: vardict_vep_somatic_only_tbi
      type: File
      outputSource: vep_annot_vardict/output_tbi
    - id: vardict_vep_somatic_only_maf
      type: File
      outputSource: vep_annot_vardict/output_maf
    - id: vardict_prepass_vcf
      type: File
      outputSource: sort_merge_vardict_vcf/merged_vcf

    steps:
    - id: vardict
      in:
      - id: input_tumor_bam
        source: input_tumor_aligned
      - id: input_tumor_name
        source: input_tumor_name
      - id: input_normal_bam
        source: input_normal_aligned
      - id: input_normal_name
        source: input_normal_name
      - id: padding
        source: padding
      - id: min_vaf
        source: min_vaf
      - id: cpus
        source: cpus
      - id: ram
        source: ram
      - id: reference
        source: indexed_reference_fasta
      - id: bed
        source: bed_invtl_split
      - id: output_basename
        source: output_basename
      scatter:
      - bed
      run:
        cwlVersion: v1.0
        class: CommandLineTool

        requirements:
        - class: ShellCommandRequirement
        - class: InlineJavascriptRequirement
        - class: ResourceRequirement
          coresMin: $(inputs.cpus)
          ramMin: ${return inputs.ram * 1000}
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/vardict:1.7.0

        inputs:
        - id: reference
          type: File
          secondaryFiles:
          - ^.dict
          - .fai
        - id: input_tumor_bam
          type: File
          secondaryFiles:
          - ^.bai
        - id: input_tumor_name
          type: string
        - id: input_normal_bam
          type: File
          secondaryFiles:
          - ^.bai
        - id: input_normal_name
          type: string
        - id: cpus
          type:
          - 'null'
          - int
          default: 9
        - id: ram
          doc: In GB
          type:
          - 'null'
          - int
          default: 18
        - id: padding
          doc: |-
            Padding to add to input intervals, recommened 0 if intervals already padded, 150 if not
          type:
          - 'null'
          - int
          default: 150
        - id: min_vaf
          doc: Recommend 0.05
          type:
          - 'null'
          - float
          default: 0.05
        - id: output_basename
          type: string
        - id: bed
          type: File

        outputs:
        - id: vardict_vcf
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
            set -eo pipefail; ${
              var ram = Math.floor(inputs.ram/1.074 - 1);
              var exp_cmd = "export VAR_DICT_OPTS='\"-Xms768m\" \"-Xmx" + ram + "g\"';";
              return exp_cmd;
            } /VarDict-1.7.0/bin/VarDict -G $(inputs.reference.path) -f $(inputs.min_vaf) -th $(inputs.cpus) --nosv -N $(inputs.output_basename) -b '$(inputs.input_tumor_bam.path)|$(inputs.input_normal_bam.path)' -z -c 1 -S 2 -E 3 -g 4 -F 0x700 -Q 10 -V 0.01 -x $(inputs.padding) $(inputs.bed.path) > vardict_results.txt && cat vardict_results.txt | /VarDict-1.7.0/bin/testsomatic.R > vardict_r_test_results.txt && cat vardict_r_test_results.txt | /VarDict-1.7.0/bin/var2vcf_paired.pl -N '$(inputs.input_tumor_name)|$(inputs.input_normal_name)' -f $(inputs.min_vaf) -M -m 4.25 > $(inputs.output_basename).result.vcf && cat $(inputs.output_basename).result.vcf | perl -e 'while(<>){if ($_ =~ /^#/){print $_;} else{@a = split /\t/,$_; if($a[3] =~ /[KMRYSWBVHDXkmryswbvhdx]/){$a[3] = "N";} if($a[4] =~ /[KMRYSWBVHDXkmryswbvhdx]/){$a[4] = "N";} if($a[3] ne $a[4]){print join("\t", @a);}}}' > $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf && bgzip  $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf && tabix  $(inputs.output_basename).canonical_base_only.$(inputs.bed.nameroot).vcf.gz
          shellQuote: false
        id: vardictjava
      out:
      - vardict_vcf
      hints:
      - class: sbg:AWSInstanceType
        value: c5.9xlarge
    - id: sort_merge_vardict_vcf
      label: GATK Sort & merge vardict
      in:
      - id: input_vcfs
        source: vardict/vardict_vcf
      - id: output_basename
        source: output_basename
      - id: reference_dict
        source: reference_dict
      - id: tool_name
        valueFrom: ${return "vardict"}
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        label: GATK Merge VCF
        doc: Merge input vcfs

        requirements:
        - class: ShellCommandRequirement
        - class: InlineJavascriptRequirement
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
        - class: ResourceRequirement
          coresMin: 4
          ramMin: 6000

        inputs:
        - id: input_vcfs
          type:
            type: array
            inputBinding:
              prefix: -I
            items: File
          secondaryFiles:
          - .tbi
          inputBinding:
            position: 1
        - id: reference_dict
          type: File
        - id: tool_name
          type: string
        - id: output_basename
          type: string

        outputs:
        - id: merged_vcf
          type: File
          secondaryFiles:
          - .tbi
          outputBinding:
            glob: '*.merged.vcf.gz'

        baseCommand:
        - /gatk
        - SortVcf
        arguments:
        - position: 0
          valueFrom: |-
            --java-options "-Xmx6g" -O $(inputs.output_basename).$(inputs.tool_name).merged.vcf --SEQUENCE_DICTIONARY $(inputs.reference_dict.path) --CREATE_INDEX false
          shellQuote: false
        - position: 2
          valueFrom: |2-

            && cat $(inputs.output_basename).$(inputs.tool_name).merged.vcf | uniq | bgzip > $(inputs.output_basename).$(inputs.tool_name).merged.vcf.gz && tabix $(inputs.output_basename).$(inputs.tool_name).merged.vcf.gz
          shellQuote: false
        id: gatk4_mergevcfs
      out:
      - merged_vcf
    - id: bcbio_filter_fp_somatic
      in:
      - id: input_vcf
        source: sort_merge_vardict_vcf/merged_vcf
      - id: output_basename
        source: output_basename
      run:
        cwlVersion: v1.0
        class: CommandLineTool

        requirements:
        - class: ShellCommandRequirement
        - class: InlineJavascriptRequirement
        - class: ResourceRequirement
          coresMin: 4
          ramMin: 8000
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/bcbio_vardict_filter

        inputs:
        - id: input_vcf
          type: File
        - id: output_basename
          type: string

        outputs:
        - id: filtered_vcf
          type: File
          secondaryFiles:
          - .tbi
          outputBinding:
            glob: '*.vcf.gz'

        baseCommand:
        - python
        arguments:
        - position: 0
          valueFrom: |-
            /bcbio_vardict_filter.py $(inputs.input_vcf.path) | grep -E "^#|STATUS=StrongSomatic" > $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf
            bgzip $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf && tabix $(inputs.output_basename).bcbio_vardict_fp_somatic_filter.vcf.gz
          shellQuote: false
        id: bcbio_vardict_fp_somatic_filter
      out:
      - filtered_vcf
    - id: gatk_selectvariants_vardict
      label: GATK Select Vardict PASS
      in:
      - id: input_vcf
        source: bcbio_filter_fp_somatic/filtered_vcf
      - id: output_basename
        source: output_basename
      - id: tool_name
        valueFrom: ${return "vardict"}
      - id: mode
        source: select_vars_mode
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        label: GATK Select PASS

        requirements:
        - class: ShellCommandRequirement
        - class: InlineJavascriptRequirement
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0
        - class: ResourceRequirement
          coresMin: 4
          ramMin: 8000

        inputs:
        - id: input_vcf
          type: File
          secondaryFiles:
          - .tbi
        - id: output_basename
          type: string
        - id: tool_name
          type: string
        - id: mode
          doc: Choose 'gatk' for SelectVariants tool, or 'grep' for grep expression
          type:
          - 'null'
          - name: select_vars_mode
            type: enum
            symbols:
            - gatk
            - grep
          default: gatk

        outputs:
        - id: pass_vcf
          type: File
          secondaryFiles:
          - .tbi
          outputBinding:
            glob: '*.PASS.vcf.gz'

        baseCommand:
        - /bin/bash
        - -c
        arguments:
        - position: 0
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
          shellQuote: false
        id: gatk4_selectvariants
      out:
      - pass_vcf
    - id: vep_annot_vardict
      in:
      - id: input_vcf
        source: gatk_selectvariants_vardict/pass_vcf
      - id: output_basename
        source: output_basename
      - id: tumor_id
        source: input_tumor_name
      - id: normal_id
        source: input_normal_name
      - id: tool_name
        valueFrom: ${return "vardict_somatic"}
      - id: reference
        source: indexed_reference_fasta
      - id: cache
        source: vep_cache
      - id: ref_build
        source: vep_ref_build
      run:
        cwlVersion: v1.0
        class: CommandLineTool

        requirements:
        - class: ShellCommandRequirement
        - class: InlineJavascriptRequirement
        - class: ResourceRequirement
          coresMin: 16
          ramMin: 24000
        - class: DockerRequirement
          dockerPull: pgc-images.sbgenomics.com/d3b-bixu/vep:r93_v2

        inputs:
        - id: reference
          label: Fasta genome assembly with index
          type: File
          secondaryFiles:
          - .fai
        - id: input_vcf
          type: File
          secondaryFiles:
          - .tbi
        - id: output_basename
          type: string
        - id: tumor_id
          type: string
        - id: normal_id
          type: string
        - id: tool_name
          type: string
        - id: cache
          label: tar gzipped cache from ensembl/local converted cache
          type: File
        - id: cache_version
          doc: Version being used, should match build version
          type:
          - 'null'
          - int
          default: 93
        - id: ref_build
          doc: Genome ref build used, should line up with cache.
          type:
          - 'null'
          - string
          default: GRCh38

        outputs:
        - id: output_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
        - id: output_tbi
          type: File
          outputBinding:
            glob: '*.vcf.gz.tbi'
        - id: output_maf
          type: File
          outputBinding:
            glob: '*.maf'
        - id: warn_txt
          type:
          - 'null'
          - File
          outputBinding:
            glob: '*.txt'

        baseCommand:
        - tar
        - -xzf
        arguments:
        - position: 1
          valueFrom: |-
            $(inputs.cache.path) && gunzip -c $(inputs.input_vcf.path) > input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf input_file.vcf --output-maf $(inputs.output_basename).$(inputs.tool_name).vep.maf --filter-vcf 0 --vep-path /ensembl-vep/ --vep-data $PWD --vep-forks 16 --ncbi-build $(inputs.ref_build) --cache-version $(inputs.cache_version) --ref-fasta $(inputs.reference.path) --tumor-id $(inputs.tumor_id) --normal-id $(inputs.normal_id) && mv input_file.vep.vcf $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf && /ensembl-vep/htslib/bgzip $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf && /ensembl-vep/htslib/tabix $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf.gz
          shellQuote: false
        id: kfdrc-vep-somatic-annotate-maf
      out:
      - output_vcf
      - output_tbi
      - output_maf
      - warn_txt
    id: kfdrc_vardict_1_7_sub_wf
  out:
  - vardict_vep_somatic_only_vcf
  - vardict_vep_somatic_only_tbi
  - vardict_vep_somatic_only_maf
  - vardict_prepass_vcf
  hints:
  - class: sbg:AWSInstanceType
    value: c5.9xlarge

hints:
- class: sbg:maxNumberOfParallelInstances
  value: 4
id: |-
  https://cavatica-api.sbgenomics.com/v2/apps/kfdrc-harmonization/kf-reference-pipeline/kfdrc_production_vardict_wf/2/raw/
sbg:appVersion:
- v1.0
sbg:content_hash: a792f7de8319d3a6e47c9f010f679af0c3f11e831bc7432b0f0bafabf37a71d87
sbg:contributors:
- brownm28
- danmiller
sbg:createdBy: danmiller
sbg:createdOn: 1641845305
sbg:id: kfdrc-harmonization/kf-reference-pipeline/kfdrc_production_vardict_wf/2
sbg:image_url: |-
  https://cavatica.sbgenomics.com/ns/brood/images/kfdrc-harmonization/kf-reference-pipeline/kfdrc_production_vardict_wf/2.png
sbg:latestRevision: 2
sbg:modifiedBy: brownm28
sbg:modifiedOn: 1655925964
sbg:project: kfdrc-harmonization/kf-reference-pipeline
sbg:projectName: kf-reference-pipeline
sbg:publisher: sbg
sbg:revision: 2
sbg:revisionNotes: |-
  Uploaded using sbpack v2022.03.16. 
  Source: 
  repo: https://github.com/kids-first/kf-somatic-workflow
  file: workflow/kfdrc_production_vardict_wf.cwl
  commit: 1.1.1-190-gc221996
sbg:revisionsInfo:
- sbg:modifiedBy: danmiller
  sbg:modifiedOn: 1641845305
  sbg:revision: 0
  sbg:revisionNotes: |-
    Uploaded using sbpack v2020.10.05. 
    Source: 
    repo: git@github.com:kids-first/kf-somatic-workflow.git
    file: workflow/kfdrc_production_vardict_wf.cwl
    commit: (uncommitted file)
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1655925712
  sbg:revision: 1
  sbg:revisionNotes: |-
    Uploaded using sbpack v2022.03.16. 
    Source: 
    repo: https://github.com/kids-first/kf-somatic-workflow
    file: workflow/kfdrc_production_vardict_wf.cwl
    commit: 1.1.1-281-gf1c453c
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1655925964
  sbg:revision: 2
  sbg:revisionNotes: |-
    Uploaded using sbpack v2022.03.16. 
    Source: 
    repo: https://github.com/kids-first/kf-somatic-workflow
    file: workflow/kfdrc_production_vardict_wf.cwl
    commit: 1.1.1-190-gc221996
sbg:sbgMaintained: false
sbg:validationErrors: []
sbg:workflowLanguage: CWL
