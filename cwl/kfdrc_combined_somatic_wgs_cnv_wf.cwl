class: Workflow
cwlVersion: v1.0
id: kfdrc_combined_somatic_wgs_cnv_wf
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: annotation_file
    type: File
    doc: refFlat.txt file
    'sbg:x': 296.890625
    'sbg:y': 1165.5
  - id: b_allele
    type: File?
    doc: >-
      germline calls, needed for BAF.  GATK HC VQSR input recommended.  Tool
      will prefilter for germline and pass if expression given
    secondaryFiles:
      - .tbi
    'sbg:x': 0
    'sbg:y': 1712
  - id: capture_regions
    type: File?
    doc: 'If not WGS, provide this bed file'
    'sbg:x': 296.890625
    'sbg:y': 909.5
  - id: cfree_sex
    type:
      - 'null'
      - type: enum
        symbols:
          - XX
          - XY
        name: sex
    doc: 'If known, XX for female, XY for male'
    'sbg:x': 596.5790405273438
    'sbg:y': 1268
  - id: chr_len
    type: File
    doc: file with chromosome lengths
    'sbg:x': 596.5790405273438
    'sbg:y': 1161
  - id: cnvkit_cnn_input
    type: File?
    doc: 'If running using an existing .cnn, supply here'
    'sbg:x': 0
    'sbg:y': 1605
  - id: cnvkit_sex
    type: string?
    doc: 'If known, choices are m,y,male,Male,f,x,female,Female'
    'sbg:x': 0
    'sbg:y': 1498
  - id: coeff_var
    type: float?
    doc: Coefficient of variation to set window size.  Default 0.05 recommended
    default: 0.05
    'sbg:x': 596.5790405273438
    'sbg:y': 793
  - id: combined_exclude_expression
    type: string?
    doc: 'Filter expression if vcf has non-PASS combined calls, use as-needed'
    'sbg:x': 0
    'sbg:y': 1391
  - id: combined_include_expression
    type: string?
    doc: 'Filter expression if vcf has non-PASS combined calls, use as-needed'
    'sbg:x': 0
    'sbg:y': 1284
  - id: contamination_adjustment
    type: boolean?
    doc: TRUE or FALSE to have ControlFreec estimate normal contam
    'sbg:x': 596.5790405273438
    'sbg:y': 686
  - id: indexed_reference_fasta
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    'sbg:x': 0
    'sbg:y': 1177
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
    'sbg:y': 1070
  - id: input_normal_name
    type: string
    'sbg:x': 0
    'sbg:y': 963
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
    'sbg:y': 856
  - id: input_tumor_name
    type: string
    'sbg:x': 0
    'sbg:y': 749
  - id: mate_orientation_control
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - FR
          - RF
          - FF
        name: mate_orientation_control
    doc: >-
      0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends),
      FF (SOLiD mate-pairs)
    default: FR
    'sbg:x': 0
    'sbg:y': 642
  - id: mate_orientation_sample
    type:
      - 'null'
      - type: enum
        symbols:
          - '0'
          - FR
          - RF
          - FF
        name: mate_orientation_sample
    doc: >-
      0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends),
      FF (SOLiD mate-pairs)
    default: FR
    'sbg:x': 0
    'sbg:y': 535
  - id: min_theta2_frac
    type: float?
    doc: >-
      Minimum fraction of genome with copy umber alterations.  Default is 0.05,
      recommend 0.01
    default: 0.01
    'sbg:x': 1018.6522216796875
    'sbg:y': 355.5
  - id: output_basename
    type: string
    'sbg:x': 1018.6522216796875
    'sbg:y': 248.5
  - id: paired_vcf
    type: File
    doc: Combined somatic and germline call file. VarDict input recommended.
    'sbg:x': 0
    'sbg:y': 428
  - id: ploidy
    type: 'int[]'
    doc: Array of ploidy possibilities for ControlFreeC to try
    'sbg:x': 0
    'sbg:y': 321
  - id: reference_fai
    type: File
    'sbg:x': 0
    'sbg:y': 214
  - id: threads
    type: int?
    doc: >-
      For ControlFreeC.  Recommend 16 max, as I/O gets saturated after that
      losing any advantage.
    default: 16
    'sbg:x': 0
    'sbg:y': 107
  - id: wgs_mode
    type: string?
    doc: 'for WGS mode, input Y. leave blank for hybrid mode'
    default: 'Y'
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: cnvkit_calls
    outputSource:
      - cnvkit/output_calls
    type: File
    'sbg:x': 1018.6522216796875
    'sbg:y': 1463.5
  - id: cnvkit_cnn_output
    outputSource:
      - cnvkit/output_cnn
    type: File?
    'sbg:x': 1018.6522216796875
    'sbg:y': 1356.5
  - id: cnvkit_cnr
    outputSource:
      - cnvkit/output_cnr
    type: File
    'sbg:x': 1018.6522216796875
    'sbg:y': 1249.5
  - id: cnvkit_gainloss
    outputSource:
      - cnvkit/output_gainloss
    type: File
    'sbg:x': 1018.6522216796875
    'sbg:y': 979.5
  - id: cnvkit_metrics
    outputSource:
      - cnvkit/output_metrics
    type: File
    'sbg:x': 1018.6522216796875
    'sbg:y': 872.5
  - id: cnvkit_seg
    outputSource:
      - cnvkit/output_seg
    type: File
    'sbg:x': 1018.6522216796875
    'sbg:y': 765.5
  - id: ctrlfreec_baf
    outputSource:
      - rename_outputs/ctrlfreec_baf
    type: File
    'sbg:x': 1979.84619140625
    'sbg:y': 1149
  - id: ctrlfreec_bam_ratio
    outputSource:
      - rename_outputs/ctrlfreec_bam_ratio
    type: File
    'sbg:x': 1979.84619140625
    'sbg:y': 1042
  - id: ctrlfreec_bam_seg
    outputSource:
      - convert_ratio_to_seg/ctrlfreec_ratio2seg
    type: File
    'sbg:x': 1979.84619140625
    'sbg:y': 935
  - id: ctrlfreec_config
    outputSource:
      - rename_outputs/ctrlfreec_config
    type: File
    'sbg:x': 1979.84619140625
    'sbg:y': 828
  - id: ctrlfreec_info
    outputSource:
      - rename_outputs/ctrlfreec_info
    type: File
    'sbg:x': 1979.84619140625
    'sbg:y': 721
  - id: ctrlfreec_pngs
    outputSource:
      - rename_outputs/ctrlfreec_pngs
    type: 'File[]'
    'sbg:x': 1979.84619140625
    'sbg:y': 614
  - id: ctrlfreec_pval
    outputSource:
      - rename_outputs/ctrlfreec_pval
    type: File
    'sbg:x': 1979.84619140625
    'sbg:y': 507
  - id: theta2_calls
    outputSource:
      - cnvkit_import_theta2/theta2_adjusted_cns
    type: File
    'sbg:x': 2335.287109375
    'sbg:y': 1016.5
  - id: theta2_seg
    outputSource:
      - cnvkit_import_theta2/theta2_adjusted_seg
    type: File
    'sbg:x': 2335.287109375
    'sbg:y': 909.5
  - id: theta2_subclonal_cns
    outputSource:
      - cnvkit_import_theta2/theta2_subclone_cns
    type: 'File[]'
    'sbg:x': 2335.287109375
    'sbg:y': 802.5
  - id: theta2_subclonal_results
    outputSource:
      - run_theta2/n3_graph
      - run_theta2/n2_results
      - run_theta2/best_results
    type: 'File[]'
    'sbg:x': 1979.84619140625
    'sbg:y': 400
  - id: theta2_subclone_seg
    outputSource:
      - cnvkit_import_theta2/theta2_subclone_seg
    type: 'File[]'
    'sbg:x': 2335.287109375
    'sbg:y': 695.5
steps:
  - id: bcftools_filter_combined_vcf
    in:
      - id: exclude_expression
        source: combined_exclude_expression
      - id: include_expression
        source: combined_include_expression
      - id: input_vcf
        source: paired_vcf
      - id: output_basename
        source: output_basename
    out:
      - id: filtered_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: bcftools_filter_vcf
      baseCommand: []
      inputs:
        - id: exclude_expression
          type: string?
        - id: include_expression
          type: string?
        - id: input_vcf
          type: File
        - id: output_basename
          type: string?
      outputs:
        - id: filtered_vcf
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
      doc: >-
        More generic tool to take in an include expression and optionally an
        exclude expresssion to filter a vcf
      arguments:
        - position: 1
          shellQuote: false
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
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 1000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'kfdrc/bvcftools:latest'
        - class: InlineJavascriptRequirement
    'sbg:x': 296.890625
    'sbg:y': 1037.5
  - id: cnvkit
    in:
      - id: annotation_file
        source: annotation_file
      - id: b_allele_vcf
        source: gatk_filter_germline/filtered_pass_vcf
      - id: capture_regions
        source: capture_regions
      - id: cnvkit_cnn
        source: cnvkit_cnn_input
      - id: input_control
        source: samtools_normal_cram2bam/bam_file
      - id: input_sample
        source: samtools_tumor_cram2bam/bam_file
      - id: output_basename
        source: output_basename
      - id: reference
        source: indexed_reference_fasta
      - id: sex
        source: cnvkit_sex
      - id: threads
        source: threads
      - id: tumor_sample_name
        source: input_tumor_name
      - id: wgs_mode
        source: wgs_mode
    out:
      - id: output_calls
      - id: output_cnn
      - id: output_cnr
      - id: output_diagram
      - id: output_gainloss
      - id: output_metrics
      - id: output_scatter
      - id: output_seg
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: cnvkit_batch
      baseCommand: []
      inputs:
        - id: annotation_file
          type: File?
          doc: 'refFlat.txt file,  needed if cnv kit cnn not already built'
        - id: b_allele_vcf
          type: File?
          doc: 'b allele germline vcf, if available'
        - id: capture_regions
          type: File?
          doc: target regions for WES
        - id: cnvkit_cnn
          type: File?
          doc: 'If running using an existing .cnn, supply here'
        - id: input_control
          type: File?
          doc: normal bam file
          secondaryFiles:
            - ^.bai
        - id: input_sample
          type: File
          doc: tumor bam file
          secondaryFiles:
            - ^.bai
        - id: output_basename
          type: string
        - id: reference
          type: File?
          doc: 'fasta file, needed if cnv kit cnn not already built'
          secondaryFiles:
            - .fai
        - id: sex
          type: string
          doc: Set sample sex.  CNVkit isn't always great at guessing it
        - default: 16
          id: threads
          type: int?
        - id: tumor_sample_name
          type: string
        - id: wgs_mode
          type: string?
          doc: 'for WGS mode, input Y. leave blank for hybrid mode'
      outputs:
        - id: output_calls
          type: File
          outputBinding:
            glob: '*.call.cns'
        - id: output_cnn
          doc: >-
            Output if starting from cnn scratch.  Should not appear if an
            existing .cnn was given as input.
          type: File?
          outputBinding:
            glob: '*_reference.cnn'
        - id: output_cnr
          type: File
          outputBinding:
            glob: '*.cnr'
        - id: output_diagram
          type: File
          outputBinding:
            glob: '*.diagram.pdf'
        - id: output_gainloss
          type: File
          outputBinding:
            glob: '*.gainloss.txt'
        - id: output_metrics
          type: File
          outputBinding:
            glob: '*.metrics.txt'
        - id: output_scatter
          type: File
          outputBinding:
            glob: '*.scatter.pdf'
        - id: output_seg
          type: File
          outputBinding:
            glob: '*.seg'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            ln -s $(inputs.input_sample.path) .; ln -s
            $(inputs.input_sample.secondaryFiles[0].path)
            ./$(inputs.input_sample.basename).bai

            ${ 
                var cmd = "";
                if (inputs.input_control != null) {
                    cmd = "ln -s " + inputs.input_control.path + " .; ln -s " + inputs.input_control.secondaryFiles[0].path + " ./" + inputs.input_control.basename + ".bai"
                }
                return cmd;
            }

            cnvkit.py batch -p $(inputs.threads) ${
                var cmd = "";
                if (inputs.wgs_mode == 'Y') {
                    cmd = " -m wgs ";
                }
                return cmd;
            } $(inputs.input_sample.path) ${
                var cmd = "";
                if (inputs.capture_regions != null) {
                    cmd = "--targets " + inputs.capture_regions.path;
                }
                return cmd;
            } ${
              if (inputs.cnv_kit_cnn == null){
                var arg = "--output-reference " + inputs.output_basename + "_cnvkit_reference.cnn --fasta " + inputs.reference.path + " --annotate " + inputs.annotation_file.path;
                if (inputs.input_control != null) {
                    arg += " --normal " + inputs.input_control.path;
                }
              }
              else{
                var arg = "--reference " + inputs.cnv_kit_cnn.path;
                var msex = ['m','y','male','Male']
                if (msex.indexOf(inputs.sex) >= 0){
                  arg += " --male-reference";
                }
              }
              return arg;
            } --diagram  --scatter

            cnvkit.py call $(inputs.input_sample.nameroot).cns ${
              var arg = "";
              if (inputs.b_allele_vcf != null){
                arg = "--vcf " + inputs.b_allele_vcf.path;
              }
              return arg;
            } ${
              var arg = "--sample-sex " + inputs.sex;
              var msex = ['m','y','male','Male']
              if (msex.indexOf(inputs.sex) >= 0){
                arg += " --male-reference";
              }
              return arg;
            } -o $(inputs.output_basename).call.cns
                  
            ln -s $(inputs.output_basename).call.cns
            $(inputs.tumor_sample_name).cns

            cnvkit.py export seg $(inputs.tumor_sample_name).cns -o
            $(inputs.output_basename).call.seg

            rm $(inputs.tumor_sample_name).cns

            cnvkit.py metrics $(inputs.input_sample.nameroot).cnr -s
            $(inputs.input_sample.nameroot).cns -o
            $(inputs.output_basename).metrics.txt

            cnvkit.py gainloss $(inputs.input_sample.nameroot).cnr -o
            $(inputs.output_basename).gainloss.txt

            mv $(inputs.input_sample.nameroot).cnr $(inputs.output_basename).cnr

            mv $(inputs.input_sample.nameroot)-diagram.pdf
            $(inputs.output_basename).diagram.pdf

            mv $(inputs.input_sample.nameroot)-scatter.pdf
            $(inputs.output_basename).scatter.pdf
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 32000
          coresMin: $(inputs.threads)
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
        - class: InlineJavascriptRequirement
    'sbg:x': 596.5790405273438
    'sbg:y': 977
  - id: cnvkit_export_theta2
    in:
      - id: normal_ID
        source: input_normal_name
      - id: paired_vcf
        source: bcftools_filter_combined_vcf/filtered_vcf
      - id: reference_cnn
        source: cnvkit/output_cnn
      - id: tumor_ID
        source: input_tumor_name
      - id: tumor_cns
        source: cnvkit/output_calls
    out:
      - id: call_interval_count
      - id: call_normal_snp
      - id: call_tumor_snp
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: cnvkit_export_theta2
      baseCommand:
        - cnvkit.py
        - export
        - theta
      inputs:
        - id: normal_ID
          type: string
        - id: paired_vcf
          type: File
        - id: reference_cnn
          type: File
        - id: tumor_ID
          type: string
        - id: tumor_cns
          type: File
      outputs:
        - id: call_interval_count
          type: File
          outputBinding:
            glob: '*.call.interval_count'
        - id: call_normal_snp
          type: File
          outputBinding:
            glob: '*.call.normal.snp_formatted.txt'
        - id: call_tumor_snp
          type: File
          outputBinding:
            glob: '*.call.tumor.snp_formatted.txt'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            -r $(inputs.reference_cnn.path) -v $(inputs.paired_vcf.path) -i
            $(inputs.tumor_ID) -n $(inputs.normal_ID) $(inputs.tumor_cns.path)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 16000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
        - class: InlineJavascriptRequirement
    'sbg:x': 1018.6522216796875
    'sbg:y': 1114.5
  - id: cnvkit_import_theta2
    in:
      - id: output_basename
        source: output_basename
      - id: theta2_best_results
        source: run_theta2/best_results
      - id: theta2_n2_results
        source: run_theta2/n2_results
      - id: tumor_cns
        source: cnvkit/output_calls
      - id: tumor_sample_name
        source: input_tumor_name
    out:
      - id: theta2_adjusted_cns
      - id: theta2_adjusted_seg
      - id: theta2_subclone_cns
      - id: theta2_subclone_seg
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: cnvkit_import_theta2
      baseCommand:
        - cnvkit.py
        - import-theta
      inputs:
        - id: output_basename
          type: string
        - id: theta2_best_results
          type: File
        - id: theta2_n2_results
          type: File
        - id: tumor_cns
          type: File
        - id: tumor_sample_name
          type: string
      outputs:
        - id: theta2_adjusted_cns
          type: File
          outputBinding:
            glob: '*.theta2.total.cns'
        - id: theta2_adjusted_seg
          type: File
          outputBinding:
            glob: '*.theta2.total.seg'
        - id: theta2_subclone_cns
          type: 'File[]'
          outputBinding:
            glob: '*.theta2.subclone*.cns'
        - id: theta2_subclone_seg
          type: 'File[]'
          outputBinding:
            glob: '*.theta2.subclone*.seg'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            $(inputs.tumor_cns.path) $(inputs.theta2_n2_results.path) -d ./

            mv $(inputs.output_basename).call-1.cns
            $(inputs.output_basename).theta2.total.cns

            ln -s $(inputs.output_basename).theta2.total.cns
            $(inputs.tumor_sample_name).cns

            cnvkit.py export seg $(inputs.tumor_sample_name).cns -o
            $(inputs.output_basename).theta2.total.seg

            rm $(inputs.tumor_sample_name).cns

            cnvkit.py import-theta $(inputs.tumor_cns.path)
            $(inputs.theta2_best_results.path) -d ./

            mv $(inputs.output_basename).call-1.cns
            $(inputs.output_basename).theta2.subclone1.cns

            ln -s $(inputs.output_basename).theta2.subclone1.cns
            $(inputs.tumor_sample_name).cns

            cnvkit.py export seg $(inputs.tumor_sample_name).cns -o
            $(inputs.output_basename).theta2.subclone1.seg

            rm $(inputs.tumor_sample_name).cns

            SC2=$(inputs.output_basename).call-2.cns;

            if [ -f "$SC2" ]; then
                mv $(inputs.output_basename).call-2.cns $(inputs.output_basename).theta2.subclone2.cns;
                ln -s $(inputs.output_basename).theta2.subclone2.cns $(inputs.tumor_sample_name).cns;
                cnvkit.py export seg $(inputs.tumor_sample_name).cns -o $(inputs.output_basename).theta2.subclone2.seg;
                rm $(inputs.tumor_sample_name).cns;
            else
              echo "second subclone file not found.  skipping!";
            fi
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 4
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/milos_nikolic/cnvkit:0.9.3'
        - class: InlineJavascriptRequirement
    'sbg:x': 1979.84619140625
    'sbg:y': 1284
  - id: control_free_c
    in:
      - id: capture_regions
        source: capture_regions
      - id: chr_len
        source: chr_len
      - id: coeff_var
        source: coeff_var
      - id: contamination_adjustment
        source: contamination_adjustment
      - id: mate_file_control
        source: samtools_normal_cram2bam/bam_file
      - id: mate_file_sample
        source: samtools_tumor_cram2bam/bam_file
      - id: mate_orientation_control
        source: mate_orientation_control
      - id: mate_orientation_sample
        source: mate_orientation_sample
      - id: max_threads
        source: threads
      - id: mini_pileup_control
        source: controlfreec_normal_mini_pileup/pileup
      - id: mini_pileup_sample
        source: controlfreec_tumor_mini_pileup/pileup
      - id: ploidy
        source:
          - ploidy
      - id: reference
        source: indexed_reference_fasta
      - id: sex
        source: cfree_sex
      - id: snp_file
        source: gatk_filter_germline/filtered_pass_vcf
    out:
      - id: GC_profile
      - id: cnvs
      - id: cnvs_pvalue
      - id: config_script
      - id: control_cpn
      - id: control_pileup
      - id: info_txt
      - id: pngs
      - id: ratio
      - id: ratio_BedGraph
      - id: sample_BAF
      - id: sample_cpn
      - id: sample_pileup
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        sbg: 'https://sevenbridges.com'
      id: brownm28/mb-controlfreec-troubleshoot/control-freec-11-6-sbg/0
      baseCommand: []
      inputs:
        - 'sbg:category': General
          id: GC_content_profile
          type: File?
          label: GC content profile
          doc: GC-content profile for a given window-size.
          'sbg:fileTypes': CNP
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: bed_graph_output
          type: boolean?
          label: Bed Graph output
          doc: >-
            Set "BedGraphOutput=TRUE" if you want an additional output in
            BedGraph format for the UCSC genome browser.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '0.8'
          id: break_point_threshold
          type: float?
          label: Break point threshold
          doc: >-
            Positive value of threshold for segmentation of normalized profiles.
            The closer it is to zero, the more breakpoints will be called. Its
            recommended value is between 0.1 and 1.2.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '2'
          id: break_point_type
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - '1'
                - '2'
                - '3'
                - '4'
              name: break_point_type
          label: Break point type
          doc: >-
            Desired behavior in the ambiguous regions (poly-N or low mappability
            regions between two different copy number values). 0: the "unknown"
            region is attached to the "known" region on the right  1: make a
            separate fragment of this “unknown” region and then attaches it to
            the left or to the right region choosing the longer one  2: make a
            separate fragment of this “unknown” region and then attaches it to
            the left or to the right region but the “ploidy” copy number has a
            priority  3: make a separate fragment of this “unknown” region and
            then attaches it to the left or to the right region choosing the
            longer one but this “known” region should make at least half-size of
            the “unknown” region  4: make a separate fragment of this “unknown”
            region and do not assign any copy number to this region at all
        - 'sbg:category': Target
          id: capture_regions
          type: File?
          label: Capture regions
          doc: >-
            Capture regions in .bed format; sorted .bed file should contain the
            following colomns: chr   0-based start   1-based end.
          'sbg:fileTypes': BED
        - 'sbg:category': File Input
          id: chr_len
          type: File
          label: Chromosomes length file
          doc: Chromosome length file in a tab-delimited format.
          'sbg:fileTypes': 'TXT, LEN, SIZES'
        - 'sbg:toolDefaultValue': '0.05'
          id: coeff_var
          type: float?
          label: Coefficient of variation
          doc: Coefficient of variation to evaluate necessary window size.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '0'
          id: contamination
          type: float?
          label: Contamination
          doc: >-
            A priori known value of tumor sample contamiantion by normal cells.
            Set "contaminationAdjustment=TRUE" to correct for the contamination.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: contamination_adjustment
          type: boolean?
          label: Contamination adjustment
          doc: >-
            Set TRUE to correct for contamination by normal cells. If
            "contamination" is not provided, it will automatically evaluate the
            level of contamination.
        - 'sbg:category': General
          'sbg:toolDefaultValue': >-
            3&4 (GC-content based normalization, WGS) or 1
            (control-read-count-based normalization, WES)
          id: degree
          type: float?
          label: Degree
          doc: Degree of polynomial.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'WGS: 0, WES: 1'
          id: force_GC_content_normalization
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - '1'
                - '2'
              name: force_GC_content_normalization
          label: Force GC content normalization
          doc: >-
            Set 1 or 2 to correct the Read Count (RC) for GC-content bias and
            low mappability even when you have a control sample. 0: simply model
            "sample RC ~ Control RC"  1: normalize the sample and the control RC
            using GC-content and then calculate the ratio "Sample RC/contol RC"
            2: model "sample RC ~ Control RC" bias, and then normalize for
            GC-content.
        - 'sbg:category': File Input
          id: gem_mappability_file
          type: File?
          label: GEM mappability file
          doc: Mappability file in GEM format.
          'sbg:fileTypes': GEM
        - 'sbg:category': Execution
          'sbg:toolDefaultValue': '1 - with GC-content, 0 - with a control dataset'
          id: intercept
          type: float?
          label: Intercept
          doc: Intercept of polynomial.
        - 'sbg:category': Control
          id: mate_copynumber_file_control
          type: File?
          label: Mate copy number file control
          doc: >-
            Raw copy number profile for a given window-size (higher priority
            than mateFile)  (don't need to provide a mateFile if
            mateCopyNumberFile is provided).
          'sbg:fileTypes': CPN
        - 'sbg:category': Sample
          id: mate_copynumber_file_sample
          type: File?
          label: Mate copy number file sample
          doc: >-
            Raw copy number profile for a given window-size (higher priority
            than mateFile)  (don't need to provide a mateFile if
            mateCopyNumberFile is provided).
          'sbg:fileTypes': CPN
        - 'sbg:category': Control
          id: mate_file_control
          type: File?
          label: Mate file Control
          doc: >-
            Mapped reads (can be single end reads, mate-pairs or paired-end
            reads).
          'sbg:fileTypes': 'SAM, BAM, PILEUP, PILEUP.GZ'
        - 'sbg:category': Sample
          id: mate_file_sample
          type: File?
          label: Mate file Sample
          doc: >-
            Mapped reads (can be single end reads, mate-pairs or paired-end
            reads).
          'sbg:fileTypes': 'SAM, BAM, PILEUP, PILEUP.GZ'
        - 'sbg:category': Control
          id: mate_orientation_control
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - RF
                - FR
                - FF
              name: mate_orientation_control
          label: Mate orientation control
          doc: >-
            Format of reads (in mateFile). 0 (for single ends), RF (Illumina
            mate-pairs),  FR (Illumina paired-ends), FF (SOLiD mate-pairs).
        - 'sbg:category': Sample
          id: mate_orientation_sample
          type:
            - 'null'
            - type: enum
              symbols:
                - '0'
                - RF
                - FR
                - FF
              name: mate_orientation_sample
          label: Mate orientation sample
          doc: >-
            Format of reads (in mateFile). 0 (for single ends), RF (Illumina
            mate-pairs),  FR (Illumina paired-ends), FF (SOLiD mate-pairs).
        - 'sbg:category': General
          'sbg:toolDefaultValue': 0.55 (change only if you run Control-FREEC on a bacterial genome)
          id: max_expected_GC
          type: float?
          label: Maximum expected GC
          doc: >-
            Maximal exptected value of the GC-content for the prior evaluation
            of "Read Count ~ GC-content" dependancy.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '8'
          id: max_threads
          type: int?
          label: Maximum threads
          doc: Number of threads (multi-threading mode).
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'WES: 3, WGS: 1'
          id: min_CNA_length
          type: int?
          label: Minimum CNA length
          doc: Minimal number of consecutive windows to call a CNA.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 0.35 (change only if you run Control-FREEC on a bacterial genome)
          id: min_expected_GC
          type: float?
          label: Minimal exptected GC
          doc: >-
            Minimal expected value of the GC-content for the prior evaluation of
            "Read Count ~ GC-content" dependency.
        - 'sbg:category': General
          'sbg:toolDefaultValue': '0.85'
          id: min_map_per_w
          type: float?
          label: Minimum mappability per window
          doc: >-
            Only windows with fraction of mappable positions higher than or
            equal to this threshold will be considered  (if "gemMappabilityFile"
            is not provided, one uses the percentage of non-N letters per
            window).
        - 'sbg:category': General
          'sbg:toolDefaultValue': >-
            100 (meaning "do not look for subclones"). Suggested: 20 (or 0.2)
            for WGS and 30 (or 0.3) for WES.
          id: min_subclone_presence
          type: int?
          label: Minimal subclone presence
          doc: Detects subclones present in x% of cell population
        - 'sbg:category': Control
          id: mini_pileup_control
          type: File?
          label: Mini pileup Control
          doc: >-
            Mini pileup file created from the corresponding BAM file dring a
            previous run of Control-FREEC - providing this file will
            significantly speed up the whole process.
          'sbg:fileTypes': PILEUP
        - 'sbg:category': Sample
          id: mini_pileup_sample
          type: File?
          label: Mini pileup Sample
          doc: >-
            Mini pileup file created from the corresponding BAM file dring a
            previous run of Control-FREEC - providing this file will
            significantly speed up the whole process.
          'sbg:fileTypes': PILEUP
        - 'sbg:category': BAF
          'sbg:toolDefaultValue': '0'
          id: minimal_coverage_per_position
          type: int?
          label: Minimal coverage per position
          doc: >-
            Minimal read coverage for a position to be considered in BAF
            analysis.
        - 'sbg:category': BAF
          'sbg:toolDefaultValue': '0'
          id: minimal_quality_per_position
          type: int?
          label: Minimal quality per position
          doc: >-
            Minimal sequencing quality for a position to be considered in BAF
            analysis. Default: 0; using this option can slow down reading of
            pileup files.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: noisy_data
          type: boolean?
          label: Noisy data
          doc: >-
            Set TRUE for target resequencing data (e.g., exome-seq) to avoid
            false positive predictions due to nonuniform capture.
        - 'sbg:category': General
          id: ploidy
          type: 'int[]'
          label: Ploidy
          doc: >-
            Genome ploidy. In case of doubt, you can set different values and
            Control-FREEC will select the one that explains most observed CNAs
            (eg. 2,3,4)
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'TRUE'
          id: print_NA
          type: boolean?
          label: Print NA
          doc: >-
            Set FALSE to avoid printing "-1" to the _ratio.txt files Useful for
            exome-seq or targeted sequencing data.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 10; recommended value >=50 for for exome data
          id: read_cnt_threshold
          type: int?
          label: Read count threshold
          doc: >-
            Threshold on the minimal number of reads per window in the control
            sample Useful for exome-seq or targeted sequencing data.
        - 'sbg:category': BAF
          id: reference
          type: File
          label: Reference file
          doc: >-
            Reference file that will be divided as needed to separate
            chromosomes and contigs.
          'sbg:fileTypes': 'FA, FASTA'
        - 'sbg:category': General
          id: sex
          type:
            - 'null'
            - type: enum
              symbols:
                - XX
                - XY
              name: sex
          label: Sex
          doc: Sample sex.
        - 'sbg:category': BAF
          'sbg:toolDefaultValue': '0'
          id: shift_in_quality
          type: int?
          label: Shift in quality
          doc: >-
            Basis for Phred quality. Default: 0; usually 33 or 64; see fastq
            quality.
        - 'sbg:category': BAF
          id: snp_file
          type: File?
          label: Known SNPs
          doc: Known SNPs.
          'sbg:fileTypes': 'TXT, VCF'
        - 'sbg:category': General
          id: step
          type: int?
          label: Step
          doc: 'Step (used only when "window" is specified - Ex: 10000).'
        - 'sbg:category': General
          'sbg:toolDefaultValue': '50000'
          id: telocentromeric
          type: int?
          label: Telocentromeric
          doc: >-
            Length of pre-telomeric and pre-centromeric regions: Control-FREEC
            will not output small CNAs and LOH found within these regions (they
            are likely to be false because of mappability/genome assembly
            issues) 50000 is OK for human/mouse genomes. Use smaller values for
            yeasts and flies.
        - 'sbg:category': Execution
          id: total_memory
          type: int?
          label: 'Total memory [MB]'
          doc: Total amount of memory in MB reserved on the instance.
        - 'sbg:category': General
          'sbg:toolDefaultValue': 'FALSE'
          id: unique_match
          type: boolean?
          label: Unique match
          doc: >-
            Use a mappability profile to correct read counts (in this case a
            mappability file must be provided with "gemMappabilityFile").
        - id: window
          type: int?
          label: Window
          doc: explicit window size (higher priority than coefficientOfVariation)
      outputs:
        - id: GC_profile
          doc: GC profile output file.
          label: GC profile
          type: File?
          outputBinding:
            glob: |-
              ${
                  return "GC_profile.targetedRegions.cnp"
              }
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': CNP
        - id: cnvs
          doc: File with coordinates of predicted copy number alterations.
          label: CNVs output
          type: File?
          outputBinding:
            glob: '*_CNVs'
            outputEval: |+
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }

          'sbg:fileTypes': TXT
        - id: cnvs_pvalue
          doc: >-
            File with coordinates of predicted copy number alterations with
            p-values.
          label: CNVs with p-value
          type: File?
          outputBinding:
            glob: '*.p.value.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: config_script
          doc: Configuration script used for running.
          label: Configuration script used for running
          type: File?
          outputBinding:
            glob: config.txt
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: control_cpn
          doc: >-
            Control CPN file to be used in the future for more efficient
            computation.
          label: Control CPN
          type: File?
          outputBinding:
            glob: '*_control.cpn'
            outputEval: |-
              ${  if (inputs.mate_file_control){
                  return inheritMetadata(self, inputs.mate_file_control)
              }
              }
          'sbg:fileTypes': CPN
        - id: control_pileup
          doc: >-
            Mini Pileup created for BAF calculation. It is used to speed up
            consequent runs with the same samples.
          label: Control Pileup
          type: File?
          outputBinding:
            glob: |-
              ${
                  if (inputs.mate_file_control) {

                      return inputs.mate_file_control.path.split('/').pop() + '_minipileup.pileup'

                  }
              }
            outputEval: |-
              ${
                  if (inputs.mate_file_control){
                  return inheritMetadata(self, inputs.mate_file_control)
              }
              }
          'sbg:fileTypes': PILEUP
        - id: info_txt
          doc: Parsable file with information about FREEC run.
          label: Info TXT
          type: File?
          outputBinding:
            glob: '*_info.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: pngs
          doc: >-
            Visalized normalized copy number profile with predicted CNAs as well
            as BAF profile (if dbSNP is provided)
          label: Copy number profile
          type: 'File[]?'
          outputBinding:
            glob: '*png'
          'sbg:fileTypes': PNG
        - id: ratio
          doc: >-
            File with ratios and predicted copy number alterations for each
            window.
          label: Ratio
          type: File?
          outputBinding:
            glob: '*_ratio.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: ratio_BedGraph
          doc: >-
            File with ratios in BedGraph format for visualization in the UCSC
            genome browser
          label: Ratio BedGraph
          type: File?
          outputBinding:
            glob: '*.BedGraph'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': BEDGRAPH
        - id: sample_BAF
          doc: >-
            File with B-allele frequencies for each possibly heterozygous SNP
            position
          label: BAF sample file
          type: File?
          outputBinding:
            glob: '*_BAF.txt'
            outputEval: |
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': TXT
        - id: sample_cpn
          doc: >-
            Sample CPN file to be used in the future for more efficient
            computation.
          label: Sample CPN
          type: File?
          outputBinding:
            glob: '*_sample.cpn'
            outputEval: |-
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': CPN
        - id: sample_pileup
          doc: >-
            Mini Pileup created for BAF calculation. It is used to speed up
            consequent runs with the same samples.
          label: Sample Pileup
          type: File?
          outputBinding:
            glob: |-
              ${
                  if (inputs.mate_file_sample) {

                      return inputs.mate_file_sample.path.split('/').pop() + '_minipileup.pileup'
                  }
              }
            outputEval: |-
              ${
                  return inheritMetadata(self, inputs.mate_file_sample)

              }
          'sbg:fileTypes': PILEUP
      doc: >-
        Control-FREEC analyzes copy-number variants and allelic imbalances in
        exome and whole-genome DNA sequencing.


        This tool automatically computes, normalizes and segments copy number
        and beta allele frequency (BAF) profiles, then calls copy number
        alterations and LOH. [1]


        *A list of **all inputs and parameters** with corresponding descriptions
        can be found at the bottom of the page.*

        ### Common Use Cases


        * The **chrLenFile** input is required and can be found in Public
        Reference Files as **Homo\_sapiens\_assembly38.fasta.sizes**,
        **ucsc.hg19.fasta.sizes** and **human\_g1k\_v37\_decoy.fasta.sizes**.


        * The **ploidy** parameter is required. In case of doubt, different
        values can be set and Control-FREEC will select the one that explains
        the most observed CNVs.


        * Normal and control sample can be provided through two possible inputs:
             * **mateFile**, a file with mapped reads
             * **mateCopyNumberFile**, a raw copy number file created for both normal and control sample, provided through **mateFile** in a first run, and can be reused in the future runs for more efficient computation.



        * **A control (matched normal) sample is optional for whole genome
        sequencing data but mandatory for whole exome or targeted sequencing
        data.**


        * Similar to **mateCopyNumberFile**, a **Mini pileup Sample** and **Mini
        pileup Control** files can be created in the first run, if the **Known
        SNPs** file is provided. Consequently, by providing these files as
        inputs in future tasks, execution time will decrease significantly.


        * If a **mateFile** is specified, the **mateOrientation** parameter must
        be set.


        * In order to create a **BAF profile**, one of the following options
        must be implemented:
            * **mateFile** + **Known SNPs** 
            * **mateCopyNumberFile** + **mateFile** + **KnownSNPs**
            * **mateCopyNumberFile** + **miniPileup** + **KnownSNPs**

        ### Changes Introduced by Seven Bridges


        * Based on the input parameters, a config file is created in order to
        properly run Control-FREEC.


        ### Common Issues and Important Notes


        * **A control (matched normal) sample is optional for whole genome
        sequencing data but mandatory for whole exome or targeted sequencing
        data.**


        * A **gemMappabilityFile** can be used only in the mode without a
        control sample.


        * If a **mateFile** is specified, the **mateOrientation** parameter must
        be set.


        * Currently, there is an issue with creating a **BAF sample file** with
        the b37 notation. The genotypes for CNV regions are, however, created.



        ### Performance Benchmarking


        The instance set for this tool is the AWS c4.2xlarge instance with 8
        vCPUs, 15 GiB of RAM and 1 TB of EBS (disk space).

        |     BAM size in GB    | Type |  Instance  | Duration | Cost ($) |

        |:--------------------:|:----:|:----------:|:--------:|:--------:|

        |         2x12 (Normal-Tumor)         |  WES | c4.2xlarge |  1h 52m 
        |    0.8   |

        |   100 (Tumor-only)   |  WGS | c4.2xlarge |  17h 43m |     7    |

        | 2x100 (Normal-Tumor) |  WGS | c4.2xlarge |   1d 8h  |    13    |

        |   100 (Tumor-only)   |  WGS | c4.8xlarge |  6h 30m |     10    |

        | 2x100 (Normal-Tumor) |  WGS | c4.8xlarge |   11h  |    18    |


        An instance with more resources can be obtained by providing inputs for
        **Maximum threads** and **Total memory [MB]**.


        *Cost can be significantly reduced by using **spot instances**. Visit
        the [Knowledge
        Center](https://docs.sevenbridges.com/docs/about-spot-instances) for
        more details.*  


        ###References

        [1] [Control-FREEC: Prediction of copy number alterations and loss of
        heterozygosity using deep-sequencing
        data](http://boevalab.com/FREEC/tutorial.html#install)
      label: Control-FREEC 11.6
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: |-
            ${
                // script for splitting the genome fasta into chromosomes
                return 'python split_fasta.py ' + inputs.reference.path

            }
        - position: 1
          shellQuote: false
          valueFrom: '&&'
        - position: 2
          prefix: ''
          shellQuote: false
          valueFrom: /opt/controlfreec/FREEC/src/freec
        - position: 3
          shellQuote: false
          valueFrom: '-conf'
        - position: 4
          shellQuote: false
          valueFrom: config.txt
        - position: 5
          shellQuote: false
          valueFrom: '&&'
        - position: 6
          shellQuote: false
          valueFrom: |-
            ${




                if (inputs.mate_file_sample) {
                    filepath = inputs.mate_file_sample.path
                    filename = filepath.split("/").pop()
                } else {
                    filepath = inputs.mate_copynumber_file_sample.path
                    filename = filepath.split("/").pop()
                }

                CNVs = filename + "_CNVs"
                ratio = filename + "_ratio" + ".txt"


                return "cat assess_significance.R | R --slave --args " + CNVs + " " + ratio
            }
        - position: 7
          shellQuote: false
          valueFrom: '&&'
        - position: 8
          shellQuote: false
          valueFrom: |-
            ${
                return "line=$(cat *info.txt | grep Output_Ploidy | sed -E 's/.+([0-9]+)/\\1/')"
            }
        - position: 9
          shellQuote: false
          valueFrom: '&&'
        - position: 10
          shellQuote: false
          valueFrom: |-
            ${
                return "cat makeGraph.R | R --slave --args"
            }
        - position: 11
          shellQuote: false
          valueFrom: $line
        - position: 12
          shellQuote: false
          valueFrom: |-
            ${
                if (inputs.mate_file_sample) {
                    filepath = inputs.mate_file_sample.path
                    filename = filepath.split("/").pop()
                } else {
                    filepath = inputs.mate_copynumber_file_sample.path
                    filename = filepath.split("/").pop()
                }

                ratio = filename + "_ratio" + ".txt"

                return ratio
            }
        - position: 13
          shellQuote: false
          valueFrom: |-
            ${

                if (inputs.snp_file) {

                    sufix = "_BAF"
                    sufix_ext = ".txt"

                    if (inputs.mate_file_sample) {
                        filepath = inputs.mate_file_sample.path
                        filename = filepath.split("/").pop()
                    } else {
                        filepath = inputs.mate_copynumber_file_sample.path
                        filename = filepath.split("/").pop()
                    }


                    new_filename = filename + sufix + sufix_ext

                    return new_filename
                }
            }
        - position: 114
          shellQuote: false
          valueFrom: |-
            ${ //conversion of file names

                if (inputs.mate_file_control) {
                    if (inputs.mate_file_control.path.split('.').pop() != 'pileup') {
                        com = ''
                        com += '&& mv sample.pileup '

                    }
                }


            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: |-
            ${
                if (inputs.total_memory) {
                    return inputs.total_memory
                } else {
                    return 15000
                }
            }
          coresMin: |-
            ${
                if (inputs.max_threads) {
                    return inputs.max_threads
                } else {
                    return 8
                }
            }
        - class: DockerRequirement
          dockerImageId: caf6947244fa
          dockerPull: 'images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1'
        - class: InitialWorkDirRequirement
          listing:
            - entryname: config.txt
              entry: |-
                ${

                    // The function returns the concatenated line for config file
                    function makeline(content, p1, p2) {
                        if (p2 != null) {
                            if (p2.path != null) {
                                p2 = p2.path
                            }
                            content = content.concat(p1)
                            content = content.concat(" = ")
                            content = content.concat(p2)
                            content = content.concat("\n")
                        }
                        return content

                    }

                    // General section
                    content = "[general]\n\n"
                    content = makeline(content, "BedGraphOutput", inputs.bed_graph_output)
                    content = content.concat("bedtools = /opt/bedtools2/bin/bedtools\n")
                    content = makeline(content, "chrLenFile", inputs.chr_len)
                    content = makeline(content, "breakPointThreshold", inputs.break_point_threshold)
                    content = makeline(content, "breakPointType", inputs.break_point_type)
                    content = makeline(content, "chrFiles", ".")
                    content = makeline(content, "coefficientOfVariation", inputs.coeff_var)

                    if (inputs.capture_regions) {
                        content = content.concat("window = 0\n")
                    } else {
                        content = makeline(content, "window", inputs.window)
                    }

                    content = makeline(content, "contamination", inputs.contamination)
                    content = makeline(content, "contaminationAdjustment", inputs.contamination_adjustment)
                    content = makeline(content, "degree", inputs.degree)
                    content = makeline(content, "forceGCcontentNormalization", inputs.force_GC_content_normalization)
                    content = makeline(content, "GCcontentProfile", inputs.GC_content_profile)
                    content = makeline(content, "gemMappabilityFile", inputs.gem_mappability_file)
                    content = makeline(content, "intercept", inputs.intercept)
                    content = makeline(content, "minCNAlength", inputs.min_CNA_length)
                    content = makeline(content, "minMappabilityPerWindow", inputs.min_map_per_w)
                    content = makeline(content, "minExpectedGC", inputs.min_expected_GC)
                    content = makeline(content, "maxExpectedGC", inputs.max_expected_GC)
                    content = makeline(content, "minimalSubclonePresence", inputs.min_subclone_presence)
                    if (inputs.max_threads) {
                        content = makeline(content, "maxThreads", inputs.max_threads)
                        content = makeline(content, "SambambaThreads", inputs.max_threads)
                    } else {
                        content = content.concat("maxThreads = 8\n")
                        content = content.concat("SambambaThreads = 8\n")
                    }
                    content = makeline(content, "noisyData", inputs.noisy_data)
                    content = makeline(content, "ploidy", inputs.ploidy.toString())
                    content = makeline(content, "printNA", inputs.print_NA)
                    content = makeline(content, "readCountThreshold", inputs.read_cnt_threshold)
                    content = content.concat("sambamba = /opt/sambamba_0.5.9/sambamba_v0.5.9\n")

                    content = content.concat("samtools = /opt/samtools-1.3.1/samtools\n")
                    content = makeline(content, "sex", inputs.sex)
                    content = makeline(content, "step", inputs.step)
                    content = makeline(content, "telocentromeric", inputs.telocentromeric)
                    content = makeline(content, "uniqueMatch", inputs.unique_match)


                    // Sample section

                    content = content.concat("\n[sample]\n\n")
                    content = makeline(content, "mateFile", inputs.mate_file_sample)
                    content = makeline(content, "mateCopyNumberFile", inputs.mate_copynumber_file_sample)
                    content = makeline(content, "miniPileup", inputs.mini_pileup_sample)
                    if (inputs.mate_file_sample) {
                        if (inputs.mate_file_sample.path.split('.').pop() == "gz") {
                            content = makeline(content, "inputFormat", inputs.mate_file_sample.path.split('.').slice(-2, -1)[0])
                        } else {
                            content = makeline(content, "inputFormat", inputs.mate_file_sample.path.split('.').pop())
                        }
                        content = makeline(content, "mateOrientation", inputs.mate_orientation_sample)
                    }


                    // Control section

                    content = content.concat("\n[control]\n\n")
                    content = makeline(content, "mateFile", inputs.mate_file_control)
                    content = makeline(content, "mateCopyNumberFile", inputs.mate_copynumber_file_control)
                    content = makeline(content, "miniPileup", inputs.mini_pileup_control)
                    if (inputs.mate_file_control) {
                        content = makeline(content, "inputFormat", inputs.mate_file_control.path.split('.').pop())
                        content = makeline(content, "mateOrientation", inputs.mate_orientation_sample)
                    }




                    // BAF section

                    content = content.concat("\n[BAF]\n\n")
                    content = makeline(content, "minimalCoveragePerPosition", inputs.minimal_coverage_per_position)
                    content = makeline(content, "minimalQualityPerPosition", inputs.minimal_quality_per_position)
                    content = makeline(content, "shiftInQuality", inputs.shift_in_quality)
                    if (inputs.snp_file) {
                        content = makeline(content, "SNPfile", inputs.snp_file)
                        if (inputs.mate_file_sample) {
                            if ((inputs.mate_file_sample.path.split('.').pop().toUpperCase() != 'PILEUP') &&
                                (inputs.mate_file_sample.path.split('.').slice(-2, -1)[0].toUpperCase() != 'PILEUP')) {
                                content = makeline(content, "makePileup", inputs.snp_file)
                                content = makeline(content, "fastaFile", inputs.reference)
                            }
                        }
                    }

                    // Target section

                    content = content.concat("\n[target]\n\n")
                    content = makeline(content, "captureRegions", inputs.capture_regions)

                    return content
                }
              writable: false
            - entryname: split_fasta.py
              entry: |-
                import sys

                with open(sys.argv[1], "r") as f:
                    fasta = f.readlines()

                ref_lines = {}
                for i in range(0, len(fasta)):
                    if fasta[i][0] == ">":
                        chrom = fasta[i].split()[0].split(">")[1]
                        print("Reading chromosome: " + chrom)
                        ref_lines[chrom] = [fasta[i]]
                    else:
                        ref_lines[chrom].append(fasta[i])

                for chromosome, lines in ref_lines.items():
                    print("Creating " + chromosome + ".fasta")
                    with open(chromosome + ".fasta", "w") as chr_fasta:
                        for line in lines:
                            chr_fasta.write(line)
              writable: false
            - entryname: assess_significance.R
              entry: "#!/usr/bin/env Rscript\n\nlibrary(rtracklayer)\n\nargs <- commandArgs()\n\ndataTable <-read.table(args[5], header=TRUE);\nratio<-data.frame(dataTable)\n\ndataTable <- read.table(args[4], header=FALSE)\ncnvs<- data.frame(dataTable) \n\nratio$Ratio[which(ratio$Ratio==-1)]=NA\n\ncnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))  \nratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio)\n\noverlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)\nnormals <- setdiff(ratio.bed,cnvs.bed)\nnormals <- subsetByOverlaps(ratio.bed,normals)\n\n#mu <- mean(score(normals),na.rm=TRUE)\n#sigma<- sd(score(normals),na.rm=TRUE)\n\n#hist(score(normals),n=500,xlim=c(0,2))\n#hist(log(score(normals)),n=500,xlim=c(-1,1))\n\n#shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])\n#qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))\n#qqline(score(normals)[which(!is.na(score(normals)))], col = 2)\n\n#shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])\n#qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))\n#qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)\n\nnumberOfCol=length(cnvs)\n\nfor (i in c(1:length(cnvs[,1]))) {\n  values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))\n  #wilcox.test(values,mu=mu)\n  W <- function(values,normals){resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)\n\tif(class(resultw)==\"try-error\") return(list(\"statistic\"=NA,\"parameter\"=NA,\"p.value\"=NA,\"null.value\"=NA,\"alternative\"=NA,\"method\"=NA,\"data.name\"=NA)) else resultw}\n  KS <- function(values,normals){resultks <- try(ks.test(values,score(normals)), silent = TRUE)\n\tif(class(resultks)==\"try-error\") return(list(\"statistic\"=NA,\"p.value\"=NA,\"alternative\"=NA,\"method\"=NA,\"data.name\"=NA)) else resultks}\n  #resultks <- try(KS <- ks.test(values,score(normals)), silent = TRUE)\n  #\tif(class(resultks)==\"try-error\") NA) else resultks\n  cnvs[i,numberOfCol+1]=W(values,normals)$p.value\n  cnvs[i,numberOfCol+2]=KS(values,normals)$p.value\n  }\n\nif (numberOfCol==5) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nif (numberOfCol==7) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"genotype\",\"uncertainty\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nif (numberOfCol==9) {\n  names(cnvs)=c(\"chr\",\"start\",\"end\",\"copy number\",\"status\",\"genotype\",\"uncertainty\",\"somatic/germline\",\"precentageOfGermline\",\"WilcoxonRankSumTestPvalue\",\"KolmogorovSmirnovPvalue\")  \n}\nwrite.table(cnvs, file=paste(args[4],\".p.value.txt\",sep=\"\"),sep=\"\\t\",quote=F,row.names=F)"
              writable: false
            - entryname: makeGraph.R
              entry: "#!/usr/bin/env Rscript\n\nargs <- commandArgs()\n\ndataTable <-read.table(args[5], header=TRUE);\n\nratio<-data.frame(dataTable)\nploidy <- type.convert(args[4])\n\n\npng(filename = paste(args[5],\".log2.png\",sep = \"\"), width = 1180, height = 1180,\n    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\nplot(1:10)\nop <- par(mfrow = c(5,5))\n\nfor (i in c(1:22,'X','Y')) {\n\ttt <- which(ratio$Chromosome==i)\n\tif (length(tt)>0) {\n\t plot(ratio$Start[tt],log2(ratio$Ratio[tt]),xlab = paste (\"position, chr\",i),ylab = \"normalized copy number profile (log2)\",pch = \".\",col = colors()[88])\n\t tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )\n\t points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = \".\",col = colors()[136])\n\t\n\t\n\ttt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)\n\t points(ratio$Start[tt],log2(ratio$Ratio[tt]),pch = \".\",col = colors()[461])\n\t tt <- which(ratio$Chromosome==i)\n\t \n\t #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:\n\t #points(ratio$Start[tt],log2(ratio$CopyNumber[tt]/ploidy), pch = \".\", col = colors()[24],cex=4)\n\t \n\t}\n\ttt <- which(ratio$Chromosome==i)\n\t\n\t#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:\n\tpoints(ratio$Start[tt],log2(ratio$MedianRatio[tt]), pch = \".\", col = colors()[463],cex=4)\n\t\n}\n\ndev.off()\n\n\npng(filename = paste(args[5],\".png\",sep = \"\"), width = 1180, height = 1180,\n    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\nplot(1:10)\nop <- par(mfrow = c(5,5))\n\nmaxLevelToPlot <- 3\nfor (i in c(1:length(ratio$Ratio))) {\n\tif (ratio$Ratio[i]>maxLevelToPlot) {\n\t\tratio$Ratio[i]=maxLevelToPlot;\n\t}\n}\n\n\nfor (i in c(1:22,'X','Y')) {\n\ttt <- which(ratio$Chromosome==i)\n\tif (length(tt)>0) {\n\t plot(ratio$Start[tt],ratio$Ratio[tt]*ploidy,ylim = c(0,maxLevelToPlot*ploidy),xlab = paste (\"position, chr\",i),ylab = \"normalized copy number profile\",pch = \".\",col = colors()[88])\n\t tt <- which(ratio$Chromosome==i  & ratio$CopyNumber>ploidy )\n\t points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[136])\n\t\n\ttt <- which(ratio$Chromosome==i  & ratio$Ratio==maxLevelToPlot & ratio$CopyNumber>ploidy)\t\n\tpoints(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[136],cex=4)\n\t \n\ttt <- which(ratio$Chromosome==i  & ratio$CopyNumber<ploidy & ratio$CopyNumber!= -1)\n\t points(ratio$Start[tt],ratio$Ratio[tt]*ploidy,pch = \".\",col = colors()[461])\n\t tt <- which(ratio$Chromosome==i)\n\t \n\t #UNCOMMENT HERE TO SEE THE PREDICTED COPY NUMBER LEVEL:\n\t #points(ratio$Start[tt],ratio$CopyNumber[tt], pch = \".\", col = colors()[24],cex=4)\n\t \n\t}\n\ttt <- which(ratio$Chromosome==i)\n\t\n\t#UNCOMMENT HERE TO SEE THE EVALUATED MEDIAN LEVEL PER SEGMENT:\n\tpoints(ratio$Start[tt],ratio$MedianRatio[tt]*ploidy, pch = \".\", col = colors()[463],cex=4)\n\t\n}\n\ndev.off()\n\n\n\n\nif (length(args)>5) {\n\tdataTable <-read.table(args[6], header=TRUE);\n\tBAF<-data.frame(dataTable)\n\n\tpng(filename = paste(args[6],\".png\",sep = \"\"), width = 1180, height = 1180,\n\t    units = \"px\", pointsize = 20, bg = \"white\", res = NA)\n\tplot(1:10)\n\top <- par(mfrow = c(5,5))\n\n\tfor (i in c(1:22,'X','Y')) {\n\t    tt <- which(BAF$Chromosome==i)\n\t    if (length(tt)>0){\n\t\tlBAF <-BAF[tt,]\n\t\tplot(lBAF$Position,lBAF$BAF,ylim = c(-0.1,1.1),xlab = paste (\"position, chr\",i),ylab = \"BAF\",pch = \".\",col = colors()[1])\n\n\t\ttt <- which(lBAF$A==0.5)\t\t\n\t\tpoints(lBAF$Position[tt],lBAF$BAF[tt],pch = \".\",col = colors()[92])\n\t\ttt <- which(lBAF$A!=0.5 & lBAF$A>=0)\n\t\tpoints(lBAF$Position[tt],lBAF$BAF[tt],pch = \".\",col = colors()[62])\n\t\ttt <- 1\n\t\tpres <- 1\n\n\t\tif (length(lBAF$A)>4) {\n\t\t\tfor (j in c(2:(length(lBAF$A)-pres-1))) {\n\t\t\t\tif (lBAF$A[j]==lBAF$A[j+pres]) {\t\n\t\t\t\t\ttt[length(tt)+1] <- j \n\t\t\t\t}\n\t\t\t}\n\t\t\tpoints(lBAF$Position[tt],lBAF$A[tt],pch = \".\",col = colors()[24],cex=4)\n\t\t\tpoints(lBAF$Position[tt],lBAF$B[tt],pch = \".\",col = colors()[24],cex=4)\t\n\t\t}\n\n\t\ttt <- 1\n\t\tpres <- 1\n\t\tif (length(lBAF$FittedA)>4) {\n\t\t\tfor (j in c(2:(length(lBAF$FittedA)-pres-1))) {\n\t\t\t\tif (lBAF$FittedA[j]==lBAF$FittedA[j+pres]) {\t\n\t\t\t\t\ttt[length(tt)+1] <- j \n\t\t\t\t}\n\t\t\t}\n\t\t\tpoints(lBAF$Position[tt],lBAF$FittedA[tt],pch = \".\",col = colors()[463],cex=4)\n\t\t\tpoints(lBAF$Position[tt],lBAF$FittedB[tt],pch = \".\",col = colors()[463],cex=4)\t\n\t\t}\n\n\t   }\n\n\t}\n\tdev.off()\n\n}"
              writable: false
        - class: InlineJavascriptRequirement
          expressionLib:
            - |-
              var updateMetadata = function(file, key, value) {
                  file['metadata'][key] = value;
                  return file;
              };


              var setMetadata = function(file, metadata) {
                  if (!('metadata' in file)) {
                      file['metadata'] = {}
                  }
                  for (var key in metadata) {
                      file['metadata'][key] = metadata[key];
                  }
                  return file
              };

              var inheritMetadata = function(o1, o2) {
                  var commonMetadata = {};
                  if (!Array.isArray(o2)) {
                      o2 = [o2]
                  }
                  for (var i = 0; i < o2.length; i++) {
                      var example = o2[i]['metadata'];
                      for (var key in example) {
                          if (i == 0)
                              commonMetadata[key] = example[key];
                          else {
                              if (!(commonMetadata[key] == example[key])) {
                                  delete commonMetadata[key]
                              }
                          }
                      }
                  }
                  if (!Array.isArray(o1)) {
                      o1 = setMetadata(o1, commonMetadata)
                  } else {
                      for (var i = 0; i < o1.length; i++) {
                          o1[i] = setMetadata(o1[i], commonMetadata)
                      }
                  }
                  return o1;
              };

              var toArray = function(file) {
                  return [].concat(file);
              };

              var groupBy = function(files, key) {
                  var groupedFiles = [];
                  var tempDict = {};
                  for (var i = 0; i < files.length; i++) {
                      var value = files[i]['metadata'][key];
                      if (value in tempDict)
                          tempDict[value].push(files[i]);
                      else tempDict[value] = [files[i]];
                  }
                  for (var key in tempDict) {
                      groupedFiles.push(tempDict[key]);
                  }
                  return groupedFiles;
              };

              var orderBy = function(files, key, order) {
                  var compareFunction = function(a, b) {
                      if (a['metadata'][key].constructor === Number) {
                          return a['metadata'][key] - b['metadata'][key];
                      } else {
                          var nameA = a['metadata'][key].toUpperCase();
                          var nameB = b['metadata'][key].toUpperCase();
                          if (nameA < nameB) {
                              return -1;
                          }
                          if (nameA > nameB) {
                              return 1;
                          }
                          return 0;
                      }
                  };

                  files = files.sort(compareFunction);
                  if (order == undefined || order == "asc")
                      return files;
                  else
                      return files.reverse();
              };
      successCodes:
        - 0
      temporaryFailCodes:
        - 1
      'sbg:appVersion':
        - v1.0
      'sbg:categories':
        - Copy-Number-Analysis
      'sbg:cmdPreview': >-
        python split_fasta.py pat/sbgpro/path/to/reference.ext &&
        /opt/controlfreec/FREEC-11.5/src/freec -conf config.txt && cat
        assess_significance.R | R --slave --args chr_19.noDup0.pileup.gz_CNVs
        chr_19.noDup0.pileup.gz_ratio.txt && line=$(cat *info.txt | grep
        Output_Ploidy | sed -E 's/.+([0-9]+)/\1/') && cat makeGraph.R | R
        --slave --args $line chr_19.noDup0.pileup.gz_ratio.txt
        chr_19.noDup0.pileup.gz_BAF.txt
      'sbg:content_hash': a8c609371987909ebbcad810349a09839fd48c2a7d4481c277e28cce897ed68cd
      'sbg:contributors':
        - brownm28
      'sbg:copyOf': bogdang/controlfreec/control-freec-11-6/1
      'sbg:createdBy': brownm28
      'sbg:createdOn': 1569186036
      'sbg:id': brownm28/mb-controlfreec-troubleshoot/control-freec-11-6-sbg/0
      'sbg:image_url': null
      'sbg:latestRevision': 0
      'sbg:license': GNU General Public License v3.0 only
      'sbg:links':
        - id: 'http://bioinfo-out.curie.fr/projects/freec/'
          label: Homepage
        - id: 'https://github.com/BoevaLab/FREEC'
          label: GitHub
        - id: >-
            http://bioinformatics.oxfordjournals.org/content/early/2011/12/05/bioinformatics.btr670
          label: Paper
      'sbg:modifiedBy': brownm28
      'sbg:modifiedOn': 1569186036
      'sbg:project': brownm28/mb-controlfreec-troubleshoot
      'sbg:projectName': MB ControlFREEC Troubleshoot
      'sbg:publisher': sbg
      'sbg:revision': 0
      'sbg:revisionNotes': Copy of bogdang/controlfreec/control-freec-11-6/1
      'sbg:revisionsInfo':
        - 'sbg:modifiedBy': brownm28
          'sbg:modifiedOn': 1569186036
          'sbg:revision': 0
          'sbg:revisionNotes': Copy of bogdang/controlfreec/control-freec-11-6/1
      'sbg:sbgMaintained': false
      'sbg:toolAuthor': Bioinformatics Laboratory of Institut Curie
      'sbg:toolkit': Control-FREEC
      'sbg:toolkitVersion': '11.6'
      'sbg:validationErrors': []
    label: Control-FREEC 11.6
    'sbg:x': 1018.6522216796875
    'sbg:y': 560.5
  - id: controlfreec_normal_mini_pileup
    in:
      - id: input_reads
        source: samtools_normal_cram2bam/bam_file
      - id: reference
        source: indexed_reference_fasta
      - id: snp_vcf
        source: gatk_filter_germline/filtered_pass_vcf
      - id: threads
        valueFrom: '${return 16}'
    out:
      - id: pileup
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: controlfreec_mini_pileup
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: input_reads
          type: File
          secondaryFiles:
            - ^.bai
        - id: reference
          type: File
          secondaryFiles:
            - .fai
        - id: snp_vcf
          type: File
          doc: Germline vcf with sites to filter pielup on
        - default: 16
          id: threads
          type: int?
      outputs:
        - id: pileup
          type: File
          outputBinding:
            glob: $(inputs.input_reads.nameroot).miniPileup
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            set -eo pipefail

            zcat $(inputs.snp_vcf.path) | grep -v "#" | awk {'printf
            ("%s\t%s\t%s\t%s\t%s\n", $1,$2-1,$2,$4,$5)'} > snps.bed

            /opt/sambamba_0.5.9/sambamba_v0.5.9 mpileup -t $(inputs.threads) -o
            $(inputs.input_reads.nameroot).miniPileup $(inputs.input_reads.path)
            --samtools -f $(inputs.reference.path) -d 8000 -Q 0 -q 1 -l snps.bed
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 10000
          coresMin: $(inputs.threads)
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1'
        - class: InlineJavascriptRequirement
    'sbg:x': 596.5790405273438
    'sbg:y': 565
  - id: controlfreec_tumor_mini_pileup
    in:
      - id: input_reads
        source: samtools_tumor_cram2bam/bam_file
      - id: reference
        source: indexed_reference_fasta
      - id: snp_vcf
        source: gatk_filter_germline/filtered_pass_vcf
      - id: threads
        valueFrom: '${return 16}'
    out:
      - id: pileup
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: controlfreec_mini_pileup
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: input_reads
          type: File
          secondaryFiles:
            - ^.bai
        - id: reference
          type: File
          secondaryFiles:
            - .fai
        - id: snp_vcf
          type: File
          doc: Germline vcf with sites to filter pielup on
        - default: 16
          id: threads
          type: int?
      outputs:
        - id: pileup
          type: File
          outputBinding:
            glob: $(inputs.input_reads.nameroot).miniPileup
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            set -eo pipefail

            zcat $(inputs.snp_vcf.path) | grep -v "#" | awk {'printf
            ("%s\t%s\t%s\t%s\t%s\n", $1,$2-1,$2,$4,$5)'} > snps.bed

            /opt/sambamba_0.5.9/sambamba_v0.5.9 mpileup -t $(inputs.threads) -o
            $(inputs.input_reads.nameroot).miniPileup $(inputs.input_reads.path)
            --samtools -f $(inputs.reference.path) -d 8000 -Q 0 -q 1 -l snps.bed
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 10000
          coresMin: $(inputs.threads)
        - class: DockerRequirement
          dockerPull: 'images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1'
        - class: InlineJavascriptRequirement
    'sbg:x': 596.5790405273438
    'sbg:y': 430
  - id: convert_ratio_to_seg
    in:
      - id: ctrlfreec_ratio
        source: control_free_c/ratio
      - id: output_basename
        source: output_basename
      - id: reference_fai
        source: reference_fai
      - id: sample_name
        source: input_tumor_name
    out:
      - id: ctrlfreec_ratio2seg
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: ubuntu_ratio2seg
      baseCommand:
        - python
        - '-c'
      inputs:
        - id: ctrlfreec_ratio
          type: File
        - id: output_basename
          type: string
        - id: reference_fai
          type: File
        - id: sample_name
          type: string
      outputs:
        - id: ctrlfreec_ratio2seg
          type: File
          outputBinding:
            glob: '*.seg'
      arguments:
        - position: 0
          valueFrom: |

            import math

            fai = open("$(inputs.reference_fai.path)")
            h = {}
            for line in fai:
                f = line.split("\t")
                if f == "chrM":
                    break
                f[0] = f[0].replace("chr","")
                h[f[0]] = f[1]
            fai.close()

            smp = "$(inputs.sample_name)"

            ratio_file = open("$(inputs.ctrlfreec_ratio.path)")
            out = open("$(inputs.output_basename).controlfreec.seg", "w")
            out.write("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n") 
            head = next(ratio_file)
            count = 0
            for line in ratio_file:
                data = line.rstrip("\n").split("\t")
                (chrom, pos, ratio, meanRatio) = (data[0], data[1], data[2], data[3])
                if meanRatio == "-1":
                    continue
                count += 1
                if count == 1:
                    start = pos
                    seg_ratio = meanRatio
                    on_chr = chrom
                else:
                    if chrom != on_chr:
                        out.write("\t".join((smp, "chr" + on_chr, start, h[on_chr], str(count))) + "\t")
                        if seg_ratio != "0":
                            out.write(str(math.log(float(seg_ratio), 2)) + "\n")
                        else:
                            out.write(str(math.log(float(seg_ratio) + 1, 2)) + "\n")
                        start = pos
                        seg_ratio = meanRatio
                        on_chr = chrom
                        count = 1
                    elif meanRatio != seg_ratio:
                        out.write("\t".join((smp, "chr" + chrom, start, str(int(pos)-1), str(count))) + "\t")
                        if seg_ratio != "0":
                            out.write(str(math.log(float(seg_ratio), 2)) + "\n")
                        else:
                            out.write(str(math.log(float(seg_ratio) + 1, 2)) + "\n")
                        start = pos
                        seg_ratio = meanRatio
                        count = 1
            ratio_file.close()
            out.close()
      requirements:
        - class: ResourceRequirement
          ramMin: 1000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'kfdrc/python:2.7.13'
        - class: InlineJavascriptRequirement
    'sbg:x': 1590.5213623046875
    'sbg:y': 1019
  - id: gatk_filter_germline
    in:
      - id: input_vcf
        source: b_allele
      - id: output_basename
        source: output_basename
      - id: reference_fasta
        source: indexed_reference_fasta
    out:
      - id: filtered_pass_vcf
      - id: filtered_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: kfdrc_gatk_variantfiltration
      baseCommand:
        - /gatk
      inputs:
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: output_basename
          type: string
        - id: reference_fasta
          type: File
      outputs:
        - id: filtered_pass_vcf
          type: File
          outputBinding:
            glob: '*.gatk.hardfiltered.PASS.vcf.gz'
          secondaryFiles:
            - .tbi
        - id: filtered_vcf
          type: File
          outputBinding:
            glob: '*.gatk.hardfiltered.vcf.gz'
          secondaryFiles:
            - .tbi
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --java-options "-Xmx7000m" SelectVariants --exclude-filtered TRUE
            -select-type SNP -V $(inputs.input_vcf.path) -O snp_pass.vcf.gz

            /gatk VariantFiltration --java-options "-Xmx7000m" --filter-name
            GATK_QD --filter-expression "QD < 2.0" --filter-name GATK_FS
            --filter-expression "FS > 60.0" --filter-name GATK_MQ
            --filter-expression "MQ < 40.0" --filter-name GATK_MQRankSum
            --filter-expression "MQRankSum < -12.5" --filter-name
            GATK_ReadPosRankSum --filter-expression "ReadPosRankSum < -8.0"
            --filter-name KFDRC_DP10 --filter-expression "DP < 10" -V
            snp_pass.vcf.gz -O
            $(inputs.output_basename).gatk.hardfiltered.vcf.gz

            /gatk SelectVariants --exclude-filtered TRUE -V
            $(inputs.output_basename).gatk.hardfiltered.vcf.gz -O
            $(inputs.output_basename).gatk.hardfiltered.PASS.vcf.gz
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 2
          coresMax: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.1.1.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 296.890625
    'sbg:y': 788.5
  - id: rename_outputs
    in:
      - id: input_files
        source:
          - control_free_c/cnvs_pvalue
          - control_free_c/config_script
          - control_free_c/ratio
          - control_free_c/sample_BAF
          - control_free_c/info_txt
      - id: input_pngs
        source:
          - control_free_c/pngs
      - id: output_basename
        source: output_basename
    out:
      - id: ctrlfreec_baf
      - id: ctrlfreec_bam_ratio
      - id: ctrlfreec_cnv
      - id: ctrlfreec_config
      - id: ctrlfreec_info
      - id: ctrlfreec_normal_cpn
      - id: ctrlfreec_pngs
      - id: ctrlfreec_pval
      - id: ctrlfreec_tumor_cpn
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: ubuntu_rename_cf_outputs
      baseCommand: []
      inputs:
        - id: input_files
          type: 'File[]'
        - id: input_pngs
          type: 'File[]'
        - id: output_basename
          type: string
      outputs:
        - id: ctrlfreec_baf
          type: File
          outputBinding:
            glob: '*.BAF.txt'
        - id: ctrlfreec_bam_ratio
          type: File
          outputBinding:
            glob: '*.ratio.txt'
        - id: ctrlfreec_cnv
          type: File
          outputBinding:
            glob: '*.CNVs'
        - id: ctrlfreec_config
          type: File
          outputBinding:
            glob: '*.config.txt'
        - id: ctrlfreec_info
          type: File
          outputBinding:
            glob: '*.info.txt'
        - id: ctrlfreec_normal_cpn
          type: File
          outputBinding:
            glob: '*.control.cpn'
        - id: ctrlfreec_pngs
          type: 'File[]'
          outputBinding:
            glob: '*.png'
        - id: ctrlfreec_pval
          type: File
          outputBinding:
            glob: '*.CNVs.p.value.txt'
        - id: ctrlfreec_tumor_cpn
          type: File
          outputBinding:
            glob: '*.sample.cpn'
      doc: Rename contrfreeec outputs
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: |-
            ${
                var cmd = "";
                for (var i=0; i < inputs.input_files.length; i++){
                  var basename = inputs.input_files[i].basename;
                  var fname = basename.replace("bam_", "");
                  var parts = fname.split(".");
                  parts.shift();
                  var check = fname.substr(fname.length - 10);
                  if (check == "config.txt") {
                      cmd += "cp " + inputs.input_files[i].path + " " + inputs.output_basename + ".controlfreec.config.txt;";
                  } else {
                  fname = inputs.output_basename + ".controlfreec." + parts.join(".");
                  cmd += " cp " + inputs.input_files[i].path + " " + fname + ";";
                  
                  }
              for (var j=0; j < inputs.input_pngs.length; j++){
                  var nameroot = inputs.input_pngs[j].nameroot;
                  var fname = nameroot.replace("bam_", "");
                  fname = fname.replace(".txt", "");
                  var parts = fname.split(".");
                  parts.shift();
                  fname = inputs.output_basename + ".controlfreec." + parts.join(".") + ".png";
                  cmd += " cp " + inputs.input_pngs[j].path + " " + fname + ";";
                  }
                }
                return cmd;
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 1000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'ubuntu:18.04'
        - class: InlineJavascriptRequirement
    'sbg:x': 1590.5213623046875
    'sbg:y': 835
  - id: run_theta2
    in:
      - id: interval_count
        source: cnvkit_export_theta2/call_interval_count
      - id: min_frac
        source: min_theta2_frac
      - id: normal_snp
        source: cnvkit_export_theta2/call_normal_snp
      - id: output_basename
        source: output_basename
      - id: tumor_snp
        source: cnvkit_export_theta2/call_tumor_snp
    out:
      - id: best_results
      - id: n2_graph
      - id: n2_results
      - id: n2_withBounds
      - id: n3_graph
      - id: n3_results
      - id: n3_withBounds
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: theta2
      baseCommand:
        - /THetA/bin/RunTHetA
      inputs:
        - id: interval_count
          type: File
        - default: 0.05
          id: min_frac
          type: float?
          doc: >-
            Minimum fraction of genome with copy umber alterations.  Default is
            0.05
        - id: normal_snp
          type: File
        - id: output_basename
          type: string
        - id: tumor_snp
          type: File
      outputs:
        - id: best_results
          type: File
          outputBinding:
            glob: '*.BEST.results'
        - id: n2_graph
          type: File
          outputBinding:
            glob: '*.n2.graph.pdf'
        - id: n2_results
          type: File
          outputBinding:
            glob: '*.n2.results'
        - id: n2_withBounds
          type: File
          outputBinding:
            glob: '*.n2.withBounds'
        - id: n3_graph
          type: File
          outputBinding:
            glob: '*.n3.graph.pdf'
        - id: n3_results
          type: File
          outputBinding:
            glob: '*.n3.results'
        - id: n3_withBounds
          type: File
          outputBinding:
            glob: '*.n3.withBounds'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --TUMOR_FILE $(inputs.tumor_snp.path) --NORMAL_FILE
            $(inputs.normal_snp.path) --OUTPUT_PREFIX $(inputs.output_basename)
            --NUM_PROCESSES 8 --MIN_FRAC $(inputs.min_frac)
            $(inputs.interval_count.path)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 32000
          coresMin: 8
        - class: DockerRequirement
          dockerPull: 'kfdrc/theta2:0.7'
        - class: InlineJavascriptRequirement
    'sbg:x': 1590.5213623046875
    'sbg:y': 630
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
    'sbg:x': 296.890625
    'sbg:y': 660.5
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
    'sbg:x': 296.890625
    'sbg:y': 539.5
hints:
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 4
requirements:
  - class: MultipleInputFeatureRequirement
