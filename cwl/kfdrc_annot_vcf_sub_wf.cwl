cwlVersion: v1.0
class: Workflow
label: kfdrc_annot_sub_wf
doc: This subworkflow normalizes, annotates, and adds soft filters to vcf inputs
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: InlineJavascriptRequirement

inputs:
- id: indexed_reference_fasta
  type: File
  secondaryFiles:
  - .fai
  - ^.dict
- id: input_vcf
  doc: Input vcf to annotate and soft filter
  type: File
  secondaryFiles:
  - .tbi
- id: input_tumor_name
  type: string
- id: input_normal_name
  type: string
- id: add_common_fields
  doc: Set to true if input is a strelka2 vcf that hasn't had common fields added
  type: boolean
  default: false
- id: bcftools_annot_columns
  doc: csv string of columns from annotation to port into the input vcf, i.e INFO/AF
  type: string
- id: bcftools_annot_vcf
  doc: bgzipped annotation vcf file
  type: File
  secondaryFiles:
  - .tbi
- id: bcftools_public_filter
  doc: Will hard filter final result to create a public version
  type:
  - 'null'
  - string
  default: FILTER="PASS"|INFO/HotSpotAllele=1
- id: gatk_filter_name
  doc: Array of names for each filter tag to add
  type:
    type: array
    items: string
- id: gatk_filter_expression
  doc: |-
    Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues
  type:
    type: array
    items: string
- id: vep_cache
  doc: tar gzipped cache from ensembl/local converted cache
  type: File
- id: vep_ref_build
  doc: Genome ref build used, should line up with cache.
  type:
  - 'null'
  - string
  default: GRCh38
- id: disable_hotspot_annotation
  doc: Disable Hotspot Annotation and skip this task.
  type:
  - 'null'
  - boolean
- id: genomic_hotspots
  doc: |-
    Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots
  type:
  - 'null'
  - type: array
    items: File
- id: protein_snv_hotspots
  doc: |-
    Column-name-containing, tab-delimited file(s) containing protein names and amino acid positions corresponding to hotspots
  type:
  - 'null'
  - type: array
    items: File
- id: protein_indel_hotspots
  doc: |-
    Column-name-containing, tab-delimited file(s) containing protein names and amino acid position ranges corresponding to hotspots
  type:
  - 'null'
  - type: array
    items: File
- id: output_basename
  type: string
- id: tool_name
  type: string
- id: retain_info
  doc: |-
    csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`
  type:
  - 'null'
  - string
- id: retain_fmt
  doc: csv string with FORMAT fields that you want to keep
  type:
  - 'null'
  - string
- id: maf_center
  doc: Sequencing center of variant called
  type:
  - 'null'
  - string
  default: .

outputs:
- id: annotated_protected_vcf
  type: File
  outputSource: hotspots_annotation/hotspots_vcf
- id: annotated_protected_maf
  type: File
  outputSource: kfdrc_vcf2maf_protected/output_maf
- id: annotated_public_vcf
  type: File
  outputSource: hard_filter_vcf/filtered_vcf
- id: annotated_public_maf
  type: File
  outputSource: kfdrc_vcf2maf_public/output_maf

steps:
- id: normalize_vcf
  in:
  - id: indexed_reference_fasta
    source: indexed_reference_fasta
  - id: input_vcf
    source: input_vcf
  - id: output_basename
    source: output_basename
  - id: tool_name
    source: tool_name
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: |-
      Before consensus calling, left align indel calls, break up multi-allelic calls, but leave mnps intact

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 4
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest

    inputs:
    - id: input_vcf
      type: File
      secondaryFiles:
      - .tbi
    - id: indexed_reference_fasta
      type: File
      secondaryFiles:
      - .fai
    - id: output_basename
      type: string
    - id: tool_name
      type: string
    - id: strip_info
      doc: |-
        If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/ANN or INFO/CSQ - check vcf
      type:
      - 'null'
      - string

    outputs:
    - id: normalized_vcf
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.bcf_vt_norm.vcf.gz'

    baseCommand:
    - /bin/bash
    - -c
    arguments:
    - position: 0
      valueFrom: |-
        set -eo pipefail
        VCF=$(inputs.input_vcf.path)
        ${
            var cmd = " >&2 echo checking if strip flag given;";
            if (inputs.strip_info != null){
              cmd += ">&2 echo strip flag given; VCF=stripped.vcf;"
              cmd += "bcftools annotate -x " + inputs.strip_info + " " + inputs.input_vcf.path + " -o $VCF"
              cmd += " || VCF=" + inputs.input_vcf.path
            }else{
              cmd += " >&2 echo no strip flag given"
            }
            return cmd;
        } && bcftools norm -m '-any' $VCF > $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcf && /vt/vt normalize -n -r $(inputs.indexed_reference_fasta.path) $(inputs.output_basename).$(inputs.tool_name).bcf_norm.vcf > $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf && bgzip $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf && tabix $(inputs.output_basename).$(inputs.tool_name).bcf_vt_norm.vcf.gz
      shellQuote: false
    id: normalize_vcf
  out:
  - normalized_vcf
- id: bcftools_strip_info
  in:
  - id: input_vcf
    source: normalize_vcf/normalized_vcf
  - id: output_basename
    source: output_basename
  - id: tool_name
    source: tool_name
  - id: strip_info
    source: bcftools_annot_columns
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: Quick tool to strip info from vcf file before re-annotation

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 4
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest

    inputs:
    - id: input_vcf
      type: File
      secondaryFiles:
      - .tbi
    - id: output_basename
      type: string
    - id: tool_name
      type: string
    - id: strip_info
      doc: |-
        If given, remove previous annotation information based on INFO file, i.e. to strip VEP info, use INFO/ANN
      type:
      - 'null'
      - string

    outputs:
    - id: stripped_vcf
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.vcf.gz'

    baseCommand:
    - /bin/bash
    - -c
    arguments:
    - position: 0
      valueFrom: |-
        set -eo pipefail
        (bcftools annotate -x $(inputs.strip_info) $(inputs.input_vcf.path) -O z  -o $(inputs.output_basename).$(inputs.tool_name).INFO_stripped.vcf.gz && tabix $(inputs.output_basename).$(inputs.tool_name).INFO_stripped.vcf.gz) || (echo "Check errors, likely does not have INFO, trying to pass input instead" >&2; cp $(inputs.input_vcf.path) .; cp $(inputs.input_vcf.secondaryFiles[0].path)  .;)
      shellQuote: false
    id: bcftools_strip_info
  out:
  - stripped_vcf
- id: add_standard_fields
  in:
  - id: strelka2_vcf
    source: bcftools_strip_info/stripped_vcf
  - id: run_tool_flag
    source: add_common_fields
  - id: tumor_name
    source: input_tumor_name
  - id: normal_name
    source: input_normal_name
  - id: output_basename
    source: output_basename
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: |-
      This tool will run a Python script that takes a Strelka2 tumor-normal VCF and adds
          canonical fields like GT and AD

    requirements:
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/add-strelka2-fields:1.0.0
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: $(inputs.cores)
      ramMin: ${ return inputs.ram * 1000 }
    - class: ShellCommandRequirement

    inputs:
    - id: strelka2_vcf
      type: File
      secondaryFiles:
      - .tbi
      inputBinding:
        prefix: --strelka2_vcf
        position: 1
    - id: run_tool_flag
      doc: When run as part of a workflow, can skip
      type:
      - 'null'
      - boolean
      default: false
    - id: tumor_name
      type: string
      inputBinding:
        prefix: --tumor_name
        position: 2
    - id: normal_name
      type: string
      inputBinding:
        prefix: --normal_name
        position: 3
    - id: output_basename
      type: string
      inputBinding:
        prefix: --output_basename
        position: 4
    - id: cores
      type:
      - 'null'
      - int
      default: 4
    - id: ram
      doc: RAM requirement in GB
      type:
      - 'null'
      - int
      default: 3

    outputs:
    - id: output
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.gz'
        outputEval: '$( inputs.run_tool_flag ? self : inputs.strelka2_vcf)'

    baseCommand: []
    arguments:
    - position: 0
      valueFrom: |-
        $(inputs.run_tool_flag ? ">&2 /usr/bin/add_strelka2_fields.py" : ">&2 echo 'User opted to skip adding fields to vcf' && exit 0;")
      shellQuote: false
    id: add_strelka2_fields
  out:
  - output
- id: vep_annotate_vcf
  in:
  - id: reference
    source: indexed_reference_fasta
  - id: input_vcf
    source: add_standard_fields/output
  - id: output_basename
    source: output_basename
  - id: tool_name
    source: tool_name
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
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/vep:r93.7

    inputs:
    - id: reference
      doc: Fasta genome assembly with index
      type: File
      secondaryFiles:
      - .fai
    - id: input_vcf
      type: File
      secondaryFiles:
      - .tbi
    - id: species
      doc: Refer to the cache dir structure to set this
      type:
      - 'null'
      - string
      default: homo_sapiens
    - id: merged
      doc: Set to true if a merged VEP cache is being used
      type:
      - 'null'
      - boolean
      default: false
    - id: use_reg
      doc: Not all caches have the regulatory feature. Set to false for dog especially
      type:
      - 'null'
      - boolean
      default: true
    - id: output_basename
      type: string
    - id: tool_name
      type: string
    - id: cache
      doc: tar gzipped cache from ensembl/local converted cache
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
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.vcf.gz'

    baseCommand:
    - mkdir
    arguments:
    - position: 1
      valueFrom: |-
        $(inputs.species) && tar --use-compress-program="pigz -p 8" -xf $(inputs.cache.path) -C $(inputs.species) && perl /ensembl-vep-release-93.7/vep --af --af_1kg --af_esp --af_gnomad --allele_number --assembly $(inputs.ref_build) --biotype --buffer_size 10000 --cache --cache_version $(inputs.cache_version) --canonical --ccds --check_existing --dir_cache $(inputs.species) --domains --failed 1 --fasta $(inputs.reference.path) --flag_pick_allele --fork 16 --format vcf --gene_phenotype --hgvs --input_file $(inputs.input_vcf.path) --no_escape --no_progress --no_stats ${
          if(inputs.merged){
            return "--merged";
          }
          else{
            return "";
          }
        } --numbers --offline --output_file $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf --pick_order canonical,tsl,biotype,rank,ccds,length --polyphen b --protein --pubmed ${
          if (!inputs.use_reg){
            return "--regulatory"
          }
          else{
            return "";
          }
        } --shift_hgvs 1 --sift b --species $(inputs.species) --symbol --total_length --tsl --uniprot --variant_class --vcf --xref_refseq && /ensembl-vep-release-93.7/htslib/bgzip $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf && /ensembl-vep-release-93.7/htslib/tabix $(inputs.output_basename).$(inputs.tool_name).PASS.vep.vcf.gz
      shellQuote: false
    id: kfdrc-vep-somatic-annotate
  out:
  - output_vcf
- id: bcftools_gnomad_annotate
  in:
  - id: input_vcf
    source: vep_annotate_vcf/output_vcf
  - id: annotation_vcf
    source: bcftools_annot_vcf
  - id: columns
    source: bcftools_annot_columns
  - id: output_basename
    source: output_basename
  - id: tool_name
    source: tool_name
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: Simple tool to annotate a vcf using bcftools and an annotation vcf

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 4
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/vcfutils:latest

    inputs:
    - id: input_vcf
      type: File
      secondaryFiles:
      - .tbi
    - id: annotation_vcf
      doc: bgzipped annotation vcf file
      type: File
      secondaryFiles:
      - .tbi
    - id: columns
      doc: csv string of columns from annotation to port into the input vcf, i.e INFO/AF
      type: string
    - id: threads
      doc: Number of compression/decompression threads
      type:
      - 'null'
      - int
      default: 4
    - id: output_basename
      type: string
    - id: tool_name
      type: string

    outputs:
    - id: bcftools_annotated_vcf
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.vcf.gz'

    baseCommand:
    - bcftools
    - annotate
    arguments:
    - position: 0
      valueFrom: |-
        --annotations $(inputs.annotation_vcf.path) --columns $(inputs.columns) -o $(inputs.output_basename).$(inputs.tool_name).bcf_annotated.vcf.gz -O z --threads $(inputs.threads) $(inputs.input_vcf.path) && tabix $(inputs.output_basename).$(inputs.tool_name).bcf_annotated.vcf.gz
      shellQuote: false
    id: bcftools_annotate_vcf
  out:
  - bcftools_annotated_vcf
- id: gatk_add_soft_filter
  in:
  - id: input_vcf
    source: bcftools_gnomad_annotate/bcftools_annotated_vcf
  - id: reference
    source: indexed_reference_fasta
  - id: filter_name
    source: gatk_filter_name
  - id: filter_expression
    source: gatk_filter_expression
  - id: output_basename
    source: output_basename
  - id: tool_name
    source: tool_name
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: Simple tool add a FILTER tag based on criteria

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 4
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/gatk:4.1.1.0

    inputs:
    - id: input_vcf
      type: File
      secondaryFiles:
      - .tbi
    - id: reference
      type: File
      secondaryFiles:
      - ^.dict
      - .fai
    - id: filter_name
      doc: Array of names for each filter tag to add
      type:
        type: array
        items: string
    - id: filter_expression
      doc: |-
        Array of filter expressions to establish criteria to tag variants with. See https://gatk.broadinstitute.org/hc/en-us/articles/360036730071-VariantFiltration for clues
      type:
        type: array
        items: string
    - id: threads
      doc: Number of compression/decompression threads
      type:
      - 'null'
      - int
      default: 4
    - id: output_basename
      type: string
    - id: tool_name
      type: string

    outputs:
    - id: gatk_soft_filtered_vcf
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.vcf.gz'

    baseCommand:
    - /gatk
    - VariantFiltration
    arguments:
    - position: 0
      valueFrom: |-
        -R $(inputs.reference.path) -V $(inputs.input_vcf.path) -O $(inputs.output_basename).$(inputs.tool_name).gatk.soft_filtered.vcf.gz ${
          var args = "";
          for (var i = 0; i < inputs.filter_name.length; i++){
            args += "--filter-name \"" + inputs.filter_name[i] + "\" --filter-expression \"" + inputs.filter_expression[i] + "\" ";
          }
          return args
        }
      shellQuote: false
    id: gatk_variantfilter_vcf
  out:
  - gatk_soft_filtered_vcf
- id: hotspots_annotation
  in:
  - id: input_vcf
    source: gatk_add_soft_filter/gatk_soft_filtered_vcf
  - id: disable_hotspot_annotation
    source: disable_hotspot_annotation
  - id: genomic_hotspots
    source: genomic_hotspots
  - id: protein_snvs
    source: protein_snv_hotspots
  - id: protein_indels
    source: protein_indel_hotspots
  - id: output_basename
    source: output_basename
  run:
    cwlVersion: v1.0
    class: CommandLineTool

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: $(inputs.cores)
      ramMin: ${ return inputs.ram * 1000 }
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/hotspots:0.1.0

    inputs:
    - id: disable_hotspot_annotation
      doc: Disable Hotspot Annotation and skip this task.
      type:
      - 'null'
      - boolean
    - id: input_vcf
      doc: VCF file to annotate hotspots.
      type: File
      secondaryFiles:
      - .tbi
      inputBinding:
        prefix: --vcf
        position: 99
        shellQuote: false
    - id: genomic_hotspots
      doc: |-
        Tab-delimited BED formatted file(s) containing hg38 genomic positions corresponding to hotspots
      type:
      - 'null'
      - type: array
        items: File
      inputBinding:
        prefix: --genomic_hotspots
        position: 1
        shellQuote: false
    - id: protein_indels
      doc: |-
        Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos>-<end_aa_pos> for INDEL hotspots
      type:
      - 'null'
      - type: array
        items: File
      inputBinding:
        prefix: --protein_indels
        position: 1
        shellQuote: false
    - id: protein_snvs
      doc: |-
        Column name-labeled, tab-delimited file(s) containing HUGO-formatted protein names and VEP-formatted positions <start_aa_pos> for SNV hotspots
      type:
      - 'null'
      - type: array
        items: File
      inputBinding:
        prefix: --protein_snvs
        position: 1
        shellQuote: false
    - id: output_basename
      doc: String to use as basename for output file
      type:
      - 'null'
      - string
      inputBinding:
        prefix: --output_basename
        position: 10
        shellQuote: false
    - id: ram
      doc: GB of RAM to allocate to this task.
      type:
      - 'null'
      - int
      default: 2
    - id: cores
      doc: CPU cores to allocate to this task.
      type:
      - 'null'
      - int
      default: 1

    outputs:
    - id: hotspots_vcf
      type: File
      secondaryFiles:
      - .tbi
      outputBinding:
        glob: '*.gz'
        outputEval: '$(inputs.disable_hotspot_annotation ? inputs.input_vcf : self)'

    baseCommand: []
    arguments:
    - position: 0
      valueFrom: |-
        $(inputs.disable_hotspot_annotation ? ">&2 echo 'User elected to skip hotspot annotation' && exit 0;" : ">&2 /hotspot.py")
      shellQuote: false
    id: hotspots_annotation
  out:
  - hotspots_vcf
- id: kfdrc_vcf2maf_protected
  in:
  - id: reference
    source: indexed_reference_fasta
  - id: input_vcf
    source: hotspots_annotation/hotspots_vcf
  - id: output_basename
    source: output_basename
  - id: tumor_id
    source: input_tumor_name
  - id: normal_id
    source: input_normal_name
  - id: tool_name
    source: tool_name
  - id: retain_info
    source: retain_info
  - id: retain_fmt
    source: retain_fmt
  - id: maf_center
    source: maf_center
  run:
    cwlVersion: v1.0
    class: CommandLineTool

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 2
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3

    inputs:
    - id: reference
      doc: Fasta genome assembly with index
      type: File
      secondaryFiles:
      - .fai
    - id: input_vcf
      doc: VEP annotated vcf file.
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
    - id: ref_build
      doc: Genome ref build used, should line up with cache.
      type:
      - 'null'
      - string
      default: GRCh38
    - id: retain_info
      doc: |-
        csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`
      type:
      - 'null'
      - string
    - id: retain_fmt
      doc: csv string with FORMAT fields that you want to keep
      type:
      - 'null'
      - string
    - id: custom_enst
      doc: Use a file with ens tx IDs for each gene to override VEP PICK
      type:
      - 'null'
      - File
    - id: maf_center
      doc: Sequencing center of variant called
      type:
      - 'null'
      - string
      default: .

    outputs:
    - id: output_maf
      type: File
      outputBinding:
        glob: '*.maf'

    baseCommand:
    - gunzip
    - -c
    arguments:
    - position: 1
      valueFrom: |-
        $(inputs.input_vcf.path) > input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf input_file.vcf --output-maf $(inputs.output_basename).$(inputs.tool_name).vep.maf --tumor-id $(inputs.tumor_id) --normal-id $(inputs.normal_id) --ncbi-build $(inputs.ref_build) --ref-fasta $(inputs.reference.path) ${
          if(inputs.maf_center){
            return "--maf-center \"" + inputs.maf_center + "\""
          }
          else{
            return "";
          }
        } ${
          if(inputs.retain_info){
            return "--retain-info " + inputs.retain_info;
          }
          else{
            return "";
          }
        } ${
          if(inputs.retain_fmt){
            return "--retain-fmt " + inputs.retain_fmt;
          }
          else{
            return "";
          }
        } ${
          if(inputs.custom_enst){
            return "--custom-enst " + inputs.custom_enst.path;
          }
          else{
            return "";
          }
        }
      shellQuote: false
    id: kf-mskcc-vcf2maf
  out:
  - output_maf
- id: hard_filter_vcf
  in:
  - id: input_vcf
    source: hotspots_annotation/hotspots_vcf
  - id: include_expression
    source: bcftools_public_filter
  - id: output_basename
    source: output_basename
  run:
    cwlVersion: v1.0
    class: CommandLineTool
    doc: |-
      More generic tool to take in an include expression and optionally an exclude expresssion to filter a vcf

    requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/bvcftools:latest
    - class: ResourceRequirement
      coresMin: 1
      ramMin: 1000
    - class: InlineJavascriptRequirement

    inputs:
    - id: input_vcf
      type: File
    - id: include_expression
      type:
      - 'null'
      - string
    - id: exclude_expression
      type:
      - 'null'
      - string
    - id: output_basename
      type:
      - 'null'
      - string

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
    id: bcftools_filter_vcf
  out:
  - filtered_vcf
- id: kfdrc_vcf2maf_public
  in:
  - id: reference
    source: indexed_reference_fasta
  - id: input_vcf
    source: hard_filter_vcf/filtered_vcf
  - id: output_basename
    source: output_basename
  - id: tumor_id
    source: input_tumor_name
  - id: normal_id
    source: input_normal_name
  - id: tool_name
    source: tool_name
  - id: retain_info
    source: retain_info
  - id: retain_fmt
    source: retain_fmt
  - id: maf_center
    source: maf_center
  run:
    cwlVersion: v1.0
    class: CommandLineTool

    requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      coresMin: 2
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: pgc-images.sbgenomics.com/d3b-bixu/kf_vcf2maf:v1.0.3

    inputs:
    - id: reference
      doc: Fasta genome assembly with index
      type: File
      secondaryFiles:
      - .fai
    - id: input_vcf
      doc: VEP annotated vcf file.
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
    - id: ref_build
      doc: Genome ref build used, should line up with cache.
      type:
      - 'null'
      - string
      default: GRCh38
    - id: retain_info
      doc: |-
        csv string with INFO fields that you want to keep, i.e. for consensus `MQ,MQ0,CAL,Hotspot`
      type:
      - 'null'
      - string
    - id: retain_fmt
      doc: csv string with FORMAT fields that you want to keep
      type:
      - 'null'
      - string
    - id: custom_enst
      doc: Use a file with ens tx IDs for each gene to override VEP PICK
      type:
      - 'null'
      - File
    - id: maf_center
      doc: Sequencing center of variant called
      type:
      - 'null'
      - string
      default: .

    outputs:
    - id: output_maf
      type: File
      outputBinding:
        glob: '*.maf'

    baseCommand:
    - gunzip
    - -c
    arguments:
    - position: 1
      valueFrom: |-
        $(inputs.input_vcf.path) > input_file.vcf && perl /vcf2maf/vcf2maf.pl --input-vcf input_file.vcf --output-maf $(inputs.output_basename).$(inputs.tool_name).vep.maf --tumor-id $(inputs.tumor_id) --normal-id $(inputs.normal_id) --ncbi-build $(inputs.ref_build) --ref-fasta $(inputs.reference.path) ${
          if(inputs.maf_center){
            return "--maf-center \"" + inputs.maf_center + "\""
          }
          else{
            return "";
          }
        } ${
          if(inputs.retain_info){
            return "--retain-info " + inputs.retain_info;
          }
          else{
            return "";
          }
        } ${
          if(inputs.retain_fmt){
            return "--retain-fmt " + inputs.retain_fmt;
          }
          else{
            return "";
          }
        } ${
          if(inputs.custom_enst){
            return "--custom-enst " + inputs.custom_enst.path;
          }
          else{
            return "";
          }
        }
      shellQuote: false
    id: kf-mskcc-vcf2maf
  out:
  - output_maf
id: |-
  https://cavatica-api.sbgenomics.com/v2/apps/kfdrc-harmonization/kf-reference-pipeline/kfdrc_annot_sub_wf/0/raw/
sbg:appVersion:
- v1.0
sbg:content_hash: a178233107186744432791c113a43ccc550703a567d71a0ea5f8a5eb37abb7818
sbg:contributors:
- brownm28
sbg:createdBy: brownm28
sbg:createdOn: 1650476054
sbg:id: kfdrc-harmonization/kf-reference-pipeline/kfdrc_annot_sub_wf/0
sbg:image_url: |-
  https://cavatica.sbgenomics.com/ns/brood/images/kfdrc-harmonization/kf-reference-pipeline/kfdrc_annot_sub_wf/0.png
sbg:latestRevision: 0
sbg:modifiedBy: brownm28
sbg:modifiedOn: 1650476054
sbg:project: kfdrc-harmonization/kf-reference-pipeline
sbg:projectName: kf-reference-pipeline
sbg:publisher: sbg
sbg:revision: 0
sbg:revisionNotes: |-
  Uploaded using sbpack v2021.10.07. 
  Source: 
  repo: https://github.com/kids-first/kf-somatic-workflow
  file: sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
  commit: 1.1.1-281-gf1c453c
sbg:revisionsInfo:
- sbg:modifiedBy: brownm28
  sbg:modifiedOn: 1650476054
  sbg:revision: 0
  sbg:revisionNotes: |-
    Uploaded using sbpack v2021.10.07. 
    Source: 
    repo: https://github.com/kids-first/kf-somatic-workflow
    file: sub_workflows/kfdrc_annot_vcf_sub_wf.cwl
    commit: 1.1.1-281-gf1c453c
sbg:sbgMaintained: false
sbg:validationErrors: []
sbg:workflowLanguage: CWL
