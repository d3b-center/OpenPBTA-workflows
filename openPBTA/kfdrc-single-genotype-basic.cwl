class: Workflow
cwlVersion: v1.0
id: kfdrc-harmonization/sd-bhjxbdqk-01/kfdrc-single-genotype-basic/0
label: kfdrc-single-genotype-basic
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: axiomPoly_resource_vcf
    type: File
    doc: Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
    'sbg:x': 881.713623046875
    'sbg:y': 588.5
  - id: dbsnp_vcf
    type: File
    doc: Homo_sapiens_assembly38.dbsnp138.vcf
    'sbg:x': 362.0625
    'sbg:y': 588.5
  - id: hapmap_resource_vcf
    type: File
    doc: Hapmap genotype SNP input vcf
    'sbg:x': 0
    'sbg:y': 1070
  - id: input_vcfs
    type: 'File[]'
    doc: Input array of individual sample gVCF files
    'sbg:x': 0
    'sbg:y': 963
  - id: mills_resource_vcf
    type: File
    doc: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    'sbg:x': 0
    'sbg:y': 856
  - id: omni_resource_vcf
    type: File
    doc: 1000G_omni2.5.hg38.vcf.gz
    'sbg:x': 0
    'sbg:y': 749
  - id: one_thousand_genomes_resource_vcf
    type: File
    doc: '1000G_phase1.snps.high_confidence.hg38.vcf.gz, high confidence snps'
    'sbg:x': 0
    'sbg:y': 642
  - id: output_basename
    type: string
    'sbg:x': 1957.722412109375
    'sbg:y': 421
  - id: ped
    type: File?
    'sbg:x': 0
    'sbg:y': 535
  - id: reference_dict
    type: File
    doc: Homo_sapiens_assembly38.dict
    'sbg:x': 0
    'sbg:y': 428
  - id: reference_fasta
    type: File
    doc: Homo_sapiens_assembly38.fasta
    'sbg:x': 0
    'sbg:y': 321
  - id: snp_sites
    type: File
    doc: 1000G_phase3_v4_20130502.sites.hg38.vcf
    'sbg:x': 0
    'sbg:y': 214
  - id: unpadded_intervals_file
    type: File
    doc: hg38.even.handcurated.20k.intervals
    'sbg:x': 0
    'sbg:y': 107
  - id: wgs_evaluation_interval_list
    type: File
    doc: wgs_evaluation_regions.hg38.interval_list
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: collectvariantcallingmetrics
    outputSource:
      - picard_collectvariantcallingmetrics/output
    type: 'File[]'
    doc: Variant calling summary and detailed metrics files
    'sbg:x': 2496.104736328125
    'sbg:y': 588.5
  - id: gatk_filtered_vcf
    outputSource:
      - gatk_variantfiltration/output
    type: File
    'sbg:x': 2759.557861328125
    'sbg:y': 535
steps:
  - id: dynamicallycombineintervals
    in:
      - id: input_vcfs
        source:
          - input_vcfs
      - id: interval
        source: unpadded_intervals_file
    out:
      - id: out_intervals
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: script_dynamicallycombineintervals
      baseCommand:
        - python
        - '-c'
      inputs:
        - id: input_vcfs
          type: 'File[]'
        - id: interval
          type: File
      outputs:
        - id: out_intervals
          type: 'File[]'
          outputBinding:
            glob: out-*.intervals
            outputEval: >-
              ${ var i; var name = []; var dict = {}; for (i = 0; i <
              self.length; ++i) { name[i] = self[i].nameroot;
              dict[self[i].nameroot] = self[i]; }; name = name.sort(); for (i =
              0; i < name.length; ++i) { self[i] = dict[name[i]]; }; return
              self; }
      arguments:
        - position: 0
          valueFrom: |-
            def parse_interval(interval):
                colon_split = interval.split(":")
                chromosome = colon_split[0]
                dash_split = colon_split[1].split("-")
                start = int(dash_split[0])
                end = int(dash_split[1])
                return chromosome, start, end
            def add_interval(chr, start, end, i):
                fn = "out-{:0>5d}.intervals".format(i)
                lw = chr + ":" + str(start) + "-" + str(end) + "\n"
                with open(fn, "w") as fo:
                    fo.writelines(lw)
                return chr, start, end
            def main():
                interval = "$(inputs.interval.path)"
                num_of_original_intervals = sum(1 for line in open(interval))
                num_gvcfs = $(inputs.input_vcfs.length)
                merge_count = int(num_of_original_intervals/num_gvcfs/2.5)
                count = 0
                i = 1
                chain_count = merge_count
                l_chr, l_start, l_end = "", 0, 0
                with open(interval) as f:
                    for line in f.readlines():
                        # initialization
                        if count == 0:
                            w_chr, w_start, w_end = parse_interval(line)
                            count = 1
                            continue
                        # reached number to combine, so spit out and start over
                        if count == chain_count:
                            l_char, l_start, l_end = add_interval(w_chr, w_start, w_end, i)
                            w_chr, w_start, w_end = parse_interval(line)
                            count = 1
                            i += 1
                            continue
                        c_chr, c_start, c_end = parse_interval(line)
                        # if adjacent keep the chain going
                        if c_chr == w_chr and c_start == w_end + 1:
                            w_end = c_end
                            count += 1
                            continue
                        # not adjacent, end here and start a new chain
                        else:
                            l_char, l_start, l_end = add_interval(w_chr, w_start, w_end, i)
                            w_chr, w_start, w_end = parse_interval(line)
                            count = 1
                            i += 1
                    if l_char != w_chr or l_start != w_start or l_end != w_end:
                        add_interval(w_chr, w_start, w_end, i)
            if __name__ == "__main__":
                main()
      requirements:
        - class: DockerRequirement
          dockerPull: 'kfdrc/python:2.7.13'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: r4.2xlarge;ebs-gp2;500
    label: Combine intervals
    doc: Merge interval lists based on number of gVCF inputs
    'sbg:x': 362.0625
    'sbg:y': 474.5
  - id: gatk_applyrecalibration
    in:
      - id: indels_recalibration
        source: gatk_indelsvariantrecalibrator/recalibration
      - id: indels_tranches
        source: gatk_indelsvariantrecalibrator/tranches
      - id: input_vcf
        source: gatk_import_genotype_filtergvcf_merge/variant_filtered_vcf
      - id: snps_recalibration
        source: gatk_snpsvariantrecalibratorscattered/recalibration
      - id: snps_tranches
        source: gatk_gathertranches/output
    out:
      - id: recalibrated_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_applyrecalibration
      baseCommand: []
      inputs:
        - id: indels_recalibration
          type: File
          secondaryFiles:
            - .idx
        - id: indels_tranches
          type: File
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: snps_recalibration
          type: File
          secondaryFiles:
            - .idx
        - id: snps_tranches
          type: File
      outputs:
        - id: recalibrated_vcf
          type: File
          outputBinding:
            glob: scatter.filtered.vcf.gz
          secondaryFiles:
            - .tbi
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR -O
            tmp.indel.recalibrated.vcf -V $(inputs.input_vcf.path) --recal-file
            $(inputs.indels_recalibration.path) --tranches-file
            $(inputs.indels_tranches.path) -ts-filter-level 99.7
            --create-output-bam-index true -mode INDEL

            /gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR -O
            scatter.filtered.vcf.gz -V tmp.indel.recalibrated.vcf --recal-file
            $(inputs.snps_recalibration.path) --tranches-file
            $(inputs.snps_tranches.path) -ts-filter-level 99.7
            --create-output-bam-index true -mode SNP
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 7000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: r4.2xlarge;ebs-gp2;500
    label: GATK ApplyVQSR
    doc: Apply recalibration to snps and indels
    scatter:
      - input_vcf
      - snps_recalibration
    scatterMethod: dotproduct
    'sbg:x': 1543.257080078125
    'sbg:y': 588.5
  - id: gatk_calculategenotypeposteriors
    in:
      - id: output_basename
        source: output_basename
      - id: ped
        source: ped
      - id: reference_fasta
        source: reference_fasta
      - id: snp_sites
        source: snp_sites
      - id: vqsr_vcf
        source: gatk_gatherfinalvcf/output
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: kfdrc_gatk_calculategenotypeposteriors
      baseCommand:
        - /gatk
      inputs:
        - id: output_basename
          type: string
        - id: ped
          type: File?
        - id: reference_fasta
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
        - id: snp_sites
          type: File
          secondaryFiles:
            - .idx
        - id: vqsr_vcf
          type: File
          secondaryFiles:
            - .tbi
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --java-options "-Xms8000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
            CalculateGenotypePosteriors -R $(inputs.reference_fasta.path) -O
            $(inputs.output_basename).postCGP.vcf.gz -V $(inputs.vqsr_vcf.path)
            --supporting $(inputs.snp_sites.path) ${
              var arg = "";
              if (inputs.ped != null){
                arg += " --pedigree " + inputs.ped.path;
              }
              return arg;
            }
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 2
          coresMax: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 2188.3115234375
    'sbg:y': 588.5
  - id: gatk_gatherfinalvcf
    in:
      - id: input_vcfs
        source:
          - gatk_applyrecalibration/recalibrated_vcf
      - id: output_basename
        source: output_basename
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_gathervcfs
      baseCommand: []
      inputs:
        - id: input_vcfs
          type:
            type: array
            items: File
            inputBinding:
              prefix: '-I'
          inputBinding:
            position: 1
        - id: output_basename
          type: string
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: $(inputs.output_basename + '.vcf.gz')
          secondaryFiles:
            - .tbi
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xmx6g -Xms6g" GatherVcfsCloud
            --ignore-safety-checks --gather-type BLOCK --output
            $(inputs.output_basename + '.vcf.gz')
        - position: 2
          shellQuote: false
          valueFrom: '&& /gatk IndexFeatureFile -F $(inputs.output_basename + ''.vcf.gz'')'
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 7000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'migbro/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
    label: GATK GatherVcfsCloud
    doc: Combine resultant VQSR VCFs
    'sbg:x': 1957.722412109375
    'sbg:y': 642
  - id: gatk_gathertranches
    in:
      - id: tranches
        source:
          - gatk_snpsvariantrecalibratorscattered/tranches
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_gathertranches
      baseCommand: []
      inputs:
        - id: tranches
          type:
            type: array
            items: File
            inputBinding:
              prefix: '--input'
          inputBinding:
            position: 1
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: snps.gathered.tranches
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xmx6g -Xms6g" GatherTranches --output
            snps.gathered.tranches
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 7000
          coresMin: 2
        - class: DockerRequirement
          dockerPull: 'migbro/gatk:4.0.12.0'
      hints:
        - class: 'sbg:AWSInstanceType'
          value: r4.2xlarge;ebs-gp2;500
    label: GATK GatherTranches
    doc: Gather tranches from SNP variant recalibrate scatter
    'sbg:x': 1957.722412109375
    'sbg:y': 528
  - id: gatk_gathervcfs
    in:
      - id: input_vcfs
        source:
          - gatk_import_genotype_filtergvcf_merge/sites_only_vcf
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_gathervcfs
      baseCommand: []
      inputs:
        - id: input_vcfs
          type:
            type: array
            items: File
            inputBinding:
              prefix: '-I'
          inputBinding:
            position: 1
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: sites_only.vcf.gz
          secondaryFiles:
            - .tbi
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xmx6g -Xms6g" GatherVcfsCloud
            --ignore-safety-checks --gather-type BLOCK --output
            sites_only.vcf.gz
        - position: 2
          shellQuote: false
          valueFrom: '&& /gatk IndexFeatureFile -F sites_only.vcf.gz'
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 10000
          coresMin: 5
        - class: DockerRequirement
          dockerPull: 'migbro/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: r4.2xlarge;ebs-gp2;500
    label: Gather VCFs
    doc: Merge VCFs scattered from previous step
    'sbg:x': 881.713623046875
    'sbg:y': 481.5
  - id: gatk_import_genotype_filtergvcf_merge
    in:
      - id: dbsnp_vcf
        source: dbsnp_vcf
      - id: input_vcfs
        source:
          - input_vcfs
      - id: interval
        source: dynamicallycombineintervals/out_intervals
      - id: reference_fasta
        source: reference_fasta
    out:
      - id: sites_only_vcf
      - id: variant_filtered_vcf
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_import_genotype_filtergvcf_merge
      baseCommand: []
      inputs:
        - id: dbsnp_vcf
          type: File
          secondaryFiles:
            - .idx
        - id: input_vcfs
          type:
            type: array
            items: File
            inputBinding:
              prefix: '-V'
          inputBinding:
            position: 1
          secondaryFiles:
            - .tbi
        - id: interval
          type: File
        - id: reference_fasta
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
      outputs:
        - id: sites_only_vcf
          type: File
          outputBinding:
            glob: sites_only.variant_filtered.vcf.gz
          secondaryFiles:
            - .tbi
        - id: variant_filtered_vcf
          type: File
          outputBinding:
            glob: variant_filtered.vcf.gz
          secondaryFiles:
            - .tbi
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xms4g" GenomicsDBImport
            --genomicsdb-workspace-path genomicsdb --batch-size 50 -L
            $(inputs.interval.path) --reader-threads 16 -ip 5
        - position: 2
          shellQuote: false
          valueFrom: '&& tar -cf genomicsdb.tar genomicsdb'
        - position: 3
          shellQuote: false
          valueFrom: >-
            && /gatk --java-options "-Xmx8g -Xms4g" GenotypeGVCFs -R
            $(inputs.reference_fasta.path) -O output.vcf.gz -D
            $(inputs.dbsnp_vcf.path) -G StandardAnnotation
            --only-output-calls-starting-in-intervals -new-qual -V
            gendb://genomicsdb -L $(inputs.interval.path)
        - position: 4
          shellQuote: false
          valueFrom: >-
            && /gatk --java-options "-Xmx3g -Xms3g" VariantFiltration
            --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet -O
            variant_filtered.vcf.gz -V output.vcf.gz
        - position: 5
          shellQuote: false
          valueFrom: >-
            && /gatk MakeSitesOnlyVcf -I variant_filtered.vcf.gz -O
            sites_only.variant_filtered.vcf.gz
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
    label: 'Genotype, filter, & merge'
    doc: >-
      Use GATK GenomicsDBImport, VariantFiltration GenotypeGVCFs, and picard
      MakeSitesOnlyVcf to genotype, filter and merge gVCF based on known sites
    scatter:
      - interval
    'sbg:x': 585.4639892578125
    'sbg:y': 514
  - id: gatk_indelsvariantrecalibrator
    in:
      - id: axiomPoly_resource_vcf
        source: axiomPoly_resource_vcf
      - id: dbsnp_resource_vcf
        source: dbsnp_vcf
      - id: mills_resource_vcf
        source: mills_resource_vcf
      - id: sites_only_variant_filtered_vcf
        source: gatk_gathervcfs/output
    out:
      - id: recalibration
      - id: tranches
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_indelsvariantrecalibrator
      baseCommand: []
      inputs:
        - id: axiomPoly_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: dbsnp_resource_vcf
          type: File
          secondaryFiles:
            - .idx
        - id: mills_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: sites_only_variant_filtered_vcf
          type: File
          secondaryFiles:
            - .tbi
      outputs:
        - id: recalibration
          type: File
          outputBinding:
            glob: indels.recal
          secondaryFiles:
            - .idx
        - id: tranches
          type: File
          outputBinding:
            glob: indels.tranches
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator -V
            $(inputs.sites_only_variant_filtered_vcf.path) -O indels.recal
            --tranches-file indels.tranches --trust-all-polymorphic --mode INDEL
            --max-gaussians 4 -resource
            mills,known=false,training=true,truth=true,prior=12:$(inputs.mills_resource_vcf.path)
            -resource
            axiomPoly,known=false,training=true,truth=false,prior=10:$(inputs.axiomPoly_resource_vcf.path)
            -resource
            dbsnp,known=true,training=false,truth=false,prior=2:$(inputs.dbsnp_resource_vcf.path)
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche
            99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0
            -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche
            90.0 -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 7000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: r4.2xlarge;ebs-gp2;500
    label: GATK VariantRecalibrator Indels
    doc: >-
      Create recalibration model for indels using GATK VariantRecalibrator,
      tranch values, and known site VCFs
    'sbg:x': 1127.369873046875
    'sbg:y': 588.5
  - id: gatk_snpsvariantrecalibratorcreatemodel
    in:
      - id: dbsnp_resource_vcf
        source: dbsnp_vcf
      - id: hapmap_resource_vcf
        source: hapmap_resource_vcf
      - id: omni_resource_vcf
        source: omni_resource_vcf
      - id: one_thousand_genomes_resource_vcf
        source: one_thousand_genomes_resource_vcf
      - id: sites_only_variant_filtered_vcf
        source: gatk_gathervcfs/output
    out:
      - id: model_report
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_snpsvariantrecalibratorcreatemodel
      baseCommand:
        - /bin/bash
        - '-c'
      inputs:
        - id: dbsnp_resource_vcf
          type: File
          secondaryFiles:
            - .idx
        - id: hapmap_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: omni_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: one_thousand_genomes_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: sites_only_variant_filtered_vcf
          type: File
          secondaryFiles:
            - .tbi
      outputs:
        - id: model_report
          type: File
          outputBinding:
            glob: snps.model.report
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            set -eo pipefail

            /gatk --java-options "-Xmx60g -Xms15g" VariantRecalibrator -V
            $(inputs.sites_only_variant_filtered_vcf.path) -O snps.recal
            --tranches-file snps.tranches --trust-all-polymorphic --mode SNP
            --output-model snps.model.report --max-gaussians 6 -resource
            hapmap,known=false,training=true,truth=true,prior=15:$(inputs.hapmap_resource_vcf.path)
            -resource
            omni,known=false,training=true,truth=true,prior=12:$(inputs.omni_resource_vcf.path)
            -resource
            1000G,known=false,training=true,truth=false,prior=10:$(inputs.one_thousand_genomes_resource_vcf.path)
            -resource
            dbsnp,known=true,training=false,truth=false,prior=7:$(inputs.dbsnp_resource_vcf.path)
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche
            99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0
            -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an
            ReadPosRankSum -an FS -an MQ -an SOR -an DP || (echo 'Failed with
            max gaussians 6, trying 4' && /gatk --java-options "-Xmx60g -Xms15g"
            VariantRecalibrator -V
            $(inputs.sites_only_variant_filtered_vcf.path) -O snps.recal
            --tranches-file snps.tranches --trust-all-polymorphic --mode SNP
            --output-model snps.model.report --max-gaussians 4 -resource
            hapmap,known=false,training=true,truth=true,prior=15:$(inputs.hapmap_resource_vcf.path)
            -resource
            omni,known=false,training=true,truth=true,prior=12:$(inputs.omni_resource_vcf.path)
            -resource
            1000G,known=false,training=true,truth=false,prior=10:$(inputs.one_thousand_genomes_resource_vcf.path)
            -resource
            dbsnp,known=true,training=false,truth=false,prior=7:$(inputs.dbsnp_resource_vcf.path)
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche
            99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0
            -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an
            ReadPosRankSum -an FS -an MQ -an SOR -an DP)
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 7000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
    label: GATK VariantRecalibrator SNPs
    doc: >-
      Create recalibration model for snps using GATK VariantRecalibrator, tranch
      values, and known site VCFs
    'sbg:x': 1127.369873046875
    'sbg:y': 432.5
  - id: gatk_snpsvariantrecalibratorscattered
    in:
      - id: dbsnp_resource_vcf
        source: dbsnp_vcf
      - id: hapmap_resource_vcf
        source: hapmap_resource_vcf
      - id: model_report
        source: gatk_snpsvariantrecalibratorcreatemodel/model_report
      - id: omni_resource_vcf
        source: omni_resource_vcf
      - id: one_thousand_genomes_resource_vcf
        source: one_thousand_genomes_resource_vcf
      - id: sites_only_variant_filtered_vcf
        source: gatk_import_genotype_filtergvcf_merge/sites_only_vcf
    out:
      - id: recalibration
      - id: tranches
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: gatk_snpsvariantrecalibratorscattered
      baseCommand: []
      inputs:
        - id: dbsnp_resource_vcf
          type: File
          secondaryFiles:
            - .idx
        - id: hapmap_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: model_report
          type: File
        - id: omni_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: one_thousand_genomes_resource_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: sites_only_variant_filtered_vcf
          type: File
          secondaryFiles:
            - .tbi
      outputs:
        - id: recalibration
          type: File
          outputBinding:
            glob: scatter.snps.recal
          secondaryFiles:
            - .idx
        - id: tranches
          type: File
          outputBinding:
            glob: scatter.snps.tranches
      arguments:
        - position: 0
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator -V
            $(inputs.sites_only_variant_filtered_vcf.path) -O scatter.snps.recal
            --tranches-file scatter.snps.tranches --output-tranches-for-scatter
            --trust-all-polymorphic --mode SNP --input-model
            $(inputs.model_report.path) --max-gaussians 6 -resource
            hapmap,known=false,training=true,truth=true,prior=15:$(inputs.hapmap_resource_vcf.path)
            -resource
            omni,known=false,training=true,truth=true,prior=12:$(inputs.omni_resource_vcf.path)
            -resource
            1000G,known=false,training=true,truth=false,prior=10:$(inputs.one_thousand_genomes_resource_vcf.path)
            -resource
            dbsnp,known=true,training=false,truth=false,prior=7:$(inputs.dbsnp_resource_vcf.path)
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche
            99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0
            -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an
            ReadPosRankSum -an FS -an MQ -an SOR -an DP
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 7000
          coresMin: 1
        - class: DockerRequirement
          dockerPull: 'migbro/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: r4.2xlarge;ebs-gp2;500
    label: GATK VariantRecalibrator Scatter
    doc: >-
      Create recalibration model for known sites from input data using GATK
      VariantRecalibrator, tranch values, and known site VCFs
    scatter:
      - sites_only_variant_filtered_vcf
    'sbg:x': 1543.257080078125
    'sbg:y': 418.5
  - id: gatk_variantfiltration
    in:
      - id: cgp_vcf
        source: gatk_calculategenotypeposteriors/output
      - id: output_basename
        source: output_basename
      - id: reference_fasta
        source: reference_fasta
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: kfdrc_gatk_variantfiltration
      baseCommand:
        - /gatk
      inputs:
        - id: cgp_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: output_basename
          type: string
        - id: reference_fasta
          type: File
          secondaryFiles:
            - ^.dict
            - .fai
      outputs:
        - id: output
          type: File
          outputBinding:
            glob: '*.vcf.gz'
          secondaryFiles:
            - .tbi
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            --java-options "-Xms8000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
            VariantFiltration -R $(inputs.reference_fasta.path) -O
            $(inputs.output_basename).postCGP.Gfiltered.vcf.gz -V
            $(inputs.cgp_vcf.path) -G-filter "GQ < 20.0" -G-filter-name lowGQ
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 8000
          coresMin: 2
          coresMax: 4
        - class: DockerRequirement
          dockerPull: 'kfdrc/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
    'sbg:x': 2496.104736328125
    'sbg:y': 467.5
  - id: picard_collectvariantcallingmetrics
    in:
      - id: dbsnp_vcf
        source: dbsnp_vcf
      - id: input_vcf
        source: gatk_gatherfinalvcf/output
      - id: output_basename
        source: output_basename
      - id: reference_dict
        source: reference_dict
      - id: wgs_evaluation_interval_list
        source: wgs_evaluation_interval_list
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      id: picard_collectgvcfcallingmetrics
      baseCommand: []
      inputs:
        - id: dbsnp_vcf
          type: File
          secondaryFiles:
            - .idx
        - id: input_vcf
          type: File
          secondaryFiles:
            - .tbi
        - id: output_basename
          type: string
        - id: reference_dict
          type: File
        - id: wgs_evaluation_interval_list
          type: File
      outputs:
        - id: output
          type: 'File[]'
          outputBinding:
            glob: '*_metrics'
      arguments:
        - position: 1
          shellQuote: false
          valueFrom: >-
            /gatk --java-options "-Xmx6g -Xms6g" CollectVariantCallingMetrics -I
            $(inputs.input_vcf.path) -O $(inputs.output_basename) --DBSNP
            $(inputs.dbsnp_vcf.path) -SD $(inputs.reference_dict.path) -TI
            $(inputs.wgs_evaluation_interval_list.path) --THREAD_COUNT 8
      requirements:
        - class: ShellCommandRequirement
        - class: ResourceRequirement
          ramMin: 7000
          coresMin: 8
        - class: DockerRequirement
          dockerPull: 'migbro/gatk:4.0.12.0'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'sbg:AWSInstanceType'
          value: r4.2xlarge;ebs-gp2;500
    label: CollectVariantCallingMetrics
    doc: picard calculate variant calling metrics
    'sbg:x': 2188.3115234375
    'sbg:y': 425.5
hints:
  - class: 'sbg:AWSInstanceType'
    value: r4.4xlarge;ebs-gp2;500
  - class: 'sbg:maxNumberOfParallelInstances'
    value: 2
requirements:
  - class: ScatterFeatureRequirement
'sbg:appVersion':
  - v1.0
'sbg:id': kfdrc-harmonization/sd-bhjxbdqk-01/kfdrc-single-genotype-basic/0
'sbg:revision': 0
'sbg:revisionNotes': null
'sbg:modifiedOn': 1571880156
'sbg:modifiedBy': ennisb
'sbg:createdOn': 1571880156
'sbg:createdBy': ennisb
'sbg:project': kfdrc-harmonization/sd-bhjxbdqk-01
'sbg:projectName': GH-CBTTC-Germline
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - ennisb
'sbg:latestRevision': 0
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': ennisb
    'sbg:modifiedOn': 1571880156
    'sbg:revisionNotes': null
'sbg:image_url': >-
  https://cavatica.sbgenomics.com/ns/brood/images/kfdrc-harmonization/sd-bhjxbdqk-01/kfdrc-single-genotype-basic/0.png
'sbg:publisher': sbg
'sbg:content_hash': aa57c8b2baef5cf9793633b7ea2e43c9f335d3df1f2069b1ce2e750e7ed169468
