cwlVersion: v1.0
class: CommandLineTool
label: annovar-yiran
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 16
  ramMin: 36000
- class: DockerRequirement
  dockerPull: kfdrc/annovar:latest
- class: InlineJavascriptRequirement

inputs:
- id: annovar_ref
  type: File
- id: input_vcf
  type: File

outputs:
- id: anno_vcf
  type: File
  outputBinding:
    glob: '*.hg38_multianno.vcf.gz'
- id: anno_tbi
  type: File
  outputBinding:
    glob: '*.hg38_multianno.vcf.gz.tbi'

baseCommand:
- tar
- -xzf
arguments:
- position: 1
  valueFrom: |-
    $(inputs.annovar_ref.path) && perl /home/TOOLS/tools/annovar/current/bin/table_annovar.pl $(inputs.input_vcf.path) ./humandb_20190319 -buildver hg38 -out $(inputs.input_vcf.basename) -remove -protocol refGene,intervar_20180118,abraom,nci60,kaviar_20150923,hrcr1,gnomad_exome,gnomad_genome,dbnsfp35a,avsnp150,clinvar_20190305,dbscsnv11,gme,regsnpintron,ljb26_all,cosmic88_coding,cosmic88_noncoding -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -argument '--hgvs --splicing_threshold 10',,,,,,,,,,,,,,,, -thread 12 && bgzip -i $(inputs.input_vcf.basename).hg38_multianno.vcf
  shellQuote: false
id: |-
  https://cavatica-api.sbgenomics.com/v2/apps/cavatica/qtx6-ewhc/annovar-20190319-yiran/0/raw/
sbg:appVersion:
- v1.0
sbg:content_hash: aa3c17ed5473b2a9198e9276b39a2538300b7ee99fa4708c7484c9097960e25fa
sbg:contributors:
- yiran
sbg:copyOf: cavatica/22q11-deletion-syndrome-project/annovar-20190319-yiran/2
sbg:createdBy: yiran
sbg:createdOn: 1553104743
sbg:id: cavatica/qtx6-ewhc/annovar-20190319-yiran/0
sbg:image_url:
sbg:latestRevision: 0
sbg:modifiedBy: yiran
sbg:modifiedOn: 1553104743
sbg:project: cavatica/qtx6-ewhc
sbg:projectName: PBTA Internal Project
sbg:publisher: sbg
sbg:revision: 0
sbg:revisionNotes: Copy of cavatica/22q11-deletion-syndrome-project/annovar-20190319-yiran/2
sbg:revisionsInfo:
- sbg:modifiedBy: yiran
  sbg:modifiedOn: 1553104743
  sbg:revision: 0
  sbg:revisionNotes: Copy of cavatica/22q11-deletion-syndrome-project/annovar-20190319-yiran/2
sbg:sbgMaintained: false
sbg:validationErrors: []
sbg:workflowLanguage: CWL
