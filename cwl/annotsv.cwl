cwlVersion: v1.0
class: CommandLineTool
label: annotsv
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: gaonkark/annotsv:latest
- class: InlineJavascriptRequirement
  expressionLib:
  - |2-

    var setMetadata = function(file, metadata) {
        if (!('metadata' in file))
            file['metadata'] = metadata;
        else {
            for (var key in metadata) {
                file['metadata'][key] = metadata[key];
            }
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

inputs:
- id: inputSV
  type: File
  inputBinding:
    prefix: -SVinputFile
    position: 0
    shellQuote: false

outputs:
- id: output
  type: File?
  outputBinding:
    glob: '*.tsv'
    outputEval: $(inheritMetadata(self, inputs.inputSV))

baseCommand:
- export
arguments:
- prefix: ''
  position: 0
  valueFrom: "${\n    \n    return \"ANNOTSV=/opt/AnnotSV_2.1 && $ANNOTSV/bin/AnnotSV\"\
    \n}"
  shellQuote: false
- prefix: ''
  position: 0
  valueFrom: -bedtools /usr/local/bin/bedtools
  shellQuote: false
- prefix: ''
  position: 99
  valueFrom: -outputDir ./
  shellQuote: false
- prefix: ''
  position: 0
  valueFrom: -genomeBuild GRCh38
  shellQuote: false
- prefix: ''
  position: 0
  valueFrom: -promoterSize 2000
  shellQuote: false
- prefix: ''
  position: 0
  valueFrom: -SVminSize 200
  shellQuote: false
- position: 0
  shellQuote: false

hints:
- class: sbg:AWSInstanceType
  value: c5.9xlarge;ebs-gp2;850
- class: sbg:maxNumberOfParallelInstances
  value: '4'
id: https://cavatica-api.sbgenomics.com/v2/apps/gaonkark/cbttc-dev/annotsv/7/raw/
sbg:appVersion:
- v1.0
sbg:content_hash: a71b866f16af33a887133cd0663552d13c43ee24cee36431570a475181afb0c8d
sbg:contributors:
- gaonkark
- danmiller
- wongj4
sbg:createdBy: gaonkark
sbg:createdOn: 1559761314
sbg:id: gaonkark/cbttc-dev/annotsv/7
sbg:image_url:
sbg:latestRevision: 8
sbg:modifiedBy: gaonkark
sbg:modifiedOn: 1571247979
sbg:project: gaonkark/cbttc-dev
sbg:projectName: PBTA-Dev
sbg:publisher: sbg
sbg:revision: 7
sbg:revisionNotes: bedtools path update
sbg:revisionsInfo:
- sbg:modifiedBy: gaonkark
  sbg:modifiedOn: 1559761314
  sbg:revision: 0
  sbg:revisionNotes: Copy of gaonkark/cnmc-test-space/annotsv/8
- sbg:modifiedBy: gaonkark
  sbg:modifiedOn: 1559831737
  sbg:revision: 1
  sbg:revisionNotes: Copy of gaonkark/cnmc-test-space/annotsv/9
- sbg:modifiedBy: gaonkark
  sbg:modifiedOn: 1565631794
  sbg:revision: 2
  sbg:revisionNotes: added instance type
- sbg:modifiedBy: wongj4
  sbg:modifiedOn: 1566580770
  sbg:revision: 3
  sbg:revisionNotes: removed zcat
- sbg:modifiedBy: wongj4
  sbg:modifiedOn: 1566581688
  sbg:revision: 4
  sbg:revisionNotes: removed zcat
- sbg:modifiedBy: gaonkark
  sbg:modifiedOn: 1566835339
  sbg:revision: 5
  sbg:revisionNotes: removed zcat ;generalized
- sbg:modifiedBy: gaonkark
  sbg:modifiedOn: 1566835417
  sbg:revision: 6
  sbg:revisionNotes: ''
- sbg:modifiedBy: gaonkark
  sbg:modifiedOn: 1571247979
  sbg:revision: 7
  sbg:revisionNotes: bedtools path update
- sbg:modifiedBy: danmiller
  sbg:modifiedOn: 1643219114
  sbg:revision: 8
  sbg:revisionNotes: stdout to stderr
sbg:sbgMaintained: false
sbg:validationErrors: []
sbg:workflowLanguage: CWL
