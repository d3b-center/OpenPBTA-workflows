class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://sevenbridges.com'
id: gaonkark/cbttc-dev/annotsv/1
baseCommand:
  - zcat
inputs:
  - id: inputSV
    type: File
    inputBinding:
      position: 0
      prefix: '-SVinputFile'
      shellQuote: false
  - id: inputvcf
    type: File?
    inputBinding:
      position: 0
      prefix: '-vcfFiles'
outputs:
  - id: output
    type: File?
    outputBinding:
      glob: '*.tsv'
      outputEval: '$(inheritMetadata(self, inputs.inputSV))'
label: annotsv
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: |-
      ${
          
          return inputs.inputSV.path +" && export ANNOTSV=/opt/AnnotSV_2.1 && $ANNOTSV/bin/AnnotSV"
      }
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: '-bedtools /tmp/bedtools2/bin/bedtools'
  - position: 99
    prefix: ''
    shellQuote: false
    valueFrom: '-outputDir ./'
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: '-genomeBuild GRCh38'
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: '-promoterSize 2000'
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: '-SVminSize 200'
  - position: 0
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'gaonkark/annotsv:latest'
  - class: InlineJavascriptRequirement
    expressionLib:
      - |-

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
'sbg:projectName': PBTA-Dev
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': gaonkark
    'sbg:modifiedOn': 1559761314
    'sbg:revisionNotes': Copy of gaonkark/cnmc-test-space/annotsv/8
  - 'sbg:revision': 1
    'sbg:modifiedBy': gaonkark
    'sbg:modifiedOn': 1559831737
    'sbg:revisionNotes': Copy of gaonkark/cnmc-test-space/annotsv/9
'sbg:image_url': null
'sbg:appVersion':
  - v1.0
'sbg:id': gaonkark/cbttc-dev/annotsv/1
'sbg:revision': 1
'sbg:revisionNotes': Copy of gaonkark/cnmc-test-space/annotsv/9
'sbg:modifiedOn': 1559831737
'sbg:modifiedBy': gaonkark
'sbg:createdOn': 1559761314
'sbg:createdBy': gaonkark
'sbg:project': gaonkark/cbttc-dev
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - gaonkark
'sbg:latestRevision': 1
'sbg:publisher': sbg
'sbg:content_hash': a2f24845abac780332194b9a20f9ed735a0dead1fcf0f963916d968139f1bf0f0
'sbg:copyOf': gaonkark/cnmc-test-space/annotsv/9
