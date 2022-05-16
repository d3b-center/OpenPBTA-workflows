{
    "class": "Workflow",
    "cwlVersion": "v1.0",
    "id": "gaonkark/cbttc-dev/annotsv-wf/1",
    "label": "annotSV_wf",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "inputs": [
        {
            "id": "inputSV",
            "type": "File[]",
            "sbg:x": -563.3968505859375,
            "sbg:y": -95
        }
    ],
    "outputs": [
        {
            "id": "output",
            "outputSource": [
                "annotsv/output"
            ],
            "type": "File?",
            "sbg:x": 149.6031494140625,
            "sbg:y": -90
        }
    ],
    "steps": [
        {
            "id": "annotsv",
            "in": [
                {
                    "id": "inputSV",
                    "source": "inputSV"
                }
            ],
            "out": [
                {
                    "id": "output"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.0",
                "$namespaces": {
                    "sbg": "https://sevenbridges.com"
                },
                "id": "gaonkark/cbttc-dev/annotsv/2",
                "baseCommand": [
                    "zcat"
                ],
                "inputs": [
                    {
                        "id": "inputSV",
                        "type": "File",
                        "inputBinding": {
                            "position": 0,
                            "prefix": "-SVinputFile",
                            "shellQuote": false
                        }
                    },
                    {
                        "id": "inputvcf",
                        "type": "File?",
                        "inputBinding": {
                            "position": 0,
                            "prefix": "-vcfFiles"
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "output",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "*.tsv",
                            "outputEval": "$(inheritMetadata(self, inputs.inputSV))"
                        }
                    }
                ],
                "label": "annotsv",
                "arguments": [
                    {
                        "position": 0,
                        "prefix": "",
                        "shellQuote": false,
                        "valueFrom": "${\n    \n    return inputs.inputSV.path +\" && export ANNOTSV=/opt/AnnotSV_2.1 && $ANNOTSV/bin/AnnotSV\"\n}"
                    },
                    {
                        "position": 0,
                        "prefix": "",
                        "shellQuote": false,
                        "valueFrom": "-bedtools /tmp/bedtools2/bin/bedtools"
                    },
                    {
                        "position": 99,
                        "prefix": "",
                        "shellQuote": false,
                        "valueFrom": "-outputDir ./"
                    },
                    {
                        "position": 0,
                        "prefix": "",
                        "shellQuote": false,
                        "valueFrom": "-genomeBuild GRCh38"
                    },
                    {
                        "position": 0,
                        "prefix": "",
                        "shellQuote": false,
                        "valueFrom": "-promoterSize 2000"
                    },
                    {
                        "position": 0,
                        "prefix": "",
                        "shellQuote": false,
                        "valueFrom": "-SVminSize 200"
                    },
                    {
                        "position": 0
                    }
                ],
                "requirements": [
                    {
                        "class": "ShellCommandRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "gaonkark/annotsv:latest"
                    },
                    {
                        "class": "InlineJavascriptRequirement",
                        "expressionLib": [
                            "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file))\n        file['metadata'] = metadata;\n    else {\n        for (var key in metadata) {\n            file['metadata'][key] = metadata[key];\n        }\n    }\n    return file\n};\n\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n        }\n    }\n    return o1;\n};"
                        ]
                    }
                ],
                "hints": [
                    {
                        "class": "sbg:AWSInstanceType",
                        "value": "c5.9xlarge;ebs-gp2;850"
                    },
                    {
                        "class": "sbg:maxNumberOfParallelInstances",
                        "value": "4"
                    }
                ],
                "sbg:projectName": "PBTA-Dev",
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "gaonkark",
                        "sbg:modifiedOn": 1559761314,
                        "sbg:revisionNotes": "Copy of gaonkark/cnmc-test-space/annotsv/8"
                    },
                    {
                        "sbg:revision": 1,
                        "sbg:modifiedBy": "gaonkark",
                        "sbg:modifiedOn": 1559831737,
                        "sbg:revisionNotes": "Copy of gaonkark/cnmc-test-space/annotsv/9"
                    },
                    {
                        "sbg:revision": 2,
                        "sbg:modifiedBy": "gaonkark",
                        "sbg:modifiedOn": 1565631794,
                        "sbg:revisionNotes": "added instance type"
                    }
                ],
                "sbg:image_url": null,
                "sbg:appVersion": [
                    "v1.0"
                ],
                "sbg:id": "gaonkark/cbttc-dev/annotsv/2",
                "sbg:revision": 2,
                "sbg:revisionNotes": "added instance type",
                "sbg:modifiedOn": 1565631794,
                "sbg:modifiedBy": "gaonkark",
                "sbg:createdOn": 1559761314,
                "sbg:createdBy": "gaonkark",
                "sbg:project": "gaonkark/cbttc-dev",
                "sbg:sbgMaintained": false,
                "sbg:validationErrors": [],
                "sbg:contributors": [
                    "gaonkark"
                ],
                "sbg:latestRevision": 2,
                "sbg:publisher": "sbg",
                "sbg:content_hash": "a4177cc48496917db314490675d46e575ac0423b1777de2e9fb44c0bcc4fb7cca"
            },
            "label": "annotsv",
            "scatter": [
                "inputSV"
            ],
            "scatterMethod": "dotproduct",
            "sbg:x": -213,
            "sbg:y": -113
        }
    ],
    "requirements": [
        {
            "class": "ScatterFeatureRequirement"
        }
    ],
    "sbg:image_url": "https://cavatica.sbgenomics.com/ns/brood/images/gaonkark/cbttc-dev/annotsv-wf/1.png",
    "sbg:projectName": "PBTA-Dev",
    "sbg:revisionsInfo": [
        {
            "sbg:revision": 0,
            "sbg:modifiedBy": "gaonkark",
            "sbg:modifiedOn": 1565621617,
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 1,
            "sbg:modifiedBy": "gaonkark",
            "sbg:modifiedOn": 1565631832,
            "sbg:revisionNotes": "updated instance to c5.9"
        }
    ],
    "sbg:appVersion": [
        "v1.0"
    ],
    "sbg:id": "gaonkark/cbttc-dev/annotsv-wf/1",
    "sbg:revision": 1,
    "sbg:revisionNotes": "updated instance to c5.9",
    "sbg:modifiedOn": 1565631832,
    "sbg:modifiedBy": "gaonkark",
    "sbg:createdOn": 1565621617,
    "sbg:createdBy": "gaonkark",
    "sbg:project": "gaonkark/cbttc-dev",
    "sbg:sbgMaintained": false,
    "sbg:validationErrors": [],
    "sbg:contributors": [
        "gaonkark"
    ],
    "sbg:latestRevision": 1,
    "sbg:publisher": "sbg",
    "sbg:content_hash": "a97a2cb52693e203030430cdf7d28c27c552257aff19c7ca94c217d3848a9f928"
}