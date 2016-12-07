class: Workflow
cwlVersion: v1.0

doc: "Epitope prediction workflow"

inputs: 

  DNA: File

  RNA: File

  algorithmDependencies: Directory

outputs:

  finalOutput:
    type: File
    outputSource: bindingAffinity/output

steps:

  snvDetection:
    run: snv.cwl
    doc: "SNC detection"
    in:
      input: DNA
      dependencies: algorithmDependencies
      parameters: { default: string }
    out: [output]

  indelDetection:
    run: indel.cwl
    doc: "Indel detection"
    in:
      input: DNA
      dependencies: algorithmDependencies
    out: [output]

  fusionsDetection:
    run: fusions.cwl
    doc: "Fusion detection"
    in:
      input: RNA
      dependencies: algorithmDependencies
    out: [output]

  RNAExp:
    run: RNAExp.cwl
    doc: "RNA expression checker"
    in:
      input: [snvDetection/output,indelDetection/output,fusionsDetection/output]
    out: [output]

  netMHC:
    run: netMHC.cwl
    doc: "netMHC binding"
    in:
      input: RNAExp/output
    out: [output]

  bindingAffinity:
    run: bindingAffinity.cwl
    doc: "Binding affinity comparison to wild type"
    in:
      input: netMHC/output
    out: [output]