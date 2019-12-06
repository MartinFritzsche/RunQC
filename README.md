# RunQC
An R Markdown script to automatically generate interactive post-run QC reports for Illumina MiSeq and NextSeq.

## Primary Goals
1. Generate plots similar to what is possible within [Illumina Sequence Analysis Viewer](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/sav/sequencing-analysis-viewer-v-2-4-software-guide-15066069-04.pdf)
2. Combine plots with general info about sequencing run (e.g. reagent info, included projects, multiplex-info)

## Secondary Goals
1. Interactive plots
2. Multi-run comparisons

## Inputs
- InterOp Folder
- runParameters.xml
- SampleSheet.csv

## Outputs
- html report
