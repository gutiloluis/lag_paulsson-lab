# Pipeline Overview

This pipeline is designed to process and analyze sequencing data from isolates using a combination of graph-based alignment and downstream analysis with the goal of estimating the error rate of the nanopore analysis of the barcoded libraries.

## Prerrequisites

- Output files from nanopore analysis (with the same format as the files usually named `output.tsv`, or `*_final_df*.pkl`), stored in a folder called `nanopore`.
- Reads from isolate sequencing (e.g. Plasmidsaurus), concatenated into a single `fasta` file for each library. Note: The code assumes a specific format for file naming. If it is different, the file `scripts/get_grnas_from_gfa_and_reads.py` needs to be modified.
- Reference `gfa` file with the sequence graph.

## Workflow Summary

1. **Graph Alignment (`graph_align` rule)**:
    - Aligns sequencing reads of isolates to a graph-based reference genome using `GraphAligner`.
    - Outputs alignment files in `.gaf` format.

2. **Barcode Extraction (`get_barcodes` rule)**:
    - Extracts barcodes from the `.gaf` alignment files of the isolates and generates a string of 1s and 0s for the barcode.
    - Saves the extracted data in `.pkl` format for downstream processing.

3. **gRNA Extraction (`get_grnas` rule)**:
    - Combines the data with the barcodes and the original sequencing reads to identify gRNA sequences from the isolates.
    - Outputs data with gRNA sequences in `.pkl` format.

4. **Mapping to Nanopore Data (`map_to_nanopore` rule)**:
    - Maps the isolate barcodes into the nanopore data, and compares if the gRNAs between the isolates and the nanopore data match.
    - Produces a summary of the mapping results and an error report.

5. **Error Report Combination (`combine_error_reports` rule)**:
    - Aggregates individual error reports into a single combined report.

## Outputs

- **Graph Alignment Files**: `.gaf`
- **Barcode Data**: `.pkl`
- **gRNA Data**: `.pkl`
- **Mapping Summary**: `.pkl`
- **Error Reports**: `.txt`

## TODO - Potential improvements

- The code 
- Show as a stacked bar plot instead of just text.
- Extract library names from the filename and save the error reports as tables.
- Currently it analyzes all isolates vs all the nanopore files. This is not necesary.