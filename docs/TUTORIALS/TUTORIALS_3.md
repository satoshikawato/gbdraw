# Tutorials

These tutorials are step-by-step lessons designed to teach core concepts.

## Comparative Genomics with BLAST


**Goal**: Learn how to use the linear mode to visualize sequence similarity between two or more genomes.

### 1. Required Inputs
For comparative plots, you need:

Two or more genome files (e.g., genome1.gbk, genome2.gbk).

A BLAST result file showing comparisons between them (e.g., g1_vs_g2.blast).

### 2. BLAST Output Format
The BLAST result file must be in tabular format, specifically outfmt 7. This format includes commented header lines that 

gbdraw uses to identify the query and subject sequences.

### 3. The gbdraw Command
The command structure involves listing all GenBank files and then all BLAST result files using the -b flag.

This will generate a linear plot showing genome1 and genome2 with ribbons connecting regions of sequence similarity defined in the BLAST file.

