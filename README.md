# Metagenomics Orthopoxvirus Pipeline X (MOPX)

This pipelines assembles genomes for orthopoxviruses from nanopore metagenomic reads.
It is build and maintained by the Outbreak Preparedness and Response (OPR) team at the Bernhard Nocht Institute for Tropical Medicine (BNITM).
Currently, it only has mpox reference genomes, but arbitrary orthopoxvirus species should work.

## Workflow

The pipelines performs the following steps
- de-novo assembly (with Canu and Flye)
- scaffolding to a reference genome
- read mapping and polishing of the draft
- correction of the terminal, non-complementary hairpin loops
- annotation of variants with respect to the reference genome
- write a html report

