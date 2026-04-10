# Metagenomics Orthopoxvirus Pipeline X (MOPX)

This pipelines assembles genomes for orthopoxviruses from nanopore metagenomic reads.
It is build and maintained by the  [Outbreak Preparedness and Response team](https://www.bnitm.de/forschung/forschungsgruppen/pathogen/abt-virologie/laborgruppe-duraffour-pahlmann/team) at the Bernhard Nocht Institute for Tropical Medicine (BNITM).
Currently, it only has mpox reference genomes, but arbitrary orthopoxvirus species should work.

## Workflow

The pipelines performs the following steps
- de-novo assembly (with Canu and Flye)
- scaffolding to a reference genome
- read mapping and polishing of the draft
- correction of the terminal, non-complementary hairpin loops
- annotation of variants with respect to the reference genome
- write a html report

## Setup

### Software dependencies

You need [nextflow](https://www.nextflow.io/) and Docker installed to run this pipeline.
It was tested on macOS Sequoia 15.7.4 with nextflow 25.10.2 and Docker version 25.0.2.

### Running the pipeline on a test case

You can download some test data [here](https://opr.bnitm.de/mopx_testcase/mpxv_cladeIIa_testcase.tar.gz) or with the command 

```wget https://opr.bnitm.de/mopx_testcase/mpxv_cladeIIa_testcase.tar.gz```

followed by the extraction of the downloaded archive with 

```tar -xzf mpxv_cladeIIa_testcase.tar.gz```

Then run the pipeline with

```
nextflow run opr-group-bnitm/mopx \
    --reads ./mpxv_cladeIIa_testcase/ \
    --out_dir output_cladeIIa_testcase \
    --reference cladeIIa_GIN \
    --min_coverage 20
```

### EPI2ME desktop

You can also use EPI2ME desktop provided by Oxford Nanopore Technologies to run this pipeline in a graphical interface.
Download and follow the instructions for setup [here](https://nanoporetech.com/software/other/epi2me-desktop-application).

## System requirements

This pipeline was tested on a MacBook Pro with 2,4 GHz 8-core Intel Core i9 processor and 64 GB RAM with macOS Sequoia 15.7.4 operating system.
Analyses on this machine took less than 10 minutes.
It should, however, also run on a much weaker laptop.

