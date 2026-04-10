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
    --reads $PWD/mpxv_cladeIIa_testcase/ \
    --out_dir output_cladeIIa_testcase \
    --reference cladeIIa_GIN \
    --min_coverage 20
```

This will download the pipeline and run it on the test data. Necessary Docker images will be downloaded automatically during the run, so you need to have an internet connection throughout your first analysis. Afterwards, the pipeline works offline. 
The option `--reads` lets you pass a directory with a fastq file (for one sample). You have to pass a full path so in this case we used $PWD to resolve to the current working directory. `--out_dir` defines where the output is written to. The test data is clade IIa, so we chose the reference with `--reference cladeIIa_GIN` which uses a full length mpox genome from guinea as the reference. Other options here are `cladeIIb` which uses the clade IIb RefSeq genome NC_063383.1 and `cladeIIb_GIN_G1` which letsyou use a clade IIb lineage G1 full length genome from the outbreak in Guinea 2025. 

On the setup described below, installation and running this example required approximately 10 minutes. Further runs (pipeline and docker containers already installed) took 3 to 4 minutes.

### EPI2ME desktop

You can also use EPI2ME desktop provided by Oxford Nanopore Technologies to run this pipeline in a graphical interface.
Download and follow the instructions for setup [here](https://nanoporetech.com/software/other/epi2me-desktop-application).

## System requirements

This pipeline was tested on a MacBook Pro with 2,4 GHz 8-core Intel Core i9 processor and 64 GB RAM with macOS Sequoia 15.7.4 operating system.
Analyses on this machine took less than 10 minutes.
It should, however, also run on a much weaker laptop.

### Output

The pipeline write the following output files:

| Directory          | File                          | Description                         |
| ------------------ | ----------------------------- | ----------------------------------- |
| `report/`          | `sample_report.html`          | Final summary report                |
| `draft/`           | `polished.fasta`              | Consensus sequence with low coverage regions NOT yet masked (INTERMEDIATE result for investigations) |
| `draft/`           | `draft_final.fasta`           | Assembly with low coverage regions masked - Final output  |
| `draft/`           | `basecounts_final.tsv`        | Base composition statistics         |
| `alignments/`      | `reads_vs_ref.bam`            | Reads mapped to reference           |
| `alignments/`      | `reads_vs_ref.bam.bai`        | BAM index                           |
| `alignments/`      | `reads_vs_polished.bam`       | Reads mapped to polished assembly   |
| `alignments/`      | `reads_vs_polished.bam.bai`   | BAM index                           |
| `alignments/`      | `contigs_vs_ref.bam`          | Contigs mapped to reference         |
| `alignments/`      | `contigs_vs_ref.bam.bai`      | BAM index                           |
| `alignments/`      | `contigs_vs_polished.bam`     | Contigs mapped to polished assembly |
| `alignments/`      | `contigs_vs_polished.bam.bai` | BAM index                           |
| `annotation/`      | `variants.vcf.gz`             | Called variants                     |
| `annotation/`      | `variants.vcf.gz.tbi`         | Index for VCF file                  |
| `annotation/`      | `alignment.fasta`             | Alignment of consensus genome to reference genome used to get variants      |
| `annotation/`      | `gene_variant_summary.txt`    | Per-gene variant summary            |
| `annotation/`      | `variant_summary.html`        | Interactive variant report from SNPEff          |
| `./`               | `nf-report-*.html`            | Nextflow execution reports          |
