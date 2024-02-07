
[![GitHub Actions CI Status](https://github.com/DOH-JDJ0303/VAPER/workflows/nf-core%20CI/badge.svg)](https://github.com/DOH-JDJ0303/VAPER/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/DOH-JDJ0303/VAPER/workflows/nf-core%20linting/badge.svg)](https://github.com/DOH-JDJ0303/VAPER/actions?query=workflow%3A%22nf-core+linting%22)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/waphlviral/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/DOH-JDJ0303/VAPER)

## Introduction

**VAPER (Viral Assembly from Probe-based EnRichment)** creates consensus-based assemblies from probe enrichment (a.k.a hybrid capture/enrichment) sequence data. One strength is that it can handle samples containing multiple viral species and/or variants. When multiple viruses are present, VAPER will generate a consensus assembly for each, so long as an appropriate reference genome is supplied and the estimated genome fraction exceeds the user defined threshold (default: 70%). To ensure all relevant species are captured, VAPER also supplies a summary of all viral sequences in the sample using [Sourmash](https://github.com/sourmash-bio/sourmash).

## Usage

:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

### Step 1: Preparing your reference genomes
VAPER creates assemblies using a consensus (i.e., reference-based) approach. As such, it is necessary to provide VAPER with appropriate references for each species/variant you intend to assemble. References are provided in a samplesheet. An example of how to create this samplesheet is shown below.

`ref-list.csv`
```csv
taxa,assembly
Influenza_A_virus_H1N1,GCF_001343785.1_ViralMultiSegProj274766_genomic.fna
Influenza_A_virus_H2N2,GCF_000866645.1_ViralMultiSegProj15620_genomic.fna
Influenza_A_virus_H3N2,GCF_000865085.1_ViralMultiSegProj15622_genomic.fna
Influenza_A_virus_H5N1,GCF_000864105.1_ViralMultiSegProj15617_genomic.fna
Influenza_A_virus_H7N9,GCF_000928555.1_ViralMultiSegProj274585_genomic.fna
Influenza_A_virus_H9N2,GCF_000851145.1_ViralMultiSegProj14892_genomic.fna
Influenza_B_virus,GCF_000820495.2_ViralMultiSegProj14656_genomic.fna
Lyssavirus_rabies,GCF_000859625.1_ViralProj15144_genomic.fna
Measles_Morbillivirus,GCF_000854845.1_ViralProj15025_genomic.fna
Mumps_orthorubulavirus,GCF_000856685.1_ViralProj15059_genomic.fna
Severe_acute_respiratory_syndrome_coronavirus_2,GCF_009858895.2_ASM985889v3_genomic.fna
West_Nile_virus,GCF_000875385.1_ViralProj30293_genomic.fna
```

### Step 2: Prepare your samplesheet
VAPER takes a standard Nextflow samplesheet as input (see example below).
`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
sample01,sample01_R1_001.fastq.gz,sample01_R2_001.fastq.gz
sample02,sample02_R1_001.fastq.gz,sample02_R2_001.fastq.gz
```

### Step 3: Run VAPER
Run VAPER using the command below, making adjustments where necessary.
```bash
nextflow run DOH-JDJ0303/VAPER \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --refs ref-list.csv \
   --outdir <OUTDIR>
```
### Step 4: Fine tuning your assembly
Adjust one or more of the options below to fine-tune your assembly.
```
options:
--mode            Reference selection mode ('fast' or 'accurate'; default: 'accurate')
--avg_depth       Minimum average depth of coverage for an assembly to be created (default: 100). Only used in 'fast' mode.
--gen_frac        Minimum genome fraction for an assembly to be created (default: 0.7). Used in 'fast' and 'accurate' mode.
--assembler       Assembler to use for Shovill (skesa, spades, velvet, or megahit) (Default: spades)
--min_contig_cov  Minimum contig coverage for Shovill (Default: 10)
--min_contig_len  Minimum contig length for Shovill (Default: 100)
--gsize           Approx. genome size for Shovill (Default: 1.0M)
--ivar_q          Minimum quality score threshold to count base for ivar (default: 20)
--ivar_t          Minimum frequency threshold(0 - 1) to call consensus for ivar (default: 0.5)
--ivar_n          (N/-) Character to print in regions with less than minimum coverage for ivar (default: N)
--ivar_m          Minimum depth to call consensus for ivar (default: 10)
```

## Pipeline output


## Credits

DOH-JDJ0303/VAPER was originally written by Jared Johnson.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#waphlviral` channel](https://nfcore.slack.com/channels/waphlviral) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  DOH-JDJ0303/VAPER for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
