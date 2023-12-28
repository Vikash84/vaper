
[![GitHub Actions CI Status](https://github.com/DOH-JDJ0303/VAPER/workflows/nf-core%20CI/badge.svg)](https://github.com/DOH-JDJ0303/VAPER/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/DOH-JDJ0303/VAPER/workflows/nf-core%20linting/badge.svg)](https://github.com/DOH-JDJ0303/VAPER/actions?query=workflow%3A%22nf-core+linting%22)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/waphlviral/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/DOH-JDJ0303/VAPER)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23waphlviral-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/waphlviral)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**VAPER (Viral Assembly from Probe-based EnRichment)** creates consensus-based assemblies from probe enrichment (a.k.a hybrid capture/enrichment) sequence data. One strength is that it can handle samples containing multiple viral species and/or variants. In the case that multiple viruses are present, VAPER will generate a consensus assembly for each, so long as an appropriate reference genome is supplied and the estimated genome fraction exceeds the user defined threshold (default: 80%). To ensure all relevant species are captured, VAPER also supplies a summary of all viral sequences in the sample using Kraken2.

## Usage

:::note
If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
with `-profile test` before running the workflow on actual data.
:::

### Step 1: Preparing your reference genomes
VAPER creates assemblies using a consensus (i.e., reference-based) approach. As such, it is necessary to provide VAPER with appropriate references for each species/variant you intend to assemble. References are provided as individual FASTA files within a tar compressed directory. Assemblies created from references containing multiple contigs will be concatenated into a single contig. See instructions below for how prepare the reference directory.

#### Gather all your reference genomes and place them into a single directory
📂refs\
 ┣ 📜sars-cov-2.fasta\
 ┣ 📜mumps.fasta\
 ┣ 📜measles.fasta\
 ┣ 📜flu-a-h1n1.fasta\
 ┣ 📜flu-a-h3n2.fasta\
 ┗ 📜flu-b.fasta

#### Compress the directory
 ```
 tar czvf refs.tar.gz refs/
 ```
#### Prepapre your reference metadata file (Optional)
Metadata can be provided for each reference assembly. This data will be incorporated into the final report and is intended to aid interpretation. The `REFERENCE` column is the only required field. Otherwise, you can provide whatever fields/information you want. See an example below.
`refs-meta.csv`:
```csv
REFERENCE,SPECIES,VARIANT
sars-cov-2.fasta,Severe acute respiratory syndrome coronavirus 2,NA
mumps.fasta,Mumps orthorubulavirus,NA
measles.fasta,Measles Morbillivirus,NA
flu-a-h1n1.fasta,Influenza A virus,H1N1
flu-a-h3n2.fasta,Influenza A virus,H3N2
flu-b.fasta,Influenza B virus,NA
```

### Step 2: Download the Kraken2 RefSeq viral database
VAPER gives you a summary of all viral species in your sample, as determined via Kraken2 and the RefSeq viral database. This step is completely independent of consensus assembly generation and is only meant to ensure that you are capturing all relevant species in your sample. You can download the most recent version of the RefSeq viral database [here](https://benlangmead.github.io/aws-indexes/k2). An example of how to do this from the command-line is shown below:
```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20231009.tar.gz
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
   --refs $PWD/refs.tar.gz \
   --refs_meta $PWD/refs-meta.csv \
   --k2db $PWD/k2_viral_20231009.tar.gz \
   --outdir <OUTDIR>
```
### Step 4: Fine tuning your assembly
Adjust one or more of the options below to fine-tune your assembly.
```
options:
--gen_frac        Minimum genome fraction for an assembly to be created (Default: 0.8)
--assembler       Assembler to use for Shovill (skesa, spades, velvet, or megahit) (Default: spades)
--min_contig_cov  Minimum contig coverage for Shovill (Default: 2)
--min_contig_len  Minimum contig length for Shovill (Default: 100)
--gsize           Approx. genome size for Shovill (Default: 1.0M)
```

:::warning
Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).
:::

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
