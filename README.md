# VAPER: Viral Assembly from Probe-based EnRichment
VAPER creates consensus-based assemblies from probe enrichment (a.k.a hybrid capture/enrichment) and shotgun metagenomic sequence data. One strength is it can handle samples containing multiple viral species and/or subtypes. When multiple viruses are present, VAPER will generate a consensus assembly for each, so long as an appropriate reference genome is supplied and the estimated genome fraction exceeds the minimum threshold (default: 10%). To ensure all relevant species are captured, VAPER supplies a summary of all viral sequences in the sample using [Sourmash](https://github.com/sourmash-bio/sourmash).

See the [wiki](https://github.com/DOH-JDJ0303/VAPER/wiki) for more information.

## Quick Start
### Step 1: Prepare your references
> Note: Nextflow requires absolute file paths.

`ref-list.csv`
```csv
taxa,assembly
Influenza-A_HA,flu-a-HA.fa
Influenza-A_NA,flu-a-NA.fa
Measles,measles.fa
```

### Step 2: Prepare your samplesheet
> Note: Nextflow requires absolute file paths.
`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
sample01,sample01_R1_001.fastq.gz,sample01_R2_001.fastq.gz
sample02,sample02_R1_001.fastq.gz,sample02_R2_001.fastq.gz
```

### Step 3: Run VAPER
```bash
nextflow run DOH-JDJ0303/VAPER \
    -r main \
    -profile <docker/singularity/.../institute> \
    --input samplesheet.csv \
    --refs ref-list.csv \
    --outdir results \
    --max_cpus 8 \
    --max_memory 16.GB
```
