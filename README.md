# VAPER: Viral Assembly from Probe-based EnRichment
VAPER is a viral assembly pipeline. Key features include:
-  Builds assemblies from probe enrichment (a.k.a hybrid capture/enrichment), shotgun metagenomic, and tiled-amplicon sequence data
-  Automated reference selection (Detects what is in your sample)
-  Can an generate multiple assemblies per sample (Useful for co-infections)
-  Can utilize iVar or IRMA assemblers (IRMA modules built on the fly!)
-  Reads associated with each assembly are exported for downstream use
-  Provides a metagenomic summary of dominant viral taxa

VAPER comes with comprehensive reference sets for the following viral species (created using [EPITOME](https://github.com/DOH-JDJ0303/epitome)):
- Influenza A (segments 1-8)
- Influenza B (segments 1-8)
- Measles morbillivirus
- Mumps orthorubulavirus
- Lyssavirus rabies
- Norovirus
- Respiratory Syncytial Virus (RSV)
- West Nile virus
- Enterovirus D68
- Hepacivirus (HCV)
- Hepatovirus (HAV)
  
See the [wiki](https://github.com/DOH-JDJ0303/VAPER/wiki) for more information.

## Quick Start
### Step 1: Prepare your samplesheet
> Note: Nextflow requires absolute file paths.
`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
sample01,sample01_R1_001.fastq.gz,sample01_R2_001.fastq.gz
sample02,sample02_R1_001.fastq.gz,sample02_R2_001.fastq.gz
```

### Step 2: Run VAPER
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
