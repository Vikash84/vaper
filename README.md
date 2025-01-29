# VAPER: Viral Assembly from Probe-based EnRichment
## Overview
VAPER is a pipeline for virus genome assembly.
### Key Features:
-  Builds assemblies from probe enrichment (a.k.a hybrid capture/enrichment), shotgun metagenomic, and tiled-amplicon sequence data
-  Automated reference selection (Detects what is in your sample)
-  Can generate multiple assemblies per sample (Useful for co-infections)
-  Predicts the taxonomy of each assembly (Also includes a *metagenomic* summary!)
-  Reads associated with each assembly are exported for downstream use
-  Can utilize iVar or IRMA assemblers (IRMA modules built on the fly!)

### Automated Reference Selection:
VAPER comes with comprehensive reference sets for the following viral taxa (created using [EPITOME](https://github.com/DOH-JDJ0303/epitome)):
|Taxon               |Segments  | No. References|No. Species | No. Input Sequences|
|:-------------------|:---------|--------------:|:-----------|-------------------:|
|Alphacoronavirus    |wg        |            201|67          |                1040|
|Alphainfluenzavirus |1 - 8     |           2394|1           |              484876|
|Betacoronavirus     |wg        |            140|53          |                3566|
|Betainfluenzavirus  |1 - 8     |             94|1           |               96657|
|Bocaparvovirus      |wg        |            104|48          |                 608|
|Enterovirus         |wg        |           1967|26          |                6740|
|Hantavirus          |l, m, s   |            244|32          |                 825|
|Hepacivirus         |wg        |            888|42          |                1060|
|Hepatovirus         |wg        |             60|14          |                 114|
|Lyssavirus          |wg        |            198|23          |                2126|
|Mastadenovirus      |wg        |            128|73          |                1507|
|Metapneumovirus     |wg        |             28|4           |                 420|
|Morbillivirus       |wg        |            100|14          |                1176|
|Norovirus           |wg        |            284|5           |                1148|
|Orthoflavivirus     |wg        |            428|83          |                7960|
|Orthopneumovirus    |wg        |             22|4           |               18344|
|Orthopoxvirus       |wg        |             15|13          |                7885|
|Orthorubulavirus    |wg        |             38|9           |                1053|
|Respirovirus        |wg        |             51|17          |                 531|

Click [here](https://github.com/DOH-JDJ0303/vaper/blob/main/assets/reference_sets/readme_EPITOME_2025-01-23.md) to learn more about these references.

### More Information:
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
> Note: This will use the default reference set. You can provide your own reference set using the `--refs` parameter.
```bash
nextflow run DOH-JDJ0303/VAPER \
    -r main \
    -profile <docker/singularity/.../institute> \
    --input samplesheet.csv \
    --outdir results \
    --max_cpus 8 \
    --max_memory 16.GB
```
## Acknowledgements
**VAPER would not be possible without the following people:**
- Holly Halstead (WA PHL) (check out her tool [varcraft](https://github.com/DOH-HNH0303/varcraft)!)
- Pauline Trinh (WA MEP)
- Allison Black (WA MEP)
- Stephanie Lunn (WA MEP)
- Kristen Waterman (WA PHL)
- Brandi Torrevillas (WSU)

VAPER was originally written by Jared Johnson for the Washington State Department of Health.

