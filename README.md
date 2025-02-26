# VAPER: Viral Assembly from Probe-based EnRichment
## Overview
VAPER is a viral (meta-)genome assembly pipeline.
### Key Features:
-  Builds assemblies from probe enrichment (a.k.a hybrid capture/enrichment), shotgun metagenomic, and tiled-amplicon sequence data
-  Automated reference selection (Detects what is in your sample)
-  Can generate multiple assemblies per sample (Useful for co-infections)
-  Predicts the taxonomy of each assembly (Also includes a viral *metagenomic* summary!)
-  Reads associated with each assembly are exported for downstream use
-  Can utilize iVar or IRMA assemblers (IRMA modules built on the fly!)

### Automated Reference Selection:
VAPER comes with comprehensive reference sets for the following viral taxa (created using [EPITOME](https://github.com/DOH-JDJ0303/epitome)):
|Taxon               |Segments  | No. References|No. Species | No. Input Sequences|
|:-------------------|:---------|--------------:|:-----------|-------------------:|
|Alphacoronavirus    |wg        |            204|66          |                1104|
|Alphainfluenzavirus |1 - 8     |           1939|1           |              385051|
|Betacoronavirus     |wg        |            135|53          |                3436|
|Betainfluenzavirus  |1 - 8     |             45|1           |              110948|
|Bocaparvovirus      |wg        |             94|48          |                 574|
|Enterovirus         |wg        |           1782|26          |                5870|
|Hantavirus          |[l; m; s] |            227|31          |                 762|
|Hepacivirus         |wg        |            790|42          |                 959|
|Hepatovirus         |wg        |             56|14          |                 110|
|Lyssavirus          |wg        |            197|23          |                2532|
|Mastadenovirus      |wg        |            115|69          |                1339|
|Metapneumovirus     |wg        |             24|3           |                 380|
|Morbillivirus       |wg        |             87|14          |                1075|
|Norovirus           |wg        |            245|5           |                 947|
|Orthoflavivirus     |wg        |            371|83          |                6657|
|Orthopneumovirus    |wg        |             18|4           |               17121|
|Orthopoxvirus       |wg        |             14|11          |                7247|
|Orthorubulavirus    |wg        |             28|9           |                1052|
|Respirovirus        |wg        |             43|16          |                 549|

Click [here](https://github.com/DOH-JDJ0303/vaper/blob/main/assets/reference_sets/EPITOME_2025-02-14.md) to learn more about these references.

### More Information:
See the [wiki](https://github.com/DOH-JDJ0303/VAPER/wiki) for more information.

## Quick Start
### Step 1: Prepare your samplesheet
> [!NOTE]
> Nextflow requires absolute file paths.
`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
sample01,sample01_R1_001.fastq.gz,sample01_R2_001.fastq.gz
sample02,sample02_R1_001.fastq.gz,sample02_R2_001.fastq.gz
```
### Step 2: Run VAPER
> [!NOTE]
> This will use the default reference set. You can provide your own refererences using the `reference` column in the samplesheet or using the `--refs` parameter.
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

