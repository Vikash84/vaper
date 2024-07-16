# VAPER: Viral Assembly from Probe-based EnRichment
VAPER is a viral assembly pipeline. Key features include:
-  Builds assemblies from probe enrichment (a.k.a hybrid capture/enrichment), shotgun metagenomic, and tiled-amplicon sequence data
-  Automated reference selection (Detects what is in your sample)
-  Can generate multiple assemblies per sample (Useful for co-infections)
-  Can utilize iVar or IRMA assemblers (IRMA modules built on the fly!)
-  Reads associated with each assembly are exported for downstream use
-  Provides a metagenomic summary of dominant viral taxa

VAPER comes with comprehensive reference sets for the following viral species (created using [EPITOME](https://github.com/DOH-JDJ0303/epitome)):
**Taxon**|**Segments**|**Input Sequences**|**Data Source**
-----|-----|-----|-----
Influenza A|1-8|78703 (per segment)|GISAID
Influenza B|1-8|17401 (per segment)|GISAID
Measles morbillivirus|wg|890|NCBI
Mumps orthorubulavirus|wg|1343|NCBI
Lyssavirus rabies|wg|2607|NCBI
Norovirus|wg|1662|NCBI
Respiratory Syncytial Virus|wg|15273|GISAID
West Nile virus|wg|1993|NCBI
Enterovirus D68|wg|590|NCBI
Hepacivirus|wg|1245|NCBI
Hepatovirus|wg|131|NCBI
Monkeypox virus|wg|2129|NCBI
Severe acute respiratory syndrome coronavirus|wg|2000 (random)|NCBI


  
See the [wiki](https://github.com/DOH-JDJ0303/VAPER/wiki) for more information.

VAPER was originally written by Jared Johnson for the Washington State Department of Health.

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
VAPER would not be possible without the following people:\
Holly Halstead (WA PHL) (check out her tool [varcraft](https://github.com/DOH-HNH0303/varcraft)!)\
Pauline Trinh (WA MEP)\
Allison Black (WA MEP)\
Stephanie Lunn (WA MEP)\
Kristen Waterman (WA PHL)
