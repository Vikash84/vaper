# QUICK INSTALL

## PACKAGE

- LABEL, Lineage Assignment by Extended Learning: <https://wonder.cdc.gov/amd/flu/label/>
- IRMA, Iterative Refinement Meta-Assembler: <https://wonder.cdc.gov/amd/flu/irma/>

## HARDWARE

We recommend a single multi-core machine with no fewer than 2 cores (8 to 12 cores work best) and at least 4 GB of RAM.  Runtime is impacted by the number of cores available on a machine. Use with MacOS and Apple ARM silicon requires installation of Rosetta.

## LABEL SOFTWARE PRE-REQUISITES

- Linux (64-bit) or Mac OS X (64-bit)
- BASH version 3 or later
- Standard utilities: sleep, cut, paste, jobs, zip, env, cat, cp, getopts.
  - License: GPL (any)
- Perl version 5.16 or later
  - Standard includes: Getopt::Long, File::Basename
- License: GPL (any)

## ADDITIONAL IRMA REQUIREMENTS

- Linux (64-bit only) or Mac OS X (64-bit)
- Perl 5.16 modules (standard includes): Storable, File::Path, POSIX
- R version 3.6+
- Standard *nix utilities like gzip

## INSTALLATION

1. Unzip the archive "flu-amd.zip".
2. Move the contents of the directory to a folder in your shell PATH, otherwise, add the directory flu-amd directory path to your path.
3. Restart your terminal emulator. LABEL and IRMA are now installed.
4. To test LABEL and IRMA, cd to the folder "tests" and run the script "test_run.sh":

```{bash}
cd tests
./test_run.sh
```

## Alternative Dockerhub Image Installation

Get the latest official IRMA version from Dockerhub: <https://hub.docker.com/r/cdcgov/irma/tags>

## LICENSE INFORMATION

Please visit <http://wonder.cdc.gov/amd/flu/irma/disclaimer.html> for a list of components, their licenses, and for the GPL for IRMA & LABEL.
