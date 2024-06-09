#!/bin/bash

# nc-qc-config.sh v1.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#------ REQUIREMENTS -----#
# check that jq is installed
if ! command -v jq &> /dev/null
then
    echo "This script requires 'jq'. Please make sure it is installed."
    exit 1
fi

#----- HELP -----#
if [[ -z "$1" ]] || [[ "$1" == "-h" ]] || [[ "$1" == "-help" ]] || [[ "$1" == "--help" ]]
then
    echo "nc-qc-config.sh [input file] [output file] [missingDataThreshold] [missingData.scoreBias] [mixedSitesThreshold]"
fi

# inputs
input=$1
output=$2
miss_th=$3
miss_bias=$4
mix_th=$5

set -euxo pipefail

cat ${input} \
    | jq ".qc.missingData.missingDataThreshold = ${miss_th}" \
    | jq ".qc.missingData.scoreBias = ${miss_bias}" \
    | jq ".qc.mixedSites.mixedSitesThreshold = ${mix_th}" \
    > ${output}