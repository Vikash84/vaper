#!/bin/bash

# val_gather.sh v1.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#----- HELP -----#
[ -z $1 ] && echo "val_gather.sh [path/to/alignment/file] [accuracy|precision] [prefix] [seq1_name] [seq2_name]" && exit

#----- INPUTS -----#
ALN=$1
METRIC=$2
S1=$3
S2=$4

# Check that the alignment exists and only has two sequences
! [ -f $ALN ] && echo "Error: Please supply a valid alignment file path." && exit 1
[ $(cat $ALN | grep '>' | wc -l) != 2 ] && echo "Error: Alignment file should contain only two sequences!" && exit 1
[[ "$METRIC" != "accuracy" ]] && [[ "$METRIC" != "precision" ]] && echo "Error: Second field must be 'accuracy' or 'precision'." && exit 1

set -euxo pipefail
#----- FORMAT ALIGNMENT -----#
# Transpose alignment so that each sequences is represented as a column
cat ${ALN} \
    | sed 's/>.*$/@&@/g' \
    | tr -d '\n' \
    | tr '@' '\n' \
    | grep -v '>' \
    | tail -n +2 \
    | sed 's/./&\t/g' \
    | sed 's/\t$//g' \
    | awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $i}} NF>p { p = NF } END { for(j=1; j<=p; j++) {str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j];} print toupper(str) }}' \
    > transposed.txt

# get only sites that were called in both sequences
cat transposed.txt | awk '$1 != "N" && $2 != "N" {print}' > comp.txt

#----- CALCULATE ACCURACY -----#
# First sequence (S1) in the alignment is treated as the "truth".
# Second sequence (S2) in the alignment is treated as the "sample".
# Only lines containing called bases in both S1 and S2 are considered.
# Total sites = sites in S1 and S2
# Compared sites = sites compared between S1 and S2
# Correct sites = sites in S2 that match S1
# Incorrect sites = sites in S2 that do not match S1
# Accuracy = 100 * Correct sites / Compared sites

if [[ "${METRIC}" == "accuracy" ]]
then
    TOTAL=$(cat transposed.txt | wc -l)
    COMPARED=$(cat comp.txt | wc -l)
    CORRECT=$(cat comp.txt | awk '$2 == $1 {print}' | wc -l)

    echo "Sample,Truth,Total,Compared,Correct,Incorrect,Accuracy" > "${S2}_result.csv"
    echo -e "${S2}\t${S1}\t${TOTAL}\t${COMPARED}\t${CORRECT}" | awk -v OFS=',' '{print $1,$2,$3,$4,$5,$4-$5,100 * $5 / ($4)}' >> "${S2}_result.csv"
fi

#----- CALCULATE PRECISION -----#
# Sequences in alignment (S1 & S2) are treated as replicates
# Only lines containing called bases in both S1 and S2 are considered.
# Total sites = sites in S1 and S2
# Compared sites = sites compared between S1 and S2
# Agreements = compared sites that are the same in S1 and S2
# Disagreements = compared sites that differ in S1 and S2
# Precision = 100 * Agreements / Compared sites

if [[ "${METRIC}" == "precision" ]]
then
    TOTAL=$(cat transposed.txt | wc -l)
    COMPARED=$(cat comp.txt | wc -l)
    AGREEMENTS=$(cat comp.txt | awk '$1 == $2 {print}' | wc -l)


    echo "Sample1,Sample2,Total,Agreements,Disagreements,Precision" > "${S1}-${S2}_result.csv"
    echo -e "${S1}\t${S2}\t${TOTAL}\t${COMPARED}\t${AGREEMENTS}" | awk -v OFS=',' '{print $1,$2,$3,$4,$5,$4-$5,100 * $5 / $4}' >> "${S1}-${S2}_result.csv"
fi

#-----CLEAN UP-----#
rm transposed.txt comp.txt








