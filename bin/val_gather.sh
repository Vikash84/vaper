#!/bin/bash

# val_gather.sh v1.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#----- INPUTS -----#
ALN=$1
METRIC=$2

# Help
[ -z $ALN ] && echo "validate.sh [path/to/alignment/file] [accuracy|precision] [prefix]" && exit

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
    | awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $i}} NF>p { p = NF } END { for(j=1; j<=p; j++) {str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j];} print str }}' \
    > transposed.txt

# extract sequence names
S1=$(cat $ALN | grep '>' | sed -n 1p | tr -d '>\n')
S2=$(cat $ALN | grep '>' | sed -n 2p | tr -d '>\n')


#----- CALCULATE ACCURACY -----#
# First sequence (S1) in the alignment is treated as the "truth"
# Second sequence (S2) in the alignment is treated as the "sample"
# True Positives (TP) = all sites in S2 that match S1
# True Negative (TN) = does not apply
# False Positives (FP) = sites in S2 that are absent in S1 (i.e., "A", "T", "C", or "G" in S2 and "-" in S1) or sites in S2 that differ from S1, excluding absent sites in S2
# False Negative (FN) = sites in S1 that are absent in S2 (i.e., "A", "T", "C", or "G" in S1 and "-" in S2)
# Accuracy = 100 * TP / (TP + FP + FN)

if [[ "${METRIC}" == "accuracy" ]]
then
    TP=$(cat transposed.txt | awk '$1 == $2 {print}' | wc -l)
    FP=$(cat transposed.txt | awk '$2 != "-" && $1 != $2 {print}' | wc -l)
    FN=$(cat transposed.txt | awk '$2 == "-" {print}' | wc -l)

    echo "Sample,Truth,TP,FP,FN,Result" > "${S2}_result.csv"
    echo -e "${S2}\t${S1}\t${TP}\t${FP}\t${FN}" | awk -v OFS=',' '{print $1,$2,$3,$4,$5, 100 * $3 / ($3+$4+$5)}' >> "${S2}_result.csv"
fi

#----- CALCULATE PRECISION -----#
# Sequences in alignment are treated as replicates
# True Positives (TP) = all sites that are the same in both sequences
# Predicted Positives (PP) = all sites in both sequences
# Precision = 100 * TP / PP

if [[ "${METRIC}" == "precision" ]]
then
    TP=$(cat transposed.txt | awk '$1 == $2 {print}' | wc -l)
    PP=$(cat transposed.txt | wc -l)

    echo "Sample1,Sample2,TP,PP,Result" > "${S1}-${S2}_result.csv"
    echo -e "${S1}\t${S2}\t${TP}\t${PP}" | awk -v OFS=',' '{print $1,$2,$3,$4, 100 * $3 / $4}' >> "${S1}-${S2}_result.csv"
fi

#-----CLEAN UP-----#
rm transposed.txt








