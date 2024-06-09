#!/bin/bash

# val_pair.sh v1.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#----- HELP -----#
[ -z $1 ] && echo "val_pait.sh [path/to/seqs1.fa] [path/to/seqs2.fa] [prefix]" && exit

#----- INPUTS -----#
F1=$1
F2=$2
PREFIX=$3

#----- CLEAN-UP -----#
# Remove any extras in the fasta header
cat $F1 | awk '{print $1}' > f1.fa
cat $F2 | awk '{print $1}' > f2.fa

#----- COMPARE -----#
# Compare sequnces in each fasta using Mash
mash sketch -i -o f1 f1.fa
mash sketch -i -o f2 f2.fa
mash dist f1.msh f2.msh > ava.txt

#------ PAIR FASTA_1 TO FASTA_2 -----#
# Find sequences in fasta_1 within 10% identity with fasta_2 
[ -f $PREFIX.pairs.txt ] && rm $PREFIX.pairs.txt
for C1 in $(cat f1.fa | grep '>' | tr -d '>')
do
    cat ava.txt | awk -v c1=$C1 '$1 == c1 {print}' | sort -gk 3 | sed -n 1p | awk -v OFS='\t' '{if($3 < 0.1){ print $1,$2 }else{ print $1,"null" }}' >> $PREFIX.pairs.txt
done

#----- GET REMAINDER -----#
# Find any sequences in fasta_2 that are missing in fasta_1
for C2 in $(cat f2.fa | grep '>' | tr -d '>')
do
    status=$(cat $PREFIX.pairs.txt | awk -v c2=$C2 '$2 == c2 {print "found"}')
    if [[ "$status" != "found"  ]]
    then
        echo -e "null\t$C2" >> $PREFIX.pairs.txt
    fi
done

#----- CREATE PAIRS -----#
# Get all sequences into a tabulated format [ <sequence header> <sequence> ]
cat f1.fa f2.fa | sed 's/>.*$/@&@/g' \
    | tr -d '\n' \
    | tr '@' '\n' \
    | tail -n +2  \
    | tr -d '>' \
    | paste - - > all.txt

# Create multi-fasta files containing identified pairs.
for PAIR in $(cat $PREFIX.pairs.txt | tr '\t' ',')
do
    C1=$(echo $PAIR | cut -f 1 -d ',')
    C2=$(echo $PAIR | cut -f 2 -d ',')

    if [[ "$C1" != "null" ]] && [[ "$C2" != "null" ]]
    then
        cat all.txt | awk -v c1=$C1 -v c2=$C2 '$1 == c1 || $1 == c2 {print}' | sed 's/^/>/g' | tr '\t' '\n' > $PREFIX."${C1}_${C2}.pair.fa"
    fi
done

