#!/bin/bash

assembly=$1

tot_len=$(cat ${assembly} | grep -v ">" | tr -d '\t\n\r ' | wc -c)
atcg_count=$(cat ${assembly} | grep -v ">" | grep -o [ATCG] | wc -l)
n_count=$(cat ${assembly} | grep -v ">" | grep -o N | wc -l)
amb_count=$(cat ${assembly} | grep -v ">" | grep -o [RYSWKMBDHV] | wc -l)

echo "assembly_length,assembly_atcg_count,assembly_amb_count,assembly_n_count,assembly_gen_frac"
echo "${tot_len},${atcg_count},${amb_count},${n_count}" | awk 'BEGIN{FS = ","; OFS=","} {gen_frac=($2+$3)/$1} END{print $0,gen_frac}'