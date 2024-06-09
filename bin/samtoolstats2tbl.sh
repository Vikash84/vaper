#!/bin/bash

# samtoolstats2tbl.sh v1.0
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#----- HELP -----#
[ -z $1 ] && echo "samtoolstats2tbl.sh [path/to/samtoolstatsfile]" && exit

#----- INPUTS -----#
stats_file=$1

#----- GATHER METRICS -----#
reads_mapped=$(cat ${stats_file} | grep "SN" | grep "reads mapped:" | sed 's/.*://g' | sed 's/#.*//g' | tr -d '\t\n\r ')
bases_mapped=$(cat ${stats_file} | grep "SN" | grep "bases mapped:" | sed 's/.*://g' | sed 's/#.*//g' |  tr -d '\t\n\r ')
average_length=$(cat ${stats_file} | grep "SN" | grep "average length:" | sed 's/.*://g' | sed 's/#.*//g' |  tr -d '\t\n\r ')
average_quality=$(cat ${stats_file} | grep "SN" | grep "average quality:" | sed 's/.*://g' | sed 's/#.*//g' |  tr -d '\t\n\r ')

#----- TABULATE -----#
echo "reads_mapped,bases_mapped,mean_mapped_read_length,mean_mapped_read_quality"
echo ${reads_mapped},${bases_mapped},${average_length},${average_quality}