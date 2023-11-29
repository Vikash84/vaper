#!/bin/bash

fastp_json=$1

total_reads_before=$(jq '.summary.before_filtering.total_reads' ${fastp_json})
total_bases_before=$(jq '.summary.before_filtering.total_bases' ${fastp_json})
read1_mean_length_before=$(jq '.summary.before_filtering.read1_mean_length' ${fastp_json})
read2_mean_length_before=$(jq '.summary.before_filtering.read2_mean_length' ${fastp_json})
q30_rate_before=$(jq '.summary.before_filtering.q30_rate' ${fastp_json})

total_reads_after=$(jq '.summary.after_filtering.total_reads' ${fastp_json})
total_bases_after=$(jq '.summary.after_filtering.total_bases' ${fastp_json})
read1_mean_length_after=$(jq '.summary.after_filtering.read1_mean_length' ${fastp_json})
read2_mean_length_after=$(jq '.summary.after_filtering.read2_mean_length' ${fastp_json})
q30_rate_after=$(jq '.summary.after_filtering.q30_rate' ${fastp_json})

echo "total_reads_before,total_bases_before,read1_mean_length_before,read2_mean_length_before,q30_rate_before,total_reads_after,total_bases_after,read1_mean_length_after,read2_mean_length_after,q30_rate_after"
echo "${total_reads_before},${total_bases_before},${read1_mean_length_before},${read2_mean_length_before},${q30_rate_before},${total_reads_after},${total_bases_after},${read1_mean_length_after},${read2_mean_length_after},${q30_rate_after}"
