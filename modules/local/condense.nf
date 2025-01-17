process CONDENSE {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(meta), path(assemblies), path(read_stats)

    output:
    tuple val(meta), path("*.fa.gz", includeInputs: true),   emit: assembly
    tuple val(meta), path("${prefix}.condense_summary.csv"), emit: summary
    tuple val(meta), path("${prefix}.condense_dist.csv"),    emit: dist
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${meta.id}"
    """
    # combine read stats into single file
    cat ${read_stats} | grep '#rname' | sort | uniq > read_stats.tsv # head -n 1 fails
    cat ${read_stats} | grep -v '#rname' >> read_stats.tsv
    # combine sequences into single file
    cat ${assemblies} > seqs.fa
    # run script
    vaper-condense.py \\
        --fasta seqs.fa \\
        --stats read_stats.tsv \\
        --dist_threshold ${params.cons_cond_dist} \\
        --prefix "${prefix}"
    
    # rename output
    mv condensed.csv ${prefix}.condense_summary.csv
    mv dists.csv ${prefix}.condense_dist.csv
    # remove condensed sequences
    rm seqs.fa
    for s in \$(cat ${prefix}.condensed.csv | tr -d '"' | tail -n +2 | awk -v FS=',' '\$2 != "" {print \$1}' | uniq)
    do
        rm \${s}.fa
    done
    gzip *.fa

    # version info
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        vaper-condense: \$(vaper-condense.py --version)
    END_VERSIONS
    """
}