process NEXTCLADE_CONFIG {
    tag "${prefix}"
    label 'process_low'

    container 'docker.io/jdj0303/vaper-base:beta'

    input:
    tuple val(ref_id), path(reference)
    path db_template

    output:
    tuple val(ref_id), path('ref.fa'), path("${prefix}-nc-config.json"), emit: config
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${ref_id}"
    """
    # concat reference contigs
    echo ">${ref_id}" > ref.fa
    zcat ${reference} | grep -v ">" | tr -d '\n\t\r ' >> ref.fa
    # Calculate QC metrics based on the reference length
    len=\$(cat ref.fa | grep -v '>' | tr -d '\t\n\r' | wc -c)
    miss_th=\$(echo \${len} | awk -v var=$params.nc_miss '{print int(var*\$0)}')
    miss_bias=\$(echo \${len} | awk -v var=$params.nc_miss_bias '{print int(var*\$0)}')
    mix_th=\$(echo \${len} | awk -v var=$params.nc_mixed '{print int(var*\$0)}')

    # configure QC thresholds for Nextclade
    nc-qc-config.sh ${db_template} ${prefix}-nc-config.json \$miss_th \$miss_bias \$mix_th
    """
}
