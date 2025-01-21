process SM2REFS {
    label 'process_low'

    input:
    tuple val(meta), path(sm_species)

    output:
    tuple val(meta), path("${prefix}.sm-refs.csv"), emit: refs, optional: true

    when:
    task.ext.when == null || task.ext.when

    version = "1.0"
    script:
    prefix = "${meta.id}"
    """
    if [[ -s ${sm_species} ]]
    then
        # get list of NCBI accessions & species names
        cat ${sm_species} | tr ';' '\n' | cut -f 1 -n -d ' ' > ACCESSIONS
        cat ${sm_species} | tr ';' '\n' | cut -c 17- | tr ' ' '_' | tr -d '%' | tr '/()[]' '_' > SPECIES
        # combine into CSV
        echo "taxa,accession" > ${prefix}.sm-refs.csv
        paste -d ',' SPECIES ACCESSIONS >> ${prefix}.sm-refs.csv
    fi
    """
}
