process NCBI_DATASETS {
    tag "${ref_id}"
    label 'process_low'

    input:
    tuple val(ref_id), val(accession)

    output:
    tuple val(ref_id), path("*.fna"), emit: assembly
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args   = task.ext.args ?: ''
    '''    
    # download assembly, unzip, and move for publishing
    datasets download genome accession !{accession} && unzip ncbi_dataset.zip && mv ncbi_dataset/data/*/*.fna ./

    # version info
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ncbi-datasets: \$(datasets --version | sed 's/.*: //g')
    END_VERSIONS
    '''
}
