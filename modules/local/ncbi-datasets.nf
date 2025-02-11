process NCBI_DATASETS {
    tag "${accession}"
    label 'process_low'
    container 'public.ecr.aws/o8h2f0o1/ncbi-datasets:16.15.0'

    input:
    val accession

    output:
    tuple val(accession), path("*.fa.gz"), emit: assembly, optional: true
    path "refinfo.csv",                    emit: refsheet, optional: true
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    """
    # Download assembly, unzip, and move for publishing
    datasets download genome accession ${accession} && unzip ncbi_dataset.zip
    # Check if it is a complete genome
    if grep -q 'Complete Genome' ncbi_dataset/data/assembly_data_report.jsonl
    then
        # Split multi-fastas
        cat ncbi_dataset/data/*/*.fna | awk '/^>/ { id=substr(\$1, 2); print ">"id > id".fa" } !/^>/ { print >> id".fa" }'
        gzip *.fa

        # Extract reference data
        echo "assembly,ref_info" > refinfo.csv
        cat ncbi_dataset/data/*/*.fna | grep '>' | tr -d '>' | tr ',' ';' | sed 's/ /.fa.gz,/1' >> refinfo.csv
    fi
    
    # version info
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ncbi-datasets: \$(datasets --version | sed 's/.*: //g')
    END_VERSIONS
    """
}
