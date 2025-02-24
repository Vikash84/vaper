process FINALIZE_REFS {
    tag "${ref_id}"
    label 'process_low'

    input:
    tuple val(ref_id), path(ref)

    output:
    tuple val(ref_id), path("*.fa.gz", includeInputs: true), emit: refs

    when:
    task.ext.when == null || task.ext.when

    script:
    def cat = ref.name.endsWith('.gz') ? 'zcat' : 'cat'
    """
    # Set header to match the file name and concatenate multi-fasta files
    echo ">${ref_id}" > tmp.fa
    ${cat} ${ref} | grep -v '>' | tr -d '\n\r\t ' >> tmp.fa
    cat tmp.fa | gzip > ${ref_id}.fa.gz
    """
}
