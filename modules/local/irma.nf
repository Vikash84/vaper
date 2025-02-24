process IRMA {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'
    
    container "public.ecr.aws/o8h2f0o1/irma:1.1.4"

    input:
    tuple val(meta), val(ref_id), path(ref), path(reads)
    path module_template

    output:
    tuple val(meta), val(ref_id), path('*.fa.gz'),          emit: consensus, optional: true
    tuple val(meta), val(ref_id), path('results/*.bam'),    emit: bam, optional: true
    tuple val(meta), val(ref_id), path('results/logs/'),    emit: logs
    tuple val(meta), val(ref_id), path('results/figures/'), emit: figures
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}_${ref_id}"
    assembly_path = params.cons_assembly_type == 'plurality' ?  "results/${ref_id}.fasta" : ( params.cons_assembly_type == 'padded' ? "results/amended_consensus/results-${ref_id}.pad.fa" : "results/amended_consensus/results-${ref_id}.fa" )
    """
    #----- PREPARE IRMA -----#
    # Copy IRMA to workdir
    irma_path=\$(which IRMA)
    cp -r \${irma_path%IRMA} ./
    irma_path="\$(dirname \${irma_path#/})/IRMA"

    # Create module
    mod=\$(shuf -er -n20  {A..Z} {a..z} {0..9} | tr -d '\n')
    mv ${module_template} \${irma_path}_RES/modules/\${mod}

    #----- CONFIG OPTIONS -----#
    # Custom parameters
    echo 'MIN_AMBIG=${params.cons_allele_ratio}' > irma.config
    echo 'MIN_CONS_SUPPORT=${params.cons_allele_depth}' >> irma.config
    echo 'MIN_CONS_QUALITY=${params.cons_allele_qual}' >> irma.config
    echo -e 'SKIP_E=${ params.cons_assembly_elong ? '0' : '1' }' >> irma.config
    # Padded Assembly
    if [ "${params.cons_assembly_type}" == 'padded' ]
    then
        echo 'ASSEM_REF=1' >> irma.config
        echo 'ALIGN_AMENDED=1' >> irma.config
        echo 'PADDED_CONSENSUS=1' >> irma.config
    fi

    #----- PREPARE REFERENCE -----#
    # Create HMM profiles
    \${irma_path%IRMA}LABEL_RES/scripts/modelfromalign_Linux ${ref_id}_hmm -alignfile ${ref}
    mv ${ref_id}_hmm.mod \${irma_path}_RES/modules/\${mod}/profiles/
    # Copy the reference to the module
    zcat ${ref} > \${irma_path}_RES/modules/\${mod}/reference/consensus.fasta

    #----- RUN IRMA -----#
    # run IRMA
    \$irma_path \${mod} --external-config irma.config ${reads[0]} ${reads[1]} results

    #----- RENAME ASSEMBLY -----#
    if [ -f "${assembly_path}" ]
    then
        echo ">${prefix}" > ${prefix}.fa
        cat ${assembly_path} | grep -v '>' >> ${prefix}.fa
        gzip ${prefix}.fa
    fi

    # nextflow doesn't like me :(
    echo -e "${task.process}:\\n    IRMA: \$(IRMA | grep "(IRMA)" | cut -f 2 -d ',' | cut -f 2 -d ' ')" > versions.yml
    """
}
