process IRMA {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'
    
    container "public.ecr.aws/o8h2f0o1/irma:1.1.4"

    input:
    tuple val(meta), path(refs), path(reads)
    path module_template

    output:
    tuple val(meta), path('*.fa'),     emit: consensus
    tuple val(meta), path('results/*.bam'),    emit: bam
    tuple val(meta), path('results/logs/'),    emit: logs
    tuple val(meta), path('results/figures/'), emit: figures
    //path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.id}"

    """
    # copy IRMA to workdir
    irma_path=\$(which IRMA)
    cp -r \${irma_path%IRMA} ./
    irma_path="\$(dirname \${irma_path#/})/IRMA"

    # create module
    mod=\$(shuf -er -n20  {A..Z} {a..z} {0..9} | tr -d '\n')
    mv ${module_template} \${irma_path}_RES/modules/\${mod}

    # apply config options
    ## choose an assembly type
    case "${params.cons_type}" in
        'plurality')
            echo 'Returning plurality assembly'
            assembly_path='results/'
            assembly_ext=".fasta" 
            ;;
        'amended')
            echo 'Returning amended assembly'
            assembly_path='results/amended_consensus/results-'
            assembly_ext=".fa" 
            ;;
        'padded')
            echo -e 'ASSEM_REF=1\nALIGN_AMENDED=1\nPADDED_CONSENSUS=1' >> irma.config
            assembly_path='results/amended_consensus/results-'
            assembly_ext=".pad.fa"
            ;;
        *)
            echo "ERROR: '--cons_type' must be 'plurality', 'amended', or 'padded'" && exit 1
            ;;
    esac
    ## set QC thresholds
    echo -e 'MIN_AMBIG=${params.cons_ratio}\nMIN_CONS_SUPPORT=${params.cons_depth}\nMIN_CONS_QUALITY=${params.cons_qual}\nDEL_TYPE=${ params.cons_amb == 'N' ? 'NNN' : params.cons_amb }' >> irma.config
    ## set elongation option
    echo -e 'SKIP_E=${ params.cons_elong ? '1' : '0' }' >> irma.config
         
    # prepare references
    for REF in \$(echo ${refs.join(',')} | tr ',' '\n')
    do
        ## concatenate multi-fasta references and change fasta header to match the fasta name
        ## this mimics the concat step performed by iVar while also allowing for grouping of data in the main workflow
        echo ">\${REF%%.*}" > \${REF%.gz}
        zcat \${REF} | grep -v '>' | tr -d '\n\r\t ' >> \${REF%.gz}
        
        ## create HMM profiles for each reference
        ls \${irma_path%IRMA}LABEL_RES/scripts/modelfromalign_Linux
        \${irma_path%IRMA}LABEL_RES/scripts/modelfromalign_Linux \${REF%%.*}_hmm -alignfile \${REF%.gz}
        mv \${REF%%.*}_hmm.mod \${irma_path}_RES/modules/\${mod}/profiles/

        # add the concatenated reference to a multi-fasta in the IRMA module
        cat \${REF%.gz} >> \${irma_path}_RES/modules/\${mod}/reference/consensus.fasta

        # remove the concatenated reference for easy reporting
        rm \${REF%.gz}
    done

    # run IRMA
    \$irma_path \${mod} --external-config irma.config ${reads[0]} ${reads[1]} results

    # rename assemblies
    for F in \$(ls results/*.fasta)
    do
        # get file name
        FILE=\${F##*/}
        REF_ID=\${FILE%.fasta}
        PREFIX="${prefix}_\${REF_ID}"

        # return fasta 
        assembly="\${assembly_path}\${REF_ID}\${assembly_ext}"
        cat \${assembly} > \${PREFIX}.fa
    done

    # clean up
    #rm -r \${irma_path}_RES/modules/\${mod}
    """
}
