process TAR2REFS {
    tag "${refs}"
    label 'process_low'
    stageInMode 'copy'

    input:
    path refs

    output:
    path "*/*.csv",  emit: refsheet
    path "*/references/*", emit: refs


    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/waphlviral/bin/
    """
    tar -xvzf ${refs}
    """
}
