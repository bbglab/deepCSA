process ONCODRIVE3D_PLOT_CHIMERAX {
    tag "$meta.id"
    label 'process_medium'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // TODO pending to push the container somewhere and be able to retrieve it
    container 'docker.io/ferriolcalvet/o3d_chimerax:latest'


    input:
    tuple val(meta), path(genes_csv)
    tuple val(meta), path(pos_csv)
    tuple val(meta), path(miss_prob_json)
    tuple val(meta), path(seq_df_tsv)
    path(datasets)

    output:
    tuple val(meta), path("chimerax/attributes/**.defattr")       , emit: chimerax_defattr
    tuple val(meta), path("chimerax/plots/**.png")                , emit: chimerax_plot
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Increase pixel_size to decrease resolution and speed up png generation 
    """
    python3 /o3d_chimera_plot.py \
        -o $prefix \\
        -g $genes_csv \\
        -p $pos_csv \\
        -d $datasets \\
        -s $seq_df_tsv \\
        -c $prefix \\
        --fragmented_proteins \\
        --pixel_size 0.1   


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 2.0
    END_VERSIONS
    """ 

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 2.0
    END_VERSIONS
    """
}

