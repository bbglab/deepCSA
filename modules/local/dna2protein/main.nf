process DNA_2_PROTEIN_MAPPING {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"


    input:
    tuple val(meta) , path(mutations_file)
    tuple val(meta2), path(panel_file)
    tuple val(meta3), path(all_samples_depths)


    output:
    tuple val(meta2), path("depths_per_position_exon_gene.tsv") , emit: depths_exons_positions
    tuple val(meta2), path("panel_exons.bed4.bed")              , emit: panel_exons_bed
    path  "versions.yml"                                        , topic: versions


    script:
    """
    cut -f 1,2,6 ${panel_file} | uniq > ${meta2.id}.panel.unique.tsv
    panels_computedna2protein.py \\
                --mutations-file ${mutations_file} \\
                --consensus-file ${meta2.id}.panel.unique.tsv \\
                --depths-file ${all_samples_depths}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
