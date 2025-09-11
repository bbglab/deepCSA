process CONCAT_PROFILES {
    tag "${meta.id}"
    label 'process_medium'
    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(mutation_profiles_list)

    output:
    path("cosine_similarity_heatmap.png")      , emit: heatmap
    path("cosine_similarity_clustermap.png")   , emit: clustermap
    path("*.compiled_profiles.tsv")            , emit: compiled_profiles
    path "versions.yml"                        , topic: versions

    script:
    """
    ls ${mutation_profiles_list} > mutation_profiles_list.txt
    concat_profiles.py \\
            --mutation-profiles-list mutation_profiles_list.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: $(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch cosine_similarity_heatmap.png
    touch cosine_similarity_clustermap.png
    touch dummy.compiled_profiles.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: $(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
