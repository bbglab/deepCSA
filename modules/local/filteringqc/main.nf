process FILTERING_REPORT {

    tag "$meta.id"
    label 'process_high'

    label 'deepcsa_core'

    input:
    tuple val(meta), path(maf)
    tuple val(meta2), path(somatic_maf)


    output:
    tuple val(meta), path("*.tsv")  ,                   emit: comparison_tables
    tuple val(meta2), path("*.pdf") ,   optional:true,  emit: comparison_plots
    path "versions.yml" ,                               topic: versions

    script:
    """
    compare_pre_post_filtering.py \\
                        --maf_path ${maf} \\
                        --somatic_maf ${somatic_maf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch contaminated_samples.tsv
    touch contaminated_samples.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
