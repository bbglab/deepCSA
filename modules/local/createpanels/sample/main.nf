process CREATESAMPLEPANELS {
    tag "$meta.id"
    label 'process_single'
    label 'time_low'

    conda "bioconda::pybedtools=0.9.1--py38he0f268d_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pybedtools:0.9.1--py38he0f268d_0' :
            'biocontainers/pybedtools:0.9.1--py38he0f268d_0' }"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //     'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta) , path(compact_captured_panel_annotation)
    tuple val(meta2), path(depths)
    val(min_depth)

    output:
    path("*.tsv")           , emit: sample_specific_panel
    path("*.bed")           , emit: sample_specific_panel_bed
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_panel4sample.py \\
                    ${compact_captured_panel_annotation} \\
                    all_samples.depths.tsv.gz \\
                    ${prefix} \\
                    ${min_depth};

    for sample_panel in \$(ls *${prefix}.tsv ); do
        bedtools merge \\
            -i <(
                tail -n +2 \$sample_panel | \\
                awk -F'\\t' '{print \$1, \$2, \$2}' OFS='\\t' | uniq
            ) | \\
            awk 'BEGIN { FS = "\\t"; OFS = "\\t" } { if (\$2 != \$3) { \$2 += 1 ; \$3 -= 1 } else { \$2 = \$2 ; \$3 = \$3 }; print }' > \${sample_panel%.tsv}.bed;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "TargetRegions"
    """
    touch *.tsv
    touch *.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
