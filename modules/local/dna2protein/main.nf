process DNA_2_PROTEIN_MAPPING {
    tag "$meta.id"

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"


    input:
    tuple val(meta) , path(mutations_file)
    tuple val(meta2), path(panel_file)

    output:
    tuple val(meta2), path("*.mapping.tsv") , emit: mapping
    path  "versions.yml"                   , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    """
    cut -f 1,2,6,9 ${panel_file} | cut -f-4 | uniq > ${meta2.id}.panel.unique.tsv
    panels_computedna2protein.py \\
                --maf ${mutations_file} \\
                --consensus-file ${meta2.id}.panel.unique.tsv \\
                --output ${prefix}.mapping.tsv;

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
