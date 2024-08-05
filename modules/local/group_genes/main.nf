process GROUP_GENES {
    tag "groups"
    label 'process_low'

    // TODO revise this conda package name
    conda "pandas:1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(mutations_table)
    path(features_table)
    tuple val(meta2), path(hotspots_file)

    output:
    path("genes2group_out.json")                        , emit: json_genes
    path("pathway_groups_out.json") , optional : true   , emit: json_pathways
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "groups"
    def separator = task.ext.separator ?: "tab"
    def hotspots = task.ext.hotspots ? "True" : "False"
    def cmd_custom = task.ext.custom ? "${features_table} ${separator} pathway_groups_out.json" : ""
    """
    awk 'NR>1 {print \$1}' ${mutations_table} | sort -u  > gene_list.txt

    features_2group_genes.py gene_list.txt ${hotspots} ${hotspots_file} genes2group_out.json ${cmd_custom}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "groups"
    """
    touch genes2group_out.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
    // features_1table2groups.py ${features_table} ${separator} features_table_information.json
