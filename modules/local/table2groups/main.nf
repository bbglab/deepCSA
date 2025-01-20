process TABLE_2_GROUP {

    tag "groups"
    label 'process_low'

    // TODO revise this conda package name
    conda "pandas:1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path(features_table)

    output:
    path("samples.json")                        , emit: json_samples
    path("groups.json")       , optional : true , emit: json_groups
    path("all_groups.json")                     , emit: json_allgroups
    path "versions.yml"                         , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "groups"
    def separator = task.ext.separator ?: "comma"
    def features = task.ext.features ?: ""
    """
    cat > features_table_information.json << EOF
    {
        ${features}
    }
    EOF

    features_1table2groups.py ${features_table} ${separator} features_table_information.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "groups"
    """
    touch all_samples.cohort.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
