process TABLE_2_GROUP {

    tag "groups"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    path(features_table)

    output:
    path("samples.json")                        , emit: json_samples
    path("groups.json")       , optional : true , emit: json_groups
    path("all_groups.json")                     , emit: json_allgroups
    path "versions.yml"                         , topic: versions


    script:
    def separator = task.ext.separator ?: "comma"
    def features = task.ext.features ?: ""
    // TODO reimplement with click
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
    """
    touch samples.json groups.json all_groups.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
