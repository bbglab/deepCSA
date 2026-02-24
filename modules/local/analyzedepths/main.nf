process ANALYZE_DEPTHS_GROUPS {

    tag "groups"
    label 'process_low'

    label 'deepcsa_core'

    input:
    path(features_table)
    // add another path with the depths per gene per sample
    // optionally add another one with the depth per sample

    output:
    // the main outputs will be the PDFs
    // path("samples.json")      , emit: json_samples
    path "versions.yml"       , topic: versions

    script:
    def separator = task.ext.separator ?: "comma"
    def custom_groups = task.ext.feature_groups ? "--groups \"${task.ext.feature_groups}\" " : ""
    def unique_identifier = task.ext.unique_identifier ? "--unique-identifier ${task.ext.unique_identifier}" : ""

    // <your_script>.py should be in bin/ and you should make sure it can be executed (check permissions of the file and add shebang if needed)
    """
    <your_script>.py \\ 
                --table-filename ${features_table} \\
                --separator ${separator} \\
                ${unique_identifier} \\
                ${custom_groups}
                // TODO add the missing parameters
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
