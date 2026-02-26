
process SIGNATURES_HDP_TO_SIGPROFILER {
    tag "$meta.id"
    label 'process_medium'
    label 'deepcsa_core'


    input:
    tuple val(meta), path (signatures)

    output:
    tuple val(meta), path("*.tsv")  , emit: signatures_for_sp
    path "versions.yml"             , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    python <<CODE
import pandas as pd
data = pd.read_table("${signatures}")
data = data.sort_index().reset_index()
data.columns = ["Type"] + [f"HDP_{x}" for x in data.columns[1:]]
data.to_csv("${signatures.baseName}.sp_formatted.tsv", sep='\\t', header=True, index=False)
CODE
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}