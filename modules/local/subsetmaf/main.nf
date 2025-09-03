process SUBSET_MAF {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(mut_files)

    output:
    tuple val(meta), path("*.mutations.tsv")  , optional : true , emit: mutations
    path "versions.yml"                                         , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def output_prefix = task.ext.output_prefix ?: ""
    def filters = task.ext.filters ?: ""
    def output_format = task.ext.output_fmt ?: ""
    def min_muts = task.ext.minimum_mutations ? "--min_mutations ${task.ext.minimum_mutations}" : ""
    """
    cat > mutations_subset.conf << EOF
    {
        ${filters}
    }
    EOF

    cat > output_formats.conf << EOF
    {
        ${output_format}
    }
    EOF

    subset_maf.py \\
                    --sample_name ${prefix} \\
                    --mut_file ${mut_files} \\
                    --out_maf ${prefix}${output_prefix}.mutations.tsv \\
                    --json_filters mutations_subset.conf \\
                    --req_fields output_formats.conf \\
                    ${min_muts} \\
                    ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.mutations.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
