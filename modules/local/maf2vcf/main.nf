process MAF_2_VCF {

    tag "${meta.id}"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"


    input:
    tuple val(meta), path (maf_file)

    output:
    path "*.vcf"        , emit: vcf_files
    path "versions.yml" , topic: versions


    script:
    def prefix = task.ext.prefix ?: ''
    prefix = "${meta.id}${prefix}"
    def args = task.ext.args ?: ""
    """
    deepcsa_maf2samplevcfs.py \\
                    --mutations-file ${maf_file} \\
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
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
