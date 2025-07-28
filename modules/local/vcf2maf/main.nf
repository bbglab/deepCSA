process VCF2MAF {
    tag "$meta.id"

    label 'cpu_low'
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta) , path(vcf)
    tuple val(meta2), path(annotation)

    output:
    tuple val(meta), path("*.tsv.gz")  , emit: maf
    path "versions.yml"                , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    vcf2maf.py 
        --vcf ${vcf} \\
        --sampleid ${prefix} \\
        --annotation_file ${annotation} \\
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
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
