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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def batch = task.ext.batch ?: "${meta.batch}"
    def level = task.ext.level ?: "high"
    def all_molecules_dp = task.ext.all_molecules_dp ?: "false"
    """
    vcf2maf.py ${vcf} ${prefix} ${batch} ${level} ${annotation} ${all_molecules_dp};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
