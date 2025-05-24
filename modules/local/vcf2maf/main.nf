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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def batch = task.ext.batch ?: "${meta.batch}"
    def level = task.ext.level ?: "high"
    def all_molecules_dp = task.ext.all_molecules_dp ?: "false"
    // TODO reimplement it with click
    // TODO level and all_molecules can be defined in the modules.config file think about making all molecules mandatory
    """
    vcf2maf.py ${vcf} ${prefix} ${batch} ${level} ${annotation} ${all_molecules_dp};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
