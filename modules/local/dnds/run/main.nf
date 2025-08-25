process RUN_DNDS {
    tag "${meta.id}"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/dnds:latest'

    input:
    tuple val(meta), path(mutations_table), path(depths)
    tuple val(meta2), path(ref_cds)
    path covariates

    output:
    tuple val(meta), path("*.out.tsv*"), emit: results
    path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    dNdS_run.R --inputfile ${mutations_table} \\
                --outputfile ${prefix}.out.tsv \\
                --samplename ${prefix} \\
                --covariates ${covariates} \\
                --referencetranscripts ${ref_cds} \\
                --genedepth ${depths}
    # --cores ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv: 0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv: 0.1.0
    END_VERSIONS
    """
}
