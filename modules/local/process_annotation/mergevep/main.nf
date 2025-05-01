process MERGE_VEP {
    label 'process_low'

    conda "bioconda::bcftools=1.15.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0' : 'biocontainers/bcftools:1.15.1--h0ea216a_0' }"

    input:
    tuple val(meta), path(vep_files)

    output:
    tuple val(meta), path("*.merged.tab.gz"), emit: tab
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Sort input files by chromosome number
    sorted_files=\$(ls ${vep_files} | sort -V)

    # Concatenate sorted files
    bcftools concat -a \$sorted_files | bgzip > ${prefix}.merged.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.merged.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tools: gzip \$(gzip --version | head -n1)
    END_VERSIONS
    """
}