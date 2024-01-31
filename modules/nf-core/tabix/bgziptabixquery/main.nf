process TABIX_BGZIPTABIX_QUERY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta)     , path(input)
    tuple val(meta2)    , path(bedfile)

    output:
    tuple val(meta), path("*.gz")   , emit: subset
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = task.ext.extension ?: "${input.getExtension()}"
    def header = task.ext.header ?: "1"
    // TODO add a variable that tells if the input file will be compressed / will contain header and act accordingly
    """
    if [[ ${input} == *.gz ]]; then
        if [[ "${header}" == "1" ]]; then
            zcat ${input} | tail -n +2 | bgzip --threads ${task.cpus} -c $args > ${prefix}.tmp.${extension}.gz;
            tabix $args2 ${prefix}.tmp.${extension}.gz;
            cat <(zcat ${input} | head -n 1) <(tabix --regions ${bedfile} $args3 ${prefix}.tmp.${extension}.gz) | gzip > ${prefix}.${extension}.gz;
        else
            zcat ${input} | bgzip --threads ${task.cpus} -c $args > ${prefix}.tmp.${extension}.gz;
            tabix $args2 ${prefix}.tmp.${extension}.gz;
            tabix --regions ${bedfile} $args3 ${prefix}.tmp.${extension}.gz | gzip > ${prefix}.${extension}.gz;
        fi
    else
        if [[ "${header}" == "1" ]]; then
            tail -n +2 ${input} | bgzip --threads ${task.cpus} -c $args > ${prefix}.tmp.${extension}.gz;
            tabix $args2 ${prefix}.tmp.${extension}.gz;
            cat <(head -n 1 ${input}) <(tabix --regions ${bedfile} $args3 ${prefix}.tmp.${extension}.gz) | gzip > ${prefix}.${extension}.gz;
        else
            bgzip --threads ${task.cpus} -c $args ${input} > ${prefix}.tmp.${extension}.gz;
            tabix $args2 ${prefix}.tmp.${extension}.gz;
            tabix --regions ${bedfile} $args3 ${prefix}.tmp.${extension}.gz | gzip > ${prefix}.${extension}.gz;
        fi
    fi
    rm ${prefix}.tmp.${extension}.gz*;
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${input.getExtension()}.gz
    touch ${prefix}.${input.getExtension()}.gz.tbi
    touch ${prefix}.${input.getExtension()}.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}