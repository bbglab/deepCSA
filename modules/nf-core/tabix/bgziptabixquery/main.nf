process TABIX_BGZIPTABIX_QUERY {
    cache true
    
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta)   , path(input)
    tuple val(meta2)  , path(bedfile)

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
    // TODO fix the issue with post generation indexing, only possible in the non-header files
    """
    if [[ ${input} == *.gz ]]; then
        if [[ "${header}" == "1" ]]; then
            echo "mode compressed and header";
            zcat ${input} | tail -n +2 | sort -k1,1 -k2,2n | bgzip --threads ${task.cpus} -c ${args} > ${prefix}.tmp.${extension}.gz;
            tabix ${args2} ${prefix}.tmp.${extension}.gz;
            cat <(zcat ${input} | head -n 1) <(tabix --regions ${bedfile} ${args3} ${prefix}.tmp.${extension}.gz) | bgzip --threads ${task.cpus} -c > ${prefix}.${extension}.gz;
        elif [[ "${header}" == "pile" ]]; then
            echo "mode pileup";
            tabix ${args2} ${input};
            tabix --regions ${bedfile} ${args3}  ${input} | bgzip --threads ${task.cpus} -c > ${prefix}.${extension}.gz;
        else
            echo "mode compressed without header";
            zcat ${input} | sort -k1,1 -k2,2n | bgzip --threads ${task.cpus} -c ${args} > ${prefix}.tmp.${extension}.gz;
            tabix ${args2} ${prefix}.tmp.${extension}.gz;
            tabix --regions ${bedfile} ${args3} ${prefix}.tmp.${extension}.gz | bgzip --threads ${task.cpus} -c > ${prefix}.${extension}.gz;
        fi
    else
        if [[ "${header}" == "1" ]]; then
            echo "mode uncompressed with header";
            tail -n +2 ${input} | sort -k1,1 -k2,2n | bgzip --threads ${task.cpus} -c ${args} > ${prefix}.tmp.${extension}.gz;
            tabix ${args2} ${prefix}.tmp.${extension}.gz;
            cat <(head -n 1 ${input}) <(tabix --regions ${bedfile} ${args3} ${prefix}.tmp.${extension}.gz) | bgzip --threads ${task.cpus} -c > ${prefix}.${extension}.gz;
        else
            echo "mode uncompressed without header";
            sort -k1,1 -k2,2n ${input} | bgzip --threads ${task.cpus} -c ${args} > ${prefix}.tmp.${extension}.gz;
            tabix ${args2} ${prefix}.tmp.${extension}.gz;
            tabix --regions ${bedfile} ${args3} ${prefix}.tmp.${extension}.gz | bgzip --threads ${task.cpus} -c > ${prefix}.${extension}.gz;
        fi
    fi
    rm -f ${prefix}.tmp.${extension}.gz*;
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