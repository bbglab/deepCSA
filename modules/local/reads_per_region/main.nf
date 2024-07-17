process READS_PER_REGION {

    tag "$meta.id"
    label 'cpu_medium_high'
    label 'process_high_memory'
    label 'time_low'


    // TODO revise if we want a custom container
    container 'docker.io/ferriolcalvet/oncodrivefml:latest'

    input:
    tuple val(meta) , path(pileup), path(pileupindex), path(fragments_coords)
    tuple val(meta2), path(regions)

    output:
    tuple val(meta), path("*.reads_x_region.tsv.gz")    , emit: read_counts
    tuple val(meta), path("*.custom.bed")               , emit: redesigned_bed
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO decide whether any of these would be useful in here
    // zcat ${pileup} | bgzip -c > ${prefix}.indexed.tsv.tmp.gz; tabix -s 1 -b 2 -e 2 ${prefix}.indexed.tsv.tmp.gz
    // tabix -s 1 -b 2 -e 2 ${prefix}.indexed.tsv.tmp.gz
    """
    reads_x_region.py ${pileup} ${fragments_coords} ${regions} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reads_x_region.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
