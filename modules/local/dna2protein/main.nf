process DNA_2_PROTEIN_MAPPING {
    tag "$meta.id"
    label 'process_single'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //     'biocontainers/pandas:1.5.2' }"
    container 'docker.io/ferriolcalvet/bgreference'


    input:
    tuple val(meta) , path(mutations_file)
    tuple val(meta2), path(panel_file)

    output:
    tuple val(meta2), path("*.mapping.tsv") , emit: mapping
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    """
    cut -f 1,2,6,9 ${panel_file} | cut -f-4 | uniq > ${meta2.id}.panel.unique.tsv
    panels_compute_dna_2_protein.py \\
                --maf mutations_file \\
                --consensus-file ${meta2.id}.panel.unique.tsv \\
                --output ${prefix}.mapping.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
