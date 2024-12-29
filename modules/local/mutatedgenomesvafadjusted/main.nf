process MUTATED_GENOMES_FROM_VAF {
    tag "$meta.id"
    label 'process_single'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //     'biocontainers/pandas:1.5.2' }"

    container 'docker.io/ferranmuinos/test_mutated_genomes'

    input:
    tuple val(meta), path(mutations)

    output:
    tuple val(meta), path("*.sample.mutated_genomes_from_vaf.tsv") , emit: mutated_epi_sample
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    [params]
memory=12G

[pre]
 . "/home/$USER/miniconda3/etc/profile.d/conda.sh"
conda activate probabilistic

    python mutgenomes_driver_priority.py
            snv-am
            --sample P19_0001_BTR_01
            --outfolder /workspace/nobackup/bladder_ts/results/2024-11-25_deepCSA/coveredurothelium/covered_genomes_snv_AM

    mutated_genomes_from_vaf.py \\
                --sample ${prefix} \\
                --filename ${mutations} ;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    // def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.${panel_version}.mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
