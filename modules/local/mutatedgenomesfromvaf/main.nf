process MUTATED_GENOMES_FROM_VAF {
    tag "$meta.id"
    label 'process_single'
    label 'memory_medium'

    container 'docker.io/ferranmuinos/test_mutated_genomes'
    // conda activate probabilistic

    input:
    tuple val(meta) , path(mutations), path(omegas)
    // fn = os.path.join(deepcsa_run_dir, 'omegagloballoc', f'output_mle.{sample}.multi.global_loc.tsv')

    output:
    tuple val(meta), path("*.sample.mutated_genomes_from_vaf.tsv") , emit: mutated_epi_sample
    path  "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.covered_genomes_snv_AM

    mutgenomes_driver_priority.py \\
                snv-am \\
                --sample ${prefix} \\
                --outfolder ${prefix}.covered_genomes_snv_AM \\
                --somatic-mutations-file ${mutations} \\
                --omega-file ${omegas} ;

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
