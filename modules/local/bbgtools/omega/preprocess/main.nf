process OMEGA_PREPROCESS {
    tag "${meta.id}"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/omega:20250716'

    input:
    tuple val(meta), path(mutations), path(depths), path(mutation_profile)
    tuple val(meta2), path(annotated_panel)
    tuple val(meta3), path(syn_muts_global)
    tuple val(meta4), path(mut_profile_global, stageAs: 'global_mutprofile.tsv')

    output:
    tuple val(meta), path("mutability_per_sample_gene_context.*.tsv"), path("mutations_per_sample_gene_impact_context.*.tsv"), emit: mutabs_n_mutations_tsv
    tuple val(meta), path("syn_muts.*.tsv"), emit: syn_muts_tsv
    path "versions.yml", topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"

    // TODO revise this fix
    def sample_name = prefix.tokenize('.')[0]
    def global_loc = task.ext.global_loc ? "--absent-synonymous infer_global_custom  --mutational-profile-global-file global_mutprofile.tsv --synonymous-mutrates-file ${syn_muts_global}" : "--absent-synonymous ignore"
    prefix = task.ext.global_loc ? "${prefix}.gLoc" : "${prefix}"
    """
    omega preprocessing --preprocessing-mode compute_mutabilities \\
                        --depths-file ${depths} \\
                        --mutations-file ${mutations} \\
                        --input-vep-postprocessed-file ${annotated_panel} \\
                        --table-observed-muts mutations_per_sample_gene_impact_context.${prefix}.tsv \\
                        --mutabilities-table mutability_per_sample_gene_context.${prefix}.tsv \\
                        --synonymous-muts-table syn_muts.${prefix}.tsv \\
                        --mutational-profile-file ${mutation_profile} \\
                        --single-sample ${sample_name} \\
                        ${global_loc}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:

    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}
