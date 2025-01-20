process OMEGA_PREPROCESS {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/omega:20250113'

    input:
    tuple val(meta) , path(mutations), path(depths), path(mutation_profile)
    tuple val(meta2), path (annotated_panel)
    tuple val(meta3), path (syn_muts_global)
    tuple val(meta4), path (mut_profile_global, stageAs: 'global_mutprofile.tsv')


    output:
    tuple val(meta), path("mutability_per_sample_gene_context.*.tsv"), path("mutations_per_sample_gene_impact_context.*.tsv") , emit: mutabs_n_mutations_tsv
    tuple val(meta), path("syn_muts.*.tsv")     , emit: syn_muts_tsv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    # $args -c $task.cpus
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

// omega preprocessing --preprocessing-mode compute_mutabilities
//                      --depths-file all_samples.subset_depths.tsv.gz
//                      --mutations-file all_samples.mutations.tsv
//                      --input-vep-postprocessed-file consensus.exons_splice_sites.tsv
//                      --table-observed-muts mutations_per_sample_gene_impact_context.all_samples2.global_loc.gLoc.tsv
//                      --mutabilities-table mutability_per_sample_gene_context.all_samples2.global_loc.gLoc.tsv
//                      --synonymous-muts-table syn_muts.all_samples2.global_loc.gLoc.tsv
//                      --mutational-profile-file all_samples.all.profile.tsv
//                      --mutational-profile-global-file P19_0033_BTR_01.all.profile.tsv
//                      --single-sample all_samples
//                      --absent-synonymous infer_global_custom
//                      --synonymous-mutrates-file mutrates_per_gene.tsv
