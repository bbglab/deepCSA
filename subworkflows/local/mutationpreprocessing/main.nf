// Annotation
include { VCF_ANNOTATE_ENSEMBLVEP       as VCFANNOTATE      } from '../../nf-core/vcf_annotate_ensemblvep/main'


include { SUMMARIZE_ANNOTATION          as SUMANNOTATION    } from '../../../modules/local/process_annotation/mutations/main'
include { CUSTOM_MUTATION_PROCESSING    as CUSTOMANNOTATION } from '../../../modules/local/process_annotation/mutations_custom/main'
include { VCF2MAF                       as VCF2MAF          } from '../../../modules/local/vcf2maf/main'
include { FILTERBED                     as FILTERPANEL      } from '../../../modules/local/filterbed/main'
include { FILTERBED                     as FILTEREXONS      } from '../../../modules/local/filterbed/main'
include { MERGE_BATCH                   as MERGEBATCH       } from '../../../modules/local/mergemafs/main'
include { FILTER_BATCH                  as FILTERBATCH      } from '../../../modules/local/filtermaf/main'
include { WRITE_MAFS                    as WRITEMAF         } from '../../../modules/local/writemaf/main'
include { SUBSET_MAF                    as SOMATICMUTATIONS } from '../../../modules/local/subsetmaf/main'
include { SUBSET_MAF                    as CLEANMUTATIONS   } from '../../../modules/local/subsetmaf/main'
include { BLACKLIST_MUTATIONS           as BLACKLISTMUTS    } from '../../../modules/local/blacklistmuts/main'
include { PLOT_MUTATIONS                as PLOTMAF          } from '../../../modules/local/plot/mutations_summary/main'
include { PLOT_MUTATIONS                as PLOTSOMATICMAF   } from '../../../modules/local/plot/mutations_summary/main'
include { PLOT_NEEDLES                  as PLOTNEEDLES      } from '../../../modules/local/plot/needles/main'
include { DOWNSAMPLE_MUTATIONS          as DOWNSAMPLEMUTS   } from '../../../modules/local/downsample/mutations/main'
include { COMPUTE_CONTAMINATION         as CONTAMINATION    } from '../../../modules/local/contamination/main'


workflow MUTATION_PREPROCESSING {

    take:
    vcfs
    bedfile
    bedfile_exons
    all_groups
    groups
    sequence_information_df
    custom_annotation_tsv

    main:


    VCFANNOTATE(vcfs,
                    params.fasta,
                    params.vep_genome,
                    params.vep_species,
                    params.vep_cache_version,
                    params.vep_cache,
                    [])

    // Join all annotated samples and put them in a channel to be summarized together
    VCFANNOTATE.out.tab.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ annotated_samples }

    hotspots_definition_file = params.hotspots_annotation ? Channel.fromPath( params.hotspots_definition_file, checkIfExists: true).first() : Channel.fromPath(params.input).first()
    SUMANNOTATION(annotated_samples, hotspots_definition_file)

    if (params.customize_annotation) {
        // Update impact of mutations in specific regions based on user preferences
        CUSTOMANNOTATION(SUMANNOTATION.out.tab, custom_annotation_tsv)
        summary_of_mutations = CUSTOMANNOTATION.out.mutations.first()
    } else {
        summary_of_mutations = SUMANNOTATION.out.tab.first()
    }

    VCF2MAF(vcfs, summary_of_mutations)

    FILTEREXONS(VCF2MAF.out.maf, bedfile_exons)

    FILTERPANEL(FILTEREXONS.out.maf, bedfile)

    // Join all samples' MAFs and put them in a channel to be merged
    FILTERPANEL.out.maf.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ samples_maf }

    MERGEBATCH(samples_maf)

    FILTERBATCH(MERGEBATCH.out.cohort_maf)

    PLOTMAF(FILTERBATCH.out.cohort_maf)

    WRITEMAF(FILTERBATCH.out.cohort_maf, all_groups)

    // Here we flatten the output of the WRITEMAF module to have a channel where each item is a sample-maf pair
    WRITEMAF.out.mafs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_mafs }

    // Remove mutations that are blacklisted
    if (params.blacklist_mutations) {
        blacklist_mutations  = Channel.fromPath( params.blacklist_mutations ?: params.input, checkIfExists: true).first()
        BLACKLISTMUTS(named_mafs, blacklist_mutations)
        _all_clean_mutations = BLACKLISTMUTS.out.mutations
    } else {
        _all_clean_mutations = named_mafs
    }

    // if (params.downsample && params.downsample_proportion < 1) {
    if (params.downsample) {
        DOWNSAMPLEMUTS(_all_clean_mutations)
        all_clean_mutations = DOWNSAMPLEMUTS.out.downsampled_muts
    } else {
        all_clean_mutations = _all_clean_mutations
    }

    // Clean mutations based on artifact filtering decisions
    CLEANMUTATIONS(all_clean_mutations)

    // Keep only somatic mutations
    SOMATICMUTATIONS(CLEANMUTATIONS.out.mutations)

    

    Channel.of([["id": "all_samples"]])
    .join(named_mafs).first()
    .set{raw_muts_all_samples}

    Channel.of([["id": "all_samples"]])
    .join(SOMATICMUTATIONS.out.mutations).first()
    .set{muts_all_samples}

    CONTAMINATION(raw_muts_all_samples, muts_all_samples)

    PLOTSOMATICMAF(muts_all_samples)

    // Filter SOMATICMUTATIONS.out.mutations by meta.id in group_keys
    SOMATICMUTATIONS.out.mutations
    .map { mut -> tuple(mut[0].id, mut) }
    .join(groups)
    .map { it[1] }
    .set { muts_for_plotting }

    PLOTNEEDLES(muts_for_plotting, sequence_information_df)


    // Compile a BED file with all the mutations that are discarded due to:
    // Other sample SNP
    //     All sites with this filter should be remove from the background.
    // Repetitive variant
    //     Mutation seen in more than X samples.
    //     All sites with this filter should be remove from the background.
    // N rich
    //     - only sites with mutations discarded by N-rich criteria
    //     ?- also non mutated sites that do have lots of Ns

    // bedfiles.sample_discarded
    bedfile.set{ bedfile_updated }


    emit:
    mafs                    = named_mafs
    somatic_mafs            = SOMATICMUTATIONS.out.mutations
    clean_mafs              = CLEANMUTATIONS.out.mutations
    mutations_all_samples   = muts_all_samples
    all_raw_vep_annotation  = SUMANNOTATION.out.tab_all
    bedfile_clean           = bedfile_updated

}
