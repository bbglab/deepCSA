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
include { BLACKLIST_MUTATIONS           as BLACKLISTMUTS    } from '../../../modules/local/blacklistmuts/main'
include { PLOT_MUTATIONS                as PLOTMAF          } from '../../../modules/local/plot/mutations_summary/main'
include { PLOT_NEEDLES                  as PLOTNEEDLES      } from '../../../modules/local/plot/needles/main'


workflow MUTATION_PREPROCESSING {

    take:
    vcfs
    vep_cache
    vep_extra_files
    bedfile
    bedfile_exons
    groups
    sequence_information_df
    custom_annotation_tsv

    main:


    VCFANNOTATE(vcfs,
                    params.fasta,
                    params.vep_genome,
                    params.vep_species,
                    params.vep_cache_version,
                    vep_cache,
                    vep_extra_files)

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

    WRITEMAF(FILTERBATCH.out.cohort_maf, groups)

    // Here we flatten the output of the WRITEMAF module to have a channel where each item is a sample-maf pair
    WRITEMAF.out.mafs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_mafs }

    SOMATICMUTATIONS(named_mafs)

    if (params.blacklist_mutations) {
        blacklist_mutations  = Channel.fromPath( params.blacklist_mutations ?: params.input, checkIfExists: true).first()
        BLACKLISTMUTS(SOMATICMUTATIONS.out.mutations, blacklist_mutations)
        somatic_mutations = BLACKLISTMUTS.out.mutations
    } else {
        somatic_mutations = SOMATICMUTATIONS.out.mutations
    }




    PLOTNEEDLES(somatic_mutations, sequence_information_df)

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
    somatic_mafs            = somatic_mutations
    all_raw_vep_annotation  = SUMANNOTATION.out.tab_all
    bedfile_clean           = bedfile_updated

}
