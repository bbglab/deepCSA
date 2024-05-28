// Annotation
include { VCF_ANNOTATE_ENSEMBLVEP   as VCFANNOTATE      } from '../../nf-core/vcf_annotate_ensemblvep/main'


include { SUMMARIZE_ANNOTATION      as SUMANNOTATION    } from '../../../modules/local/summarize_annotation/main'
include { VCF2MAF                   as VCF2MAF          } from '../../../modules/local/vcf2maf/main'
include { FILTERBED                 as FILTERPANEL      } from '../../../modules/local/filterbed/main'
include { MERGE_BATCH               as MERGEBATCH       } from '../../../modules/local/mergemafs/main'
include { FILTER_BATCH              as FILTERBATCH      } from '../../../modules/local/filtermaf/main'
include { WRITE_MAFS                as WRITEMAF         } from '../../../modules/local/writemaf/main'
include { SUBSET_MAF                as SOMATICMUTATIONS } from '../../../modules/local/subsetmaf/main'
include { PLOT_MUTATIONS            as PLOTMAF          } from '../../../modules/local/plot/mutations_summary/main'


workflow MUTATION_PREPROCESSING {

    take:
    vcfs
    vep_cache
    vep_extra_files
    bedfile
    groups

    main:

    ch_versions = Channel.empty()

    VCFANNOTATE(vcfs,
                    params.fasta,
                    params.vep_genome,
                    params.vep_species,
                    params.vep_cache_version,
                    vep_cache,
                    vep_extra_files)
    ch_versions = ch_versions.mix(VCFANNOTATE.out.versions.first())

    // Join all annotated samples and put them in a channel to be summarized together
    VCFANNOTATE.out.tab.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ annotated_samples }


    SUMANNOTATION(annotated_samples)
    ch_versions = ch_versions.mix(SUMANNOTATION.out.versions)


    VCF2MAF(vcfs, SUMANNOTATION.out.tab)
    ch_versions = ch_versions.mix(VCF2MAF.out.versions.first())


    FILTERPANEL(VCF2MAF.out.maf, bedfile)
    ch_versions = ch_versions.mix(FILTERPANEL.out.versions.first())

    // Join all samples' MAFs and put them in a channel to be merged
    FILTERPANEL.out.maf.map{ it -> it[1] }.collect().map{ it -> [[ id:"all_samples" ], it]}.set{ samples_maf }


    MERGEBATCH(samples_maf)
    ch_versions = ch_versions.mix(MERGEBATCH.out.versions)


    FILTERBATCH(MERGEBATCH.out.cohort_maf)
    ch_versions = ch_versions.mix(FILTERBATCH.out.versions)

    PLOTMAF(FILTERBATCH.out.cohort_maf)
    ch_versions = ch_versions.mix(PLOTMAF.out.versions)

    WRITEMAF(FILTERBATCH.out.cohort_maf, groups)
    ch_versions = ch_versions.mix(WRITEMAF.out.versions)

    // Here we flatten the output of the WRITEMAF module to have a channel where each item is a sample-maf pair
    WRITEMAF.out.mafs.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ named_mafs }


    SOMATICMUTATIONS(named_mafs)
    ch_versions = ch_versions.mix(SOMATICMUTATIONS.out.versions)

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
    all_raw_vep_annotation  = SUMANNOTATION.out.tab_all
    bedfile_clean           = bedfile_updated
    versions                = ch_versions

}
