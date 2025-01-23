/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; paramsSummaryMap; samplesheetToList } from 'plugin/nf-schema'

// def summary_params = paramsSummaryMap(workflow)

// def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
// def citation = '\n' + WorkflowMain.citation(workflow) + '\n'

// // Print parameter summary log to screen
// log.info logo + paramsSummaryLog(workflow) + citation

// WorkflowDeepcsa.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Consisting of a mix of local and nf-core/modules.
*/

// SUBWORKFLOW
include { INPUT_CHECK                                       } from '../subworkflows/local/input_check'

include { DEPTH_ANALYSIS            as DEPTHANALYSIS        } from '../subworkflows/local/depthanalysis/main'
include { CREATE_PANELS             as CREATEPANELS         } from '../subworkflows/local/createpanels/main'

include { PLOT_DEPTHS               as PLOTDEPTHSALLCONS    } from '../subworkflows/local/plotdepths/main'
include { PLOT_DEPTHS               as PLOTDEPTHSEXONS      } from '../subworkflows/local/plotdepths/main'
include { PLOT_DEPTHS               as PLOTDEPTHSEXONSCONS  } from '../subworkflows/local/plotdepths/main'

include { MUTATION_PREPROCESSING    as MUT_PREPROCESSING    } from '../subworkflows/local/mutationpreprocessing/main'

include { MUTATION_RATE             as MUTRATEALL           } from '../subworkflows/local/mutationrate/main'
include { MUTATION_RATE             as MUTRATEPROT          } from '../subworkflows/local/mutationrate/main'
include { MUTATION_RATE             as MUTRATENONPROT       } from '../subworkflows/local/mutationrate/main'
include { MUTATION_RATE             as MUTRATESYNONYMOUS    } from '../subworkflows/local/mutationrate/main'

include { MUTATIONAL_PROFILE        as MUTPROFILEALL        } from '../subworkflows/local/mutationprofile/main'
include { MUTATIONAL_PROFILE        as MUTPROFILENONPROT    } from '../subworkflows/local/mutationprofile/main'
include { MUTATIONAL_PROFILE        as MUTPROFILEEXONS      } from '../subworkflows/local/mutationprofile/main'
include { MUTATIONAL_PROFILE        as MUTPROFILEINTRONS    } from '../subworkflows/local/mutationprofile/main'

include { MUTABILITY                as MUTABILITYALL        } from '../subworkflows/local/mutability/main'
include { MUTABILITY                as MUTABILITYNONPROT    } from '../subworkflows/local/mutability/main'

include { ONCODRIVEFML_ANALYSIS     as ONCODRIVEFMLALL      } from '../subworkflows/local/oncodrivefml/main'
include { ONCODRIVEFML_ANALYSIS     as ONCODRIVEFMLNONPROT  } from '../subworkflows/local/oncodrivefml/main'
include { ONCODRIVE3D_ANALYSIS      as ONCODRIVE3D          } from '../subworkflows/local/oncodrive3d/main'
include { ONCODRIVECLUSTL_ANALYSIS  as ONCODRIVECLUSTL      } from '../subworkflows/local/oncodriveclustl/main'

include { OMEGA_ANALYSIS            as OMEGA                } from '../subworkflows/local/omega/main'
include { OMEGA_ANALYSIS            as OMEGANONPROT         } from '../subworkflows/local/omega/main'
include { OMEGA_ANALYSIS            as OMEGAMULTI           } from '../subworkflows/local/omega/main'
include { OMEGA_ANALYSIS            as OMEGANONPROTMULTI    } from '../subworkflows/local/omega/main'

include { INDELS_SELECTION          as INDELSSELECTION      } from '../subworkflows/local/indels/main'

include { MUTATED_EPITHELIUM        as MUTATEDEPITHELIUM    } from '../subworkflows/local/mutatedepithelium/reads/main'
include { MUTATED_EPITHELIUM_VAF    as MUTATEDEPITHELIUMVAF } from '../subworkflows/local/mutatedepithelium/vaf/main'

include { EXPECTED_MUTRATE          as EXPECTEDMUTRATE      } from '../subworkflows/local/mutatedepithelium/expected/main'

include { SIGNATURES                as SIGNATURESALL        } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESNONPROT    } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESEXONS      } from '../subworkflows/local/signatures/main'
include { SIGNATURES                as SIGNATURESINTRONS    } from '../subworkflows/local/signatures/main'

include { PLOT_SELECTION_METRICS    as PLOTSELECTION        } from '../modules/local/plot/selection_metrics/main'

include { DNDS                      as DNDS                 } from '../subworkflows/local/dnds/main'

include { TABIX_BGZIPTABIX_QUERY    as DEPTHSALLCONS        } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSEXONSCONS      } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSPROTCONS       } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSNONPROTCONS    } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSINTRONSCONS    } from '../modules/nf-core/tabix/bgziptabixquery/main'
include { TABIX_BGZIPTABIX_QUERY    as DEPTHSSYNONYMOUSCONS } from '../modules/nf-core/tabix/bgziptabixquery/main'

include { SELECT_MUTRATES           as SYNMUTRATE           } from '../modules/local/select_mutrate/main'
include { SELECT_MUTRATES           as SYNMUTREADSRATE      } from '../modules/local/select_mutrate/main'
include { DNA_2_PROTEIN_MAPPING     as DNA2PROTEINMAPPING   } from '../modules/local/dna2protein/main'


// Download annotation cache if needed
include { PREPARE_CACHE                                     } from '../subworkflows/local/prepare_cache/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Installed directly from nf-core/modules.
*/

// MODULE
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { ANNOTATE_DEPTHS           as ANNOTATEDEPTHS       } from '../modules/local/annotatedepth/main'
include { TABLE_2_GROUP             as TABLE2GROUP          } from '../modules/local/table2groups/main'
include { MUTATIONS_2_SIGNATURES    as MUTS2SIGS            } from '../modules/local/mutations2sbs/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEEPCSA{

    // // Input channel definitions
    features_table  = Channel.fromPath( params.features_table ?: params.input, checkIfExists: true)

    wgs_trinucs     = params.wgs_trinuc_counts 
                            ? Channel.fromPath( params.wgs_trinuc_counts, checkIfExists: true).first() 
                            : Channel.empty()
    cosmic_ref      = params.cosmic_ref_signatures 
                            ? Channel.fromPath( params.cosmic_ref_signatures, checkIfExists: true).first() 
                            : Channel.empty()
    datasets3d      = params.datasets3d 
                            ? Channel.fromPath( params.datasets3d, checkIfExists: true).first() 
                            : Channel.empty()
    annotations3d   = params.annotations3d 
                            ? Channel.fromPath( params.annotations3d, checkIfExists: true).first() 
                            : Channel.empty()
    seqinfo_df      = params.datasets3d 
                            ? Channel.fromPath( "${params.datasets3d}/seq_for_mut_prob.tsv", checkIfExists: true).first() 
                            : Channel.empty()
    cadd_scores     = params.cadd_scores 
                            ? Channel.of([ 
                                file(params.cadd_scores, checkIfExists : true), 
                                file(params.cadd_scores_ind, checkIfExists : true) 
                                ]).first() 
                            : Channel.empty()


    // if the user wants to use custom gene groups, import the gene groups table
    // otherwise I am using the input csv as a dummy value channel
    custom_groups_table = params.custom_groups_file 
                                ? Channel.fromPath( params.custom_groups_file, checkIfExists: true).first() 
                                : Channel.fromPath(params.input)

    // if the user wants to use custom BED file for computing the depths, import the BED file
    // otherwise I am using the input csv as a dummy value channel
    custom_bed_file     = params.custom_bedfile 
                                ? Channel.fromPath( params.custom_bedfile, checkIfExists: true).first() 
                                : Channel.fromPath(params.input)

    // if the user wants to define hotspots for omega, import the hotspots definition BED file
    // otherwise I am using the input csv as a dummy value channel
    hotspots_bed_file   = params.omega_hotspots_bedfile 
                                ? Channel.fromPath( params.omega_hotspots_bedfile, checkIfExists: true).first() 
                                : Channel.fromPath(params.input)

    // Initialize booleans based on user params
    def run_mutabilities    = (params.oncodrivefml || params.oncodriveclustl || params.oncodrive3d)
    def run_mutrate         = (params.mutationrate || params.omega)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )


    //
    // Separate input BAMs and VCFs
    //
    INPUT_CHECK.out.mutations
    .map{ it -> [it[0], it[1]]}
    .set{ meta_vcfs_alone }

    INPUT_CHECK.out.mutations
    .map{ it -> [it[0], it[2]]}
    .set{ meta_bams_alone }

    INPUT_CHECK.out.mutations
    .map{ it -> [it[0], it[3], it[4]]}
    .set{ meta_pileupbamindex_alone }


    // TODO: test if downloading VEP cache works
    // Download Ensembl VEP cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache ? [] : Channel.of([ [ id:"${params.vep_genome}.${params.vep_cache_version}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    if (params.download_cache) {
        PREPARE_CACHE(ensemblvep_info)
        vep_cache = PREPARE_CACHE.out.ensemblvep_cache.map{ _meta, cache -> [ cache ] }
    } else {
        vep_cache = params.vep_cache
    }
    vep_extra_files = []


    TABLE2GROUP(features_table)


    // Depth analysis: compute and plots
    DEPTHANALYSIS(meta_bams_alone, custom_bed_file)

    // Panels generation: all modalities
    CREATEPANELS(DEPTHANALYSIS.out.depths, vep_cache, vep_extra_files)

    ANNOTATEDEPTHS(DEPTHANALYSIS.out.depths, CREATEPANELS.out.all_panel, TABLE2GROUP.out.json_allgroups, file(params.input))
    ANNOTATEDEPTHS.out.annotated_depths.flatten().map{ it -> [ [id : it.name.tokenize('.')[0]] , it]  }.set{ annotated_depths }


    if (params.plot_depths){
        PLOTDEPTHSALLCONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.all_consensus_bed, CREATEPANELS.out.all_consensus_panel)
        PLOTDEPTHSEXONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.exons_bed, CREATEPANELS.out.exons_panel)
        PLOTDEPTHSEXONSCONS(ANNOTATEDEPTHS.out.all_samples_depths, CREATEPANELS.out.exons_consensus_bed, CREATEPANELS.out.exons_consensus_panel)
    }

    // Mutation preprocessing
    MUT_PREPROCESSING(meta_vcfs_alone, vep_cache, vep_extra_files,
                        CREATEPANELS.out.all_consensus_bed,
                        CREATEPANELS.out.exons_bed,
                        TABLE2GROUP.out.json_allgroups, seqinfo_df)
    positive_selection_results = MUT_PREPROCESSING.out.somatic_mafs

    Channel.of([["id": "all_samples"]])
    .join(MUT_PREPROCESSING.out.somatic_mafs).first()
    .set{mutations_all}

    if (params.vep_species == 'homo_sapiens'){
        DNA2PROTEINMAPPING(mutations_all, CREATEPANELS.out.exons_consensus_panel)
    }


    // Intersect BED of desired sites with samples' depths
    DEPTHSALLCONS(annotated_depths, CREATEPANELS.out.all_consensus_bed)
    DEPTHSEXONSCONS(annotated_depths, CREATEPANELS.out.exons_consensus_bed)
    DEPTHSPROTCONS(annotated_depths, CREATEPANELS.out.prot_consensus_bed)
    DEPTHSNONPROTCONS(annotated_depths, CREATEPANELS.out.nonprot_consensus_bed)
    DEPTHSSYNONYMOUSCONS(annotated_depths, CREATEPANELS.out.synonymous_consensus_bed)


    if (run_mutrate){
        // Mutation Rate
        MUTRATEALL(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSALLCONS.out.subset, CREATEPANELS.out.all_consensus_bed, CREATEPANELS.out.all_consensus_panel)
        MUTRATEPROT(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSPROTCONS.out.subset, CREATEPANELS.out.prot_consensus_bed, CREATEPANELS.out.prot_consensus_panel)
        MUTRATENONPROT(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSNONPROTCONS.out.subset, CREATEPANELS.out.nonprot_consensus_bed, CREATEPANELS.out.nonprot_consensus_panel)
        MUTRATESYNONYMOUS(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSSYNONYMOUSCONS.out.subset, CREATEPANELS.out.synonymous_consensus_bed, CREATEPANELS.out.synonymous_consensus_panel)

        Channel.of([ [ id: "all_samples" ] ])
        .join( MUTRATESYNONYMOUS.out.mutrates )
        .set{ all_samples_syn_mutrate }

        SYNMUTRATE(all_samples_syn_mutrate)

        SYNMUTREADSRATE(all_samples_syn_mutrate)


        // Concatenate all outputs into a single file
        mutrate_empty = Channel.empty()
        mutrate_empty
        .concat(MUTRATEALL.out.mutrates.map{ it -> it[1]}.flatten())
        .concat(MUTRATEPROT.out.mutrates.map{ it -> it[1]}.flatten())
        .concat(MUTRATENONPROT.out.mutrates.map{ it -> it[1]}.flatten())
        .concat(MUTRATESYNONYMOUS.out.mutrates.map{ it -> it[1]}.flatten())
        .set{ all_mutrates }
        all_mutrates.collectFile(name: "all_mutrates.tsv", storeDir:"${params.outdir}/mutrate", skip: 1, keepHeader: true)

    }


    // Mutational profile
    if (params.profileall){
        MUTPROFILEALL(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSALLCONS.out.subset, CREATEPANELS.out.all_consensus_bed, wgs_trinucs)
    }
    if (params.profilenonprot){
        MUTPROFILENONPROT(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSNONPROTCONS.out.subset, CREATEPANELS.out.nonprot_consensus_bed, wgs_trinucs)
    }
    if (params.profileexons){
        MUTPROFILEEXONS(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSEXONSCONS.out.subset, CREATEPANELS.out.exons_consensus_bed, wgs_trinucs)
    }
    if (params.profileintrons){
        DEPTHSINTRONSCONS(annotated_depths, CREATEPANELS.out.introns_consensus_bed)
        MUTPROFILEINTRONS(MUT_PREPROCESSING.out.somatic_mafs, DEPTHSINTRONSCONS.out.subset, CREATEPANELS.out.introns_consensus_bed, wgs_trinucs)
    }


    if (run_mutabilities) {
        if (params.profileall){
            MUTABILITYALL(MUT_PREPROCESSING.out.somatic_mafs,
                            annotated_depths,
                            MUTPROFILEALL.out.profile,
                            CREATEPANELS.out.exons_consensus_panel,
                            CREATEPANELS.out.exons_consensus_bed
                            )
        }
        if (params.profilenonprot){
            MUTABILITYNONPROT(MUT_PREPROCESSING.out.somatic_mafs,
                                annotated_depths,
                                MUTPROFILENONPROT.out.profile,
                                CREATEPANELS.out.exons_consensus_panel,
                                CREATEPANELS.out.exons_consensus_bed
                                )
        }
    }


    if (params.indels){
        INDELSSELECTION(MUT_PREPROCESSING.out.somatic_mafs,
                        CREATEPANELS.out.all_consensus_bed
                        )
        positive_selection_results = positive_selection_results.join(INDELSSELECTION.out.indels, remainder: true)
    }

    if (params.expected_mutation_rate){

        EXPECTEDMUTRATE(mutations_all,
                        CREATEPANELS.out.exons_consensus_bed,
                        CREATEPANELS.out.exons_consensus_panel,
                        ANNOTATEDEPTHS.out.all_samples_depths,
                        CREATEPANELS.out.full_panel_annotated
                        )
    }



    if (params.mutated_epithelium_vaf){
        MUT_PREPROCESSING.out.somatic_mafs
        .join(meta_vcfs_alone.map{it -> [ ["id" : it[0].id] , it[1] ] })
        .map{it -> [ it[0] , it[1] ]  }
        .set{ sample_mutations_only }

        MUTATEDEPITHELIUMVAF(sample_mutations_only,
                                CREATEPANELS.out.exons_consensus_bed)
    }

    if (params.mutated_epithelium){
        MUT_PREPROCESSING.out.somatic_mafs
        .join(meta_vcfs_alone.map{it -> [ ["id" : it[0].id] , it[1] ] })
        .map{it -> [ it[0] , it[1] ]  }
        .set{ sample_mutations_only }

        MUTATEDEPITHELIUM(sample_mutations_only,
                            CREATEPANELS.out.exons_consensus_bed,
                            CREATEPANELS.out.exons_consensus_panel,
                            meta_pileupbamindex_alone,
                            params.fasta
                            )

        if (params.pileup_all_duplex) {
            // Concatenate all outputs into a single file
            mut_epithelium_empty = Channel.empty()
            mut_epithelium_empty
            .concat(MUTATEDEPITHELIUM.out.mut_epi_exon.map{ it -> it[1]}.flatten())
            .set{ all_mutepiexon }
            all_mutepiexon.collectFile(name: "all_mutepithelium_exon.tsv", storeDir:"${params.outdir}/mutatedgenomesfromreadsam", skip: 1, keepHeader: true)

            mut_epithelium_empty2 = Channel.empty()
            mut_epithelium_empty2
            .concat(MUTATEDEPITHELIUM.out.mut_epi_gene.map{ it -> it[1]}.flatten())
            .set{ all_mutepigene }
            all_mutepigene.collectFile(name: "all_mutepithelium_gene.tsv", storeDir:"${params.outdir}/mutatedgenomesfromreadsam", skip: 1, keepHeader: true)

            mut_epithelium_empty3 = Channel.empty()
            mut_epithelium_empty3
            .concat(MUTATEDEPITHELIUM.out.mut_epi_sample.map{ it -> it[1]}.flatten())
            .set{ all_mutepisample }
            all_mutepisample.collectFile(name: "all_mutepithelium_sample.tsv", storeDir:"${params.outdir}/mutatedgenomesfromreadsam", skip: 1, keepHeader: true)
        } else {
            // Concatenate all outputs into a single file
            mut_epithelium_empty = Channel.empty()
            mut_epithelium_empty
            .concat(MUTATEDEPITHELIUM.out.mut_epi_exon.map{ it -> it[1]}.flatten())
            .set{ all_mutepiexon }
            all_mutepiexon.collectFile(name: "all_mutepithelium_exon.tsv", storeDir:"${params.outdir}/mutatedgenomesfromreads", skip: 1, keepHeader: true)

            mut_epithelium_empty2 = Channel.empty()
            mut_epithelium_empty2
            .concat(MUTATEDEPITHELIUM.out.mut_epi_gene.map{ it -> it[1]}.flatten())
            .set{ all_mutepigene }
            all_mutepigene.collectFile(name: "all_mutepithelium_gene.tsv", storeDir:"${params.outdir}/mutatedgenomesfromreads", skip: 1, keepHeader: true)

            mut_epithelium_empty3 = Channel.empty()
            mut_epithelium_empty3
            .concat(MUTATEDEPITHELIUM.out.mut_epi_sample.map{ it -> it[1]}.flatten())
            .set{ all_mutepisample }
            all_mutepisample.collectFile(name: "all_mutepithelium_sample.tsv", storeDir:"${params.outdir}/mutatedgenomesfromreads", skip: 1, keepHeader: true)

        }
    }



    //
    // Positive selection
    //

    // OncodriveFML
    if (params.oncodrivefml){
        if (params.profileall){
            ONCODRIVEFMLALL(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYALL.out.mutability,
                                CREATEPANELS.out.exons_consensus_panel,
                                cadd_scores
                            )
            // positive_selection_results = positive_selection_results.join(ONCODRIVEFMLALL.out.results, remainder: true)
            positive_selection_results = positive_selection_results.join(ONCODRIVEFMLALL.out.results_snvs, remainder: true)
        }
        if (params.profilenonprot){
            ONCODRIVEFMLNONPROT(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYNONPROT.out.mutability,
                                    CREATEPANELS.out.exons_consensus_panel,
                                    cadd_scores
                                )
        }
    }

    if (params.oncodrive3d){
        if (params.profileall){
            // Oncodrive3D
            ONCODRIVE3D(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_bed,
                        datasets3d, annotations3d, MUT_PREPROCESSING.out.all_raw_vep_annotation)
        }
    }

    // if (params.expected_mutation_rate & params.dnds){
    if (params.dnds){
        covariates = params.dnds_covariates ? Channel.fromPath( params.dnds_covariates, checkIfExists: true).first() : Channel.empty()
        ref_transcripts = params.dnds_ref_transcripts ? Channel.fromPath( params.dnds_ref_transcripts, checkIfExists: true).first() : Channel.empty()
        DNDS(MUT_PREPROCESSING.out.somatic_mafs,
                    DEPTHSEXONSCONS.out.subset,
                    CREATEPANELS.out.exons_consensus_bed,
                    CREATEPANELS.out.exons_consensus_panel,
                    covariates,
                    ref_transcripts
                    )
    }

    if (params.omega){
        // Omega
        if (params.profileall){
            OMEGA(MUT_PREPROCESSING.out.somatic_mafs,
                    DEPTHSEXONSCONS.out.subset,
                    MUTPROFILEALL.out.profile,
                    CREATEPANELS.out.exons_consensus_bed,
                    CREATEPANELS.out.exons_consensus_panel,
                    custom_groups_table,
                    hotspots_bed_file,
                    SYNMUTRATE.out.mutrate,
                    CREATEPANELS.out.panel_annotated_rich
                    )
            positive_selection_results = positive_selection_results.join(OMEGA.out.results, remainder: true)
            positive_selection_results = positive_selection_results.join(OMEGA.out.results_global, remainder: true)

            // Omega multi
            OMEGAMULTI(MUT_PREPROCESSING.out.somatic_mafs,
                        DEPTHSEXONSCONS.out.subset,
                        MUTPROFILEALL.out.profile,
                        CREATEPANELS.out.exons_consensus_bed,
                        CREATEPANELS.out.exons_consensus_panel,
                        custom_groups_table,
                        hotspots_bed_file,
                        SYNMUTREADSRATE.out.mutrate,
                        CREATEPANELS.out.panel_annotated_rich
                        )
            positive_selection_results = positive_selection_results.join(OMEGAMULTI.out.results, remainder: true)
            positive_selection_results = positive_selection_results.join(OMEGAMULTI.out.results_global, remainder: true)
        }
        if (params.profilenonprot){
            OMEGANONPROT(MUT_PREPROCESSING.out.somatic_mafs,
                            DEPTHSEXONSCONS.out.subset,
                            MUTPROFILENONPROT.out.profile,
                            CREATEPANELS.out.exons_consensus_bed,
                            CREATEPANELS.out.exons_consensus_panel,
                            custom_groups_table,
                            hotspots_bed_file,
                            SYNMUTRATE.out.mutrate,
                            CREATEPANELS.out.panel_annotated_rich
                            )

            OMEGANONPROTMULTI(MUT_PREPROCESSING.out.somatic_mafs,
                                DEPTHSEXONSCONS.out.subset,
                                MUTPROFILENONPROT.out.profile,
                                CREATEPANELS.out.exons_consensus_bed,
                                CREATEPANELS.out.exons_consensus_panel,
                                custom_groups_table,
                                hotspots_bed_file,
                                SYNMUTREADSRATE.out.mutrate,
                                CREATEPANELS.out.panel_annotated_rich
                                )
        }

    }

    if (params.oncodriveclustl){
        // OncodriveClustl
        if (params.profileall){
            ONCODRIVECLUSTL(MUT_PREPROCESSING.out.somatic_mafs, MUTABILITYALL.out.mutability, CREATEPANELS.out.exons_consensus_panel)
        }
    }


    if (params.signatures){

        // Signature Analysis
        if (params.profileall){
            SIGNATURESALL(MUTPROFILEALL.out.matrix_sigprof, MUTPROFILEALL.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)

            MUT_PREPROCESSING.out.somatic_mafs
            .join(SIGNATURESALL.out.mutation_probs)
            .set{mutations_n_sbs}
            MUTS2SIGS(mutations_n_sbs)


        }
        if (params.profilenonprot){
            SIGNATURESNONPROT(MUTPROFILENONPROT.out.matrix_sigprof, MUTPROFILENONPROT.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)
        }
        if (params.profileexons){
            SIGNATURESEXONS(MUTPROFILEEXONS.out.matrix_sigprof, MUTPROFILEEXONS.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)
        }
        if (params.profileintrons){
            SIGNATURESINTRONS(MUTPROFILEINTRONS.out.matrix_sigprof, MUTPROFILEINTRONS.out.wgs_sigprofiler, params.cosmic_ref_signatures, TABLE2GROUP.out.json_samples)
        }
    }

    if ( params.indels & params.oncodrivefml & params.omega ){
        positive_selection_results_ready = positive_selection_results.map { element -> [element[0], element[1..-1]] }
        PLOTSELECTION(positive_selection_results_ready, seqinfo_df)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        Channel.topic('versions').unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // Info required for completion email and summary
    def multiqc_report = []
    
    def summary_params = paramsSummaryMap(workflow)
    workflow_summary    = WorkflowDeepcsa.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDeepcsa.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.dump_parameters(workflow, params)
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
