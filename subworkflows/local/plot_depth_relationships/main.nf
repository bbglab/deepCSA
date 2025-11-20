
include { PLOT_VAF_MUTDENSITY_DEPTH    as PLOTVAFMUTDENSITY     } from '../../../modules/local/plot/vaf_mutdensity_depth/main'
include { PLOT_SELECTION_DEPTH          as PLOTSELECTIONDEPTH    } from '../../../modules/local/plot/selection_depth/main'


workflow PLOT_DEPTH_RELATIONSHIPS {

    take:
    somatic_mutations       // Channel: tuple(meta, maf)
    all_mutdensities        // Channel: tuple(meta, mutdensity) - per sample
    depth_per_gene          // Channel: tuple(meta, depth) - multiple depth files per sample
    omega_results           // Channel: tuple(meta, omega) - optional
    oncodrivefml_results    // Channel: tuple(meta, oncodrivefml) - optional


    main:

    // For VAF and mutation density plotting
    // We need to join MAF with mutation density and depth files
    // First, get one depth file per sample (pick first from depths output)
    depth_per_gene
    .map { meta, depths -> 
        // PLOT_DEPTHS outputs multiple files, we need depth_per_gene file
        def depth_file = depths instanceof List ? depths.find { it.name.contains('depth_per_gene') } : depths
        tuple(meta, depth_file)
    }
    .set { depth_single_ch }
    
    // Join MAF with mutation density and depth
    somatic_mutations
    .join( all_mutdensities, remainder: true )
    .join( depth_single_ch, remainder: true )
    .map { meta, maf, mutdens, depth ->
        // Handle missing files
        def maf_file = maf ?: file('NO_FILE')
        def mutdens_file = mutdens ?: file('NO_FILE')
        def depth_file = depth ?: file('NO_FILE')
        tuple(meta, maf_file, mutdens_file, depth_file)
    }
    .set { vaf_mutdens_input }
    
    PLOTVAFMUTDENSITY(vaf_mutdens_input)
    
    // For selection metrics (omega and OncodriveFML)
    // Join with depth and process
    omega_results
    .join( oncodrivefml_results, remainder: true )
    .join( depth_single_ch, remainder: true )
    .map { meta, omega, ofml, depth ->
        def omega_file = omega ?: file('NO_FILE')
        def ofml_file = ofml ?: file('NO_FILE')
        def depth_file = depth ?: file('NO_FILE')
        tuple(meta, omega_file, ofml_file, depth_file)
    }
    .filter { meta, omega, ofml, depth -> 
        !depth.name.equals('NO_FILE') && (!omega.name.equals('NO_FILE') || !ofml.name.equals('NO_FILE'))
    }
    .set { selection_input }
    
    PLOTSELECTIONDEPTH(selection_input)

    emit:
    vaf_mutdensity_plots    = PLOTVAFMUTDENSITY.out.plots
    selection_plots         = PLOTSELECTIONDEPTH.out.plots

}
