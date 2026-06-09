import sys

with open('subworkflows/local/omega/main.nf', 'r') as f:
    lines = f.readlines()

new_lines = []
for i, line in enumerate(lines):
    new_lines.append(line)
    if "include { OMEGA_MULTITEST           as OMEGAMULTIPLETESTGLOBALLOC       } from '../../../modules/local/omega_multipletesting/main'" in line:
        new_lines.append("include { HOTSPOTS_SELECTION        as HOTSPOTSSELECTION                } from '../../../modules/local/hotspots_selection/main'\n")
        new_lines.append("include { HOTSPOTS_SELECTION        as HOTSPOTSSELECTIONGLOBALLOC       } from '../../../modules/local/hotspots_selection/main'\n")
        
    if "site_comparison_results_flattened }" in line:
        new_lines.append("""
    if (params.hotspots_annotation && params.hotspots_definition_file) {
        hotspots_file = channel.fromPath(params.hotspots_definition_file, checkIfExists: true).first()
        
        HOTSPOTSSELECTION(
            site_comparison_results,
            QUERYPANEL.out.subset.first(),
            hotspots_file
        )
        // If needed, we can also collect or emit these results
    }
""")

with open('subworkflows/local/omega/main.nf', 'w') as f:
    f.writelines(new_lines)
