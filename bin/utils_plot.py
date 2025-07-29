import matplotlib as mpl
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt



snv_color = {
    'C>A' : '#03bcee',
    'C>G' : '#010101',
    'C>T' : '#e32926',
    'T>A' : '#cac9c9',
    'T>C' : '#a1ce63',
    'T>G' : '#ebc6c4'
}

def plot_profile(frequencies, title='Mutational profile',
                        yticks=[0, 0.1],
                        ymax=0.3,
                        yaxis_name = "% SBS",
                        output_f=None):

    # Plot params
    fig, ax = plt.subplots(
        2,
        sharex='col',
        figsize=(6, 1.5),
        gridspec_kw={'height_ratios': [0.05, 0.95]}
    )
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.01,
                        hspace=0)

    plt.title(title, fontsize=7, y = 1.2)

    #### Top
    axis_key = 0
    ax[axis_key].spines['top'].set_visible(False)
    ax[axis_key].spines['right'].set_visible(False)
    ax[axis_key].spines['left'].set_visible(False)
    ax[axis_key].spines['bottom'].set_visible(False)

    ax[axis_key].set_yticks([])

    fontsize_subs_title = 5
    subs = 'C>A'
    rect1 = Rectangle((0, 0.15), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(15.75/2, 0.75, subs, fontsize=fontsize_subs_title, weight='normal', ha='center')
    subs = 'C>G'
    rect1 = Rectangle((16,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(16 + 15.75/2 -1, 0.75, subs, fontsize=fontsize_subs_title, weight='normal', ha='center')
    subs = 'C>T'
    rect1 = Rectangle((32,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(32 + 15.75/2 -1, 0.75, subs, fontsize=fontsize_subs_title, weight='normal', ha='center')
    subs = 'T>A'
    rect1 = Rectangle((48,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(48 + 15.75/2 -1, 0.75, subs, fontsize=fontsize_subs_title, weight='normal', ha='center')
    subs = 'T>C'
    rect1 = Rectangle((64,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(64 + 15.75/2 -1, 0.75, subs, fontsize=fontsize_subs_title, weight='normal', ha='center')
    subs = 'T>G'
    rect1 = Rectangle((80,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(80 + 15.75/2 -1, 0.75, subs, fontsize=fontsize_subs_title, weight='normal', ha='center')


    #### PLOT 1
    axis_key = 1

    probabilities = []
    colors = []
    labels = []
    for snv_type in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']:
        ref, alt = snv_type.split('>')
        for nuc5 in ['A', 'C', 'G', 'T']:
            for nuc3 in ['A', 'C', 'G', 'T']:
                trinucleotide = nuc5 + ref + nuc3
                trinucleotide_change = f'{trinucleotide}>{alt}'
                probabilities.append(frequencies.get(trinucleotide_change, 0))
                colors.append(snv_color[snv_type])
                labels.append(trinucleotide)

    ax[axis_key].bar(list(range(0, 96)), probabilities, width=0.8, bottom=None, align='center', color=colors)


    ax[axis_key].set_xlim(xmin=-1, xmax=96)
    # ax[axis_key].set_ylim(ymin=0, ymax=ymax)

    ax[axis_key].set_xticks([])
    ax[axis_key].set_xticks(list(range(0, 96)))
    ax[axis_key].set_xticklabels(labels, fontsize=3.5, color='black', va="top", ha="center")
    plt.xticks(rotation=90)

    ax[axis_key].set_ylabel(yaxis_name, fontsize = 6)
    ax[axis_key].set_yticklabels(ax[axis_key].get_yticklabels(), fontsize=3.5, color='black', va="center", ha="center")


    ########## Y axis
    ylabels_dict = {
        0: '',
        1: yaxis_name
    }
    for key, value in ylabels_dict.items():
        # ax[key].set_ylabel(value, fontsize=4, rotation=90, labelpad=2)
        ax[key].set_axisbelow(True)
        for tick in ax[key].yaxis.get_major_ticks():
            tick.label1.set_fontsize(3.5)

        for location in ['top', 'bottom', 'left', 'right']:
            ax[key].spines[location].set_linewidth(0.25)

        ax[key].tick_params(axis = "x", which = "both", bottom=False)
        ax[key].tick_params(axis = "y", which = "both", left=False)
        ax[key].tick_params(axis='y', which='major', pad=2)
        ax[key].tick_params(axis='x', which='major', pad=0)

    if output_f:
        plt.savefig(output_f, bbox_inches='tight', dpi=600)
#     plt.close()


metrics_colors_dictionary = {"ofml"        : "viridis_r",
                                "ofml_score"  : "#6A33E0",
                                "omega_trunc" : "#FA5E32",
                                "omega_synon" : "#89E4A2",
                                "omega_miss"  : "#FABE4A",
                                "o3d_score"   : "#6DBDCC",
                                "o3d_cluster" : "skyblue",
                                "o3d_prob"    : "darkgray",
                                "frameshift"  : "#E4ACF4",
                                "inframe"     : "C5",
                                "hv_lines"    : "lightgray", # horizontal and vertical lines,
                                "hv_lines_needle" : "gray",
                                "needle_obs"  : "#003366",
                                "omega_miss_tert" : "#f5840c",
                                "omega_synon_tert": "#378c12",
                                "nonsense" : "#FA5E32",
                                "synonymous" : "#89E4A2",
                                "missense"  : "#FABE4A",
                                "indel"       : "#ECC4F7",
                                "splicing"    : "#A1C5DF",
                                }

mutation_type_colors = {
    "C>A": "#5abdeb",
    "C>G": "#050708",
    "C>T": "#d43c32",
    "T>A": "#cbcacb",
    "T>C": "#aacb72",
    "T>G": "#e7c9c6"
}

plots_general_config = {

                        # fonsizes
                        "ylabel_fontsize": 6,
                        "xlabel_fontsize": 6,
                        "xylabel_fontsize": 6,
                        "title_fontsize": 7,
                        "xyticks_fontsize": 5,
                        "xticks_fontsize": 5,
                        "yticks_fontsize": 5,
                        "legend_fontsize": 5,
                        "annots_fontsize": 5,

                        "dot_size_scplot": 15,
                        "dot_size_coeffplot": 5,
                        "dot_sizebelow_coeffplot": 40,
                        "dot_color_coeffplot": "#D3D3D3",
                        "dot_colorabove_coeffplot": "#D62728",
                        "dot_colorbelow_coeffplot": "#f29c9e",
                        "dot_edgethres_coeffplot": 0.2,
                        "dot_edgewidth_coeffplot": 0.5
                        }

mpl.rcParams.update({
#    'font.family': 'Arial',            # Enforce Arial
   'pdf.fonttype': 42,                # TrueType for PDF
   'ps.fonttype': 42,                 # TrueType for PS/EPS
   'svg.fonttype': 'none',            # Keep text as editable text
})


# mpl.rcParams.update({
#     'axes.titlesize'    : plots_general_config["title_fontsize"],       # Title font size
#     'axes.labelsize'    : plots_general_config["xylabel_fontsize"],     # X and Y axis labels
#     'xtick.labelsize'   : plots_general_config["xyticks_fontsize"],     # X tick labels
#     'ytick.labelsize'   : plots_general_config["xyticks_fontsize"],     # Y tick labels
#     'legend.fontsize'   : plots_general_config["legend_fontsize"],      # Legend text
#     'figure.titlesize'  : plots_general_config["title_fontsize"],       # Figure suptitle (if used)
# })
