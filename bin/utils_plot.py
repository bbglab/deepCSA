import matplotlib
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rcParams



snv_color = {
    'C>A' : '#03bcee',
    'C>G' : '#010101',
    'C>T' : '#e32926',
    'T>A' : '#cac9c9',
    'T>C' : '#a1ce63',
    'T>G' : '#ebc6c4'
}

def plot_profile(frequencies, title='title', ylabels=[0, 0.1], ymax=0.3, output_f=None):

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

    plt.title(title, fontsize=7)

    #### Top
    axis_key = 0
    ax[axis_key].spines['top'].set_visible(False)
    ax[axis_key].spines['right'].set_visible(False)
    ax[axis_key].spines['left'].set_visible(False)
    ax[axis_key].spines['bottom'].set_visible(False)

    ax[axis_key].set_yticks([])

    subs = 'C>A'
    rect1 = Rectangle((0, 0.15), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(15.75/2, 0.75, subs, fontsize=3.5, weight='normal', ha='center')
    subs = 'C>G'
    rect1 = Rectangle((16,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(16 + 15.75/2 -1, 0.75, subs, fontsize=3.5, weight='normal', ha='center')
    subs = 'C>T'
    rect1 = Rectangle((32,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(32 + 15.75/2 -1, 0.75, subs, fontsize=3.5, weight='normal', ha='center')
    subs = 'T>A'
    rect1 = Rectangle((48,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(48 + 15.75/2 -1, 0.75, subs, fontsize=3.5, weight='normal', ha='center')
    subs = 'T>C'
    rect1 = Rectangle((64,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(64 + 15.75/2 -1, 0.75, subs, fontsize=3.5, weight='normal', ha='center')
    subs = 'T>G'
    rect1 = Rectangle((80,0.1), 15.75, 0.5, color=snv_color[subs])
    ax[axis_key].add_patch(rect1)
    ax[axis_key].text(80 + 15.75/2 -1, 0.75, subs, fontsize=3.5, weight='normal', ha='center')


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


    ########## Y axis
    ylabels_dict = {
        0: '',
        1: 'Frequency'
    }

    ax[axis_key].set_yticks(ylabels)
    ax[axis_key].set_yticklabels(ylabels, fontsize=3.5, color='black', va="top", ha="center")
    ax[axis_key].set_xticks([])
    ax[axis_key].set_xticks(list(range(0, 96)))
    ax[axis_key].set_xticklabels(labels, fontsize=2, color='black', va="top", ha="center")
    plt.xticks(rotation=90)

    for key, value in ylabels_dict.items():
        ax[key].set_ylabel(value, fontsize=4, rotation=90, labelpad=2)
        ax[key].set_axisbelow(True)
        for tick in ax[key].yaxis.get_major_ticks():
            tick.label1.set_fontsize(3.5)

        for location in ['top', 'bottom', 'left', 'right']:
            ax[key].spines[location].set_linewidth(0.25)

        ax[key].tick_params(axis = "x", which = "both", bottom=False)
        ax[key].tick_params(axis = "y", which = "both", left=False)
        ax[key].tick_params(axis='y', which='major', pad=0)
        ax[key].tick_params(axis='x', which='major', pad=0)

    ax[axis_key].set_xlim(xmin=-1, xmax=96)
    ax[axis_key].set_ylim(ymin=0, ymax=ymax)

    if output_f:
        plt.savefig(output_f, bbox_inches='tight', dpi=600)
#     plt.close()
