import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler



clean_maf_file = f"all_samples.clean.mutations.tsv"
somatic_maf_file = f"all_samples.somatic.mutations.tsv"

somatic_maf = pd.read_table(somatic_maf_file)
print(somatic_maf.shape)
### Explore SNPs PCA
clean_mutations_all = pd.read_table(clean_maf_file)
clean_mutations_all = clean_mutations_all.merge(somatic_maf[["SAMPLE_ID", "MUT_ID", "VAF"]], on=["SAMPLE_ID", "MUT_ID"],
                                                suffixes = ("", "_somatic"),
                                                how = 'outer')
clean_mutations_all = clean_mutations_all[~(clean_mutations_all["VAF_somatic"].notnull())].reset_index(drop = True)
clean_mutations_all

# create bins from 0-0.25, 0.25-0.75, 0.5-0.75, 0.75-1
clean_mutations_all["genotype_bin"] = np.digitize(clean_mutations_all["VAF"], right = True, bins = [0, 0.25, 0.75, 1])
sns.boxplot(data = clean_mutations_all,
            x = "genotype_bin",
            y = "VAF",
            showfliers=False
            )
sns.stripplot(data = clean_mutations_all,
            x = "genotype_bin",
            y = "VAF",
            )
subset_clean_mutations_all = clean_mutations_all[clean_mutations_all["genotype_bin"].isin([1, 2])].reset_index(drop = True)
subset_clean_mutations_all["genotype_bin"].value_counts()
# subset_clean_mutations_all = clean_mutations_all
# Build data for PCA: samples (rows) x genes (columns)
pca_data = subset_clean_mutations_all.pivot_table(index='SAMPLE_ID', columns='MUT_ID', values='genotype_bin', aggfunc='first')
pca_data.fillna(0, inplace=True)  # Fill NaN with 0 (no mutation)
# # Build data for PCA: samples (rows) x genes (columns)
# pca_data = (
#     mutdensity[

#         # Filter table
#         (mutdensity['SAMPLE_ID'].str.contains('IC')) &
#         (mutdensity['GENE'].isin(genes_pos["all_genes"])) & # only for positively selected genes selected for regressions
#         (~mutdensity['GENE'].str.contains('--')) & # discard exons 
#         (mutdensity['REGIONS'] == 'protein_affecting') &
#         (mutdensity['MUTTYPES'] == 'all_types')
#     ]
#     .pivot_table(index='SAMPLE_ID', columns='GENE', values='MUTDENSITY_MB', aggfunc='first')
# )

# display(pca_data.head())

# # We need to plot the normalized mutation density data across genes in that sample (mutdensity in geneX /mutdensity sum of all subset of genes)
# relative_pca_data = pca_data.apply(lambda x: x / x.sum(), axis=1)
# # display(relative_pca_data)
# pca_data = relative_pca_data.copy()
# display(pca_data.head())
pca_data.describe()


# Scale and run PCA
X_scaled = StandardScaler().fit_transform(pca_data)
number_of_components = 6  # Choose the components you might want to test, arbitrary decision
pca = PCA(n_components=number_of_components)
pca_result = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(number_of_components)], index=pca_data.index)


# # Merge metadata for coloring
# pca_df = pca_df.join(
#     full_cohort_metadata.set_index('SAMPLE_ID'),
#     how='left'
# )

## ELBOW PLOT
# Do these number of components explain all the variance of the data? 
plt.plot(pca.explained_variance_ratio_, marker='o', linestyle='-')
plt.ylim(0, pca.explained_variance_ratio_.max() + 0.05)
plt.xlabel('PC')
plt.ylabel('Percent of variance explained')
plt.show()


component_x = 'PC1'
component_x_index = 0
component_y = 'PC2' # modify for any component that explains most of the sample variance from the elbow plot
component_y_index = 1

fig, ax = plt.subplots(figsize=(6, 6))
sc = ax.scatter(pca_df[component_x], pca_df[component_y],
                s=80, alpha=0.85, edgecolors='k', linewidths=0.3)
for s in pca_df.index:
    ax.annotate(s, (pca_df.loc[s, component_x], pca_df.loc[s, component_y]),
                fontsize=6, ha='left', va='bottom')
ax.set_xlabel(f'{component_x} ({pca.explained_variance_ratio_[component_x_index]*100:.1f}% variance explained)')
ax.set_ylabel(f'{component_y} ({pca.explained_variance_ratio_[component_y_index]*100:.1f}% variance explained)')
ax.set_title('PCA of SNPs within Panel Genes')
plt.tight_layout()
plt.show()

# # PCA
# for variable_to_plot in ['cohort', 'race_unif']:
#     component_x = 'PC1'
#     component_x_index = 0
#     component_y = 'PC2' # modify for any component that explains most of the sample variance from the elbow plot
#     component_y_index = 1

#     if pca_df[variable_to_plot].dtype == 'O':
#         fig, ax = plt.subplots(figsize=(6, 6))
#         groups = pca_df[variable_to_plot].dropna().unique()
#         palette = sns.color_palette('Set2', n_colors=len(groups))
#         for grp, color in zip(groups, palette):
#             mask = pca_df[variable_to_plot] == grp
#             ax.scatter(pca_df.loc[mask, component_x], pca_df.loc[mask, component_y],
#                     label=grp, color=color, s=80, alpha=0.85, edgecolors='k', linewidths=0.3)
#         for s in pca_df.index:
#             ax.annotate(s, (pca_df.loc[s, component_x], pca_df.loc[s, component_y]),
#                         fontsize=6, ha='left', va='bottom')
#         ax.set_xlabel(f'{component_x} ({pca.explained_variance_ratio_[component_x_index]*100:.1f}% variance explained)')
#         ax.set_ylabel(f'{component_y} ({pca.explained_variance_ratio_[component_y_index]*100:.1f}% variance explained)')
#         ax.set_title('PCA of Mutation Density within Panel Genes')
#         ax.legend(title='Sample location', bbox_to_anchor=(1.01, 1), loc='upper left')
#         plt.tight_layout()
#         plt.show()

#     else:
#         fig, ax = plt.subplots(figsize=(6, 6))
#         sc = ax.scatter(pca_df[component_x], pca_df[component_y], c=pca_df[variable_to_plot],
#                         cmap='viridis', s=80, alpha=0.85, edgecolors='k', linewidths=0.3)
#         for s in pca_df.index:
#             ax.annotate(s, (pca_df.loc[s, component_x], pca_df.loc[s, component_y]),
#                         fontsize=6, ha='left', va='bottom')
#         ax.set_xlabel(f'{component_x} ({pca.explained_variance_ratio_[component_x_index]*100:.1f}% variance explained)')
#         ax.set_ylabel(f'{component_y} ({pca.explained_variance_ratio_[component_y_index]*100:.1f}% variance explained)')
#         ax.set_title('PCA of SNPs within Panel Genes')
#         cbar = plt.colorbar(sc)
#         cbar.set_label(variable_to_plot)
#         plt.tight_layout()
#         plt.show()

#     print(f"Feature matrix: {pca_data.shape[0]} samples × {pca_data.shape[1]} mutations")
#     print(f"Variance explained: {component_x}={pca.explained_variance_ratio_[component_x_index]*100:.1f}%, {component_y}={pca.explained_variance_ratio_[component_y_index]*100:.1f}%")

