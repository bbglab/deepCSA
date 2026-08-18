"""
From access to deepCSA output file, this module parses the data and returns a pandas 
dataframe with the relevant information to carry out the analysis of selection using bambi. 

The output dataframe contains the following columns:

* sample_id: the sample id
* gene_id: the gene id
* n: the number of mutation counts
* is_syn: 1 if the mutation is synonymous, 0 otherwise
* is_nonsyn: 1 if the mutation is non-synonymous, 0 otherwise
* csqn_type: the consequence type of the mutation (syn or nonsyn)
* offset: the offset for the mutation counts
* trinucleotide context: the trinucleotide context of the mutations (e.g. ACA>T, ACG>G, etc.)

To get this information, the following inputs are required:

* relative mutability file: a file containing the relative mutability of each trinucleotide context for each sample
* mappable genome trinucleotide composition file: trinucleotide composition of one representative megabase in the reference mappable genome
* depth file: a file containing the depth per position for each sample
* observed mutations file: a file containing the observed mutations for each sample

"""

import os
import tqdm
import glob
import json
from functools import cached_property
from itertools import product
import pandas as pd
import numpy as np


# # reference genome trinucleotide counts
# context_counts_path = '/data/bbg/datasets/genomes/context_counts/sigprofiler/context_counts_GRCh38_96.csv'

# # deepCSA output
# deepcsa_output_folder = '/data/bbg/nobackup/prominent/SantPauCH/deepCSA/sp_vhio/custom/2026-05-31_sp_vhio_H_T0sonly_custom'
# depths_path = lambda sample_id: os.path.join(deepcsa_output_folder, f'depths/individual/{sample_id}.depths.annotated.tsv.gz')  # CHROM, POS, CONTEXT, SAMPLE1, SAMPLE2, ...
# mutability_path = lambda sample_id: os.path.join(deepcsa_output_folder, f'selection/omega/preprocessing/mutability_per_sample_gene_context.{sample_id}.tsv')
# mutations_path = lambda sample_id: os.path.join(deepcsa_output_folder, f'selection/omega/preprocessing/mutations_per_sample_gene_impact_context.{sample_id}.tsv')

# # samples
# samples = list(map(lambda x: os.path.basename(x).split('.')[0], glob.glob(os.path.join(deepcsa_output_folder, 'depths/individual/*.depths.annotated.tsv.gz'))))

# impact
csqn_types = {
    'missense': ['missense'],
    'synonymous': ['synonymous'],
    'truncating': ['nonsense', 'splice_region_variant', 'essential_splice']
    }
inv_csqn_types = {v: k for k, vs in csqn_types.items() for v in vs}

# # consensus panel
# consensus_panel_path = os.path.join(deepcsa_output_folder, 'regions/consensuspanels/consensus.all.tsv')
# consensus_panel = pd.read_csv(consensus_panel_path, sep='\t')
# consensus_panel['IMPACT'] = consensus_panel['IMPACT'].map(inv_csqn_types)
# genes = consensus_panel['GENE'].unique().tolist()

# delta is a standard neutral density satisfying the following property:
# weight * delta = expected number of mutations when the mutation density is equivalent to 1 mut/sequenced-Mb
# with weight = np.dot(relative mutability, trinucleotide context content at a given region of the genome
# e.g. missense mutations in TP53


def standard_trinucleotide_contexts():

    subs = [''.join(z) for z in product('CT', 'ACGT') if z[0] != z[1]]  
    flanks = [''.join(z) for z in product('ACGT', repeat=2)]  
    contexts = sorted([(s, f) for s, f in product(subs, flanks)], key=lambda x: (x[0], x[1]))  
    contexts_with_format = [f[0] + s[0] + f[1] + '>' + s[1] for s, f in contexts]
    return contexts_with_format


standard_contexts = standard_trinucleotide_contexts()


# def get_megabase_content():

#     """
#     Trinucleotide abundance in one representative megabase of the reference mappable genome.
#     Output format: 96-channel list of trinucleotide counts in the standard order
#     """

#     df = pd.read_csv(context_counts_path, index_col=0)
#     dg = df.sum(axis=1) / df.sum(axis=1).sum()
#     d = dict(dg.apply(lambda x: int(1e6 * x)))
    
#     standard_contexts = standard_trinucleotide_contexts()
#     res = [d[c[:3]] for c in standard_contexts]
#     return res


# megabase_content = np.array(get_megabase_content())


# def get_relative_mutability(sample_id):

#     """
#     Relative mutability of each trinucleotide context for a given sample.
#     Output format: 96-channel list of relative mutability in the standard order
#     """

#     df = pd.read_csv(mutability_path(sample_id), sep='\t')
#     dg = df.groupby('CONTEXT_MUT').agg({sample_id: 'sum'}).sort_values('CONTEXT_MUT')[sample_id] / df[sample_id].sum()
#     dg = dg.reindex(standard_contexts, fill_value=0)
#     return dg


# def get_offsets_aggregated(sample_id, consensus_panel):

#     """
#     Offset per sample, gene and consequence type
#     Format: 96-channel list of offsets for each trinucleotide context in the standard order
#     """

#     rel_mut = get_relative_mutability(sample_id)
#     one_megabase_weight = np.dot(megabase_content, rel_mut)
#     delta = 1 / one_megabase_weight


#     depths = pd.read_csv(depths_path(sample_id), sep='\t')
#     panel_with_depths = pd.merge(consensus_panel, depths, on=['CHROM', 'POS', 'CONTEXT'], how='left')
#     dg = panel_with_depths.groupby(['IMPACT', 'GENE', 'CONTEXT_MUT']).agg({sample_id: 'sum'})

#     # create a multiindex with all combinations of gene, consequence type and trinucleotide context
#     index = pd.MultiIndex.from_product([csqn_types.keys(), genes, standard_contexts], names=['IMPACT', 'GENE', 'CONTEXT_MUT'])
#     dg = dg.reindex(index, fill_value=0).reset_index()  # dataframe with columns: IMPACT, GENE, CONTEXT_MUT, sample_id

#     # Pivot the dataframe into a matrix shape
#     matrix = dg.pivot(index=['GENE', 'IMPACT'], columns='CONTEXT_MUT', values=sample_id).fillna(0)
#     # Ensure rel_mut is perfectly aligned with the matrix columns
#     rel_mut_aligned = rel_mut.reindex(matrix.columns).fillna(0)
#     # Perform the matrix-vector dot product and reset the index
#     df_result = matrix.dot(rel_mut_aligned).reset_index(name='weight')
#     df_result['offset'] = df_result['weight'] * delta

#     return df_result


# def build_table_aggregated():
#     """
#     Build a table with the following columns:
#     sample_id, gene_id, csqn_type, n, n_syn, offset
#     """

#     data = []
#     for s in tqdm.tqdm(samples):
#         dg = get_offsets_aggregated(s, consensus_panel)
#         muts = pd.read_csv(mutations_path(s), sep='\t')
#         muts['IMPACT'] = muts['IMPACT'].map(inv_csqn_types)
#         muts_grouped = muts.groupby(['GENE', 'IMPACT']).agg({s: 'sum'}).reset_index()
#         dg = pd.merge(dg, muts_grouped, on=['GENE', 'IMPACT'], how='left').fillna(0)
#         dg['sample_id'] = s
#         dg.rename(columns={s: 'n'}, inplace=True)
#         data.append(dg)

#     df = pd.concat(data, axis=0)
#     df['is_nonsyn'] = df['IMPACT'].apply(lambda x: 0 if x == 'synonymous' else 1)
#     df['log_offset'] = df['offset'].apply(lambda x: np.log(x + 1e-6))
    
#     # TODO: collapse the table so that all consequence types are represented within each row
    
#     return df


# def postprocess_aggregated(df):

#     # basic filters
#     df_nonzero = df[df['offset'] > 0]  # gene-impact must have positive offset = being covered by the sequencing experiment
    
#     # sample name filters
#     df_nonzero = df_nonzero[df_nonzero['sample_id'].str.startswith('CHa')]
    
#     # data type formatting
#     df_nonzero['n'] = df_nonzero['n'].astype(int)

#     # re-arrage data items
#     res = {'gene_id': [], 'sample_id': [], 'offset_syn': [], 'offset_trunc': [], 'offset_mis': [], 'n_syn': [], 'n_trunc': [], 'n_mis': []}
#     for gene, sample in df_nonzero.groupby(['GENE', 'sample_id']).groups.keys():
        
#         dg = df_nonzero[(df_nonzero['GENE'] == gene) & (df_nonzero['sample_id'] == sample)]
#         res['gene_id'].append(gene)
#         res['sample_id'].append(sample)

#         for impact in ['missense', 'synonymous', 'truncating']:
#             if (impact == 'synonymous'):
#                 if dg[dg['IMPACT'] == impact].shape[0] > 0:
#                     res['offset_syn'].append(dg[dg['IMPACT'] == impact]['offset'].values[0])
#                     res['n_syn'].append(dg[dg['IMPACT'] == impact]['n'].values[0])
#                 else:
#                     res['offset_syn'].append(None)
#                     res['n_syn'].append(None)
#             if (impact == 'truncating'):
#                 if dg[dg['IMPACT'] == impact].shape[0] > 0:
#                     res['offset_trunc'].append(dg[dg['IMPACT'] == impact]['offset'].values[0])
#                     res['n_trunc'].append(dg[dg['IMPACT'] == impact]['n'].values[0])
#                 else:
#                     res['offset_trunc'].append(None)
#                     res['n_trunc'].append(None)
#             if (impact == 'missense'): 
#                 if dg[dg['IMPACT'] == impact].shape[0] > 0:
#                     res['offset_mis'].append(dg[dg['IMPACT'] == impact]['offset'].values[0])
#                     res['n_mis'].append(dg[dg['IMPACT'] == impact]['n'].values[0])
#                 else:
#                     res['offset_mis'].append(None)
#                     res['n_mis'].append(None)
#     res = pd.DataFrame(res)

#     # data type formatting
#     res.dropna(inplace=True)
#     res['n_syn'] = res['n_syn'].astype(int)
#     res['n_mis'] = res['n_mis'].astype(int)
#     res['n_trunc'] = res['n_trunc'].astype(int)

#     # add the principal components roadmap of epigenomics

#     pc_roadmap = pd.read_csv('./covariates_hg19_hg38_epigenome_pcawg.tsv', sep='\t', index_col=0)
#     pc_roadmap.reset_index(inplace=True)
#     pc_roadmap.rename(columns={'index': 'gene_id'}, inplace=True)
#     res = pd.merge(res, pc_roadmap, on=['gene_id'], how='left')

#     res.dropna(inplace=True)
#     return res


class deepCSAparser:

    @cached_property
    def _from_json(self):
        
        with open(self.config_path, "r") as f:
            data = json.load(f)
        
        deepcsa_output_path = data['path']['deepcsa_output_path']
        context_counts_path = data['path']['context_counts_path']
        covariates_path = data['path']['covariates_path']
        samples_json_path = data['path'].get('samples', None)
        sample_filter = data['filter']['sample']
        
        return deepcsa_output_path, context_counts_path, covariates_path, samples_json_path, sample_filter
    
    @property
    def deepcsa_output_path(self):
        return self._from_json[0]
    
    @property
    def context_counts_path(self):
        return self._from_json[1]
    
    @property
    def covariates_path(self):
        return self._from_json[2]

    @property
    def samples_list(self):
        if self._from_json[3] is not None:
            with open(self._from_json[3], "r") as f:
                samples = list(json.load(f).keys())
            return samples
        return None

    @property
    def sample_filter(self):
        return self._from_json[4]


    def __init__(self, config_path):

        self.config_path = config_path

        # deepCSA output
        self.depths_path = lambda sample_id: os.path.join(self.deepcsa_output_path, f'depths/individual/{sample_id}.depths.annotated.tsv.gz')  # CHROM, POS, CONTEXT, SAMPLE1, SAMPLE2, ...
        self.mutability_path = lambda sample_id: os.path.join(self.deepcsa_output_path, f'selection/omega/preprocessing/mutability_per_sample_gene_context.{sample_id}.tsv')
        self.mutations_path = lambda sample_id: os.path.join(self.deepcsa_output_path, f'selection/omega/preprocessing/mutations_per_sample_gene_impact_context.{sample_id}.tsv')

        # samples
        if self.samples_list is not None:
            self.samples = self.samples_list
        else:
            self.samples = list(map(lambda x: os.path.basename(x).split('.')[0], 
                glob.glob(os.path.join(self.deepcsa_output_path, 'depths/individual/*.depths.annotated.tsv.gz'))))

        # impact
        self.csqn_types = {
            'missense': ['missense'],
            'synonymous': ['synonymous'],
            'truncating': ['nonsense', 'splice_region_variant', 'essential_splice']
            }
        self.inv_csqn_types = {v: k for k, vs in self.csqn_types.items() for v in vs}

        # consensus panel
        consensus_panel_path = os.path.join(self.deepcsa_output_path, 'regions/consensuspanels/consensus.all.tsv')
        self.consensus_panel = pd.read_csv(consensus_panel_path, sep='\t')
        self.consensus_panel['IMPACT'] = self.consensus_panel['IMPACT'].map(self.inv_csqn_types)
        self.genes = self.consensus_panel['GENE'].unique().tolist()

        # standard trinucleotide contexts
        self.standard_contexts = standard_trinucleotide_contexts()

        # mappable genome average per megabase content
        self.megabase_content = np.array(self.get_megabase_content())
    

    def get_megabase_content(self):

        """
        Trinucleotide abundance in one representative megabase of the reference mappable genome.
        Output format: 96-channel list of trinucleotide counts in the standard order
        """

        df = pd.read_csv(self.context_counts_path, index_col=0)
        dg = df.sum(axis=1) / df.sum(axis=1).sum()
        d = dict(dg.apply(lambda x: int(1e6 * x)))
        
        standard_contexts = standard_trinucleotide_contexts()
        res = [d[c[:3]] for c in standard_contexts]
        return res


    def get_relative_mutability(self, sample_id):

        """
        Relative mutability of each trinucleotide context for a given sample.
        Output format: 96-channel list of relative mutability in the standard order
        """

        df = pd.read_csv(self.mutability_path(sample_id), sep='\t')
        dg = df.groupby('CONTEXT_MUT').agg({sample_id: 'sum'}).sort_values('CONTEXT_MUT')[sample_id] / df[sample_id].sum()
        dg = dg.reindex(standard_contexts, fill_value=0)
        return dg


    def get_offsets_aggregated(self, sample_id, consensus_panel):

        """
        Offset per sample, gene and consequence type
        Format: 96-channel list of offsets for each trinucleotide context in the standard order
        """

        rel_mut = self.get_relative_mutability(sample_id)
        one_megabase_weight = np.dot(self.megabase_content, rel_mut)
        delta = 1 / one_megabase_weight

        depths = pd.read_csv(self.depths_path(sample_id), sep='\t')
        panel_with_depths = pd.merge(self.consensus_panel, depths, on=['CHROM', 'POS', 'CONTEXT'], how='left')
        dg = panel_with_depths.groupby(['IMPACT', 'GENE', 'CONTEXT_MUT']).agg({sample_id: 'sum'})

        # create a multiindex with all combinations of gene, consequence type and trinucleotide context
        index = pd.MultiIndex.from_product([csqn_types.keys(), self.genes, self.standard_contexts], names=['IMPACT', 'GENE', 'CONTEXT_MUT'])
        dg = dg.reindex(index, fill_value=0).reset_index()  # dataframe with columns: IMPACT, GENE, CONTEXT_MUT, sample_id

        # Pivot the dataframe into a matrix shape
        matrix = dg.pivot(index=['GENE', 'IMPACT'], columns='CONTEXT_MUT', values=sample_id).fillna(0)
        # Ensure rel_mut is perfectly aligned with the matrix columns
        rel_mut_aligned = rel_mut.reindex(matrix.columns).fillna(0)
        # Perform the matrix-vector dot product and reset the index
        df_result = matrix.dot(rel_mut_aligned).reset_index(name='weight')
        df_result['offset'] = df_result['weight'] * delta

        return df_result


    def build_table_aggregated(self):
        """
        Build a table with the following columns:
        sample_id, gene_id, csqn_type, n, n_syn, offset
        """

        data = []
        for s in tqdm.tqdm(self.samples):
            dg = self.get_offsets_aggregated(s, self.consensus_panel)
            muts = pd.read_csv(self.mutations_path(s), sep='\t')
            muts['IMPACT'] = muts['IMPACT'].map(inv_csqn_types)
            muts_grouped = muts.groupby(['GENE', 'IMPACT']).agg({s: 'sum'}).reset_index()
            dg = pd.merge(dg, muts_grouped, on=['GENE', 'IMPACT'], how='left').fillna(0)
            dg['sample_id'] = s
            dg.rename(columns={s: 'n'}, inplace=True)
            data.append(dg)

        df = pd.concat(data, axis=0)
        df['is_nonsyn'] = df['IMPACT'].apply(lambda x: 0 if x == 'synonymous' else 1)
        df['log_offset'] = df['offset'].apply(lambda x: np.log(x + 1e-6))
        
        return df


    def postprocess_aggregated(self, df):

        # basic filters
        df_nonzero = df[df['offset'] > 0]  # gene-impact must have positive offset = being covered by the sequencing experiment
        
        # sample name filters
        res = []
        for value in self.sample_filter:
            res.append(df_nonzero[df_nonzero['sample_id'].str.startswith(value)])
        df_nonzero = pd.concat(res, axis=0)
        
        # data type formatting
        df_nonzero['n'] = df_nonzero['n'].astype(int)

        # re-arrage data items
        res = {'gene_id': [], 'sample_id': [], 'offset_syn': [], 'offset_trunc': [], 'offset_mis': [], 'n_syn': [], 'n_trunc': [], 'n_mis': []}
        for gene, sample in df_nonzero.groupby(['GENE', 'sample_id']).groups.keys():
            
            dg = df_nonzero[(df_nonzero['GENE'] == gene) & (df_nonzero['sample_id'] == sample)]
            res['gene_id'].append(gene)
            res['sample_id'].append(sample)

            for impact in ['missense', 'synonymous', 'truncating']:
                if (impact == 'synonymous'):
                    if dg[dg['IMPACT'] == impact].shape[0] > 0:
                        res['offset_syn'].append(dg[dg['IMPACT'] == impact]['offset'].values[0])
                        res['n_syn'].append(dg[dg['IMPACT'] == impact]['n'].values[0])
                    else:
                        res['offset_syn'].append(None)
                        res['n_syn'].append(None)
                if (impact == 'truncating'):
                    if dg[dg['IMPACT'] == impact].shape[0] > 0:
                        res['offset_trunc'].append(dg[dg['IMPACT'] == impact]['offset'].values[0])
                        res['n_trunc'].append(dg[dg['IMPACT'] == impact]['n'].values[0])
                    else:
                        res['offset_trunc'].append(None)
                        res['n_trunc'].append(None)
                if (impact == 'missense'): 
                    if dg[dg['IMPACT'] == impact].shape[0] > 0:
                        res['offset_mis'].append(dg[dg['IMPACT'] == impact]['offset'].values[0])
                        res['n_mis'].append(dg[dg['IMPACT'] == impact]['n'].values[0])
                    else:
                        res['offset_mis'].append(None)
                        res['n_mis'].append(None)
        res = pd.DataFrame(res)

        # data type formatting
        res.dropna(inplace=True)
        res['n_syn'] = res['n_syn'].astype(int)
        res['n_mis'] = res['n_mis'].astype(int)
        res['n_trunc'] = res['n_trunc'].astype(int)

        # add the principal components roadmap of epigenomics
        pc_roadmap = pd.read_csv(self.covariates_path, sep='\t', index_col=0)
        pc_roadmap.reset_index(inplace=True)
        pc_roadmap.rename(columns={'index': 'gene_id'}, inplace=True)
        res = pd.merge(res, pc_roadmap, on=['gene_id'], how='left')

        res.dropna(inplace=True)
        return res

    
    def parse_dataset(self):

        df = self.build_table_aggregated()
        dg = self.postprocess_aggregated(df)
        return dg


    
