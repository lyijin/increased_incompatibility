#!/usr/bin/env python3

"""
> plot_clustermap.symbiosis_genes.py <

Plots a clustered per-gene, per-strain heatmap of beta values for
symbiosis-related genes, to visually ascertain whether SSB01 behaves
differently from CC7 and H2.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# read genes of interest
guoxin_all = open('guoxin_meta.all.txt').read().strip().split('\n')
guoxin_up = open('guoxin_meta.up.txt').read().strip().split('\n')
guoxin_down = open('guoxin_meta.down.txt').read().strip().split('\n')

# read data
data = {}
for strain in ['CC7', 'SSB01', 'H2']:
    temp = pd.read_csv(
        '{}_temp32_vs_temp25_spliced_rt.csv'.format(strain.lower()),
        index_col=0, usecols=[0, 2, 3])
    
    # subset by genes of interest
    temp = temp[temp.index.isin(guoxin_up)]
    
    # force beta to zero if gene is not significant (i.e. q_value is not < 0.05)
    temp['b'] = temp.apply(lambda x: x['b'] if x['qval'] < 0.05 else 0, axis=1)
    
    # rename columns so that they can be left-joined
    temp.columns = ['_'.join([strain, x]) for x in temp.columns]
    
    data[strain] = temp

# concat dataframes together as one
data = pd.concat([data[x] for x in data], axis=1)

# subselect betas. qvalues are not impt because non-significant genes had had
# their betas set to 0
data = data.filter(regex='_b')
data.columns = [x.replace('_b', '') for x in data.columns]

# remove rows that only consists of NAs and/or 0
data = data.replace(0, np.NaN)
data = data.dropna(how='all')
data = data.replace(np.NaN, 0)

# convert beta from ln to log2
data = data.applymap(lambda x: x / np.log(2))

# enforce CC7 / SSB01 / H2 order
data = data[['CC7', 'SSB01', 'H2']]

# get row colours
rowcolors = data.index.map(lambda x: '#d53e4f')

# plot the heatmap!
fig, ax = plt.subplots()
#sns.set(font_scale=1.2)

cg = sns.clustermap(data,
                    method='ward',
                    figsize=(3, 6),
                    vmin=-2, vmax=2,
                    col_cluster=False,
                    row_colors=rowcolors,
                    cmap='RdBu_r', annot=False, linewidth=0,
                    yticklabels=False,
                    cbar_kws={'ticks': [-2, -1, 0, 1, 2]})

# save figure
fig.tight_layout()
fig = plt.gcf()

# without bbox_inches, the saved figure has truncated axes.
output_filename = 'clustermap.symbiosis_genes.up.pdf'
fig.savefig(output_filename, bbox_inches='tight')

# also save the values of the heatmap in a plain text file
cg.data2d.to_csv('clustermap.symbiosis_genes.up.tsv', sep='\t')
