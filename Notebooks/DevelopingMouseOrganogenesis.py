"""
Created on Sun May 30 19:40:39 2021

@author: trmabdelaal
"""

import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics import f1_score
import sys
sys.path.insert(1,'SIRV/')
from main import SIRV

# load preprocessed scRNA-seq and spatial datasets
RNA = scv.read('SIRV Datasets/Mouse_organogenesis/RNA_adata.h5ad')
SeqFISH = scv.read('SIRV Datasets/Mouse_organogenesis/SeqFISH_Embryo1_adata.h5ad')  # This is embryo1, can be changed to 2 or 3

# Apply SIRV to integrate both datasets and predict the un/spliced expressions
# for the spatially measured genes, additionally transfer 'celltype' label
# annotations from scRNA-seq to spatial data
SeqFISH_imputed = SIRV(SeqFISH,RNA,50,['celltype'])

# Normalize the imputed un/spliced expressions, this will also re-normalize the
# full spatial mRNA 'X', this needs to be undone 
scv.pp.normalize_per_cell(SeqFISH_imputed, enforce=True)

# Undo the double normalization of the full mRNA 'X'
SeqFISH_imputed.X = SeqFISH.to_df()[SeqFISH_imputed.var_names]

# Getting cell type colors
SeqFISH_imputed.uns = SeqFISH.uns

# Data is already annotated, no need for clustering analysis
sc.pl.umap(SeqFISH_imputed, color='celltype_mapped_refined')

sc.pl.scatter(SeqFISH_imputed, basis='xy_loc',color='celltype_mapped_refined')
plt.gca().invert_yaxis()

# Calculating RNA velocities and projecting them on the UMAP embedding and spatial
# coordinates of the tissue
scv.pp.moments(SeqFISH_imputed, n_pcs=50, n_neighbors=30)
scv.tl.velocity(SeqFISH_imputed)
scv.tl.velocity_graph(SeqFISH_imputed)

scv.pl.velocity_embedding_stream(SeqFISH_imputed, basis='umap', color='celltype_mapped_refined',size=40,legend_loc='right')

# Fig. 2A
scv.pl.velocity_embedding_stream(SeqFISH_imputed, basis='xy_loc', color='celltype_mapped_refined',size=40,legend_loc='right')
plt.gca().invert_yaxis()

# Cell-level RNA velocities 
# # Supplementary Fig. S1
scv.pl.velocity_embedding(SeqFISH_imputed,basis='xy_loc', color='celltype_mapped_refined')
plt.gca().invert_yaxis()

# Evaluating label transfer
# loading data from all embryos
SeqFISH = scv.read('SIRV Datasets/Mouse_organogenesis/SeqFISH_AllEmbryos_adata.h5ad')

# Apply SIRV to integrate both datasets and transfer 'celltype' label
# annotations from scRNA-seq to spatial data
SeqFISH_imputed = SIRV(SeqFISH,RNA,50,['celltype'])

# Comparing SeqFISH cell labels with transferred celltype labels
def Norm(x):
    return (x/np.sum(x))

SIRV_f1 = f1_score(SeqFISH_imputed.obs.celltype_mapped_refined, SeqFISH_imputed.obs.celltype,
                   labels=np.unique(SeqFISH_imputed.obs.celltype_mapped_refined),average=None)
SIRV_f1 = pd.Series(SIRV_f1,index = np.unique(SeqFISH_imputed.obs.celltype_mapped_refined))
print("SIRV median F1-score = ", np.median(SIRV_f1))
cont_mat = contingency_matrix(SeqFISH_imputed.obs.celltype_mapped_refined,SeqFISH_imputed.obs.celltype)
df_cont_mat = pd.DataFrame(cont_mat,index = np.unique(SeqFISH_imputed.obs.celltype_mapped_refined), 
                           columns=np.unique(SeqFISH_imputed.obs.celltype))

df_cont_mat = df_cont_mat.apply(Norm,axis=1)
# Supplementary Fig. S2A
plt.figure()
sns.heatmap(df_cont_mat,annot=True,fmt='.2f')
plt.yticks(np.arange(df_cont_mat.shape[0])+0.5,df_cont_mat.index)
plt.xticks(np.arange(df_cont_mat.shape[1])+0.5,df_cont_mat.columns)
