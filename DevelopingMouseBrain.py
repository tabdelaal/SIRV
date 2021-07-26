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
import sys
sys.path.insert(1,'SIRV/')
from main import SIRV
import warnings
warnings.filterwarnings('ignore')

# load preprocessed scRNA-seq and spatial datasets
RNA = scv.read('SRV_data/RNA_adata.h5ad')
HybISS = scv.read('SRV_data/HybISS_adata.h5ad')

# Apply SRV to integrate both datasets and predict the un/spliced expressions
# for the spatially measured genes, additionally transfer 'Region', 'Class' and
# 'Subclass' label annotations from scRNA-seq to spatial data
HybISS_imputed = SIRV(HybISS,RNA,50,['Tissue','Region','Class','Subclass'])

# Normalize the imputed un/spliced expressions, this will also re-normalize the
# full spatial mRNA 'X', this needs to be undone 
scv.pp.normalize_per_cell(HybISS_imputed, enforce=True)

# Undo the double normalization of the full mRNA 'X'
HybISS_imputed.X = HybISS.to_df()[HybISS_imputed.var_names]

# Zero mean and unit variance scaling, PCA, building neibourhood graph, running
# umap and cluster the HybISS spatial data using Leiden clustering
sc.pp.scale(HybISS_imputed)
sc.tl.pca(HybISS_imputed)
sc.pl.pca_variance_ratio(HybISS_imputed, n_pcs=50, log=True)
sc.pp.neighbors(HybISS_imputed, n_neighbors=30, n_pcs=30)
sc.tl.umap(HybISS_imputed)
sc.tl.leiden(HybISS_imputed)
# Fig. 2A
sc.pl.umap(HybISS_imputed, color='leiden')
# Fig. 2B
sc.pl.scatter(HybISS_imputed, basis='xy_loc',color='leiden')

# Calculating RNA velocities and projecting them on the UMAP embedding and spatial
# coordinates of the tissue
scv.pp.moments(HybISS_imputed, n_pcs=30, n_neighbors=30)
scv.tl.velocity(HybISS_imputed)
scv.tl.velocity_graph(HybISS_imputed)
# Fig. 2C
scv.pl.velocity_embedding_stream(HybISS_imputed, basis='umap', color='leiden')
# Fig. 2D
scv.pl.velocity_embedding_stream(HybISS_imputed, basis='xy_loc', color='leiden',size=60,legend_fontsize=4,legend_loc='right')

# Cell-level RNA velocities 
# Fig. 3
scv.pl.velocity_embedding(HybISS_imputed,basis='xy_loc', color='leiden')

# Visualizing transferred label annotations on UMAP embedding and spatial coordinates
# Fig. 4A
sc.pl.umap(HybISS_imputed, color='Region')
# Fig. 4B
sc.pl.scatter(HybISS_imputed, basis='xy_loc',color='Region')
# Fig. 4C
sc.pl.umap(HybISS_imputed, color='Subclass')
# Fig. 4D
sc.pl.scatter(HybISS_imputed, basis='xy_loc',color='Subclass')
# Supplementary Fig. S3A
sc.pl.umap(HybISS_imputed, color='Class')
# Supplementary Fig. S3B
sc.pl.scatter(HybISS_imputed, basis='xy_loc',color='Class')

# Intepretation of RNA velocities using transferred label annotations
# Fig. 5
scv.pl.velocity_embedding(HybISS_imputed,basis='xy_loc', color='Subclass')
# Supplementary Fig. S3C
scv.pl.velocity_embedding(HybISS_imputed,basis='xy_loc', color='Class')

# Comparing cell clusters with transferred 'Subclass' and 'Class' annotations
def Norm(x):
    return (x/np.sum(x))

# Subclass annotation
cont_mat = contingency_matrix(HybISS_imputed.obs.leiden.astype(np.int_),HybISS_imputed.obs.Subclass)
df_cont_mat = pd.DataFrame(cont_mat,index = np.unique(HybISS_imputed.obs.leiden.astype(np.int_)), 
                           columns=np.unique(HybISS_imputed.obs.Subclass))

df_cont_mat = df_cont_mat.apply(Norm,axis=1)
# Supplementary Fig. S5A
plt.figure()
sns.heatmap(df_cont_mat,annot=True,fmt='.2f')
plt.yticks(np.arange(df_cont_mat.shape[0])+0.5,df_cont_mat.index)
plt.xticks(np.arange(df_cont_mat.shape[1])+0.5,df_cont_mat.columns)

# Class annotation
cont_mat = contingency_matrix(HybISS_imputed.obs.leiden.astype(np.int_),HybISS_imputed.obs.Class)
df_cont_mat = pd.DataFrame(cont_mat,index = np.unique(HybISS_imputed.obs.leiden.astype(np.int_)), 
                           columns=np.unique(HybISS_imputed.obs.Class))

df_cont_mat = df_cont_mat.apply(Norm,axis=1)
# Supplementary Fig. S5B
plt.figure()
sns.heatmap(df_cont_mat,annot=True,fmt='.2f')
plt.yticks(np.arange(df_cont_mat.shape[0])+0.5,df_cont_mat.index)
plt.xticks(np.arange(df_cont_mat.shape[1])+0.5,df_cont_mat.columns)
