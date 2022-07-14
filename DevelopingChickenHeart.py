"""
Created on Fri Jan 14 19:26:23 2022

@author: trmabdelaal
"""

import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
import sys
sys.path.insert(1,'SIRV/')
from main import SIRV

# load preprocessed scRNA-seq and spatial datasets
RNA = scv.read('SIRV Datasets/Chicken_heart/RNA_D14_adata.h5ad')
Visium = scv.read('SIRV Datasets/Chicken_heart/Visium_D14_adata.h5ad')

# Apply SIRV to integrate both datasets and predict the un/spliced expressions
# for the spatially measured gene
Visium_imputed = SIRV(Visium,RNA,50)

# Normalize the imputed un/spliced expressions, this will also re-normalize the
# full spatial mRNA 'X', this needs to be undone
scv.pp.normalize_per_cell(Visium_imputed, enforce=True)

# Undo the double normalization of the full mRNA 'X'
Visium_imputed.X = Visium.to_df()[Visium_imputed.var_names]

#### Measured RNA velocity

# Zero mean and unit variance scaling, PCA, building neibourhood graph, running
# umap 
Visium = Visium[:,Visium_imputed.var_names]
sc.pp.scale(Visium)
sc.tl.pca(Visium)
sc.pl.pca_variance_ratio(Visium, n_pcs=50, log=True)
sc.pp.neighbors(Visium, n_neighbors=30, n_pcs=20)
sc.tl.umap(Visium)

sc.pl.umap(Visium, color='celltype_prediction')
sc.pl.scatter(Visium, basis='xy_loc',color='celltype_prediction')
plt.gca().invert_xaxis()

# Calculating RNA velocities and projecting them on the UMAP embedding and spatial
# coordinates of the tissue
scv.pp.moments(Visium, n_pcs=20, n_neighbors=30)
scv.tl.velocity(Visium)
scv.tl.velocity_graph(Visium)

scv.pl.velocity_embedding_stream(Visium, basis='umap', color='celltype_prediction',legend_loc='right')
# Fig. 3A
scv.pl.velocity_embedding_stream(Visium, basis='xy_loc', color='celltype_prediction',legend_loc='right')
plt.gca().invert_xaxis()

# Spot-level RNA velocities
# Supplementary Fig. S3A
scv.pl.velocity_embedding(Visium,basis='xy_loc', color='celltype_prediction')
plt.gca().invert_xaxis()


#### SIRV estimated RNA velocity

# Similar umap and cell type colors to the original Visium data 
Visium_imputed.obsm['X_umap'] = Visium.obsm['X_umap']
Visium_imputed.uns['celltype_prediction_colors'] = Visium.uns['celltype_prediction_colors']

# Calculating RNA velocities and projecting them on the UMAP embedding and spatial
# coordinates of the tissue
scv.pp.moments(Visium_imputed, n_pcs=20, n_neighbors=30)
scv.tl.velocity(Visium_imputed)
scv.tl.velocity_graph(Visium_imputed)

scv.pl.velocity_embedding_stream(Visium_imputed, basis='umap', color='celltype_prediction',legend_loc='right')
# Fig. 3B
scv.pl.velocity_embedding_stream(Visium_imputed, basis='xy_loc', color='celltype_prediction',legend_loc='right')
plt.gca().invert_xaxis()

# Spot-level RNA velocities
# Supplementary Fig. S3B
scv.pl.velocity_embedding(Visium_imputed,basis='xy_loc', color='celltype_prediction')
plt.gca().invert_xaxis()

# Quantitative evaluation: Weighted similarity score

Sim_velo = np.diag(cosine_similarity(Visium.layers['velocity'],Visium_imputed.layers['velocity']))
Sim_xy = np.diag(cosine_similarity(Visium.obsm['velocity_xy_loc'],Visium_imputed.obsm['velocity_xy_loc']))

Mag_velo = np.linalg.norm(Visium.layers['velocity'],axis=1)
Mag_xy = np.linalg.norm(Visium.obsm['velocity_xy_loc'],axis=1)

# Supplementary Fig. S3C
fig,ax=plt.subplots()
scv.pl.scatter(Visium,basis='xy_loc', color=np.log2(Mag_xy),ax=ax)
ax.invert_xaxis()

# Supplementary Fig. S3D
fig,ax=plt.subplots()
scv.pl.scatter(Visium,basis='xy_loc', color=Sim_xy,ax=ax)
ax.invert_xaxis()

Mag_velo = Mag_velo/np.sum(Mag_velo)
Mag_xy = Mag_xy/np.sum(Mag_xy)

print("Weighted Similarity score for high-dimensional velocity = ", np.dot(Sim_velo,Mag_velo))
print("Weighted Similarity score for spatial velocity = ", np.dot(Sim_xy,Mag_xy))

# Fig. 3C
fig,ax=plt.subplots()
scv.pl.scatter(Visium,basis='xy_loc', color=Mag_xy*Sim_xy,ax=ax)
ax.invert_xaxis()
