""" SIRV [1]
@author: Tamim Abdelaal
This function integrates two single-cell datasets, spatial and scRNA-seq, and 
predicts the spliced and unspliced expression of the spatially measured genes
from the scRNA-seq data, in order to apply RNA velocity analysis.
The integration is performed using the domain adaption method PRECISE [2]
	
References
-------
    [1] Abdelaal T., Lelieveldt B.P.F., Reiders M.J.T., Mahfouz A. (2021)
    SIRV: Spatial inference of RNA velocity at the single-cell resolution
    [2] Mourragui S., Loog M., Reinders M.J.T., Wessels L.F.A. (2019)
    PRECISE: A domain adaptation approach to transfer predictors of drug response
    from pre-clinical models to tumors
"""

import numpy as np
import pandas as pd
import scipy.stats as st
import scvelo as scv
from sklearn.neighbors import NearestNeighbors
from principal_vectors import PVComputation

def SIRV(Spatial_data,RNA_data,n_pv,metadata_to_transfer=None):
    """
        @author: Tamim Abdelaal
        This function integrates two single-cell datasets, spatial and scRNA-seq, 
        and enhance the spatial data by predicting the expression of the spatially 
        unmeasured genes from the scRNA-seq data.
        
        Parameters
        -------
        Spatial_data : scvelo Anndata Object
            Normalized Spatial data matrix (cells X genes).
        RNA_data : scvelo Anndata Object
            Contains normalized gene expression, spliced and unspliced data 
            (cells X genes), and cellular metadata (obs).
        n_pv : int
            Number of principal vectors to find from the independently computed
            principal components, and used to align both datasets. This should
            be <= number of shared genes between the two datasets.
        metadata_to_predict : str array 
            list of metadata identifiers (obs.keys) to be transferred from the 
            scRNA-seq to the spatial data. Default is None, in such case, no 
            metadata transfer is performed.
            
        Return
        -------
        Spatial_data : scvelo Anndata Object
            Enriched Spatial data, having spliced and unspliced predicted gene
            expression, and predicted cellular annotations for metadata 
            identified (only if metadata_to_transfer is not None).
    """
    
    if metadata_to_transfer is SIRV.__defaults__[0]:
        Annotate_cells = False
    else:
        Annotate_cells = True
        
    RNA_data_scaled = pd.DataFrame(data=st.zscore(RNA_data.to_df(),axis=0),
                            index = RNA_data.to_df().index,columns=RNA_data.to_df().columns)
    Spatial_data_scaled = pd.DataFrame(data=st.zscore(Spatial_data.to_df(),axis=0),
                            index = Spatial_data.to_df().index,columns=Spatial_data.to_df().columns)
    Common_data = RNA_data_scaled[np.intersect1d(RNA_data_scaled.columns,Spatial_data_scaled.columns)]
    
    RNA_data_spliced = pd.DataFrame(RNA_data.layers['spliced'],columns=RNA_data.to_df().columns)
    RNA_data_unspliced = pd.DataFrame(RNA_data.layers['unspliced'],columns=RNA_data.to_df().columns)
    
    genes_to_predict = np.intersect1d(RNA_data.to_df().columns,Spatial_data.to_df().columns)
    
    Imp_Genes_spliced = pd.DataFrame(np.zeros((Spatial_data.shape[0],len(genes_to_predict))),
                                 index=Spatial_data_scaled.index,columns=genes_to_predict)
    Imp_Genes_unspliced = pd.DataFrame(np.zeros((Spatial_data.shape[0],len(genes_to_predict))),
                                 index=Spatial_data_scaled.index,columns=genes_to_predict)
    
    pv_Spatial_RNA = PVComputation(
            n_factors = n_pv,
            n_pv = n_pv,
            dim_reduction = 'pca',
            dim_reduction_target = 'pca'
    )
    
    pv_Spatial_RNA.fit(Common_data,Spatial_data_scaled[Common_data.columns])
    
    S = pv_Spatial_RNA.source_components_.T
        
    Effective_n_pv = sum(np.diag(pv_Spatial_RNA.cosine_similarity_matrix_) > 0.3)
    S = S[:,0:Effective_n_pv]
    
    Common_data_projected = Common_data.dot(S)
    Spatial_data_projected = Spatial_data_scaled[Common_data.columns].dot(S)
        
    nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',
                            metric = 'cosine').fit(Common_data_projected)
    distances, indices = nbrs.kneighbors(Spatial_data_projected)
    
    weights = np.zeros((Spatial_data.shape[0],50))
    
    for j in range(0,Spatial_data.shape[0]):
    
        weights[j,:] = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
        weights[j,:] = weights[j,:]/(len(weights[j,:])-1)
        Imp_Genes_spliced.iloc[j,:] = np.dot(weights[j,:],RNA_data_spliced[genes_to_predict].iloc[indices[j,:][distances[j,:] < 1]])
        Imp_Genes_unspliced.iloc[j,:] = np.dot(weights[j,:],RNA_data_unspliced[genes_to_predict].iloc[indices[j,:][distances[j,:] < 1]])
        
    dic = {}
    dic['spliced'] = Imp_Genes_spliced
    dic['unspliced'] = Imp_Genes_unspliced
    
    if (Annotate_cells):
        for i in metadata_to_transfer:
            RNA_labels=RNA_data.obs[i]
            Pred_Labels_prob = pd.DataFrame(np.zeros((Spatial_data.shape[0],len(np.unique(RNA_labels)))),
                                 index=Spatial_data_scaled.index,columns=np.unique(RNA_labels))
            Pred_Labels = pd.Series(index=Spatial_data_scaled.index)
            
            for j in range(0,Spatial_data.shape[0]):
                for k in range(50):
                    Pred_Labels_prob[RNA_labels[indices[j,k]]][j] += weights[j,k]
    
                Pred_Labels.iloc[j] = Pred_Labels_prob.columns[np.argmax(Pred_Labels_prob.iloc[j,:])]
            Spatial_data.obs[i] = np.array(Pred_Labels).astype(np.str_)
            
            
    Spatial_data = scv.datasets.AnnData(X = Spatial_data.to_df()[genes_to_predict], 
                            obs = Spatial_data.obs, var = Spatial_data.var.loc[genes_to_predict,:], 
                             obsm = Spatial_data.obsm, layers = dic)       
        
    return Spatial_data