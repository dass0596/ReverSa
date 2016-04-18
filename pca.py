'''
Created on 10 de mar. de 2016

@author: Daniela Sanchez
'''
from sklearn.decomposition import PCA
import pandas as pd

class Reduction():
    '''
    Reduce dimensions of the clustering by Principal Component Analysis   
    '''
    def __init__(self, join_matrix, n_components):
        self.joined = join_matrix
        self.components= n_components
        self.pca_object = None
    
    def perform_pca(self):
        #dimensional reduction using SVD of joined clustering
        self.pca_object = PCA(self.components).fit(self.joined)
        transform_df = pd.DataFrame(self.pca_object.transform(self.joined), index=self.joined.index)
        return transform_df
        