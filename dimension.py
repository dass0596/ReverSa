'''
Created on 24 de feb. de 2016

@author: Daniela Sánchez
'''
from sklearn.decomposition import PCA

class dimensions():
    '''
    Reduce dimensions of the data by PCA or DISTATIS
    '''


    def __init__(self, joined_data, n_components ):
        self.data = joined_data
        self.components= n_components
        
    def PCA(self):
        pca_object = PCA(self.components).fit(self.data)
        return pca_object.transform(self.data), pca_object
        
        
        