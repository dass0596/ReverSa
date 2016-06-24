'''
Created on 14 de mar. de 2016

@author: Daniela Sanchez
'''
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

class Clustering():
    def __init__(self, pcaData):
        self.pcaData = pcaData
        self.r_correlation = None
        self.clusters = None
        self.silinfo = None

    def correlation_data(self):
        #Calculate correlation data from the PCA transform data
        data_transpose = self.pcaData.transpose()
        correlation = data_transpose.corr(method='pearson')
        correlation.to_csv('pcaDat.csv')
        self.r_correlation = pandas2ri.py2ri(correlation)
        
    def calculate_k(self):
        self.correlation_data()
        #Calculate k as number of cluster needed
        base = importr('base')
        #Fetch @ form the instance's dictionary of attributes
        sign = base.__dict__['@']
        apcluster = importr('apcluster')
        num_clust = apcluster.apcluster(apcluster.negDistMat(r=2, method='minkowski'), self.r_correlation)
        self.clusters = len(sign(num_clust, 'clusters'))
        
    def plot(self):
        #cluster our windows using PAM(Partitioning Around Medoids), save the plot as pdf and info as txt
        self.calculate_k()
        cluster = importr('cluster')
        cluspam = cluster.pam(self.r_correlation, self.clusters)
        self.silinfo = cluspam[6]
        grdevices = importr('grDevices')
        grdevices.pdf('cluster_plot.pdf')
        cluster.clusplot(self.r_correlation, cluspam[2], color= 'TRUE', 
                         shade='TRUE', main="Cluster plot", labels=2, lines=0, cex=0.3)
        grdevices.dev_off()
        with open("cluster_info.txt", "w") as text_file:
            text_file.write(str(cluspam[6]))
        return self.silinfo
        