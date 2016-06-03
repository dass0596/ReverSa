'''
Created on 24 de mar. de 2016

@author: Daniela Sanchez
'''
import sys
import re
import pandas.rpy.common as com

from rpy2.robjects.packages import importr

from Bio import SeqIO
from rpy2.robjects import pandas2ri
pandas2ri.activate()

class ClusterReport():
    '''
    Filter the false positives clusters.
    '''

    def __init__(self, clust_obj):
        self.clust_obj = clust_obj
        self.silinfo = None
        self.cl_width = None
        self.avg_wid = None
        self.clust_valid = None
        
    def clust_read(self):
        base = importr('base')
        #Fetch $ form the instance's dictionary of attributes
        dolar = base.__dict__['$']
        clust = dolar(self.clust_obj, 'widths')
        clus_width = dolar(self.clust_obj, 'clus.avg.widths')
        avg_width = dolar(self.clust_obj, 'avg.width')
        #Convert to pandas object
        self.cl_width = pandas2ri.ri2py(clus_width)
        self.avg_wid = pandas2ri.ri2py(avg_width)
        #pylist1 = pandas2ri.ri2py_dataframe(clus) CAMBIAR ALTERNATIVA
        pylist = com.convert_robj(clust)
        
        #Transform the first data object of cluster information
        data = pylist.reset_index()
        df = data.set_index('cluster')
        df.rename(columns = {'index':'win_id'}, inplace = True)
        dd = df.reset_index()
        #Create a list with cluster number and win_id, transform to dataframe 
        gb = dd.groupby(('cluster'))
        result = gb['win_id'].unique()
        self.silinfo = result.to_frame()
        
    
    def filter(self):  
        #Discard false positives depend on avg_width and multispecies of cluster.
        self.clust_read()  
        pattern = re.compile(r'_')
        #limit = 0.96
        self.clust_valid = []
        for i in self.silinfo.index.tolist():
            win_cluster = self.silinfo.loc[i,'win_id']
            id1 = pattern.split(win_cluster[0])[0]
            for j in win_cluster:
                id2 = pattern.split(j)[0]
                if id1 != id2:
                    #add to the tuple the clusters that pass the condition
                    self.clust_valid = self.clust_valid + [i,] 
                    #if (self.avg_wid/self.cl_width[i-1]) <= limit:
                        #self.clust_valid = self.clust_valid + [i,]
                    break  
        len_valid = len(self.clust_valid)
        print >> sys.stderr, 'Performed filter,resulted in %i cluster multispecies' %len_valid 
            
    def output_queryseq(self):
        #Select only the windows that pass filter and save them in a fastafile
        self.filter()    
        file_id = []
        seqClust = {}
        for a in self.clust_valid:
            valid_win = self.silinfo.loc[a,'win_id']
            seqClust[a] = valid_win
            #save for each cluster a fasta file with the sequences of the cluster
            handle = 'query_sequences_%i_.fasta' %a
            file_id = file_id + [handle,]
            output_file = open(handle, 'w')
            for b in valid_win:
                win_fasta = open("windows_sequence.fasta" , "r")
                for fasta in SeqIO.parse(win_fasta, "fasta"):
                    if b == fasta.id:
                        output_file.write('>' + fasta.id +'\n' + str(fasta.seq) + '\n')                 
            output_file.close()               
        return file_id, seqClust