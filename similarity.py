'''
Created on 10 de mar. de 2016

@author: Daniela Sanchez
'''
import dendropy
import pandas as pd
import numpy as np
import re

from Bio.Emboss.Applications import WaterCommandline
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import RaxmlCommandline

from collections import OrderedDict
from mcl.mcl_clustering import mcl
from itertools import groupby


class Similarity():
    
    def __init__(self, filename, score_adj):
        self.score = score_adj
        self.fasta_seq = filename
        self.similitud_m = None
        self.norm_matrix = None 
    
    def water_alignment(self):
        #perform water alignment for each sequence versus all, store scores into a matrix
        similitud_d = OrderedDict()
        pattern= re.compile(r'Score: (.*)')
        for seq in SeqIO.parse('windows_sequence.fasta', 'fasta'):
            cline = WaterCommandline(asequence= 'asis:%s' %str(seq.seq), bsequence='windows_sequence.fasta', gapopen=10, gapextend=0.5, outfile='water.txt')
            cline()
            scores = []
            for line in open('water.txt'):
                scores += pattern.findall(line)
                similitud_d[seq.id] = scores               
        similitud = pd.DataFrame.from_dict(similitud_d, orient='index', dtype=float)
        #set the same names for index and columns
        similitud.columns = similitud.index.tolist()
        #Normalize the scores 1-(score/score max of alignment)
        for i in  similitud.index.tolist():
            for j in similitud.columns.tolist():
                if i == j:
                    max_score = similitud.loc[i,j]
                similitud.loc[i,j] = 1-(similitud.loc[i,j]/max_score)
        #similitud.to_csv('data_sim.csv')
        self.similitud = similitud
        
    def generate_dist(self):
        mafft_cline = MafftCommandline(input=self.fasta_seq, maxiterate = 1000, localpair = True, phylipout=True)
        stdout, stderr = mafft_cline()
        #Save alignments into  FASTA and PHYLIP format
        phyFile = 'testing/alignment.phy'
        outPhy = open( phyFile, 'w')
        outPhy.write(stdout)
        outPhy.close()
        fastaFile = 'testing/align.fasta'
        SeqIO.convert(phyFile, 'phylip', fastaFile, 'fasta')
        #Create phylogenetic tree of the original sequences
        raxml_cline = RaxmlCommandline(sequences=phyFile, model='GTRGAMMA', name='reversatest', working_dir='~/workspace/reversa/testing')
        raxml_cline()
        #Calculate the phylo distances between each branch of the tree
        tree = dendropy.Tree.get_from_path("testing/RAxML_result.reversatest", "newick")
        pdm = tree.phylogenetic_distance_matrix()
        pdm.write_csv('distance.csv')
        
            
    def incorporate_distance(self):
        #Add the phylo distances to the alignment matrix
        self.water_alignment()
        self.generate_dist()
        pattern = re.compile(r'_')
        df = pd.read_csv('distance.csv', index_col = 0)
        for i in self.similitud.index.tolist():
            left_row = str(pattern.split(i)[0])
            for j in self.similitud.columns.tolist():
                left_column = str(pattern.split(j)[0])
                #normalize the distance value
                dist_score =  1-df.loc[left_row,left_column]
                self.similitud.loc[i,j] = self.similitud.loc[i,j] * dist_score
        #self.similitud.to_csv('sim_dist.csv')         
            
    def calculate_adjacency(self):
        #calculate adjacency depend on the sequence(index!=column) and the score(K>Value[i,j])
        self.incorporate_distance()
        pattern = re.compile(r'_')
        k = self.score
        for i in self.similitud.index.tolist():
            left_row = pattern.split(i)[0]
            for j in self.similitud.columns.tolist():
                valor = 0
                left_column = pattern.split(j)[0]
                if left_row != left_column:
                    if k > self.similitud.loc[i,j]:
                        valor = 1
                self.similitud.loc[i,j]= valor
        
    def mcl_perform(self):
        #Cluster all values(1) using Markov Cluster Algorithm 
        self.calculate_adjacency()
        #self.similitud = pd.read_csv('dist_adjmatrix.csv', index_col = 0)
        count = 0 
        d_values = self.similitud.values
        M, clusters = mcl(d_values)
        valid_dict = {}
        #Remove duplicate clusters
        for i in clusters:
            if len(clusters[i]) > 1:
                valid_dict[i] = clusters[i]
        keys = sorted(valid_dict)
        dict1 = dict(zip(keys, (x for x, y in groupby(valid_dict[k] for k in keys))))   
        #Replace value[i,j] with the number of cluster that belongs             
        for h in dict1.values():
            count += 1
            for o in h:
                for m in h:
                    if o != m:
                        self.similitud.iloc[o,m] = count
        #self.similitud.to_csv('adjacency_modif.csv')
        adjacency = np.sum(self.similitud,axis=1)
        #Replace 0 to avoid inf values in the normalization
        ad = adjacency.replace(0,0.00001)
        #Normalize
        matrix_ad = np.log(ad) 
        self.norm_matrix = pd.DataFrame(matrix_ad, columns=['adj'])
        #self.norm_matrix.to_csv('adj_col.csv')
        return self.norm_matrix
