'''
Created on 10 de mar. de 2016

@author: Daniela Sanchez
'''
from Bio.Emboss.Applications import WaterCommandline
from Bio import SeqIO
from collections import OrderedDict

import re
import pandas as pd
import numpy as np

class Similarity():
    
    def __init__(self, score_adj):
        self.score = score_adj
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
        self.similitud = similitud
        
    def calculate_adjacency(self):
        #calculate adjacency depend on the sequence(index!=column) and the score(K>Value[i,j])
        self.water_alignment()
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
        adjacency = np.sum(self.similitud,axis=1)
        #replace 0 to avoid inf values in the normalization
        ad = adjacency.replace(0,0.00001)
        #normalize
        matrix_ad = np.log(ad) 
        self.norm_matrix = pd.DataFrame(matrix_ad, columns=['adj'])
        return self.norm_matrix
