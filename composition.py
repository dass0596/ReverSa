'''
Created on 10 de mar. de 2016

@author: Daniela Sanchez
'''
from Bio.SeqUtils import GC
from itertools import product, tee, izip

from collections import OrderedDict
from Bio import SeqIO

import numpy as np
import pandas as pd
        
class Composition:
    '''
    Calculate composition features(gc content, kmers) of the sequences  
    '''
    def __init__(self, kmer_len):
        self.kmer_len = kmer_len
        self.gc_data = None 
        self.composition_df = None
        self.join = None
        
    def openFastaFile(self):
        self.fastaFile = open("windows_sequence.fasta" , "r")
        
    def gc_content(self):
        #Calculates G+C content for each window
        gc_d = OrderedDict()
        self.openFastaFile()
        for seq in SeqIO.parse(self.fastaFile, "fasta"):
            cal = GC(str(seq.seq))
            gc_d[seq.id] = cal / 100
        gc = pd.DataFrame.from_dict(gc_d, orient='index', dtype=float)
        gc.columns = ['GC']
        self.gc_data = np.log(gc['GC'])
        self.fastaFile.close()
        
    
    def generate_feature_mapping(self):
        #Generate kmer dictionary
        BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
        kmer_hash = {}
        counter = 0
        for kmer in product("ATGC",repeat=self.kmer_len):
            if kmer not in kmer_hash:
                kmer_hash[kmer] = counter
                rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
                kmer_hash[rev_compl] = counter
                counter += 1
        return kmer_hash, counter    
        
    def calculate_kmer(self):
        kmer_hash , counter = self.generate_feature_mapping()
        composition_d = OrderedDict()
        self.openFastaFile()
        for seq in SeqIO.parse(self.fastaFile, "fasta"):
            kmers = [
                kmer_hash[kmer_tuple]
                for kmer_tuple 
                in self.window(str(seq.seq))
                if kmer_tuple in kmer_hash
                ]
            kmers.append(counter - 1)
            composition_v = np.bincount(np.array(kmers))
            composition_v[-1] -= 1
            composition_d[seq.id] = composition_v + np.ones(counter)
        composition_p = pd.DataFrame.from_dict(composition_d, orient='index', dtype=float)
        #Normalize kmer frequencies, log(p_ij) = log[(X_ij +1) / rowSum(X_ij+1)]
        self.composition_df = np.log(composition_p.divide(composition_p.sum(axis=1),axis=0))
        self.fastaFile.close()
                 
    def window(self, str_seq):
        els = tee(str_seq, self.kmer_len)
        for i,el in enumerate(els):
            for _ in xrange(i):
                next(el, None)
        return izip(*els)
        
    def joined(self):
    #Join both clustering(GC and kmers) 
        self.gc_content()
        self.calculate_kmer()
        self.join = pd.concat([self.composition_df, self.gc_data], axis= 1, join='inner')
        return self.composition_df
        