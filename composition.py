from Bio.SeqUtils import GC
from itertools import product, tee, izip

from collections import OrderedDict
from Bio import SeqIO

import numpy as np
import pandas as p

class Composition:
    def __init__(self, fasta_file, kmer_len = 4, space = 2):
        self.fasta = fasta_file
        self.kmer_len = kmer_len
        self.space = space
        
    def gc_content(self):
        GC_d = OrderedDict()
        for seq in SeqIO.parse(self.fasta, "fasta"):
            cal = GC(str(seq.seq))
            GC_d[seq.id] = cal / 100.0
        GC_data = p.DataFrame.from_dict(GC_d, orient='index', dtype=float)
        print GC_data
    
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
        for seq in SeqIO.parse(self.fasta, "fasta"):
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
        composition = p.DataFrame.from_dict(composition_d, orient='index', dtype=float)
        #Normalize kmer frequencies 
        #log(p_ij) = log[(X_ij +1) / rowSum(X_ij+1)]
        composition = np.log(composition.divide(composition.sum(axis=1),axis=0))
        return composition
            
    def window(self, str_seq):
        els = tee(str_seq, self.kmer_len)
        for i,el in enumerate(els):
            for _ in xrange(i):
                next(el, None)
        return izip(*els)

    def joined(self):
    #Join both data 
        comp = self.calculate_kmer()
        gc = self.gc_content()
        
        
#out = open("composition.csv", "w")
in_seq = open("zamora_verde.fasta" , "r")
Comp_results = Composition(in_seq, 4 , 2)
Comp_results.gc_content()
#Comp_results.calculate_kmer()
#Comp_results.joined()
#out.write()
#out.close()