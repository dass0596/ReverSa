'''
Created on 31 de jul. de 2016

@author: daniela
'''
#!/usr/bin/env python

import sys
import os
import pandas as pd

from preprocessing import Preprocessing
from similarity import Similarity
from composition import Composition
from pca import Reduction
from clustering import Clustering
from cluster import ClusterReport
from postprocessing import Validate
from jplace import ParseJplace
from stats import Profiles
from binomial import StatsBinom
from clean import cleaning



def basicMode(config,fasta_file, profilePath):
    
    #create output folders
    wDir = os.getcwd()
    folders = ['pplacer','testing']
    for folder in folders:
        os.mkdir(os.path.join(wDir,folder))
        
    
    #Instance Preprocessing class
    window = Preprocessing(fasta_file, config['win_length'], config['win_step'])
    window.output_window()
    print >> sys.stderr, "Creating windows_sequence.fasta"
    
    #Instance Similarity and Composition class
    sim = Similarity(fasta_file, config['score_adj'],wDir)
    sim_matrix = sim.mcl_perform() 
    comp_results = Composition(config['kmer_len'])
    comp_matrix = comp_results.joined()
    #Join similarity and composition matrix for PCA
    join = pd.concat([comp_matrix, sim_matrix], axis= 1, join='inner')
    print >> sys.stderr, "Calculating similarity and composition matrix"
    
    #Instance Reduction class
    pca = Reduction(join, config['pca_comp'])
    pca_data = pca.perform_pca()
    print >> sys.stderr, "Performing PCA"
    
    #Instance Clustering class
    cluster = Clustering(pca_data)
    clust_obj = cluster.plot()
    print >> sys.stderr, "Performing clustering plot"
    
    #Instance ClusterReport class
    report = ClusterReport(clust_obj)
    file_name, querySeq = report.output_queryseq()
    print >> sys.stderr, "Doing report of clusters"

    #Instance Validate class
    valid = Validate(file_name, fasta_file,wDir)
    jfileComp, jfileMinus = valid.roundTwo()
    print >> sys.stderr, "Validation of results"
    
    #Instance ParseJplace Class
    parsing = ParseJplace(jfileComp, jfileMinus)
    corrMat = parsing.correlation()
    print >> sys.stderr, "Doing profiles"
    
    #Instance Profile Class
    ttest = Profiles(corrMat, querySeq, wDir, profilePath)
    bestWin = ttest.windowsAssigment()
    print >>sys.stderr, "Doing permutations"
    
    #Instance StatsBinom
    finalResult = StatsBinom(fasta_file, config['win_length'],bestWin)
    finalResult.binomial()
    
    cleaning(file_name)
    



            