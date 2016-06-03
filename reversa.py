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
from parse import arguments


def main(args, config):
    #Instance Preprocessing class
    window = Preprocessing(args.fasta_file, config['win_length'], config['win_step'])
    window.output_window()
    print >> sys.stderr, "Created windows_sequence.fasta"
    
    #Instance Similarity and Composition class
    sim = Similarity(args.fasta_file, config['score_adj'])
    sim_matrix = sim.mcl_perform() 
    comp_results = Composition(config['kmer_len'])
    comp_matrix = comp_results.joined()
    #Join similarity and composition matrix for PCA
    join = pd.concat([comp_matrix, sim_matrix], axis= 1, join='inner')
    print >> sys.stderr, "Calculated similarity and composition matrix"
    
    #Instance Reduction class
    pca = Reduction(join, config['pca_comp'])
    pca_data = pca.perform_pca()
    print >> sys.stderr, "Performed PCA"
    
    #Instance Clustering class
    cluster = Clustering(pca_data)
    clust_obj = cluster.plot()
    print >> sys.stderr, "Performed clustering plot"
    
    #Instance ClusterReport class
    report = ClusterReport(clust_obj)
    file_name, querySeq = report.output_queryseq()
    print >> sys.stderr, "Done report of clusters"

    #Instance Validate class
    valid = Validate(file_name, args.fasta_file)
    jfileComp, jfileMinus = valid.roundTwo()
    print >> sys.stderr, "Validation of results"
    
    #Instance ParseJplace Class
    parsing = ParseJplace(jfileComp, jfileMinus)
    corrMat = parsing.correlation()
    print >> sys.stderr, "Doing profiles"
    
    #Instance Profile Class
    ttest = Profiles(corrMat, querySeq)
    bestWin = ttest.windowsAssigment()
    print >>sys.stderr, "Doing permutations"
    
    #Instance StatsBinom
    finalResult = StatsBinom(args.fasta_file, config['win_length'],bestWin)
    finalResult.binomial()
    
    cleaning()


if __name__ == '__main__':
    args, config = arguments()
    results = main(args, config)
    print >> sys.stderr, "ReverSa Finished"
