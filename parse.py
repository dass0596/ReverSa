'''
Created on 31 de jul. de 2016

@author: Daniela Sanchez
'''
from argparse import ArgumentParser
import ConfigParser
import ast

def arguments():
    parser = ArgumentParser()
    config = ConfigParser.ConfigParser()
    
    #Input files
    parser.add_argument('--mode')
    parser.add_argument('--fasta_file')
    parser.add_argument('--compare')
    parser.add_argument('--config_file')
    args = parser.parse_args()
    #Read config_file
    config.read(args.config_file)
    #Create a dictionary to store values of the config_file
    config_d = {}
    config_d['win_length'] = ast.literal_eval(config.get('Preprocessing','window_length'))
    config_d['win_step'] = ast.literal_eval(config.get('Preprocessing','window_step'))
    config_d['score_adj'] = ast.literal_eval(config.get('Similarity','score_adjacency'))
    config_d['kmer_len'] = ast.literal_eval(config.get('Composition','kmer_length'))
    config_d['pca_comp'] = ast.literal_eval(config.get('PCA','variance_components'))
    return args,config_d
    