'''
Created on 15 de mar. de 2016

@author: Daniela Sanchez
'''
import ConfigParser


def configFile(parameter):
    Config = ConfigParser.ConfigParser()
    
    cfgfile = open("args_mode.ini",'w')
    
    # add the settings to the structure of the file
    Config.add_section('Preprocessing')
    Config.set('Preprocessing','window_length',parameter['win_length'])
    Config.set('Preprocessing','window_step',parameter['win_step'])
    Config.add_section('Similarity')
    Config.set('Similarity','score_adjacency',parameter['score_adj'])
    Config.add_section('Composition')
    Config.set('Composition','kmer_length', parameter['kmer_len'])
    Config.add_section('PCA')
    Config.set('PCA','variance_components', parameter['pca_comp'])
    Config.write(cfgfile)
    cfgfile.close()

