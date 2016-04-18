'''
Created on 15 de mar. de 2016

@author: daniela
'''
import ConfigParser

Config = ConfigParser.ConfigParser()

cfgfile = open("args.ini",'w')

# add the settings to the structure of the file, and lets write it out...
Config.add_section('Preprocessing')
Config.set('Preprocessing','window_length',240)
Config.set('Preprocessing','window_step',240)
Config.add_section('Similarity')
Config.set('Similarity','score_adjacency',0.1)
Config.add_section('Composition')
Config.set('Composition','kmer_length', 4)
Config.add_section('PCA')
Config.set('PCA','variance_components', 0.9)
Config.write(cfgfile)
cfgfile.close()

