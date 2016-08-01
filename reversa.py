'''
Created on 31 de jul. de 2016

@author: daniela
'''
from principal import basicMode
from parse import arguments
from config_file import configFile
from test import compareResults
from plot import graph

import os, datetime
import sys

def main(args, config):
    
    fasta_file = os.path.join(os.getcwd(),args.fasta_file)
    profilePath = os.path.join(os.getcwd(),'profiles')
    wDir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    os.mkdir(wDir)      
    
    
    if int(args.mode) == 1:
        os.chdir(wDir)
        basicMode(config, fasta_file, profilePath)
        
    if int(args.mode) == 2:   
        newConfig = config.copy()
        dfCompare = os.path.join(os.getcwd(),args.compare)
        os.chdir(wDir)
        
        allResults = {}
        for i,j in config.items():
            if type(j) is list:
                label = i
                valStart= j[0]
                valStop= j[1]
                valStep= j[2]
                
                k = int(valStart)
                while k <= int(valStop):
                    newConfig[i] = k
                    folder = i + str(k)
                    subDir = os.path.join(wDir, folder)
                    os.mkdir(subDir)
                    os.chdir(subDir)
                    if i == 'win_length':
                        newConfig['win_step'] = k
                    
                    configFile(newConfig)
                    basicMode(newConfig, fasta_file, profilePath)
                    
                    result, trueEvents = compareResults(dfCompare)
                    allResults[k] = result
                        
                    k = k + int(valStep)     
                    os.chdir('..')
                        
                    
        graph(allResults, label, trueEvents)
            

if __name__ == '__main__':
    args, config = arguments()
    results = main(args, config)
    print >> sys.stderr, "ReverSa Finished"
