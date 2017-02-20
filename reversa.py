'''
Created on 31 de jul. de 2016

@author: Daniela Sanchez
'''
#!/usr/bin/env python

from check import Check
from principal import basicMode
from parse import arguments
from config_file import configFile
from compare import compareResults, baseReport
from plot import graph

from collections import OrderedDict
import os, datetime
import sys

def main(args, config):
    
    originalPath= os.getcwd()
    fasta_file = os.path.join(os.getcwd(),args.fasta_file)
    profilePath = os.path.join(os.getcwd(),'profiles')
    wDir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    os.mkdir(wDir)      
    os.chdir(wDir)
    
    #add a check module
    fileWork = os.path.join(os.getcwd(),'final_sequences.fasta')
    validFile = Check(fasta_file, args, config)
    validFile.lenFasta()
    validFile.reverseFile()
    validFile.makeLog()
    
    
    if int(args.mode) == 1:
        basicMode(config, fasta_file, profilePath)
        
    if int(args.mode) == 2:   
        newConfig = config.copy()
        dfCompare = os.path.join(originalPath, args.compare)
        
        
        xList = []
        yList = []
        countResults = OrderedDict()
        for i,j in config.items():
            if type(j) is list:
                label = i
                valStart= j[0]
                valStop= j[1]
                valStep= j[2]
                
                k = valStart
                while k <= float(valStop):
                    newConfig[i] = k
                    folder = i + str(k)
                    subDir = os.path.join(wDir, folder)
                    os.mkdir(subDir)
                    os.chdir(subDir)
                    if i == 'win_length':
                        newConfig['win_step'] = k
                    
                    configFile(newConfig)
                    basicMode(newConfig, fasta_file, profilePath)
                    
                    result, trueEvents, times, baseEvents, seqDict = compareResults(dfCompare)
                    xList.append(k)
                    yList.append(result)
                    countResults[k] = times
    
                    baseReport(countResults, label, baseEvents, config['win_length'], seqDict)                       
                    
                    k = k + valStep     
                    os.chdir('..')
                        
                    
        graph(xList, yList, label, trueEvents)
        #print xList
        #print yList
        #print countResults
        baseReport(countResults, label, baseEvents, config['win_length'], seqDict)
            

if __name__ == '__main__':
    args, config = arguments()
    results = main(args, config)
    print >> sys.stderr, "ReverSa Finished"


