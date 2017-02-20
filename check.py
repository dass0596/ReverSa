'''
Created on 9 de nov. de 2016

@author: daniela
'''

from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline
from collections import OrderedDict

import re
import os
import sys
import pandas as pd

class Check():
    
    def __init__(self, fastaFile, args, config):
        self.file = fastaFile
        self.args = args
        self.configFile = config
        self.dfQuery = None
        self.dfRever = None
        self.dfFinal = None
        self.record_dict = None
         
    def lenFasta(self):
        #check if number of sequences is enough to run reversa
        self.record_dict = SeqIO.to_dict(SeqIO.parse(self.file, "fasta"))
        lenFasta = len(self.record_dict.keys()) 
        if lenFasta <= 4:
            print 'Not enough sequences to run ReverSa'
            sys.exit(0)
        #print lenFasta

    
    def reverseFile(self):
        handle = open('rever_seq.fasta', 'w')
        
        for seq in SeqIO.parse(self.file, "fasta"):
            #a= str(seq.seq)
            rever = str(seq.seq.reverse_complement())
            #print seq.seq.reverse_complement()
            handle.write('>%s\n%s\n' %(seq.id,rever))
            #print a.reverse_complement()
        handle.close()
    
    def makeLog(self):
        logFile = open('log.txt', 'w')
        logFile.write('Arguments:\n%s\nConfiguration file:\n' %self.args)
        for key, value in self.configFile.iteritems():
            logFile.write('%s %s\n' %(key,value))
        logFile.close()    
                
    def localAlign(self):    
        self.reverseFile()
        queryDict = OrderedDict()
        reverDict = OrderedDict()
        patt = re.compile(r'\# Score: (.*)')
        for seq in SeqIO.parse(self.file, "fasta"):
            #local alignment of the original fasta 
            cline = WaterCommandline(asequence= 'asis:%s' %str(seq.seq), bsequence=self.file, gapopen=10, gapextend=0.5, outfile='water.txt')
            cline()
            matchObj = []
            for l in open('water.txt', 'r'):
                findScore = patt.match(l) 
                if findScore is not None:
                    matchObj.append(str(findScore.group(1))) 
            #print matchObj        
            queryDict[str(seq.id)] = matchObj
            
            
            #local alignment with the reverse sequences of the original fasta
            clineRever = WaterCommandline(asequence= 'asis:%s' %str(seq.seq), bsequence='rever_seq.fasta', gapopen=10, gapextend=0.5, outfile='waterRever.txt')
            clineRever()
            scores = []
            for line in open('waterRever.txt', 'r'):
                reverseScore = patt.match(line) 
                if reverseScore is not None:
                    scores.append(str(reverseScore.group(1))) 
            reverDict[str(seq.id)] = scores
                    
        #save the results of the alignment in a data frame
        self.dfQuery = pd.DataFrame.from_dict(queryDict, orient='index', dtype=float) 
        self.dfQuery.columns = self.dfQuery.index.tolist()
        self.dfRever = pd.DataFrame.from_dict(reverDict, orient='index', dtype=float)
        self.dfRever.columns = self.dfRever.index.tolist()
        
        # print self.dfQuery
        #print self.dfRever
    
    def compare(self):
        self.localAlign()
        dictSeq = OrderedDict()
        for i in  self.dfQuery.index.tolist():
            finalScores = []
            for j in self.dfQuery.columns.tolist():
                reScore = self.dfRever.loc[i,j]
                fScore = self.dfQuery.loc[i,j]
                if reScore > fScore:
                    finalScores.append(1)
                else:
                    finalScores.append(0)
            dictSeq[i] = finalScores
        
        self.dfFinal = pd.DataFrame.from_dict(dictSeq, orient='index', dtype=float)
        self.dfFinal.columns = self.dfFinal.index.tolist()
        
    def parseResults(self):
        self.compare()
        lenDf = len(self.dfFinal.index)
        k = 0.60 * lenDf
        fastaDict= {}
        #####
        logFile = open('log.txt', 'w')
        logFile.write('Arguments:\n%s\nConfiguration file:\n' %self.args)
        for key, value in self.configFile.iteritems():
            logFile.write('%s %s\n' %(key,value))
        logFile.write('\nSense strand:\n')
        #######
        #Open reverse sequences file as dictionary
        reverse_seq = SeqIO.to_dict(SeqIO.parse('rever_seq.fasta', "fasta"))
        
        #Row sum
        self.dfFinal['sum'] = self.dfFinal.sum(axis=1)
        
        for i in self.dfFinal.index.tolist():
            value = self.dfFinal.loc[i,'sum']
                       
            if value >= k:
                fastaDict[i] = reverse_seq[i].seq
                logFile.write('reverse %s\n' %i)
                #print 'reverse %s' %i
            else:
                fastaDict[i] = self.record_dict[i].seq
                logFile.write('forward %s\n' %i)
                #print 'forward %s' %i

        #print self.dfFinal
        #print fastaDict
        logFile.close()
        
        outFile = open('final_sequences.fasta','w')
        for i, j in fastaDict.items():
            outFile.write('>%s\n%s\n' %(i,j))
        outFile.close()
        
        #remove files created
        removeFiles= ['water.txt', 'rever_seq.fasta', 'waterRever.txt']
        for j in removeFiles:
            os.remove(j)
            
               
        


    
