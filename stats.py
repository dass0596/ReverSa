'''
Created on 16 de may. de 2016

@author: Daniela Sanchez
'''
#import itertools
#import math
import numpy as np
import scipy.stats
from statsmodels.stats.multitest import multipletests
import ast
import os
import linecache

class Profiles():
    '''
    Calculate permutation and probability of profiles
    '''


    def __init__(self, interceptsMat, querySeq, wDir, profilePath):
        self.intercepts = interceptsMat
        self.seqDict= querySeq
        self.wDir = wDir
        self.profiles = profilePath
        self.posWindows = None
        self.bestpvalue = None
        self.bestProfile = None
        


    def windowsAssigment(self):
         
        os.chdir(self.profiles)
        bestWin = {}
        
        for b, c in self.seqDict.items():
            print b, c
            intercepts =  []
            windows = c
            for a in c:
                numInter = self.intercepts.loc[a,"intercept"]
                intercepts.append(numInter)
                
            #checkedProfiles = []
            rawpvalues = []
            leng = len(intercepts)
            
            if leng ==2 or leng==3:
                bestWin[int(b)] = c
                
            else:
                compressFile= "profile_%s.bz2" %leng
                os.system("bunzip2 -k %s" %compressFile)
                inFile = "profile_%s" %leng
                
        
                for line in open(inFile):
                    
                    l = ast.literal_eval(line)
                    indexes1 = [i for i, x in enumerate(l) if x == 1]
                    newList1 = [intercepts[x] for x in indexes1]
                    indexes0 = [i for i, x in enumerate(l) if x == 0]
                    newList0 = [intercepts[x] for x in indexes0]
                    pvalue = scipy.stats.ttest_ind(newList1, newList0)[1]
                    rawpvalues.append(pvalue)
                    
                    
    
                # adjust the raw pvalues for multiple hypothesis testing using Bonferroni correction for FDR
                adjpvalues = multipletests(rawpvalues, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)[1]
               
                # return the bestpvalue (corresponding to the best class assignment)  and the ID of the positive windows 
                minPvalues = adjpvalues.argmin() + 1
                getProfile = linecache.getline(inFile, minPvalues) 
                bestProfile = ast.literal_eval(getProfile)
                enumerate1 = [i for i, x in enumerate(bestProfile) if x == 1]
                mean1 = np.mean([intercepts[x] for x in enumerate1])
                enumerate0 = [i for i, x in enumerate(bestProfile) if x == 0]
                mean0 = np.mean([intercepts[x] for x in enumerate0])
                
                if mean0 < mean1:
                    posWindows0 = [windows[x] for x in enumerate0]
                    bestWin[int(b)] = posWindows0
                
                else:
                    posWindows1 = [windows[x] for x in enumerate1]
                    bestWin[int(b)] = posWindows1
                
                os.remove(inFile)
        os.chdir(self.wDir)
        return bestWin
    
