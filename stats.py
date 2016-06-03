'''
Created on 16 de may. de 2016

@author: Daniela Sanchez
'''
import itertools
import scipy.stats
from statsmodels.stats.multitest import multipletests
import math


class Profiles():
    '''
    Calculate permutation and probability of profiles
    '''


    def __init__(self, interceptsMat, querySeq):
        self.intercepts = interceptsMat
        self.seqDict= querySeq
        self.posWindows = None
        self.bestpvalue = None
        self.bestProfile = None
        


    def windowsAssigment(self):
        bestWin = {}
        for b, c in self.seqDict.items():
            intercepts =  []
            windows = c
            for a in c:
                numInter = self.intercepts.loc[a,"intercept"]
                intercepts.append(numInter)
                
            checkedProfiles = []
            rawpvalues = []
            leng = len(intercepts)
            iterProfiles = itertools.product([0,1], repeat=leng)
            for l in iterProfiles:
                
                if (l.count(1) >= 2 and l.count(1) <= math.ceil(leng)): # Check to avoid going over REDUNTANT permutations
                    indexes1 = [i for i, x in enumerate(l) if x == 1]
                    newList1 = [intercepts[x] for x in indexes1]
                    indexes0 = [i for i, x in enumerate(l) if x == 0]
                    newList0 = [intercepts[x] for x in indexes0]
                    pvalue = scipy.stats.ttest_ind(newList1, newList0)[1]
                    checkedProfiles.append(l)
                    rawpvalues.append(pvalue)

            # adjust the raw pvalues for multiple hypothesis testing using Bonferroni correction for FDR
            adjpvalues = multipletests(rawpvalues, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)[1]
            
            # return the bestpvalue (corresponding to the best class assignment)  and the ID of the positive windows 
            inposWindows = [i for i, x in enumerate(checkedProfiles[adjpvalues.argmin()]) if x == 1]
            self.posWindows = [windows[x] for x in inposWindows] 
            self.bestpvalue = min(adjpvalues)
            self.bestProfile = checkedProfiles[adjpvalues.argmin()]
            bestWin[b] = self.posWindows
            #print self.bestpvalue, self.bestProfile, self.posWindows
        return bestWin
    