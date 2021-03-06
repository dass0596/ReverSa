'''
Created on 31 de jul. de 2016

@author: Daniela Sanchez
'''
import pandas as pd
import re


def compareResults(compareFile):
    
    df = pd.read_csv(compareFile, index_col = 0)
    dfReversa = pd.read_csv('results.csv', index_col = 0)
    
    patt = re.compile('(.*)_(\d+)-(\d+)')
    
    count = []
    baseEvent = {}
    dicEvents = {}
    
    for a in df.index.tolist():
        times = 0
        donor = df.loc[a,'Seq2']
        recombinant = df.loc[a,'Seq1']
        beginFrag = int(df.loc[a,'Begin'])
        endFrag = int(df.loc[a,'End'])
        bases = endFrag - beginFrag
        baseEvent[a] = bases
    
        for i in  dfReversa.index.tolist():
            #name each option of the results file
            sequence1 = dfReversa.loc[i,'seq1']
            win= patt.match(sequence1)
            seq1 = win.group(1)
            begin1 = int(win.group(2))
            end1= int(win.group(3))
            
            sequence2 = dfReversa.loc[i,'seq2']
            win= patt.match(sequence2)
            seq2 = win.group(1)
            begin2 = int(win.group(2))
            end2= int(win.group(3))
            
            if (donor == seq1 or donor ==seq2) and (recombinant == seq1 or recombinant == seq2):  
                if (beginFrag <= begin1 <= endFrag or beginFrag <= end1 <= endFrag) \
                and (beginFrag <= begin2 <= endFrag  or beginFrag <= end2 <= endFrag):
                    times += 1
                    if a not in count:
                        count.append(a)
        dicEvents[a] = times    
        
                    #print seq1, begin1, end1, '\t', seq2, begin2, end2
                    

    return len(count), len(df.index), dicEvents, baseEvent

#do report of all results
def baseReport(i, paramChange, bpEvents, winLength):
    handle = open('report.txt', 'w')
    for a, b in i.iteritems():
        handle.write('With a %s of %s \n' %(paramChange,a))
        for x, y in b.iteritems():
            if y > 0:
                if str(paramChange) == 'win_length':
                    numbases = y * a
                    handle.write('In event %s, ReverSa found %s bases of %s \n' %(x, numbases, bpEvents[x]))
                else:
                    nBases = winLength * y
                    handle.write('In event %s, ReverSa found %s bases of %s \n' %(x, nBases, bpEvents[x]))
                    
    handle.close()

