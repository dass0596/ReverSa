'''
Created on 21 de may. de 2016

@author: Daniela Sanchez
'''
from scipy import stats
from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO
import re
import itertools


class StatsBinom():
    
    def __init__(self,fastaFile, size, bestseq):
        self.fastaFile = fastaFile
        self.sizeWin = size
        self.bestWin = bestseq
        self.seq = None
            
                
    def binomial(self):
        patt  = re.compile(r'_')
        pattern= re.compile('\#\sIdentity\:\s+[0-9\/]+\s+\(([0-9\.]+)\%\)')
        patt1= re.compile('\#\sIdentity\:\s+(\d+)\/(\d+)\s+\([0-9\.]+\%\)')
        patt2= re.compile('\#\sGaps\:\s+(\d+)\/\d+')
        out = open("reversa_result.txt", "w")
        
        for j, i in self.bestWin.items():
            out.write('cluster %s \n' %j)
            for a, b in itertools.combinations(i,2):
                seq1 = patt.split(a)[0]
                seq2 = patt.split(b)[0]
                if seq1 != seq2:
                    for fasta in SeqIO.parse(self.fastaFile, "fasta"):
                        if fasta.id == seq1:
                            aseq = str(fasta.seq)
                        if fasta.id == seq2:
                            bseq = str(fasta.seq)
                         
                    needle_cline = NeedleCommandline(asequence="asis:%s" %str(aseq), bsequence="asis:%s" %str(bseq),gapopen=10, gapextend=0.5, outfile="needle.txt")
                    needle_cline()
                    
                    
                    for line in open('needle.txt'):
                        ident= pattern.search(line)
                        if ident is not None:
                            identity = ident.group(1)
                            
                    for seq in SeqIO.parse("windows_sequence.fasta", "fasta"):
                        if seq.id == a:
                            awin = str(seq.seq)
                        if seq.id == b:
                            bwin = str(seq.seq)
                    
                    ncline = NeedleCommandline(asequence="asis:%s" %str(awin), bsequence="asis:%s" %str(bwin),gapopen=10, gapextend=0.5, outfile="needle1.txt")
                    ncline()
                    
                    for line1 in open('needle1.txt'):
                        #print line1
                        identwin= patt1.search(line1)
                        if identwin is not None:
                            iwin = identwin.group(1) 
                            ilen = identwin.group(2)
                        gapPatt= patt2.search(line1)
                        #print gapPatt
                        if gapPatt is not None:
                            gap = gapPatt.group(1)
                    
                    difWin = int(ilen) - (int(iwin) + int(gap))
                    size = self.sizeWin
                    totalDif= 1-(float(identity)/100)
                    pvalue = stats.binom.cdf(difWin, size, totalDif)
                    
                    if pvalue <= 1.00e-20:
                        #print 'combination:', a, b, 'pvalue: ', pvalue
                        out.write('combination: %s  %s pvalue: %s \n' %(a,b,pvalue))
        
        out.close()
                