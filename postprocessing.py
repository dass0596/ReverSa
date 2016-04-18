'''
Created on 1 de abr. de 2016

@author: Daniela Sanchez
'''
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import SeqIO

import os
import re


class Validate():
    '''
    Validate results by a reference tree and a multiple alignment
    '''
    
    def __init__(self, query_names, fasta):
        #dictionary with query sequences =
        self.query_name = query_names
        self.fasta_seq = fasta
        
    def roundOne(self):
        pattern = re.compile(r'_')
         
        mafft_cline = MafftCommandline(input=self.fasta_seq, maxiterate = 1000, localpair = True, phylipout=True)
        stdout, stderr = mafft_cline()
        #Save alignments into  FASTA and PHYLIP format
        phyFile = 'testing/alignment.phy'
        outPhy = open( phyFile, 'w')
        outPhy.write(stdout)
        outPhy.close()
        fastaFile = 'testing/align.fasta'
        SeqIO.convert(phyFile, 'phylip', fastaFile, 'fasta')
        #Create reference tree
        raxml_cline = RaxmlCommandline(sequences=phyFile, model='GTRGAMMA', name='reversatest', working_dir='~/workspace/reversa/testing')
        raxml_cline()
        #Add query sequences to the previous alignment 
        alignFile = []
        for i in self.query_name:
            clust = pattern.split(i)[2]
            print i , type(i)
            file_name = 'multiple_ali%s.fasta' %clust
            os.system('mafft --add %s --quiet --reorder testing/align.fasta >testing/%s' %(i,file_name)) 
            alignFile = alignFile + [file_name,]
            for a in alignFile:
                #wrap pplacer
                os.system('~/Escritorio/software/pplacer/./pplacer --out-dir pplacer  -p -t testing/RAxML_result.reversatest -s testing/RAxML_info.reversatest testing/%s' %a)
        
    def roundTwo(self):
        self.roundOne()
        pattern = re.compile(r'_')
        for j in self.query_name:
            clust_id = pattern.split(j)[2]
            for query in SeqIO.parse(j, 'fasta'):
                seq = pattern.split(query.id)[0]
                #Create special identifier for each round of files
                number = '%s_minus%s' %(clust_id,seq)
                print number
                openFasta = open(self.fasta_seq, 'r')
                record_dict = SeqIO.to_dict(SeqIO.parse(openFasta, 'fasta'))
                for i in record_dict.keys():
                    if i == seq:
                        del record_dict[i]
                        
                #write the principal fasta
                target = 'testing/win%s.fasta' % number 
                out = open(target, 'w')
                SeqIO.write(record_dict.values(), out, 'fasta')
                out.close()
            
                #Perform multiple alignment
                mafft_line = MafftCommandline(input=target, maxiterate = 1000, localpair = True, phylipout=True)
                stdout, stderr = mafft_line()
                #Save alignments into  FASTA and PHYLIP format
                phy_name = 'testing/alignment%s.phy' % number
                outphy = open( phy_name, 'w')
                outphy.write(stdout)
                outphy.close()
                fasta_name = 'testing/align%s.fasta' % number
                SeqIO.convert(phy_name, 'phylip', fasta_name, 'fasta')

                #Create reference tree
                rax_name = 'reversatest%s' % number
                raxml_line = RaxmlCommandline(sequences=phy_name, model='GTRGAMMA', name=rax_name, working_dir= '~/workspace/reversa/testing' )
                raxml_line()
                #add loop
                
                #Add query sequences to the previous alignment
                multiali_name = 'testing/multiple_ali%s.fasta' %number
                os.system('mafft --add %s --quiet --reorder %s >%s'% (j, fasta_name, multiali_name))  
                #wrap pplacer
                os.system('~/Escritorio/software/pplacer/./pplacer --out-dir pplacer  -p -t testing/RAxML_result.%s -s testing/RAxML_info.%s %s' % (rax_name, rax_name, multiali_name)) 
