'''
Created on 1 de abr. de 2016

@author: Daniela Sanchez
'''
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

import os
import re
import sys
import os.path


class Validate():
    '''
    Validate results by a reference tree and a multiple alignment
    '''
    
    def __init__(self, query_names, fasta, wDir):
        #"query names" list with  file names of query sequences 
        self.query_name = query_names
        self.fasta_seq = fasta
        self.cwPath = "%s/testing/" %wDir
        self.jfile = None
        self.jfileMinus = None
        
    def roundOne(self):
        pattern = re.compile(r'_')
        #Add query sequences to the previous alignment 
        alignFile = []
        self.jfile =[] 
        for i in self.query_name:
            clust = pattern.split(i)[2]
            
            file_name = 'multiple_ali%s.fasta' %clust
            os.system('mafft --add %s --quiet --reorder testing/align.fasta >testing/%s' %(i,file_name)) 
            alignFile = alignFile + [file_name,]
            jason_name = 'multiple_ali%s.jplace' %clust
            self.jfile.append(jason_name)
            for a in alignFile:
                #wrap pplacer
                os.system('pplacer --out-dir pplacer  -p -t testing/RAxML_result.reversatest -s testing/RAxML_info.reversatest testing/%s' %a)
        
    def roundTwo(self):
        self.roundOne()
        self.jfileMinus = []   
        pattern = re.compile(r'_')
        for j in self.query_name:
            clust_id = pattern.split(j)[2]
            for query in SeqIO.parse(j, 'fasta'):
                seq = pattern.split(query.id)[0]
                #Create special identifier for each round of files
                number = 'minus%s' %seq
                id_jfile = '%s_minus%s' %(clust_id,seq)
                
                rax_name = 'reversatest%s' % number
                fasta_name = 'testing/align%s.fasta' % number
                
                if not os.path.isfile('testing/alignment%s.phy' % number):
                    edited = MultipleSeqAlignment([])
                    openPhy = open('testing/alignment.phy')
                    record = AlignIO.read(openPhy, 'phylip')
                    for i in record:  
                        if i.id != seq:
                            edited.append(i)
                            
                    #write the alignment minus a sequence
                    phy_name = 'testing/alignment%s.phy' % number
                    out = open(phy_name , 'w')
                    AlignIO.write(edited, out, 'phylip')
                    out.close()
                    #convert FASTA to PHYLIP format 
                    SeqIO.convert(phy_name, 'phylip', fasta_name, 'fasta', )
    
                    #Create reference tree
                    raxml_line = RaxmlCommandline(sequences=phy_name, model='GTRGAMMA', name=rax_name, working_dir=self.cwPath)
                    raxml_line()
                
                #Add query sequences to the previous alignment
                multiali_name = 'testing/multiple_ali%s.fasta' %id_jfile
                
                if not os.path.isfile('testing/alignment%s.phy' %id_jfile):
                    os.system('mafft --add %s --quiet --reorder %s >%s'% (j, fasta_name, multiali_name))  
                    
                    jason_name = 'multiple_ali%s.jplace' %id_jfile

                    #wrap pplacer
                    if not os.path.isfile('pplacer/%s' %jason_name):
                        self.jfileMinus.append(jason_name)
                        os.system('pplacer --out-dir pplacer  -p -t testing/RAxML_result.%s -s testing/RAxML_info.%s %s' % (rax_name, rax_name, multiali_name))       
        print self.jfile
        print self.jfileMinus
       
        return self.jfile, self.jfileMinus