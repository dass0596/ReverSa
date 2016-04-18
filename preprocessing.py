'''
Created on 10 de mar. de 2016

@author: Daniela Sanchez
'''
from Bio import SeqIO

class Preprocessing():
    
    def __init__(self, filename, win_length, win_step):
        self.filename = filename
        self.win_length = win_length
        self.win_step = win_step
        self.iniseq = None
               
    def window_generator(self, seq):
        #generate sliding windows depend on size and step
        for i in xrange(0, len(seq)- self.win_length + 1, self.win_step):
            self.iniseq= i
            yield "_%s-%s\n%s" %(str(i+1),str(i+self.win_length),seq[i:i+self.win_length])
     
    def output_window(self): 
        #save the windows generate in a fasta
        k =  0.5
        out = open("windows_sequence.fasta", "w")
        for fasta in SeqIO.parse(self.filename, "fasta"):
            seq , name = str(fasta.seq), fasta.id
            lenseq = len(seq)
            #Generate id for the last sequence
            firstSeq = (lenseq-self.win_length)+1
            for subseq in self.window_generator(seq):
                out.write('>' + name + subseq + '\n')
            #save the last 'win_length' characters of the  string 
            if (self.win_length * k) < (lenseq - (self.iniseq + self.win_length)):
                out.write('>'+name+'_'+str(firstSeq)+'-'+str(lenseq)+'\n' + seq[-(self.win_length):]+'\n')      
        out.close()      
