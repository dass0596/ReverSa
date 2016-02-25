from Bio import SeqIO

def window_generator(seq, win_length, win_step):
	for i in xrange(0, len(seq) - win_length + 1, win_step):
		yield "_" + str(i+1)  + ":" + str(i+win_length) + "\n" + seq[i:i+win_length]

out = open("seq_windows.fasta", "w")
for fasta in SeqIO.parse("Control.fasta", "fasta"):
	sequence , name = str(fasta.seq), fasta.id
	for subseq in window_generator(sequence, 240, 25):
		out.write(">" + name + subseq + "\n")

out.close()