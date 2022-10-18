# 3PrimSeq

R codes to detect the 3-prime end of mitochondrial RNA sequences.

This code expect a file table outputted from Blastn (-outfmt: 6), which aligned preprocessed RNAseq to mitochondrial reference sequence with J01415 accession number.

#blastn -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qframe qcovhsp qseq" -db J01415 -query RNAseq.fasta -out Output.blt -max_target_seqs 1

Example output:

#qseqid	qlen	sseqid	slen	qstart	qend	sstart	send	evalue	bitscore	score	length	pident	nident	mismatch	positive	gapopen	gaps	ppos	qframe	qcovhsp	qseq
#NS500198:409:HKTHLAFXY:1:11101:3613:16183	56	ND4_loc_10760_12137	52	1	52	1	52	5.76e-25	97.1	52	52	100.00052	0	52	0	0	100.00	1	93	CCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTT
