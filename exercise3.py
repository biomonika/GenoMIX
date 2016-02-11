import pysam


#read BAM file
bam="hESC_chr22_wg.bam"
bamfile = pysam.AlignmentFile(bam, "rb")

for pileupcolumn in bamfile.pileup('chr22', 16958160,16958170):
	if pileupcolumn.pos == 16958162:
		print ("\nBases at position %s = " % (pileupcolumn.pos))
		for pileupread in pileupcolumn.pileups:
			if not pileupread.is_del and not pileupread.is_refskip:
			# query position is None if is_del or is_refskip is set.
				print ("%s" % (pileupread.alignment.query_sequence[pileupread.query_position]))
		
	
bamfile.close()
