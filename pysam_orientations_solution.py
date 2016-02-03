import pysam as ps

samfile = ps.AlignmentFile("hESC_chr22.bam", "rb", check_sq=False)
fr=0	#convergent count
rf=0	#divergent count
ff=0	#parallel count
lineCount=0
currPercentage=0	#to track progress
for read in samfile:
	if read.next_reference_start>read.reference_start:
		mateAhead=True
	else:
		mateAhead=False
	if (not read.is_reverse and read.mate_is_reverse and mateAhead) or (read.is_reverse and not read.mate_is_reverse and not mateAhead):	
		fr+=1
	elif (read.is_reverse and not read.mate_is_reverse and mateAhead) or (not read.is_reverse and read.mate_is_reverse and not mateAhead):	
		rf+=1
	else:
		ff+=1
	lineCount+=1
	newPercentage = currPercentage+1
	if float(lineCount)/8672163 >= float(newPercentage)/100: #if at least X% of the file has been read, print percentage	
		currPercentage = newPercentage
		print str(currPercentage) + "% complete"
samfile.close()

print "{0:.1f}% convergent".format(float(fr)/lineCount*100)
print "{0:.1f}% divergent".format(float(rf)/lineCount*100)
print "{0:.1f}% parallel".format(float(ff)/lineCount*100)
