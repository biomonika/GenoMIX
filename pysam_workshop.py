#######With out root access.
###STEP 1
##Install your own version of python

wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
tar -xzvf virtualenv-1.10.1.tar.gz
cd virtualenv-1.10.1/
python virtualenv.py ~/myEnv


###STEP 2
##Install pysam

~/myEnv/bin/pip install pysam

###STEP 3
##Download the BAM file.
wget 'http://sites.psu.edu/biomonika/wp-content/uploads/sites/14384/2016/02/workshop.zip' -O workshop.zip 
unzip workshop.zip  -d destination_folder
###STEP 4
##start python and import pysam 
~/myEnv/bin/python

##in python console
import pysam


#read BAM file
bam="hESC_chr22_wg.bam"
bamfile = pysam.AlignmentFile(bam, "rb")
####IGNORE this message: [W::bam_hdr_read] EOF marker is absent. The input is probably truncated.

##parse required reads
for read in bamfile.fetch('chr22', 16958180,16958190):
	print read

	
bamfile = pysam.AlignmentFile(bam, "rb")
pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=bamfile)
for read in bamfile.fetch('chr22', 16958180,16958190):
	if read.is_paired:
		if not read.is_duplicate:
			pairedreads.write(read)

pairedreads.close()
bamfile.close()

##Create a header
header = { 'HD': {'VN': '1.0'},
			'SQ': [{'LN': 51304566, 'SN': 'chr22'},
				   {'LN': 48129895, 'SN': 'chr21'}] }


#create a new read				   
a = pysam.AlignedSegment()

#assign values to each attribute
a.query_name = "read_28833_29006_6945“
a.query_sequence="AGCTTAGCTA"
a.flag = 99
a.reference_id = 0
a.reference_start = 32
a.mapping_quality = 20
a.cigar = ((0,10), (2,1), (0,25))
a.next_reference_id = 0
a.next_reference_start=199
a.template_length=167
a.query_qualities = pysam.qualitystring_to_array("<<<<<<<AAA")
a.tags = (("NM", 1),("RG", "L1"))

		   
bamfile = pysam.AlignmentFile(bam, "rb")
for pileupcolumn in bamfile.pileup('chr22', 16958180,16958190):
    print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))

for pileupcolumn in bamfile.pileup('chr22', 16958160,16958170):
	print ("\nBases at position %s = " % (pileupcolumn.pos))
	for pileupread in pileupcolumn.pileups:
		if not pileupread.is_del and not pileupread.is_refskip:
		# query position is None if is_del or is_refskip is set.
			print ("%s" % (pileupread.alignment.query_sequence[pileupread.query_position]))
	
	
bamfile.close()

		   
				   