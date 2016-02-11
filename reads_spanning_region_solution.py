import pysam

#open SAM/BAM file with alignments
samfile = pysam.AlignmentFile("hESC_chr22.bam", "rb")

def number_of_reads_spanning(chromosome,start,end):
    counter=0
    for read in samfile.fetch():
         if (read.is_unmapped==False):
              read_id=read.reference_id
              reference_name=samfile.getrname(read_id)

              if ((read.reference_start < start) and (read.reference_end > end) and (reference_name==chromosome)):
                   counter=counter+1
    return counter

#read content of file with coordinates
with open("coordinates.txt") as coordinates:
    for coordinate in coordinates:

        region_chromosome=str(coordinate.split(' ')[0])
        region_start=int(coordinate.split(' ')[1])
        region_end=int(coordinate.split(' ')[2])

        print ("chromosome: " + region_chromosome + " start: " + str(region_start) + " end: " + str(region_end))
        print ("number of reads spanning this region: " + str(number_of_reads_spanning(region_chromosome,region_start,region_end)))

samfile.close()