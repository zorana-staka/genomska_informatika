
# coding: utf-8

# In[10]:


import pysam
samfile = pysam.AlignmentFile("/sbgenomics/project-files/merged-tumor.bam", "rb")
numOfReads = 0
numOfUnmappedReads = 0
numOfReadsMappingQualityZero = 0
cumulativeMappingQuality = 0

firstRead = pysam.AlignedSegment
i = 0

for read in samfile.fetch():
    numOfReads += 1
    cumulativeMappingQuality += read.mapping_quality
    if i == 0: 
        firstRead = read
        i += 1
    if read.is_unmapped: 
        numOfUnmappedReads += 1
    if read.mapping_quality == 0:
        numOfReadsMappingQualityZero += 1

samfile.close()

avgMappingQuality = cumulativeMappingQuality / numOfReads
avgMappingQualityWithoutZeros = cumulativeMappingQuality / (numOfReads - numOfReadsMappingQualityZero)

print("Unmapped reads: " + str(numOfUnmappedReads))
print("Reads with mapping quality zero: " + str(numOfReadsMappingQualityZero))
print("Reads: " + str(numOfReads))
print("Average mapping quality: " + str(avgMappingQuality))
print("Average mapping quality without reads with zero quality mapping: " 
      + str(avgMappingQualityWithoutZeros))

# inspecting flag fields in the first read from AlignmentFile
print("Flag property in first read: " + str(firstRead.flag))
# the line above returns number 1187; if we go to 
# https://broadinstitute.github.io/picard/explain-flags.html and put value of SAM flag to 1187 we get:
# read paired (0x1)
# read mapped in proper pair (0x2)
# mate reverse strand (0x20)
# second in pair (0x80)
# read is PCR or optical duplicate (0x400)

print("For the first read we have: ")
print("Is paired: " + str(firstRead.is_paired))
print("Is proper pair: " + str(firstRead.is_proper_pair))
print("Mate is reverse: " + str(firstRead.mate_is_reverse))
print("Is read2: " + str(firstRead.is_read2))
print("Is duplicate: " + str(firstRead.is_duplicate))

# On the other hand is_unmapped will return false 
print("Is the first read unmapped: " + str(firstRead.is_unmapped))

# Some other fields in AlignedSegment
# AlignedSegment is a class acting as a handle to the samtools (interacting with high-throughput data). 
# It represents aligned read. Not sure how deep inspect should go and how many fileds should include - 
# there are just few for learning and getting big picture.
# For the complete list of properties and methods of AlignedSegment I've checked 
# https://readthedocs.org/projects/pysam/downloads/pdf/latest/
print("Mapping quality of the first read: " + str(firstRead.mapping_quality))
print("Query length: " + str(firstRead.query_length))
print("Query name: " + str(firstRead.query_name))

