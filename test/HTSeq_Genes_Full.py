import sys, os, re
import numpy, math
import HTSeq


###########################################
def  addBin(array1, binLen):
    binNum = math.ceil( float(len(array1)) / float(binLen) )   ## NOT math.floor
    binNum = int(binNum)
    array2 = numpy.zeros( binNum, dtype='float64' ) 
    for i in range(0, binNum):
        n1 = binLen*i
        n2 = binLen*(i+1) 
        if(n2 > len(array1)):
            n2 = len(array1)
        array2[i] = sum(array1[n1:n2]) / float(n2-n1) 
    return array2        
###########################################




#####################################################################################################################################################################
## example:       python  AllGenes_HTSeq_full.py      -inputFile  5-unique/day0_rep1.bam       -fragSize  200      -outputFile  7-HTSeq/1-ReadsDensity/full/day0_rep1

if len(sys.argv)!= 7:
    print ("The parameter number is wrong!!! (#1)")
    raise SystemExit (1)

if ( ( re.search("^\S+\.py$",    sys.argv[0]) ) and \
     ( re.search("^-inputFile$", sys.argv[1]) ) and \
     ( re.search("^\S+$",        sys.argv[2]) ) and \
     ( re.search("^-fragSize$",  sys.argv[3]) ) and \
     ( re.search("^\d+$",        sys.argv[4]) ) and \
     ( re.search("^-outputFile$",sys.argv[5]) ) and \
     ( re.search("^\S+$",        sys.argv[6]) ) ): 
    pass
else:
    print ("The parameters are wrong!!! (#2)")
    raise SystemExit(1)

inputFile  = str(sys.argv[2]);
fragSize   = int(sys.argv[4]);
outputFile = str(sys.argv[6]);
#####################################################################################################################################################################


readsNum = 0
bamFile  = HTSeq.BAM_Reader(inputFile )  
coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )   
for oneRead in bamFile: 
    if  (oneRead.aligned) and (oneRead.iv.start > fragSize) :
        readsNum = readsNum + 1
        oneRead.iv.length = fragSize
        coverage[ oneRead.iv ] += 1
fOut0 = open(outputFile+".log.txt",  'w')
fOut0.write( inputFile + ":\t" + str(readsNum) + "\n" )
fOut0.close()

halfWidth = 5000
normalizeFactor = float(10**7)/float(readsNum)
binLength = 20


 
output1 = outputFile + ".txStart"
output2 = outputFile + ".txEnd"      
output3 = outputFile + ".cdsStart"   
output4 = outputFile + ".cdsEnd"     
output5 = outputFile + ".5UTR"  
output6 = outputFile + ".3UTR"    
output7 = outputFile + ".Exons"      
output8 = outputFile + ".Introns"     
output9 = outputFile + ".txGeneBody"  
fOut1 = open(output1,  'w')
fOut2 = open(output2,  'w')
fOut3 = open(output3,  'w')
fOut4 = open(output4,  'w')
fOut5 = open(output5,  'w')
fOut6 = open(output6,  'w')
fOut7 = open(output7,  'w')
fOut8 = open(output8,  'w')
fOut9 = open(output9,  'w')


for line in open( "0-Regions/isoform_exp_diff_2.bed" ):  
    fields = line.split( "\t" )
    
    ## 1. txStart
    window1 = HTSeq.GenomicInterval( fields[2], int(fields[4]) - halfWidth, int(fields[4]) + halfWidth, "." )
    if (fields[3] == '-' ):
        window1 = HTSeq.GenomicInterval( fields[2], int(fields[5]) - halfWidth, int(fields[5]) + halfWidth, "." )
    wincvg1  = numpy.fromiter( coverage[window1], dtype='float64' )   
    profile1 = numpy.zeros( len(wincvg1), dtype='float64' ) 
    if fields[3] == "+":
        profile1 += wincvg1
    else:
        profile1 += wincvg1[::-1]
    fOut1.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(fields[4])  + "\t" + str(fields[5]) )
    profile1a = addBin(profile1, binLength) * normalizeFactor
    for i in range(0,len(profile1a)):
        fOut1.write("\t" + str(profile1a[i]))
    fOut1.write("\n")
    
    ## 2. txEnd
    window2 = HTSeq.GenomicInterval( fields[2], int(fields[5]) - halfWidth, int(fields[5]) + halfWidth, "." )
    if (fields[3] == '-' ):
        window2 = HTSeq.GenomicInterval( fields[2], int(fields[4]) - halfWidth, int(fields[4]) + halfWidth, "." )
    wincvg2 = numpy.fromiter( coverage[window2], dtype='float64' )
    profile2 = numpy.zeros( len(wincvg2), dtype='float64' ) 
    if fields[3] == "+":
        profile2 += wincvg2
    else:
        profile2 += wincvg2[::-1]
    fOut2.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(fields[4])  + "\t" + str(fields[5]) )
    profile2a = addBin(profile2, binLength) * normalizeFactor
    for i in range(0,len(profile2a)):
        fOut2.write("\t" + str(profile2a[i]))
    fOut2.write("\n")
    
    ## 3. cdsStart
    window3 = HTSeq.GenomicInterval( fields[2], int(fields[6]) - halfWidth, int(fields[6]) + halfWidth, "." )
    if (fields[3] == '-' ):
        window3 = HTSeq.GenomicInterval( fields[2], int(fields[7]) - halfWidth, int(fields[7]) + halfWidth, "." )
    wincvg3 = numpy.fromiter( coverage[window3], dtype='float64' )
    profile3 = numpy.zeros( len(wincvg3), dtype='float64' ) 
    if fields[3] == "+":
        profile3 += wincvg3
    else:
        profile3 += wincvg3[::-1]
    fOut3.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(fields[6])  + "\t" + str(fields[7]) )
    profile3a = addBin(profile3, binLength) * normalizeFactor
    for i in range(0,len(profile3a)):
        fOut3.write("\t" + str(profile3a[i]))
    fOut3.write("\n")
    
    ## 4. cdsEnd
    window4 = HTSeq.GenomicInterval( fields[2], int(fields[7]) - halfWidth, int(fields[7]) + halfWidth, "." )
    if (fields[3] == '-' ):
        window4 = HTSeq.GenomicInterval( fields[2], int(fields[6]) - halfWidth, int(fields[6]) + halfWidth, "." )
    wincvg4  = numpy.fromiter( coverage[window4], dtype='float64' )
    profile4 = numpy.zeros( len(wincvg4), dtype='float64' ) 
    if fields[3] == "+":
        profile4 += wincvg4
    else:
        profile4 += wincvg4[::-1]
    fOut4.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(fields[6])  + "\t" + str(fields[7]) )
    profile4a = addBin(profile4, binLength) * normalizeFactor
    for i in range(0,len(profile4a)):
        fOut4.write("\t" + str(profile4a[i]))
    fOut4.write("\n")

    ## 5. 5UTR
    window5 = HTSeq.GenomicInterval( fields[2], int(fields[4]) - halfWidth, int(fields[6]) + halfWidth, "." )
    if (fields[3] == '-' ):
        window5 = HTSeq.GenomicInterval( fields[2], int(fields[7]) - halfWidth, int(fields[5]) + halfWidth, "." )
    wincvg5  = numpy.fromiter( coverage[window5], dtype='float64' )   
    profile5 = numpy.zeros( len(wincvg5), dtype='float64' ) 
    if fields[3] == "+":
        profile5 += wincvg5
    else:
        profile5 += wincvg5[::-1]
    fOut5.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(fields[4])  + "\t" + str(fields[5]) )
    profile5a = addBin(profile5, binLength) * normalizeFactor
    for i in range(0,len(profile5a)):
        fOut5.write("\t" + str(profile5a[i]))
    fOut5.write("\n")
    
    ## 6. 3UTR
    window6 = HTSeq.GenomicInterval( fields[2], int(fields[7]) - halfWidth, int(fields[5]) + halfWidth, "." )
    if (fields[3] == '-' ):
        window6 = HTSeq.GenomicInterval( fields[2], int(fields[4]) - halfWidth, int(fields[6]) + halfWidth, "." )
    wincvg6  = numpy.fromiter( coverage[window6], dtype='float64' )   
    profile6 = numpy.zeros( len(wincvg6), dtype='float64' ) 
    if fields[3] == "+":
        profile6 += wincvg6
    else:
        profile6 += wincvg6[::-1]
    fOut6.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(fields[4])  + "\t" + str(fields[5]) )
    profile6a = addBin(profile6, binLength) * normalizeFactor
    for i in range(0,len(profile6a)):
        fOut6.write("\t" + str(profile6a[i]))
    fOut6.write("\n")

    ## 7. Exons
    for i in range(0, int(fields[8]) ): 
        exonStarts = fields[9].split( "," )
        exonEnds   = fields[10].split( "," )
        window7 = HTSeq.GenomicInterval( fields[2], int(exonStarts[i])-halfWidth, int(exonEnds[i])+halfWidth, "." )
        wincvg7 = numpy.fromiter( coverage[window7], dtype='float64' )
        profile7 = numpy.zeros( len(wincvg7), dtype='float64' ) 
        if fields[3] == "+":
            profile7 += wincvg7
        else:
            profile7 += wincvg7[::-1]
        fOut7.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(exonStarts[i])  + "\t" + str(exonEnds[i]) )
        profile7a = addBin(profile7, binLength) * normalizeFactor
        for i in range(0,len(profile7a)):
            fOut7.write("\t" + str(profile7a[i]))
        fOut7.write("\n")

    ## 8. Introns
    for i in range(0, int(fields[8])-1 ): 
        exonStarts = fields[9].split( "," )
        exonEnds   = fields[10].split( "," )
        window8 = HTSeq.GenomicInterval( fields[2],  int(exonEnds[i])-halfWidth, int(exonStarts[i+1])+halfWidth, "." )
        wincvg8 = numpy.fromiter( coverage[window8], dtype='float64' )
        profile8 = numpy.zeros( len(wincvg8), dtype='float64' ) 
        if fields[3] == "+":
            profile8 += wincvg8
        else:
            profile8 += wincvg8[::-1]
        fOut8.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(exonStarts[i])  + "\t" + str(exonEnds[i]) )
        profile8a = addBin(profile8, binLength) * normalizeFactor
        for i in range(0,len(profile8a)):
            fOut8.write("\t" + str(profile8a[i]))
        fOut8.write("\n")

    ## 9. txGeneBody
    window9 = HTSeq.GenomicInterval( fields[2], int(fields[4])-halfWidth, int(fields[5])+halfWidth, "." )
    wincvg9 = numpy.fromiter( coverage[window9], dtype='float64' )
    profile9 = numpy.zeros( len(wincvg9), dtype='float64' ) 
    if fields[3] == "+":
        profile9 += wincvg9
    else:
        profile9 += wincvg9[::-1]
    fOut9.write( str(fields[1]) + "\t" + str(fields[2])  + "\t" + str(fields[3])  + "\t" + str(fields[4])  + "\t" + str(fields[5]) )
    profile9a = addBin(profile9, binLength) * normalizeFactor
    for i in range(0,len(profile9a)):
        fOut9.write("\t" + str(profile9a[i]))
    fOut9.write("\n")


fOut1.close()
fOut2.close()
fOut3.close()
fOut4.close()
fOut5.close()
fOut6.close()
fOut7.close()
fOut8.close()
fOut9.close()










