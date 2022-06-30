




mkdir		mm39														
rsem-prepare-reference 	     --bowtie2  	-p 8        ../../../14_Genomes/Gencode/Mouse/gencode.vM29.transcripts.fa	   mm39/mm39	      	>> mm39.buildIndex.runLog	2>&1            


mkdir		hg38														
rsem-prepare-reference 	     --bowtie2  	-p 8         ../../../14_Genomes/Gencode/Human/gencode.v40.transcripts.fa      hg38/hg38	  	    >> hg38.buildIndex.runLog	2>&1                                                                              





mkdir		danRer11														
rsem-prepare-reference 	     --bowtie2  	-p 8         ../../../14_Genomes/UCSC/Zebrafish/mrna.fa	   danRer11/danRer11	       >> danRer11.buildIndex.runLog	2>&1            




mkdir		micMur2														
rsem-prepare-reference 	     --bowtie2  	-p 8         ../../../14_Genomes/UCSC/Mouse_lemur/mrna.fa	   micMur2/micMur2	       >> micMur2.buildIndex.runLog	2>&1            




mkdir		chlSab2														
rsem-prepare-reference 	     --bowtie2  	-p 8         ../../../14_Genomes/UCSC/Green_monkey/mrna.fa	   chlSab2/chlSab2	       >> chlSab2.buildIndex.runLog	2>&1            





