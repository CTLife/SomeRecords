


mkdir		SARS-CoV-2														
snap-aligner index     ../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa	   SARS-CoV-2	 >> SARS-CoV-2.buildIndex.runLog	2>&1            


mkdir		mm39														
snap-aligner index    ../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa	   mm39	     -sm	>> mm39.buildIndex.runLog	2>&1            


mkdir		hg38														
snap-aligner index     ../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa      hg38	 -sm	>> hg38.buildIndex.runLog	2>&1                                                                              





mkdir		danRer11														
snap-aligner index     ../../../14_Genomes/UCSC/Zebrafish/danRer11.fa	   danRer11	  -sm    >> danRer11.buildIndex.runLog	2>&1            




mkdir		micMur2														
snap-aligner index     ../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa	   micMur2	  -sm    >> micMur2.buildIndex.runLog	2>&1            




mkdir		chlSab2														
snap-aligner index     ../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa	   chlSab2	  -sm    >> chlSab2.buildIndex.runLog	2>&1            




