##################################################################################################################
## All samples must be overlapped. (100% overlap)
## Suffixes of all self-defined global variables must be "_g".
## Example:  
## Rscript  MaXY_allPups_cluster.R     2-SplitXY/A-rmXY     1000A_MaXY_allPups_cluster_rmXY   >     1000A_MaXY_allPups_cluster_rmXY.runLog  2>&1     
         
args_g <- commandArgs(TRUE)
print("args: ")
print(args_g[1])   
print(args_g[2])     
print("#############")

inputDir_g = args_g[1];     ## the path of input files
outDir_g   = args_g[2];     ## the path of output files
# inputDir_g =  "2-SplitXY/D-onlyY"
# outDir_g   =  "1000A_MaXY_allPups_cluster_rmXY"
print(inputDir_g)   
print(outDir_g)

if( ! file.exists(inputDir_g) ) { print("##### Error-1 #####") }
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }

lowestCoverage_g = 3   
##################################################################################################################


 

##################################################################################################################
library(methylKit)
library(genomation)
library(ggplot2) 
library(ggfortify)
library(cluster)
library(lfda)
library(MASS)
library(factoextra)
library(magrittr)  
library(dplyr)  
library(rgl)
library(gdata)
library(ggrepel)
library(scatterplot3d)
library(car) 
library(plotly)
library(plot3D)
library(FactoMineR)
library(fpc)   
library(DSS) 
library(dendextend)
library(ggpubr)


Files_NC_g <- c(
paste( inputDir_g,   "PCOS_10_P10C-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_10_P10D-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_1_P1C-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_2_P2C-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_2_P2D-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_3_P3C-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_4_P4C-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_4_P4D-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_5_P5C-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_5_P5D-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_9_P9D-boy-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),  

paste( inputDir_g,   "PCOS_1_P1D-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_3_P3D-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_6_P6C-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_6_P6D-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_7_P7C-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_7_P7D-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_8_P8C-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_8_P8D-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "PCOS_9_P9C-girl-PCOS_Rep1.CpG.bismark.cov",  sep="/" )   
)

mySampleID_NC_g <- rep( x="PCOS", times=length(Files_NC_g) )
Tech_NC_g <- mySampleID_NC_g

for(i in c(1:length(mySampleID_NC_g)) ) {
  mySampleID_NC_g[i] = paste(mySampleID_NC_g[i], i, sep="_")
}

myTreatment_NC_g <- rep( x=0,  times=length(Files_NC_g) )

Sex_NC_g = rep( x="boy",  times=length(Files_NC_g) )  
for(i in c(1:length(Sex_NC_g)) ) {
  Sex_NC_g[i] = "boy"
  if(i>=12) { Sex_NC_g[i] = "girl" }
}



Files_IVF_fresh_g <- c(
  paste( inputDir_g,   "76_ART-E18-D-Boy_Rep2.bismark.cov",              sep="/" ),
  paste( inputDir_g,   "77_E69-ART-C_Rep3.bismark.cov",                  sep="/" ),
  paste( inputDir_g,   "78_ART-E101-D-Boy_Rep4.bismark.cov" ,            sep="/" ),
  paste( inputDir_g,   "79_E72-ART-D_Rep3.bismark.cov",                  sep="/" ), 
  paste( inputDir_g,   "70_E113C-boy-ART_Rep3.bismark.cov",              sep="/" ),
  paste( inputDir_g,   "70_E113D-boy-ART_Rep3.bismark.cov",              sep="/" ),
  paste( inputDir_g,   "71_ART-W58-C-Boy_Rep1.bismark.cov",              sep="/" ),
  paste( inputDir_g,   "71_ART-W58-D-Boy_Rep1.bismark.cov",              sep="/" ),
  paste( inputDir_g,   "72_ART-W779-C-Boy_Rep1.bismark.cov",             sep="/" ),
  paste( inputDir_g,   "72_W779D-ART-boy_Rep1.bismark.cov",              sep="/" ),
  paste( inputDir_g,   "28_E23C-boy-IVF-fresh_Rep1.bismark.cov",         sep="/" ),
  paste( inputDir_g,   "29_W76C-boy-IVF-fresh_Rep1.bismark.cov",         sep="/" ),
  paste( inputDir_g,   "30_Q22-W1452C-boy-IVF-fresh_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "chenwei_111_E18D-boy-IVFfresh_Rep9.CpG.bismark.cov",  sep="/" ),  #14
  paste( inputDir_g,   "chenwei_113_E101C-girl-IVFfresh_Rep9.CpG.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "chenwei_113_E101D-boy-IVFfresh_Rep9.CpG.bismark.cov",  sep="/" ),
  
  paste( inputDir_g,     "76_ART-E18-C-Girl_Rep4.bismark.cov",        sep="/"),
  paste( inputDir_g,     "77_E69-ART-D_Rep3.bismark.cov",             sep="/"),
  paste( inputDir_g,     "78_ART-E101-C-Girl_Rep2.bismark.cov" ,      sep="/"),
  paste( inputDir_g,     "79_E72-ART-C_Rep3.bismark.cov",             sep="/"),
  paste( inputDir_g,     "64_ART-BS18-C-Girl-merge_Rep1.bismark.cov", sep="/"),
  paste( inputDir_g,     "64_ART-BS18-D-Girl_Rep1.bismark.cov",       sep="/"),
  paste( inputDir_g,     "65_ART-BS29-C-Girl_Rep2.bismark.cov",       sep="/"),
  paste( inputDir_g,     "65_ART-BS29-D-Girl_Rep2.bismark.cov",       sep="/"),
  paste( inputDir_g,     "66_ART-E29-C-Girl_Rep2.bismark.cov",        sep="/"),
  paste( inputDir_g,     "66_ART-E29-D-Girl_Rep1.bismark.cov",        sep="/"),
  paste( inputDir_g,     "42_W53C-girl-IVF-fresh_Rep1.bismark.cov",   sep="/"),
  paste( inputDir_g,     "43_W81C-girl-IVF-fresh_Rep1.bismark.cov",   sep="/"),
  paste( inputDir_g,     "44_W1694-IVF-fresh-C_Rep1.bismark.cov",     sep="/"),  
  paste( inputDir_g,   "chenwei_11_BS18C-girl-IVFfresh_Rep9.CpG.bismark.cov",  sep="/" ),  #30
  paste( inputDir_g,   "chenwei_11_BS18D-girl-IVFfresh_Rep9.CpG.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "chenwei_12_BS29C-girl-IVFfresh_Rep9.CpG.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "chenwei_12_BS29D-girl-IVFfresh_Rep9.CpG.bismark.cov",  sep="/" )
)

mySampleID_IVF_fresh_g <- rep( x="IVF_fresh", times=length(Files_IVF_fresh_g) )
Tech_IVF_fresh_g <- mySampleID_IVF_fresh_g

for(i in c(1:length(mySampleID_IVF_fresh_g)) ) {
  mySampleID_IVF_fresh_g[i] = paste(mySampleID_IVF_fresh_g[i], i, sep="_")
}

myTreatment_IVF_fresh_g <- rep( x=1,  times=length(Files_IVF_fresh_g) )

Sex_IVF_fresh_g = rep( x="boy",      times=length(Files_IVF_fresh_g) ) 
for(i in c(1:length(Sex_IVF_fresh_g)) ) {
  Sex_IVF_fresh_g[i] = "boy"
  if(i>=17) { Sex_IVF_fresh_g[i] = "girl" }
}


Files_IVF_frozen_g <- c(
  paste( inputDir_g,   "16_W839C-IVF-frozen_Rep1.bismark.cov",         sep="/" ),
  paste( inputDir_g,   "17_W1524D-IVF-frozen_Rep1.bismark.cov",        sep="/" ),
  paste( inputDir_g,   "18_W1387C-IVF-frozen_Rep1.bismark.cov" ,       sep="/" ),
  paste( inputDir_g,   "19_W1458D-IVF-frozen_Rep1.bismark.cov" ,       sep="/" ), 
  paste( inputDir_g,   "9_W1365C-boy-IVF-frozen_Rep1.bismark.cov",     sep="/" ),
  paste( inputDir_g,   "9_Q5-W1365D-boy-IVF-frozen_Rep1.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "10_W1733C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "11_W1398C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "11_W1398D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "31_W1149C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "33_W104C-boy-IVF-frozen_Rep1.bismark.cov",     sep="/" ),
  
  paste( inputDir_g,     "16_W839D-IVF-frozen_Rep1.bismark.cov",        sep="/"),
  paste( inputDir_g,     "17_Q24-W1524C-IVF-frozen_Rep1.bismark.cov",   sep="/"),
  paste( inputDir_g,     "18_W1387D-IVF-frozen_Rep1.bismark.cov" ,      sep="/"), 
  paste( inputDir_g,     "19_W1458C-IVF-frozen_Rep1.bismark.cov",       sep="/"),
  paste( inputDir_g,     "1_A42C-girl-IVF-frozen_Rep1.bismark.cov" ,    sep="/"), 
  paste( inputDir_g,     "1_A42D-girl-IVF-frozen_Rep1.bismark.cov" ,    sep="/"),
  paste( inputDir_g,     "2_W1785C-girl-IVF-frozen_Rep1.bismark.cov" ,  sep="/"), 
  paste( inputDir_g,     "2_W1785D-girl-IVF-frozen_Rep1.bismark.cov" ,  sep="/"), 
  paste( inputDir_g,     "3_Q9-W811C-girl-IVF-frozen_Rep1.bismark.cov", sep="/"),
  paste( inputDir_g,     "3_W811D-girl-IVF-frozen_Rep1.bismark.cov" ,   sep="/"),
  paste( inputDir_g,     "4_W28C-girl-IVF-frozen_Rep1.bismark.cov" ,    sep="/"), 
  paste( inputDir_g,     "4_W28D-girl-IVF-frozen_Rep1.bismark.cov",     sep="/"), 
  paste( inputDir_g,     "45_E14C-girl-IVF-frozen_Rep1.bismark.cov",    sep="/"),
  paste( inputDir_g,     "46_W857C-girl-IVF-frozen_Rep1.bismark.cov",   sep="/"),
  paste( inputDir_g,     "47_E32C-girl-IVF-frozen_Rep1.bismark.cov",    sep="/")   
)

mySampleID_IVF_frozen_g <- rep( x="IVF_frozen", times=length(Files_IVF_frozen_g) )
Tech_IVF_frozen_g <- mySampleID_IVF_frozen_g

for(i in c(1:length(mySampleID_IVF_frozen_g)) ) {
  mySampleID_IVF_frozen_g[i] = paste(mySampleID_IVF_frozen_g[i], i, sep="_")
}

myTreatment_IVF_frozen_g <- rep( x=2,  times=length(Files_IVF_frozen_g) )

Sex_IVF_frozen_g = rep( x="boy",      times=length(Files_IVF_frozen_g) )
for(i in c(1:length(Sex_IVF_frozen_g)) ) {
  Sex_IVF_frozen_g[i] = "boy"
  if(i>=13) { Sex_IVF_frozen_g[i] = "girl" }
}


Files_ICSI_fresh_g <- c(
  paste( inputDir_g,   "20_W1191D-ICSI-fresh-merge_Rep1.bismark.cov",      sep="/" ),
  paste( inputDir_g,   "21_W1507C-ICSI-fresh_Rep1.bismark.cov",            sep="/" ),
  paste( inputDir_g,   "22_W1636C-ICSI-fresh_Rep1.bismark.cov",            sep="/" ),
  paste( inputDir_g,   "12_W1579C-boy-ICSI-fresh-merge_Rep1.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "12_Q17-W1579D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_W1647C-boy-ICSI-fresh_Rep1.bismark.cov",        sep="/" ),
  paste( inputDir_g,   "13_W1647D-boy-ICSI-fresh_Rep1.bismark.cov",        sep="/" ),
  paste( inputDir_g,   "14_W1719C-boy-ICSI-fresh_Rep1.bismark.cov",        sep="/" ),
  paste( inputDir_g,   "14_W1719D-boy-ICSI-fresh_Rep1.bismark.cov",        sep="/" ),
  paste( inputDir_g,   "34_W200C-boy-ICSI-fresh_Rep1.bismark.cov",         sep="/" ),
  paste( inputDir_g,   "35_W766C-boy-ICSI-fresh-merge_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "36_W928C-boy-ICSI-fresh-merge_Rep1.bismark.cov",   sep="/" ),
  
  paste( inputDir_g,    "20_W1191C-ICSI-fresh_Rep1.bismark.cov",        sep="/"),
  paste( inputDir_g,    "21_W1507D-ICSI-fresh_Rep1.bismark.cov",        sep="/"),
  paste( inputDir_g,    "22_W1636D-ICSI-fresh_Rep1.bismark.cov",        sep="/"),
  paste( inputDir_g,    "5_E95C-girl-ICSI-fresh_Rep1.bismark.cov" ,     sep="/"),
  paste( inputDir_g,    "5_E95D-girl-ICSI-fresh-merge_Rep1.bismark.cov",sep="/"),
  paste( inputDir_g,    "6_W655C-girl-ICSI-fresh_Rep1.bismark.cov" ,    sep="/"),
  paste( inputDir_g,    "6_W655D-girl-ICSI-fresh_Rep1.bismark.cov" ,    sep="/"),
  paste( inputDir_g,    "7_W1276C-girl-ICSI-fresh_Rep1.bismark.cov" ,   sep="/"),
  paste( inputDir_g,    "7_W1276D-girl-ICSI-fresh_Rep1.bismark.cov" ,   sep="/"), 
  paste( inputDir_g,    "48_W808C-girl-ICSI-fresh_Rep1.bismark.cov",    sep="/"),
  paste( inputDir_g,    "49_W924C-girl-ICSI-fresh_Rep1.bismark.cov",    sep="/"),
  paste( inputDir_g,    "50_W934C-girl-ICSI-fresh_Rep1.bismark.cov",    sep="/")     
)

mySampleID_ICSI_fresh_g <- rep( x="ICSI_fresh", times=length(Files_ICSI_fresh_g) )
Tech_ICSI_fresh_g <- mySampleID_ICSI_fresh_g

for(i in c(1:length(mySampleID_ICSI_fresh_g)) ) {
  mySampleID_ICSI_fresh_g[i] = paste(mySampleID_ICSI_fresh_g[i], i, sep="_")
}

myTreatment_ICSI_fresh_g <- rep( x=3,  times=length(Files_ICSI_fresh_g) )

Sex_ICSI_fresh_g = rep( x="boy",      times=length(Files_ICSI_fresh_g) )  
for(i in c(1:length(Sex_ICSI_fresh_g)) ) {
  Sex_ICSI_fresh_g[i] = "boy"
  if(i>=13) { Sex_ICSI_fresh_g[i] = "girl" }
}


Files_ICSI_frozen_g <- c(
paste( inputDir_g,   "EM_10_E10C-boy_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_10_E10D-boy_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_3_E3D-boy_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_6_E6C-boy_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_7_E7C-boy_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_9_E9D-boy_Rep1.CpG.bismark.cov",  sep="/" ),

paste( inputDir_g,   "EM_1_E1D-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_2_E2C-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_2_E2D-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_3_E3C-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_4_E4C-boy_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_4_E4D-boy_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_5_E5C-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_5_E5D-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_7_E7D-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_8_E8C-girl_Rep1.CpG.bismark.cov",  sep="/" ),
paste( inputDir_g,   "EM_9_E9C-girl_Rep1.CpG.bismark.cov",  sep="/" ) 
)

mySampleID_ICSI_frozen_g <- rep( x="EM", times=length(Files_ICSI_frozen_g) )
Tech_ICSI_frozen_g <- mySampleID_ICSI_frozen_g

for(i in c(1:length(mySampleID_ICSI_frozen_g)) ) {
  mySampleID_ICSI_frozen_g[i] = paste(mySampleID_ICSI_frozen_g[i], i, sep="_")
}

myTreatment_ICSI_frozen_g <- rep( x=4,  times=length(Files_ICSI_frozen_g) )

Sex_ICSI_frozen_g = rep( x="boy",      times=length(Files_ICSI_frozen_g) )  
for(i in c(1:length(Sex_ICSI_frozen_g)) ) {
  Sex_ICSI_frozen_g[i] = "boy"
  if(i>=10) { Sex_ICSI_frozen_g[i] = "girl" }
}



Files_All_vector_g <- c(
  Files_NC_g,
  Files_IVF_fresh_g,
  Files_IVF_frozen_g,
  Files_ICSI_fresh_g,
  Files_ICSI_frozen_g  
)
Files_All_list_g <- as.list( Files_All_vector_g )


mySampleID_All_vector_g <- c( 
  mySampleID_NC_g,
  mySampleID_IVF_fresh_g,
  mySampleID_IVF_frozen_g,
  mySampleID_ICSI_fresh_g,
  mySampleID_ICSI_frozen_g  
)
mySampleID_All_list_g <- as.list( mySampleID_All_vector_g )


myTreatment_All_vector_g <- c( 
  myTreatment_NC_g,
  myTreatment_IVF_fresh_g,
  myTreatment_IVF_frozen_g,
  myTreatment_ICSI_fresh_g,
  myTreatment_ICSI_frozen_g 
)       
myTreatment_All_list_g <- as.list( myTreatment_All_vector_g )


mySex_All_vector_g <- c( 
  Sex_NC_g,
  Sex_IVF_fresh_g,
  Sex_IVF_frozen_g,
  Sex_ICSI_fresh_g,
  Sex_ICSI_frozen_g 
)       
mySex_All_list_g <- as.list( mySex_All_vector_g )


myTech_All_vector_g <- c( 
  Tech_NC_g,
  Tech_IVF_fresh_g,
  Tech_IVF_frozen_g,
  Tech_ICSI_fresh_g,
  Tech_ICSI_frozen_g 
)       
myTech_All_list_g <- as.list( myTech_All_vector_g )


## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple
mySex_All_shape_g  = c( "boy"=16, "girl"=17, "father"=1, "mother"=2 ) 
myTech_All_color_g = c( "PCOS"="black", "IVF_fresh"="blue", "IVF_frozen"="green",  "ICSI_fresh"="red", "EM"="purple" )
mySex_All_shape_g
myTech_All_color_g


length( Files_All_vector_g )
length( Files_All_list_g )
length( mySampleID_All_vector_g )
length( mySampleID_All_list_g )
length( myTreatment_All_vector_g )
length( myTreatment_All_list_g )
length( mySex_All_vector_g )
length( mySex_All_list_g )
length( myTech_All_vector_g )
length( myTech_All_list_g )

allSamples_info_g = cbind(Files_All_vector_g, mySampleID_All_vector_g, myTreatment_All_vector_g, mySex_All_vector_g, myTech_All_vector_g)

write.table(allSamples_info_g , 
            file = paste(outDir_g,   "allSamples_info.txt",  sep="/"), 
            append =FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

##################################################################################################################






##################################################################################################################
continue_on_error <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())' " )
}

# This option is very important.
options(error=continue_on_error) 



MySex_Shape_g <- function(  mySex_vector   ) {
  mySex_shape2  = mySex_vector
  for(i in c(1:length(mySex_shape2)) ) {
    if(mySex_shape2[i] == "boy")    { mySex_shape2[i] = c("boy"=16)   }
    if(mySex_shape2[i] == "girl")   { mySex_shape2[i] = c("girl"=17)  }
    if(mySex_shape2[i] == "father") { mySex_shape2[i] = c("father"=1) }
    if(mySex_shape2[i] == "mother") { mySex_shape2[i] = c("mother"=2) }
  }
  mySex_shape2 = as.numeric(mySex_shape2)
  names(mySex_shape2) = mySex_vector
  return(mySex_shape2)
}


MyTech_color_g <- function(  myTech_vector  ) {
  myTech_color  = myTech_vector
  for(i in c(1:length(myTech_color)) ) {
    if(myTech_color[i] == "PCOS")          { myTech_color[i] = c("PCOS"="black") }
    if(myTech_color[i] == "IVF_fresh")   { myTech_color[i] = c("IVF_fresh"="blue") }
    if(myTech_color[i] == "IVF_frozen")  { myTech_color[i] = c("IVF_frozen"="green") }
    if(myTech_color[i] == "ICSI_fresh")  { myTech_color[i] = c("ICSI_fresh"="red") }
    if(myTech_color[i] == "EM") { myTech_color[i] = c("EM"="purple") }
  }
  names(myTech_color) = myTech_vector
  return(myTech_color)
}


MyTech_number_g <- function(  myTech_vector  ) {
  myTech_color  = myTech_vector
  for(i in c(1:length(myTech_color)) ) {
    if(myTech_color[i] == "PCOS")          { myTech_color[i] = 1 }
    if(myTech_color[i] == "IVF_fresh")   { myTech_color[i] = 2 }
    if(myTech_color[i] == "IVF_frozen")  { myTech_color[i] = 3 }
    if(myTech_color[i] == "ICSI_fresh")  { myTech_color[i] = 4 }
    if(myTech_color[i] == "EM") { myTech_color[i] = 5 }
  }
  return(myTech_color)
}





MyTheme_1_g <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
  theme(  
    line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                                                        ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                                                 ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all text elements.           "serif" for a serif font. 所有文本相关属性.
    title = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    ## aspect.ratio = 1,   ##高宽比
    
    axis.title    = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis label (element_text; inherits from axis.title)
    axis.title.y  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=90,      lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis label (element_text; inherits from axis.title)
    axis.text     = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1,  lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    
    axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                      ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	              ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                  ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                  ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.spacing      = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.spacing.x    = grid::unit(1, "mm", data=NULL) ,
    panel.spacing.y    = grid::unit(1, "mm", data=NULL) ,
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=0.5, vjust=0.5,   angle=NULL, lineheight=NULL),     ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL,    size=NULL, linetype=NULL, fill=NULL ), 	                                                      ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 


MySaveGgplot2_1_g <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200,   device=cairo_ps)         
}


## df contains two columns, the first column (cond_col=1) is sample type, the second column (val_col=2) is value. (must be).
whisk_1_g <- function(df, cond_col=1, val_col=2) {  
  require(reshape2)
  condname <- names(df)[cond_col]  ## save the name of the first column.
  names(df)[cond_col] <- "cond" 
  names(df)[val_col]  <- "value"
  b   <- boxplot(value~cond, data=df, plot=FALSE)   
  df2 <- cbind(as.data.frame(b$stats), c("min","lq","m","uq","max"))
  names(df2) <- c(levels(df$cond), "pos")
  df2 <- melt(df2, id="pos", variable.name="cond")
  df2 <- dcast(df2, cond~pos)   
  names(df2)[1] <- condname 
  print(df2)
  df2
}


MyBoxViolinPlot_1_g <- function(vector2,   sampleType2,  colours2,   path2,   fileName2,  title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
  vector2[vector2>Ymax2] <- Ymax2
  vector2[vector2<Ymin2] <- Ymin2
  DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
  if( ! file.exists(path2) ) { dir.create(path2, recursive = TRUE) }
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) +  
    geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1_g(DataFrame_Local),   width=0.2, size=0.5 ) +
    geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=NA,  outlier.size=0, size=0.5, fill=colours2 ) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour=colours2, geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot",    sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) +  
    geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1_g(DataFrame_Local),   width=0.2, size=0.5 ) +
    geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=NA,  outlier.size=0, size=0.5, fill=colours2 ) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour=colours2, geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot_noErrorBar",    sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray" ) +   
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 2) +  
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-2Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-3Adjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 4) +  
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-4Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)  
}  






MyCluster_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) { 
  if( nrow(mymeth1)*(1-sdThres1) > 10 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  dev.off()
  }
}


MyCluster_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {  
  myTempValue = tryCatch(
    clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE, plot=FALSE  ), 
    error = function(err){"000"}
  )

  if( length( as.character(myTempValue) ) > 1) {
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  dev.off()
  }
}


MyCluster_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
  myTempFunction <- function() {
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1P.kept0.5percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.995)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1Q.kept0.1percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.999)                    
 
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    )                     
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.18.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.22.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.22 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.28.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.28 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.30.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30 )                    
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  )
}






MyPCA_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) { 
  if( nrow(mymeth1)*(1-sdThres1) > 10 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  dev.off()
  pdf( file=paste(path1, "/", file1, ".screeplot.pdf", sep="") , width=width1/2, height=height1/2  )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  dev.off()
  }
}


MyPCA_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {
  myTempFunction <- function() {  
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=FALSE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  dev.off()
  pdf( file=paste(path1, "/", file1, ".screeplot.pdf", sep="") , width=width1, height=height1  )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, screeplot=TRUE,   scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  dev.off()
  }

  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  )
}


MyPCA_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
  myTempFunction <- function() {
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1P.kept0.5percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.995)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1Q.kept0.1percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.999)                    
  
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    )                     
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.18.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.22.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.22 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.28.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.28 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.30.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30 )                    
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  ) 
}






myHierarchicalClustering_1_g  <- function(  mat_3three,   path_temp1,  width1, height1   )  {
  class_3three <- colnames(mat_3three)
  class_3three <- gsub(pattern="_\\d+$", replacement="", x=class_3three, ignore.case = FALSE, perl = TRUE )
  MyTech_color_g( class_3three )
  
  res.dist1_3three <- get_dist( t(mat_3three) ,   method = "euclidean")
  res.dist2_3three <- get_dist( t(mat_3three) ,   method = "maximum"  )
  res.dist3_3three <- get_dist( t(mat_3three) ,   method = "manhattan")
  res.dist4_3three <- get_dist( t(mat_3three) ,   method = "canberra" )
  res.dist5_3three <- get_dist( t(mat_3three) ,   method = "binary"   )
  res.dist6_3three <- get_dist( t(mat_3three) ,   method = "minkowski")
  res.dist7_3three <- get_dist( t(mat_3three) ,   method = "pearson"  )
  res.dist8_3three <- get_dist( t(mat_3three) ,   method = "spearman" )
  
  pdf( file = paste(path_temp1, "1A_visualizing-distance-matrix.pdf",  sep="/"), width=12, height=10 )
  print( fviz_dist(res.dist1_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) )
  print( fviz_dist(res.dist2_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist3_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist4_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist5_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist6_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist7_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist8_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  dev.off() 
  
  sink( file = paste(path_temp1, "1B_all-distance-matrix.txt",  sep="/") )
  print("################### euclidean: ")
  print(res.dist1_3three ) 
  print("################### maximum: ")
  print(res.dist2_3three ) 
  print("################### manhattan: ")
  print(res.dist3_3three ) 
  print("################### canberra: ")
  print(res.dist4_3three ) 
  print("################### binary: ")
  print(res.dist5_3three ) 
  print("################### minkowski: ")
  print(res.dist6_3three ) 
  print("################### pearson: ")
  print(res.dist7_3three ) 
  print("################### spearman: ")
  print(res.dist8_3three ) 
  sink() 
  
  
  # Compute hierarchical clustering by "ward.D"
  res.hc1_3three <- hclust(res.dist1_3three, method = "ward.D"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "ward.D"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "ward.D"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "ward.D"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "ward.D"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "ward.D"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "ward.D"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "ward.D"  ) 
  
  pdf( file = paste(path_temp1, "2A_hierarchical-ward.D-rectangle.pdf",  sep="/") , width=width1, height=height1 )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="ward.D, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="ward.D, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="ward.D, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="ward.D, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="ward.D, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="ward.D, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="ward.D, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="ward.D, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "2B_hierarchical-ward.D-phylogenic.pdf",  sep="/")  )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="ward.D, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="ward.D, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="ward.D, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="ward.D, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="ward.D, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="ward.D, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="ward.D, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="ward.D, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  res.hca.color1_3three <- as.dendrogram( res.hc1_3three )   
  res.hca.color2_3three <- as.dendrogram( res.hc2_3three )   
  res.hca.color3_3three <- as.dendrogram( res.hc3_3three )  
  res.hca.color4_3three <- as.dendrogram( res.hc4_3three )   
  res.hca.color5_3three <- as.dendrogram( res.hc5_3three )  
  res.hca.color6_3three <- as.dendrogram( res.hc6_3three )  
  res.hca.color7_3three <- as.dendrogram( res.hc7_3three )  
  res.hca.color8_3three <- as.dendrogram( res.hc8_3three ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three)]
  
  labels_colors(res.hca.color1_3three) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three) <-   colors_to_use6
  labels_colors(res.hca.color7_3three) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "2C_hierarchical-ward.D-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three,   main="ward.D, euclidean" )  )
  print(  plot(res.hca.color2_3three,   main="ward.D, maximum"   )  )
  print(  plot(res.hca.color3_3three,   main="ward.D, manhattan" )  )
  print(  plot(res.hca.color4_3three,   main="ward.D, canberra"  )  )
  print(  plot(res.hca.color5_3three,   main="ward.D, binary"    )  )
  print(  plot(res.hca.color6_3three,   main="ward.D, minkowski" )  )
  print(  plot(res.hca.color7_3three,   main="ward.D, pearson"   )  )
  print(  plot(res.hca.color8_3three,   main="ward.D, spearman"  )  )          
  dev.off() 
  
  
  
  # Compute hierarchical clustering by "ward.D2"
  res.hc1_3three1 <- hclust(res.dist1_3three, method = "ward.D2"  )   
  res.hc2_3three1 <- hclust(res.dist2_3three, method = "ward.D2"  )   
  res.hc3_3three1 <- hclust(res.dist3_3three, method = "ward.D2"  )  
  res.hc4_3three1 <- hclust(res.dist4_3three, method = "ward.D2"  )   
  res.hc5_3three1 <- hclust(res.dist5_3three, method = "ward.D2"  )  
  res.hc6_3three1 <- hclust(res.dist6_3three, method = "ward.D2"  )  
  res.hc7_3three1 <- hclust(res.dist7_3three, method = "ward.D2"  )  
  res.hc8_3three1 <- hclust(res.dist8_3three, method = "ward.D2"  ) 
  
  pdf( file = paste(path_temp1, "3A_hierarchical-ward.D2-rectangle.pdf",  sep="/") , width=width1, height=height1 )
  print(  fviz_dend(res.hc1_3three1,   type="rectangle", main="ward.D2, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three1,   type="rectangle", main="ward.D2, maximum"    )  )
  print(  fviz_dend(res.hc3_3three1,   type="rectangle", main="ward.D2, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three1,   type="rectangle", main="ward.D2, canberra"   )  )
  print(  fviz_dend(res.hc5_3three1,   type="rectangle", main="ward.D2, binary"     )  )
  print(  fviz_dend(res.hc6_3three1,   type="rectangle", main="ward.D2, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three1,   type="rectangle", main="ward.D2, pearson"    )  )
  print(  fviz_dend(res.hc8_3three1,   type="rectangle", main="ward.D2, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "3B_hierarchical-ward.D2-phylogenic.pdf",  sep="/")  )
  print(  fviz_dend(res.hc1_3three1,   type="phylogenic", main="ward.D2, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three1,   type="phylogenic", main="ward.D2, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three1,   type="phylogenic", main="ward.D2, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three1,   type="phylogenic", main="ward.D2, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three1,   type="phylogenic", main="ward.D2, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three1,   type="phylogenic", main="ward.D2, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three1,   type="phylogenic", main="ward.D2, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three1,   type="phylogenic", main="ward.D2, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  
  res.hca.color1_3three1 <- as.dendrogram( res.hc1_3three1 )   
  res.hca.color2_3three1 <- as.dendrogram( res.hc2_3three1 )   
  res.hca.color3_3three1 <- as.dendrogram( res.hc3_3three1 )  
  res.hca.color4_3three1 <- as.dendrogram( res.hc4_3three1 )   
  res.hca.color5_3three1 <- as.dendrogram( res.hc5_3three1 )  
  res.hca.color6_3three1 <- as.dendrogram( res.hc6_3three1 )  
  res.hca.color7_3three1 <- as.dendrogram( res.hc7_3three1 )  
  res.hca.color8_3three1 <- as.dendrogram( res.hc8_3three1 ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three1)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three1)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three1)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three1)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three1)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three1)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three1)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three1)]
  
  labels_colors(res.hca.color1_3three1) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three1) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three1) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three1) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three1) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three1) <-   colors_to_use6
  labels_colors(res.hca.color7_3three1) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three1) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "3C_hierarchical-ward.D2-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three1,   main="ward.D2, euclidean" )  )
  print(  plot(res.hca.color2_3three1,   main="ward.D2, maximum"   )  )
  print(  plot(res.hca.color3_3three1,   main="ward.D2, manhattan" )  )
  print(  plot(res.hca.color4_3three1,   main="ward.D2, canberra"  )  )
  print(  plot(res.hca.color5_3three1,   main="ward.D2, binary"    )  )
  print(  plot(res.hca.color6_3three1,   main="ward.D2, minkowski" )  )
  print(  plot(res.hca.color7_3three1,   main="ward.D2, pearson"   )  )
  print(  plot(res.hca.color8_3three1,   main="ward.D2, spearman"  )  )          
  dev.off() 
  
  
  
  
  # Compute hierarchical clustering by "single"
  res.hc1_3three <- hclust(res.dist1_3three, method = "single"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "single"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "single"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "single"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "single"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "single"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "single"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "single"  ) 
  
  pdf( file = paste(path_temp1, "4A_hierarchical-single-rectangle.pdf",  sep="/") , width=width1, height=height1 )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="single, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="single, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="single, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="single, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="single, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="single, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="single, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="single, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "4B_hierarchical-single-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="single, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="single, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="single, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="single, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="single, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="single, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="single, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="single, spearman"   , repel = TRUE)  )         
  dev.off() 

  
  
  
  res.hca.color1_3three <- as.dendrogram( res.hc1_3three )   
  res.hca.color2_3three <- as.dendrogram( res.hc2_3three )   
  res.hca.color3_3three <- as.dendrogram( res.hc3_3three )  
  res.hca.color4_3three <- as.dendrogram( res.hc4_3three )   
  res.hca.color5_3three <- as.dendrogram( res.hc5_3three )  
  res.hca.color6_3three <- as.dendrogram( res.hc6_3three )  
  res.hca.color7_3three <- as.dendrogram( res.hc7_3three )  
  res.hca.color8_3three <- as.dendrogram( res.hc8_3three ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three)]
  
  labels_colors(res.hca.color1_3three) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three) <-   colors_to_use6
  labels_colors(res.hca.color7_3three) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "4C_hierarchical-single-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three,   main="single, euclidean" )  )
  print(  plot(res.hca.color2_3three,   main="single, maximum"   )  )
  print(  plot(res.hca.color3_3three,   main="single, manhattan" )  )
  print(  plot(res.hca.color4_3three,   main="single, canberra"  )  )
  print(  plot(res.hca.color5_3three,   main="single, binary"    )  )
  print(  plot(res.hca.color6_3three,   main="single, minkowski" )  )
  print(  plot(res.hca.color7_3three,   main="single, pearson"   )  )
  print(  plot(res.hca.color8_3three,   main="single, spearman"  )  )          
  dev.off() 
  
  
  
  
  # Compute hierarchical clustering by "complete"
  res.hc1_3three <- hclust(res.dist1_3three, method = "complete"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "complete"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "complete"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "complete"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "complete"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "complete"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "complete"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "complete"  ) 
  
  pdf( file = paste(path_temp1, "5A_hierarchical-complete-rectangle.pdf",  sep="/") , width=width1, height=height1 )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="complete, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="complete, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="complete, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="complete, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="complete, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="complete, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="complete, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="complete, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "5B_hierarchical-complete-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="complete, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="complete, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="complete, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="complete, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="complete, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="complete, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="complete, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="complete, spearman"   , repel = TRUE)  )          
  dev.off() 

  res.hca.color1_3three <- as.dendrogram( res.hc1_3three )   
  res.hca.color2_3three <- as.dendrogram( res.hc2_3three )   
  res.hca.color3_3three <- as.dendrogram( res.hc3_3three )  
  res.hca.color4_3three <- as.dendrogram( res.hc4_3three )   
  res.hca.color5_3three <- as.dendrogram( res.hc5_3three )  
  res.hca.color6_3three <- as.dendrogram( res.hc6_3three )  
  res.hca.color7_3three <- as.dendrogram( res.hc7_3three )  
  res.hca.color8_3three <- as.dendrogram( res.hc8_3three ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three)]
  
  labels_colors(res.hca.color1_3three) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three) <-   colors_to_use6
  labels_colors(res.hca.color7_3three) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "5C_hierarchical-complete-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three,   main="complete, euclidean" )  )
  print(  plot(res.hca.color2_3three,   main="complete, maximum"   )  )
  print(  plot(res.hca.color3_3three,   main="complete, manhattan" )  )
  print(  plot(res.hca.color4_3three,   main="complete, canberra"  )  )
  print(  plot(res.hca.color5_3three,   main="complete, binary"    )  )
  print(  plot(res.hca.color6_3three,   main="complete, minkowski" )  )
  print(  plot(res.hca.color7_3three,   main="complete, pearson"   )  )
  print(  plot(res.hca.color8_3three,   main="complete, spearman"  )  )          
  dev.off() 
  
  
  
  # Compute hierarchical clustering by "average"
  res.hc1_3three <- hclust(res.dist1_3three, method = "average"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "average"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "average"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "average"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "average"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "average"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "average"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "average"  ) 
  
  pdf( file = paste(path_temp1, "6A_hierarchical-average-rectangle.pdf",  sep="/") , width=width1, height=height1 )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="average, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="average, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="average, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="average, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="average, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="average, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="average, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="average, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "6B_hierarchical-average-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="average, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="average, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="average, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="average, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="average, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="average, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="average, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="average, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  
  res.hca.color1_3three <- as.dendrogram( res.hc1_3three )   
  res.hca.color2_3three <- as.dendrogram( res.hc2_3three )   
  res.hca.color3_3three <- as.dendrogram( res.hc3_3three )  
  res.hca.color4_3three <- as.dendrogram( res.hc4_3three )   
  res.hca.color5_3three <- as.dendrogram( res.hc5_3three )  
  res.hca.color6_3three <- as.dendrogram( res.hc6_3three )  
  res.hca.color7_3three <- as.dendrogram( res.hc7_3three )  
  res.hca.color8_3three <- as.dendrogram( res.hc8_3three ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three)]
  
  labels_colors(res.hca.color1_3three) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three) <-   colors_to_use6
  labels_colors(res.hca.color7_3three) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "6C_hierarchical-average-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three,   main="average, euclidean" )  )
  print(  plot(res.hca.color2_3three,   main="average, maximum"   )  )
  print(  plot(res.hca.color3_3three,   main="average, manhattan" )  )
  print(  plot(res.hca.color4_3three,   main="average, canberra"  )  )
  print(  plot(res.hca.color5_3three,   main="average, binary"    )  )
  print(  plot(res.hca.color6_3three,   main="average, minkowski" )  )
  print(  plot(res.hca.color7_3three,   main="average, pearson"   )  )
  print(  plot(res.hca.color8_3three,   main="average, spearman"  )  )          
  dev.off() 

  
  
  # Compute hierarchical clustering by "mcquitty"
  res.hc1_3three <- hclust(res.dist1_3three, method = "mcquitty"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "mcquitty"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "mcquitty"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "mcquitty"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "mcquitty"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "mcquitty"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "mcquitty"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "mcquitty"  ) 
  
  pdf( file = paste(path_temp1, "7A_hierarchical-mcquitty-rectangle.pdf",  sep="/") , width=width1, height=height1 )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="mcquitty, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="mcquitty, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="mcquitty, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="mcquitty, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="mcquitty, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="mcquitty, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="mcquitty, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="mcquitty, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "7B_hierarchical-mcquitty-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="mcquitty, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="mcquitty, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="mcquitty, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="mcquitty, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="mcquitty, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="mcquitty, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="mcquitty, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="mcquitty, spearman"   , repel = TRUE)  )          
  dev.off() 


  res.hca.color1_3three <- as.dendrogram( res.hc1_3three )   
  res.hca.color2_3three <- as.dendrogram( res.hc2_3three )   
  res.hca.color3_3three <- as.dendrogram( res.hc3_3three )  
  res.hca.color4_3three <- as.dendrogram( res.hc4_3three )   
  res.hca.color5_3three <- as.dendrogram( res.hc5_3three )  
  res.hca.color6_3three <- as.dendrogram( res.hc6_3three )  
  res.hca.color7_3three <- as.dendrogram( res.hc7_3three )  
  res.hca.color8_3three <- as.dendrogram( res.hc8_3three ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three)]
  
  labels_colors(res.hca.color1_3three) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three) <-   colors_to_use6
  labels_colors(res.hca.color7_3three) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "7C_hierarchical-mcquitty-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three,   main="mcquitty, euclidean" )  )
  print(  plot(res.hca.color2_3three,   main="mcquitty, maximum"   )  )
  print(  plot(res.hca.color3_3three,   main="mcquitty, manhattan" )  )
  print(  plot(res.hca.color4_3three,   main="mcquitty, canberra"  )  )
  print(  plot(res.hca.color5_3three,   main="mcquitty, binary"    )  )
  print(  plot(res.hca.color6_3three,   main="mcquitty, minkowski" )  )
  print(  plot(res.hca.color7_3three,   main="mcquitty, pearson"   )  )
  print(  plot(res.hca.color8_3three,   main="mcquitty, spearman"  )  )          
  dev.off() 
  
  
  
  # Compute hierarchical clustering by "median"
  res.hc1_3three <- hclust(res.dist1_3three, method = "median"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "median"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "median"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "median"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "median"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "median"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "median"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "median"  ) 
  
  pdf( file = paste(path_temp1, "8A_hierarchical-median-rectangle.pdf",  sep="/") , width=width1, height=height1 )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="median, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="median, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="median, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="median, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="median, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="median, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="median, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="median, spearman"   )  )
  dev.off() 
  
  
  res.hca.color1_3three <- as.dendrogram( res.hc1_3three )   
  res.hca.color2_3three <- as.dendrogram( res.hc2_3three )   
  res.hca.color3_3three <- as.dendrogram( res.hc3_3three )  
  res.hca.color4_3three <- as.dendrogram( res.hc4_3three )   
  res.hca.color5_3three <- as.dendrogram( res.hc5_3three )  
  res.hca.color6_3three <- as.dendrogram( res.hc6_3three )  
  res.hca.color7_3three <- as.dendrogram( res.hc7_3three )  
  res.hca.color8_3three <- as.dendrogram( res.hc8_3three ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three)]
  
  labels_colors(res.hca.color1_3three) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three) <-   colors_to_use6
  labels_colors(res.hca.color7_3three) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "8C_hierarchical-median-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three,   main="median, euclidean" )  )
  print(  plot(res.hca.color2_3three,   main="median, maximum"   )  )
  print(  plot(res.hca.color3_3three,   main="median, manhattan" )  )
  print(  plot(res.hca.color4_3three,   main="median, canberra"  )  )
  print(  plot(res.hca.color5_3three,   main="median, binary"    )  )
  print(  plot(res.hca.color6_3three,   main="median, minkowski" )  )
  print(  plot(res.hca.color7_3three,   main="median, pearson"   )  )
  print(  plot(res.hca.color8_3three,   main="median, spearman"  )  )          
  dev.off() 
  
  
  
  
  
  # Compute hierarchical clustering by "centroid"
  res.hc1_3three <- hclust(res.dist1_3three, method = "centroid"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "centroid"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "centroid"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "centroid"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "centroid"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "centroid"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "centroid"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "centroid"  ) 
  
  pdf( file = paste(path_temp1, "9A_hierarchical-centroid-rectangle.pdf",  sep="/"), width=width1, height=height1  )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="centroid, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="centroid, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="centroid, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="centroid, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="centroid, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="centroid, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="centroid, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="centroid, spearman"   )  )
  dev.off() 
  
  
  
  res.hca.color1_3three <- as.dendrogram( res.hc1_3three )   
  res.hca.color2_3three <- as.dendrogram( res.hc2_3three )   
  res.hca.color3_3three <- as.dendrogram( res.hc3_3three )  
  res.hca.color4_3three <- as.dendrogram( res.hc4_3three )   
  res.hca.color5_3three <- as.dendrogram( res.hc5_3three )  
  res.hca.color6_3three <- as.dendrogram( res.hc6_3three )  
  res.hca.color7_3three <- as.dendrogram( res.hc7_3three )  
  res.hca.color8_3three <- as.dendrogram( res.hc8_3three ) 
  
  colors_to_use <- as.numeric( MyTech_number_g( class_3three ) ) 
  length(colors_to_use)
  
  colors_to_use1 <- colors_to_use[order.dendrogram(res.hca.color1_3three)]
  colors_to_use2 <- colors_to_use[order.dendrogram(res.hca.color2_3three)]
  colors_to_use3 <- colors_to_use[order.dendrogram(res.hca.color3_3three)]
  colors_to_use4 <- colors_to_use[order.dendrogram(res.hca.color4_3three)]
  colors_to_use5 <- colors_to_use[order.dendrogram(res.hca.color5_3three)]
  colors_to_use6 <- colors_to_use[order.dendrogram(res.hca.color6_3three)]
  colors_to_use7 <- colors_to_use[order.dendrogram(res.hca.color7_3three)]
  colors_to_use8 <- colors_to_use[order.dendrogram(res.hca.color8_3three)]
  
  labels_colors(res.hca.color1_3three) <-   colors_to_use1 
  labels_colors(res.hca.color2_3three) <-   colors_to_use2 
  labels_colors(res.hca.color3_3three) <-   colors_to_use3 
  labels_colors(res.hca.color4_3three) <-   colors_to_use4 
  labels_colors(res.hca.color5_3three) <-   colors_to_use5 
  labels_colors(res.hca.color6_3three) <-   colors_to_use6
  labels_colors(res.hca.color7_3three) <-   colors_to_use7 
  labels_colors(res.hca.color8_3three) <-   colors_to_use8 
  
  pdf( file = paste(path_temp1, "9C_hierarchical-centroid-plot.pdf",  sep="/") , width=width1, height=height1 )
  print(  plot(res.hca.color1_3three,   main="centroid, euclidean" )  )
  print(  plot(res.hca.color2_3three,   main="centroid, maximum"   )  )
  print(  plot(res.hca.color3_3three,   main="centroid, manhattan" )  )
  print(  plot(res.hca.color4_3three,   main="centroid, canberra"  )  )
  print(  plot(res.hca.color5_3three,   main="centroid, binary"    )  )
  print(  plot(res.hca.color6_3three,   main="centroid, minkowski" )  )
  print(  plot(res.hca.color7_3three,   main="centroid, pearson"   )  )
  print(  plot(res.hca.color8_3three,   main="centroid, spearman"  )  )          
  dev.off() 
 
}



#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
MyPrcompObj_1_g <- function(  prcompObj2,   path2,   file2,  dataFrame_temp2   ) {
  
  sink( file = paste(path2, "/1A_PCA_results_",  file2,  ".txt",  sep="") )
  print( prcompObj2 )
  print( names(prcompObj2) )
  sink()
  
  sink( file = paste(path2, "/1B_PCA_summary_",  file2,  ".txt",   sep="") )
  print( summary(prcompObj2) )
  sink()
  
  sink( file = paste(path2, "/2_PCA_all_",  file2,  ".txt",   sep="") )
  print("####################### prcompObj2$sdev #########################")
  print(prcompObj2$sdev)
  print("####################### prcompObj2$rotation #########################")
  print(prcompObj2$rotation)
  print("####################### prcompObj2$center #########################")
  print(prcompObj2$center)
  print("####################### prcompObj2$scale #########################")
  print(prcompObj2$scale)
  print("####################### prcompObj2$x #########################")
  print(prcompObj2$x)
  sink()
  
  pdf( file = paste(path2, "/3_PCA_info_",  file2,  ".pdf",   sep="")  )
  plot(prcompObj2, type="lines")
  print( fviz_eig(prcompObj2) )
  dev.off() 
  
  my_fviz_pca_ind1 <- fviz_pca_ind(prcompObj2,
                                   col.ind = "cos2", # Color by the quality of representation
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                   repel = TRUE     # Avoid text overlapping
  )
  my_fviz_pca_ind2 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE
  )
  my_fviz_pca_ind3 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none"
  )
  my_fviz_pca_ind4 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   ggtheme = MyTheme_1_g()
  )
  my_fviz_pca_ind5 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none",
                                   ggtheme = MyTheme_1_g()
  )
  
  pdf( file=paste(path2, "/4_PCA-2D_", file2, ".pdf",  sep="") )
  print(my_fviz_pca_ind1)
  print(my_fviz_pca_ind2)
  print(my_fviz_pca_ind3)
  print(my_fviz_pca_ind4)
  print(my_fviz_pca_ind5)
  dev.off() 
  
  #############################
  prcompObj2_matrix <- prcompObj2$x
  prcompObj2_Contri  <- (prcompObj2$sdev)^2
  prcompObj2_Contri  <- prcompObj2_Contri/sum(prcompObj2_Contri)
  prcompObj2_Contri  <- prcompObj2_Contri * 100
  prcompObj2_Contri  <- round(prcompObj2_Contri, 2)
  
  label1_2two <-   paste( "PC1 ",  "(", prcompObj2_Contri[1], "%)", sep="" )
  label2_2two <-   paste( "PC2 ",  "(", prcompObj2_Contri[2], "%)", sep="" )
  label3_2two <-   paste( "PC3 ",  "(", prcompObj2_Contri[3], "%)", sep="" ) 
  
  dataframeA_2two  <- data.frame( as.data.frame(prcompObj2_matrix), mySex= as.vector(dataFrame_temp2$mysex), 
                                  myTech=as.vector(dataFrame_temp2$mytech),    myLabel=as.vector(dataFrame_temp2$mysampleID)   ) 
  
  FigureTemp1_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=3, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1_2two ,  path1=path2, fileName1="5A_PCA-PC1-PC2",  height1=2.5,  width1=4.3)
  
  FigureTemp2_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=5, alpha=0.5  )+ xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2_2two ,  path1=path2, fileName1="5B_PCA-PC1-PC2-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp4 <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=6, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +     
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4_2two ,  path1=path2, fileName1="5D_PCA-PC1-PC2-big",   height1=2.5,  width1=4.3)
  
  FigureTemp5_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=2, alpha=0.7  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp5_2two ,  path1=path2, fileName1="5E_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  
  FigureTemp6_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=5, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp6_2two ,  path1=path2, fileName1="5F_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  

  ## PC1, 2, and 3     dev.off() 
  pdf( file = paste(path2, "6_PCA-3d-by-scatter3D.pdf",  sep="/") )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),     
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),          
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ##############
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ######
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
             col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  dev.off()
}


#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
MyPCAobj_FactoMineR_g <- function(  PCAobj2,   path2,   file2,  dataFrame_temp2   ) {
  
  sink( file = paste(path2, "/1A_PCA_results_",  file2,  ".txt",  sep="") )
  print( PCAobj2 )
  print( names(PCAobj2) )
  sink()
  
  sink( file = paste(path2, "/1B_PCA_summary_",  file2,  ".txt",   sep="") )
  print( summary(PCAobj2) )
  print( "#################################"  )
  myeig.val_PCAobj2 <- get_eigenvalue( PCAobj2 )
  print(myeig.val_PCAobj2)
  sink()
  
  sink( file = paste(path2, "/2_PCA_all_",  file2,  ".txt",   sep="") )
  print("####################### PCAobj2$ind$coord #########################")
  print(PCAobj2$ind$coord)
  print("####################### PCAobj2$ind #########################")
  print(PCAobj2$ind)
  sink()

  pdf( file = paste(path2, "/3_PCA_info_",  file2,  ".pdf",   sep="")  )
  print( fviz_eig(PCAobj2, addlabels = TRUE ) )
  print( fviz_screeplot(X=PCAobj2, choice = "variance", geom = "line",
                 barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                 ncp = 10, addlabels = TRUE) )
  print( fviz_screeplot(X=PCAobj2, choice = "eigenvalue", geom = "line",
                 barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                 ncp = 10, addlabels = TRUE) )
  print( fviz_screeplot(X=PCAobj2, choice = "variance", geom = "bar",
                 barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                 ncp = 10, addlabels = TRUE) )
  print( fviz_screeplot(X=PCAobj2, choice = "eigenvalue", geom = "bar",
                 barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                 ncp = 10, addlabels = TRUE) )
  dev.off() 

  
  my_fviz_pca_ind1 <- fviz_pca_ind(PCAobj2,
                                   col.ind = "cos2", # Color by the quality of representation
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                   repel = TRUE     # Avoid text overlapping
  )
  my_fviz_pca_ind2 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE
  )
  my_fviz_pca_ind3 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none"
  )
  my_fviz_pca_ind4 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   ggtheme = MyTheme_1_g()
  )
  my_fviz_pca_ind5 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none",
                                   ggtheme = MyTheme_1_g()
  )
  
  pdf( file=paste(path2, "/4_PCA-2D_", file2, ".pdf",  sep="") )
  print(my_fviz_pca_ind1)
  print(my_fviz_pca_ind2)
  print(my_fviz_pca_ind3)
  print(my_fviz_pca_ind4)
  print(my_fviz_pca_ind5)
  dev.off() 
  

  #############################
  PCAobj2_matrix <- PCAobj2$ind$coord 
  PCAobj2_Contri  <- (PCAobj2$eig)[,2]
  PCAobj2_Contri  <- PCAobj2_Contri/sum(PCAobj2_Contri)
  PCAobj2_Contri  <- PCAobj2_Contri * 100
  PCAobj2_Contri  <- round(PCAobj2_Contri, 2)
  
  label1_2two <-   paste( "PC1 ",  "(", PCAobj2_Contri[1], "%)", sep="" )
  label2_2two <-   paste( "PC2 ",  "(", PCAobj2_Contri[2], "%)", sep="" )
  label3_2two <-   paste( "PC3 ",  "(", PCAobj2_Contri[3], "%)", sep="" ) 
  
  dataframeA_2two  <- data.frame( as.data.frame(PCAobj2_matrix), mySex= as.vector(dataFrame_temp2$mysex), 
                                  myTech=as.vector(dataFrame_temp2$mytech),    myLabel=as.vector(dataFrame_temp2$mysampleID)   )  

  FigureTemp1_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=3, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1_2two ,  path1=path2, fileName1="5A_PCA-PC1-PC2",  height1=2.5,  width1=4.3)
  
  FigureTemp2_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=5, alpha=0.5  )+ xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2_2two ,  path1=path2, fileName1="5B_PCA-PC1-PC2-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp4 <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=6, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +     
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4_2two ,  path1=path2, fileName1="5D_PCA-PC1-PC2-big",   height1=2.5,  width1=4.3)
  
  FigureTemp5_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=2, alpha=0.7  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp5_2two ,  path1=path2, fileName1="5E_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  
  FigureTemp6_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=5, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp6_2two ,  path1=path2, fileName1="5F_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  
  
  ## PC1, 2, and 3     dev.off() 
  pdf( file = paste(path2, "6_PCA-3d-by-scatter3D.pdf",  sep="/") )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),     
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),          
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ##############
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ######
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  dev.off() 
}






MyPCA_1A_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1,  dataFrame_temp1    ) {   
  if( nrow(mymeth1)*(1-sdThres1) > 10 ) {   
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, file1, "filterByQuantile_1", sep="/")
  path1_2 = paste(path1, file1, "filterByQuantile_2", sep="/")
  path1_3 = paste(path1, file1, "filterByQuantile_3", sep="/")
  path1_4 = paste(path1, file1, "filterByQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  }
}


MyPCA_2A_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1,  dataFrame_temp1   ) { 
  myTempFunction <- function() {
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }   
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, file1, "filterNoQuantile_1", sep="/")
  path1_2 = paste(path1, file1, "filterNoQuantile_2", sep="/")
  path1_3 = paste(path1, file1, "filterNoQuantile_3", sep="/")
  path1_4 = paste(path1, file1, "filterNoQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )  
  }
  
  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  )
}


MyPCA_3A_g <- function(  mymeth2 ,  path2,   file2, width2, height2,  dataFrame_temp2   )  {
  myTempFunction <- function() {
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1P.kept0.5percent", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.995 ,dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1Q.kept0.1percent", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.999 ,dataFrame_temp1=dataFrame_temp2   )                     
  
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.18",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.22",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.22 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.28",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.28 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.30",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30 ,  dataFrame_temp1=dataFrame_temp2   )                     
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"000"}
  ) 
}





myDMRs_annotation_g <- function(  myDiffDMR_5,   path2_5  ) {
  if( ! file.exists(path2_5) ) { dir.create(path2_5, recursive = TRUE) }
  
  ########################## all genes.
  Ann_gene = annotateWithGeneParts( as(myDiffDMR_5,  "GRanges"),  gene.obj_g)
  sink( file=paste(path2_5, "1000A-distribution-onGenes.txt", sep="/")   )
  getFeatsWithTargetsStats( Ann_gene,  percentage=TRUE)
  print(Ann_gene)
  sink()
  pdf( file=paste(path2_5, "1000B-distribution-onGenes.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_gene,  precedence=TRUE, main="Annotation on all genes") )
  dev.off()
  
  ########################## Repeats
  Ann_Repeats = annotateWithFeatureFlank( as(myDiffDMR_5,"GRanges"), myrepeat.obj_g$Repeats,  myrepeat.obj_g$shores,
                                          feature.name="Repeats",  flank.name="shores" )
  sink( file=paste(path2_5, "1001A-distribution-onRepeats.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_Repeats,  percentage=TRUE)
  print(Ann_Repeats)
  sink()
  pdf( file=paste(path2_5, "1001B-distribution-onRepeats.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_Repeats, precedence=TRUE, main="Annotation on Repeats"))
  dev.off()
  
  ########################## ImprintedRegions 1
  diffImprinted1Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint1.obj_g$ImprintedRegions,  flank=imprint1.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002A-distribution-onImprintedRegions1.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted1Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002B-distribution-onImprintedRegions1.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted1Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 2
  diffImprinted2Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint2.obj_g$ImprintedRegions,  flank=imprint2.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002C-distribution-onImprintedRegions2.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted2Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002D-distribution-onImprintedRegions2.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted2Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 3
  diffImprinted3Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint3.obj_g$ImprintedRegions,  flank=imprint3.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002E-distribution-onImprintedRegions3.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted3Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002F-distribution-onImprintedRegions3.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted3Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 4
  diffImprinted4Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint4.obj_g$ImprintedRegions,  flank=imprint4.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002G-distribution-onImprintedRegions4.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted4Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002H-distribution-onImprintedRegions4.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted4Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  ########################## ImprintedRegions 5
  diffImprinted5Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiffDMR_5, "GRanges"),
                                                              feature=imprint5.obj_g$ImprintedRegions,  flank=imprint5.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  sink( file=paste(path2_5, "1002I-distribution-onImprintedRegions5.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted5Ann_2two_sub2_hypo)
  sink()
  pdf( file=paste(path2_5, "1002J-distribution-onImprintedRegions5.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted5Ann_2two_sub2_hypo, precedence=TRUE, main="Annotation on Imprinted Regions"))
  dev.off()
  
  
  
  ########################## H3K4me1 peaks
  Ann_H3K4me1 = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                         feature = region_1_H3K4me1.obj_g$H3K4me1,  flank = region_1_H3K4me1.obj_g$shores,
                                         feature.name = "H3K4me1",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "1A-distribution-on-H3K4me1.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K4me1,  percentage=TRUE)
  print(Ann_H3K4me1)
  sink()
  pdf( file=paste(path2_5, "1B-distribution-on-H3K4me1.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K4me1, precedence=TRUE, main="Annotation on H3K4me1 peaks"))
  dev.off()
  
  
  ########################## H3K4me3 peaks
  Ann_H3K4me3 = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                         feature = region_2_H3K4me3.obj_g$H3K4me3,  flank = region_2_H3K4me3.obj_g$shores,
                                         feature.name = "H3K4me3",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "2A-distribution-on-H3K4me3.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K4me3,  percentage=TRUE)
  print(Ann_H3K4me3)
  sink()
  pdf( file=paste(path2_5, "2B-distribution-on-H3K4me3.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K4me3, precedence=TRUE, main="Annotation on H3K4me3 peaks"))
  dev.off()
  
  
  
  ########################## H3K27ac peaks
  Ann_H3K27ac = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                         feature = region_3_H3K27ac.obj_g$H3K27ac,  flank = region_3_H3K27ac.obj_g$shores,
                                         feature.name = "H3K27ac",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "3A-distribution-on-H3K27ac.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K27ac,  percentage=TRUE)
  print(Ann_H3K27ac)
  sink()
  pdf( file=paste(path2_5, "3B-distribution-on-H3K27ac.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K27ac, precedence=TRUE, main="Annotation on H3K27ac peaks"))
  dev.off()
  
  
  ########################## H3K27me3 peaks
  Ann_H3K27me3 = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                          feature = region_4_H3K27me3.obj_g$H3K27me3,  flank = region_4_H3K27me3.obj_g$shores,
                                          feature.name = "H3K27me3",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "4A-distribution-on-H3K27me3.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_H3K27me3,  percentage=TRUE)
  print(Ann_H3K27me3)
  sink()
  pdf( file=paste(path2_5, "4B-distribution-on-H3K27me3.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_H3K27me3, precedence=TRUE, main="Annotation on H3K27me3 peaks"))
  dev.off()
  
  
  
  ########################## 380imprint regions
  Ann_380imprint = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_5_380imprint.obj_g$imprint,  flank = region_5_380imprint.obj_g$shores,
                                            feature.name = "imprint",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "5A-distribution-on-380imprint.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_380imprint,  percentage=TRUE)
  print(Ann_380imprint)
  sink()
  pdf( file=paste(path2_5, "5B-distribution-on-380imprint.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_380imprint, precedence=TRUE, main="Annotation on 380imprint regions"))
  dev.off()
  
  
  
  ########################## centromere regions
  Ann_centromere = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_6_centromeres.obj_g$centromere,  flank = region_6_centromeres.obj_g$shores,
                                            feature.name = "centromere",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "6A-distribution-on-centromere.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_centromere,  percentage=TRUE)
  print(Ann_centromere)
  sink()
  pdf( file=paste(path2_5, "6B-distribution-on-centromere.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_centromere, precedence=TRUE, main="Annotation on centromere regions"))
  dev.off()
  
  
  
  ########################## CpGislands regions
  Ann_CpGislands = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_7_CpGislands.obj_g$CpGislands,  flank = region_7_CpGislands.obj_g$shores,
                                            feature.name = "CpGislands",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "7A-distribution-on-CpGislands.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_CpGislands,  percentage=TRUE)
  print(Ann_CpGislands)
  sink()
  pdf( file=paste(path2_5, "7B-distribution-on-CpGislands.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_CpGislands, precedence=TRUE, main="Annotation on CpGislands regions"))
  dev.off()
  
  
  ########################## activeEnhancers 
  Ann_activeEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_8A_activeEn.obj_g$activeEnhancers,  flank = region_8A_activeEn.obj_g$shores,
                                                 feature.name = "activeEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "8A-distribution-on-activeEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_activeEnhancers,  percentage=TRUE)
  print(Ann_activeEnhancers)
  sink()
  pdf( file=paste(path2_5, "8A-distribution-on-activeEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_activeEnhancers, precedence=TRUE, main="Annotation on activeEnhancers"))
  dev.off()
  
  
  
  ########################## otherEnhancers 
  Ann_otherEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                feature = region_8B_otherEn.obj_g$otherEnhancers,  flank = region_8B_otherEn.obj_g$shores,
                                                feature.name = "otherEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "8B-distribution-on-otherEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_otherEnhancers,  percentage=TRUE)
  print(Ann_otherEnhancers)
  sink()
  pdf( file=paste(path2_5, "8B-distribution-on-otherEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_otherEnhancers, precedence=TRUE, main="Annotation on otherEnhancers"))
  dev.off()
  
  ########################## posiedEnhancers 
  Ann_posiedEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_8C_posiedEn.obj_g$posiedEnhancers,  flank = region_8C_posiedEn.obj_g$shores,
                                                 feature.name = "posiedEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "8C-distribution-on-posiedEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_posiedEnhancers,  percentage=TRUE)
  print(Ann_posiedEnhancers)
  sink()
  pdf( file=paste(path2_5, "8C-distribution-on-posiedEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_posiedEnhancers, precedence=TRUE, main="Annotation on posiedEnhancers"))
  dev.off()
  
  
  
  ########################## primedEnhancers 
  Ann_primedEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_8D_primedEn.obj_g$primedEnhancers,  flank = region_8D_primedEn.obj_g$shores,
                                                 feature.name = "primedEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "8D-distribution-on-primedEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_primedEnhancers,  percentage=TRUE)
  print(Ann_primedEnhancers)
  sink()
  pdf( file=paste(path2_5, "8D-distribution-on-primedEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_primedEnhancers, precedence=TRUE, main="Annotation on primedEnhancers"))
  dev.off()
  
  
  ########################## HouseKeeping genes
  Ann_HouseKeeping = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                              feature = region_9_HouseKeeping.obj_g$HouseKeeping,  flank = region_9_HouseKeeping.obj_g$shores,
                                              feature.name = "HouseKeeping",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "9A-distribution-on-HouseKeeping.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_HouseKeeping,  percentage=TRUE)
  print(Ann_HouseKeeping)
  sink()
  pdf( file=paste(path2_5, "9B-distribution-on-HouseKeeping.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_HouseKeeping, precedence=TRUE, main="Annotation on HouseKeeping genes"))
  dev.off()
  
  
  
  ########################## rRNA_Genes genes
  Ann_rRNA_Genes = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                            feature = region_10_rRNA.obj_g$rRNA_Genes,  flank = region_10_rRNA.obj_g$shores,
                                            feature.name = "rRNA_Genes",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "10A-distribution-on-rRNA_Genes.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_rRNA_Genes,  percentage=TRUE)
  print(Ann_rRNA_Genes)
  sink()
  pdf( file=paste(path2_5, "10B-distribution-on-rRNA_Genes.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_rRNA_Genes, precedence=TRUE, main="Annotation on rRNA_Genes genes"))
  dev.off()
  


  ########################## activeEnhancers 
  Ann_activeEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_11A_activeEn.obj_g$activeEnhancers,  flank = region_11A_activeEn.obj_g$shores,
                                                 feature.name = "activeEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11A-distribution-on-activeEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_activeEnhancers,  percentage=TRUE)
  print(Ann_activeEnhancers)
  sink()
  pdf( file=paste(path2_5, "11A-distribution-on-activeEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_activeEnhancers, precedence=TRUE, main="Annotation on activeEnhancers"))
  dev.off()
  
  
  
  ########################## otherEnhancers 
  Ann_otherEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                feature = region_11B_otherEn.obj_g$otherEnhancers,  flank = region_11B_otherEn.obj_g$shores,
                                                feature.name = "otherEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11B-distribution-on-otherEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_otherEnhancers,  percentage=TRUE)
  print(Ann_otherEnhancers)
  sink()
  pdf( file=paste(path2_5, "11B-distribution-on-otherEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_otherEnhancers, precedence=TRUE, main="Annotation on otherEnhancers"))
  dev.off()
  
  ########################## posiedEnhancers 
  Ann_posiedEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_11C_posiedEn.obj_g$posiedEnhancers,  flank = region_11C_posiedEn.obj_g$shores,
                                                 feature.name = "posiedEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11C-distribution-on-posiedEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_posiedEnhancers,  percentage=TRUE)
  print(Ann_posiedEnhancers)
  sink()
  pdf( file=paste(path2_5, "11C-distribution-on-posiedEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_posiedEnhancers, precedence=TRUE, main="Annotation on posiedEnhancers"))
  dev.off()
  
  
  
  ########################## primedEnhancers 
  Ann_primedEnhancers = annotateWithFeatureFlank(target = as(myDiffDMR_5, "GRanges"),
                                                 feature = region_11D_primedEn.obj_g$primedEnhancers,  flank = region_11D_primedEn.obj_g$shores,
                                                 feature.name = "primedEnhancers",  flank.name = "shores"
  )
  sink( file=paste(path2_5, "11D-distribution-on-primedEnhancers.txt", sep="/")   )
  getFeatsWithTargetsStats(Ann_primedEnhancers,  percentage=TRUE)
  print(Ann_primedEnhancers)
  sink()
  pdf( file=paste(path2_5, "11D-distribution-on-primedEnhancers.pdf", sep="/")   )
  print(plotTargetAnnotation(Ann_primedEnhancers, precedence=TRUE, main="Annotation on primedEnhancers"))
  dev.off()
  
 
}  




#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myDiff_DMC_DMR_g <- function(  methobj2,   path2  ) {
  sink(paste(path2,  "0_methobj2-for-diffMe.txt",  sep="/"))
  print(methobj2)
  sink() 
  print("################## start calculateDiffMeth:")
  myDiff_2two_sub2 = calculateDiffMeth(methobj2,  num.cores=16 )   ##  If you have replicates, the function will automatically use logistic regression.
  print("################## End calculateDiffMeth") 
  dim(myDiff_2two_sub2)
  names(myDiff_2two_sub2)
  head(myDiff_2two_sub2)
  
  myQvalue_2two_sub2 = myDiff_2two_sub2$qvalue
  myMethDi_2two_sub2 = myDiff_2two_sub2$meth.diff
  
  pdf(paste(path2,  "1A_qvalue_distribution.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=FALSE) )   
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=TRUE ) )   
  dev.off() 

  myQvalue_2two_sub2_log10 <- -log10(myQvalue_2two_sub2)
  pdf(paste(path2,  "1B_qvalue_distribution_log10.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=TRUE ) )
  dev.off() 
  
  pdf(paste(path2,  "1C_methDiff_distribution.pdf",  sep="/"))
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=TRUE ) )
  dev.off() 
             
  qvalue_cutoff = 0.001
  methDiff_cutoff = 10
  
  myColor1_2two_sub2 <- rep( "no",   times= length(myQvalue_2two_sub2) )
  myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
  number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.01
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.05
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.1
    methDiff_cutoff = 5
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.5
    methDiff_cutoff = 1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  
  print("##############################")
  print("##############################")
  print("The final parameters:")
  print(number_yes)
  print(qvalue_cutoff)
  print(methDiff_cutoff)
  print("##############################")
  print("##############################")

  
  DataFrame2_2two_sub2 <- data.frame(myx1 = myMethDi_2two_sub2,   
                                     myy1 =  -log10(myQvalue_2two_sub2),  
                                     mycolor1 = myColor1_2two_sub2 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + ylim(0, 200)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  + ylim(0, 100)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)   + ylim(0, 30)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )

  
  
  myDiff25p.hypo_2two_sub2  = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hypo" )  ## less enrich in ART
  myDiff25p.hyper_2two_sub2 = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hyper")  ## more enrich in ART
  myDiff25p_2two_sub2       = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff)
  myDiffTemp_2two_sub2      = getMethylDiff(myDiff_2two_sub2, difference=0,  qvalue=0.05)
  
  write.table(myDiff_2two_sub2 , file = paste(path2,"2A_diffMe-allsites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hypo_2two_sub2 , file = paste(path2,"2B_diffMe-hypo.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hyper_2two_sub2 , file = paste(path2,"2C_diffMe-hyper.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p_2two_sub2 , file = paste(path2,"2D_AlldiffMesites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiffTemp_2two_sub2 , file = paste(path2,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  sink( file=paste(path2, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=FALSE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  sink()
  
  pdf( file=paste(path2, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=TRUE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  dev.off()
  

  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hypo_2two_sub2,    path2_5 = paste(path2, "Annotation_Hypo",  sep="/") ),
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hyper_2two_sub2,   path2_5 = paste(path2, "Annotation_Hyper", sep="/") ), 
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p_2two_sub2,         path2_5 = paste(path2, "Annotation_All",   sep="/") ),
        error = function(err){"000"}
  )
 
}


myDiff_overdispersion_Ftest_g <- function(  methobj2,   path2  ) {
  sink(paste(path2,  "0_methobj2-for-diffMe.txt",  sep="/"))
  print(methobj2)
  sink() 
  print("################## start calculateDiffMeth:")
  myDiff_2two_sub2 = calculateDiffMeth(methobj2,  overdispersion = "MN",  test = "F", num.cores=16 ) 
  print("################## End calculateDiffMeth") 
  dim(myDiff_2two_sub2)
  names(myDiff_2two_sub2)
  head(myDiff_2two_sub2)
  
  myQvalue_2two_sub2 = myDiff_2two_sub2$qvalue
  myMethDi_2two_sub2 = myDiff_2two_sub2$meth.diff
    
  pdf(paste(path2,  "1A_qvalue_distribution.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=FALSE) )   
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=TRUE ) )   
  dev.off() 

  myQvalue_2two_sub2_log10 <- -log10(myQvalue_2two_sub2)
  pdf(paste(path2,  "1B_qvalue_distribution_log10.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=TRUE ) )
  dev.off() 
  
  pdf(paste(path2,  "1C_methDiff_distribution.pdf",  sep="/"))
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=TRUE ) )
  dev.off() 
             
  qvalue_cutoff = 0.01
  methDiff_cutoff = 10
  
  myColor1_2two_sub2 <- rep( "no",   times= length(myQvalue_2two_sub2) )
  myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
  number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.05
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.2
    methDiff_cutoff = 5
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.5
    methDiff_cutoff = 1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  
  print("##############################")
  print("##############################")
  print("The final parameters:")
  print(number_yes)
  print(qvalue_cutoff)
  print(methDiff_cutoff)
  print("##############################")
  print("##############################")
  
  
  DataFrame2_2two_sub2 <- data.frame(myx1 = myMethDi_2two_sub2,   
                                     myy1 =  -log10(myQvalue_2two_sub2),  
                                     mycolor1 = myColor1_2two_sub2 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + ylim(0, 200)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  + ylim(0, 10)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)   + ylim(0, 3)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  
  myDiff25p.hypo_2two_sub2  = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hypo" )  ## less enrich in ART
  myDiff25p.hyper_2two_sub2 = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hyper")  ## more enrich in ART
  myDiff25p_2two_sub2       = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff)
  myDiffTemp_2two_sub2      = getMethylDiff(myDiff_2two_sub2, difference=0,  qvalue=0.05)
  
  write.table(myDiff_2two_sub2 , file = paste(path2,"2A_diffMe-allsites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hypo_2two_sub2 , file = paste(path2,"2B_diffMe-hypo.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hyper_2two_sub2 , file = paste(path2,"2C_diffMe-hyper.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p_2two_sub2 , file = paste(path2,"2D_AlldiffMesites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiffTemp_2two_sub2 , file = paste(path2,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  sink( file=paste(path2, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=FALSE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  sink()
  
  pdf( file=paste(path2, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=TRUE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  dev.off()
  

  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hypo_2two_sub2,    path2_5 = paste(path2, "Annotation_Hypo",  sep="/") ),
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hyper_2two_sub2,   path2_5 = paste(path2, "Annotation_Hyper", sep="/") ), 
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p_2two_sub2,         path2_5 = paste(path2, "Annotation_All",   sep="/") ),
        error = function(err){"000"}
  )
 
}


myDiff_overdispersion_Chisq_g <- function(  methobj2,   path2  ) {
  sink(paste(path2,  "0_methobj2-for-diffMe.txt",  sep="/"))
  print(methobj2)
  sink() 
  print("################## start calculateDiffMeth:")
  myDiff_2two_sub2 = calculateDiffMeth(methobj2,  overdispersion = "MN",  test = "Chisq", num.cores=16 )   ##  If you have replicates, the function will automatically use logistic regression.
  print("################## End calculateDiffMeth") 
  dim(myDiff_2two_sub2)
  names(myDiff_2two_sub2)
  head(myDiff_2two_sub2)
  
  myQvalue_2two_sub2 = myDiff_2two_sub2$qvalue
  myMethDi_2two_sub2 = myDiff_2two_sub2$meth.diff
  
  pdf(paste(path2,  "1A_qvalue_distribution.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=FALSE) )   
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=TRUE ) )   
  dev.off() 

  myQvalue_2two_sub2_log10 <- -log10(myQvalue_2two_sub2)
  pdf(paste(path2,  "1B_qvalue_distribution_log10.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=TRUE ) )
  dev.off() 
  
  pdf(paste(path2,  "1C_methDiff_distribution.pdf",  sep="/"))
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=TRUE ) )
  dev.off() 
             
  qvalue_cutoff = 0.01
  methDiff_cutoff = 10
  
  myColor1_2two_sub2 <- rep( "no",   times= length(myQvalue_2two_sub2) )
  myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
  number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.05
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.2
    methDiff_cutoff = 5
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.5
    methDiff_cutoff = 1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  
  print("##############################")
  print("##############################")
  print("The final parameters:")
  print(number_yes)
  print(qvalue_cutoff)
  print(methDiff_cutoff)
  print("##############################")
  print("##############################")
  
  
  DataFrame2_2two_sub2 <- data.frame(myx1 = myMethDi_2two_sub2,   
                                     myy1 =  -log10(myQvalue_2two_sub2),  
                                     mycolor1 = myColor1_2two_sub2 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + ylim(0, 200)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  + ylim(0, 10)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)   + ylim(0, 3)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  
  myDiff25p.hypo_2two_sub2  = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hypo" )  ## less enrich in ART
  myDiff25p.hyper_2two_sub2 = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hyper")  ## more enrich in ART
  myDiff25p_2two_sub2       = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff)
  myDiffTemp_2two_sub2      = getMethylDiff(myDiff_2two_sub2, difference=0,  qvalue=0.05)
  
  write.table(myDiff_2two_sub2 , file = paste(path2,"2A_diffMe-allsites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hypo_2two_sub2 , file = paste(path2,"2B_diffMe-hypo.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hyper_2two_sub2 , file = paste(path2,"2C_diffMe-hyper.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p_2two_sub2 , file = paste(path2,"2D_AlldiffMesites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiffTemp_2two_sub2 , file = paste(path2,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  sink( file=paste(path2, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=FALSE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  sink()
  
  pdf( file=paste(path2, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=TRUE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  dev.off()
  

  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hypo_2two_sub2,    path2_5 = paste(path2, "Annotation_Hypo",  sep="/") ),
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hyper_2two_sub2,   path2_5 = paste(path2, "Annotation_Hyper", sep="/") ), 
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p_2two_sub2,         path2_5 = paste(path2, "Annotation_All",   sep="/") ),
        error = function(err){"000"}
  )
 
}


myDiff_calculateDiffMethDSS_g <- function(  methobj2,   path2  ) {
  sink(paste(path2,  "0_methobj2-for-diffMe.txt",  sep="/"))
  print(methobj2)
  sink() 
  print("################## start calculateDiffMethDSS:")
  myDiff_2two_sub2 = calculateDiffMethDSS(methobj2,  mc.cores=16 )   ##  It calculates the differential methylation statistics using a beta-binomial model with parameter shrink-age.
  print("################## End calculateDiffMethDSS") 
  
  dim(myDiff_2two_sub2)
  names(myDiff_2two_sub2)
  head(myDiff_2two_sub2)
  
  myQvalue_2two_sub2 = myDiff_2two_sub2$qvalue
  myMethDi_2two_sub2 = myDiff_2two_sub2$meth.diff
  
  
  pdf(paste(path2,  "1A_qvalue_distribution.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=FALSE) )   
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1),    freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.2),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.1),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.05), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2,  nclass=20,  xlim=c(0, 0.01), freq=TRUE ) )   
  dev.off() 

  myQvalue_2two_sub2_log10 <- -log10(myQvalue_2two_sub2)
  pdf(paste(path2,  "1B_qvalue_distribution_log10.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=FALSE) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=100, xlim=c(0, 100), freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 10),  freq=TRUE ) )
  print( hist( myQvalue_2two_sub2_log10,  nclass=20,  xlim=c(0, 2),   freq=TRUE ) )
  dev.off() 
  
  pdf(paste(path2,  "1C_methDiff_distribution.pdf",  sep="/"))
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100),  freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=50,  xlim=c(0, 50),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=30,  xlim=c(0, 30),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=20,  xlim=c(0, 20),   freq=TRUE ) )
  print( hist(myMethDi_2two_sub2, nclass=10,  xlim=c(0, 10),   freq=TRUE ) )
  dev.off() 
             
  qvalue_cutoff = 0.01
  methDiff_cutoff = 10
  
  myColor1_2two_sub2 <- rep( "no",   times= length(myQvalue_2two_sub2) )
  myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
  number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.05
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.2
    methDiff_cutoff = 5
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  if( number_yes < 30 ) {
    qvalue_cutoff = 0.5
    methDiff_cutoff = 1
    myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>methDiff_cutoff) & (myQvalue_2two_sub2<qvalue_cutoff) ]  <- "yes"
    number_yes = length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  }
  
  print("##############################")
  print("##############################")
  print("The final parameters:")
  print(number_yes)
  print(qvalue_cutoff)
  print(methDiff_cutoff)
  print("##############################")
  print("##############################")
  
  
  DataFrame2_2two_sub2 <- data.frame(myx1 = myMethDi_2two_sub2,   
                                     myy1 =  -log10(myQvalue_2two_sub2),  
                                     mycolor1 = myColor1_2two_sub2 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + ylim(0, 200)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  + ylim(0, 5)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.5 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)   + ylim(0, 3)  
  ggsave( filename = paste(path2,  "1E_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  
  myDiff25p.hypo_2two_sub2  = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hypo" )  ## less enrich in ART
  myDiff25p.hyper_2two_sub2 = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff, type="hyper")  ## more enrich in ART
  myDiff25p_2two_sub2       = getMethylDiff(myDiff_2two_sub2, difference=methDiff_cutoff, qvalue=qvalue_cutoff)
  myDiffTemp_2two_sub2      = getMethylDiff(myDiff_2two_sub2, difference=0,  qvalue=0.05)
  
  write.table(myDiff_2two_sub2 , file = paste(path2,"2A_diffMe-allsites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hypo_2two_sub2 , file = paste(path2,"2B_diffMe-hypo.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hyper_2two_sub2 , file = paste(path2,"2C_diffMe-hyper.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p_2two_sub2 , file = paste(path2,"2D_AlldiffMesites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiffTemp_2two_sub2 , file = paste(path2,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  sink( file=paste(path2, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=FALSE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  sink()
  
  pdf( file=paste(path2, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=TRUE,  qvalue.cutoff=qvalue_cutoff, meth.cutoff=methDiff_cutoff) )
  dev.off()
  

  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hypo_2two_sub2,    path2_5 = paste(path2, "Annotation_Hypo",  sep="/") ),
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p.hyper_2two_sub2,   path2_5 = paste(path2, "Annotation_Hyper", sep="/") ), 
        error = function(err){"000"}
  )
  tryCatch(
        myDMRs_annotation_g(  myDiffDMR_5 = myDiff25p_2two_sub2,         path2_5 = paste(path2, "Annotation_All",   sep="/") ),
        error = function(err){"000"}
  )
 
}





#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myMainFunction_1_g  <- function(  myobj_temp1,   path_temp1,   binSize_temp1, binBases_temp1, dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }

  myobj_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  write.table(dataFrame_temp1 , 
              file = paste(path_temp1_sub1,   "0_dataFrame.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(myobj_2two) )
  print( getSampleID(myobj_2two) )
  print( getTreatment(myobj_2two) )
  sink()
  
  tiles_2two = tileMethylCounts( myobj_2two,   win.size=binSize_temp1,   step.size=binSize_temp1,   cov.bases = binBases_temp1  )   
  meth_2two  = unite( tiles_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2A_dimensions-tiles-eachSample.txt", sep="/")  )
  print( tiles_2two )
  sink()

  sink( file=paste(path_temp1_sub1, "2B_only-dimensions-tiles-eachSample.txt", sep="/")  )
  for( i in c(1:length(tiles_2two)) ) {
    print(   dim(tiles_2two[[i]])  )
  }
  sink()
  
  sink( file=paste(path_temp1_sub1 , "2C_dimensions-tiles-merged.txt", sep="/")  )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
  
  path_temp1_sub1_1 = paste(path_temp1_sub1, "tile_methylation", sep="/")
  if( ! file.exists(path_temp1_sub1_1) ) { dir.create(path_temp1_sub1_1, recursive = TRUE) }
  for( i in c(1:length(tiles_2two)) ) {
    file_name = Files_All_vector_g[i]
    file_name = gsub(inputDir_g, i, file_name)
    file_name = gsub("/", "__", file_name)
    write.table( tiles_2two[[i]] , 
                file = paste(path_temp1_sub1_1,   file_name,  sep="/"), 
                append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
                row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  }

  
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(getData(meth_2two)[,1:4] , 
              file = paste(path_temp1_sub1,   "3C_regions-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf( file=paste(path_temp1_sub1, "4A_MethylationStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getMethylationStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4B_MethylationStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste(as.vector(dataFrame_temp1$mysampleID)[i],  ":", sep="") )
    print( getMethylationStats( tiles_2two[[i]] )  )
  }
  sink()
  
  pdf( file=paste(path_temp1_sub1, "5A_CoverageStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getCoverageStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  sink( file=paste(path_temp1_sub1, "5B_CoverageStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste( as.vector(dataFrame_temp1$mysampleID)[i],   ":", sep="") )
    print( getCoverageStats( tiles_2two[[i]] )  )
  }
  sink()
  
  #sink( file=paste(path_temp1_sub1, "6A_Correlation-tiles.txt", sep="/")  )
  #pdf( file=paste(path_temp1_sub1, "6A_Correlation-tiles.pdf", sep="/")  )
  #getCorrelation(meth_2two, method = "pearson",   plot=TRUE  )
  #getCorrelation(meth_2two, method = "spearman",  plot=TRUE  )
  #dev.off()
  #sink()
  
  sink( file=paste(path_temp1_sub1, "6B_pearsonCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "pearson",   plot=FALSE  )
  sink()
  
  sink( file=paste(path_temp1_sub1, "6C_spearmanCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "spearman",  plot=FALSE  )
  sink()
  
  
  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  tryCatch(
      MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=30,   height2=10 ),
      error = function(err){"000"}
  )
  
  path_temp1_sub3 = paste(path_temp1, "3_PCA_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub3) ) { dir.create(path_temp1_sub3, recursive = TRUE) }
  #tryCatch(
  #    MyPCA_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub3,     file2="PCA_byMethylKit_",   width2=10,   height2=10 ),
  #    error = function(err){"000"}
  #)
  
  path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  #tryCatch(
  #    MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=10,   height2=10,  dataFrame_temp2=dataFrame_temp1   ),
  #    error = function(err){"000"}
  #)
  
  path_temp1_sub5 = paste(path_temp1, "5_PCA_DirectlyByPrcomp", sep="/")
  if( ! file.exists(path_temp1_sub5) ) { dir.create(path_temp1_sub5, recursive = TRUE) }
  PCA_2two_sub5 <- prcomp( t(mat_2two) )
  tryCatch(
      MyPrcompObj_1_g(  prcompObj2=PCA_2two_sub5,   path2=path_temp1_sub5,   file2="PCA_DirectlyByPrcomp",  dataFrame_temp2=dataFrame_temp1   ) ,
      error = function(err){"000"}
  )
  
  path_temp1_sub6 = paste(path_temp1, "6_PCA_byFactoMineR", sep="/")
  if( ! file.exists(path_temp1_sub6) ) { dir.create(path_temp1_sub6, recursive = TRUE) }
  PCA_2two_sub6 <- PCA( t(mat_2two) , graph=FALSE)
  #tryCatch(
  #    MyPCAobj_FactoMineR_g(  PCAobj2=PCA_2two_sub6,   path2=path_temp1_sub6,   file2="PCA_byFactoMineR",  dataFrame_temp2=dataFrame_temp1   ),
  #    error = function(err){"000"}
  #)
                                                                     
  path_temp1_sub7 = paste(path_temp1, "7_HierarchicalClustering", sep="/")
  if( ! file.exists(path_temp1_sub7) ) { dir.create(path_temp1_sub7, recursive = TRUE) }  
  tryCatch(
        myHierarchicalClustering_1_g(  mat_3three=mat_2two,   path_temp1=path_temp1_sub7,   width1=30, height1=10  ) ,
        error = function(err){"000"}
  )
                                                                    
  path_temp1_sub8 = paste(path_temp1, "8_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  tryCatch(
        myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  ) ,
        error = function(err){"000"}
  )
 
  
  path_temp1_sub9 = paste(path_temp1, "9_DMR_overdispersion_Ftest", sep="/")
  if( ! file.exists(path_temp1_sub9) ) { dir.create(path_temp1_sub9, recursive = TRUE) }
  #myDiff_overdispersion_Ftest_g(  methobj2=meth_2two,   path2=path_temp1_sub9  )

  
  path_temp1_sub10 = paste(path_temp1, "10_DMR_overdispersion_Chisq", sep="/")
  if( ! file.exists(path_temp1_sub10) ) { dir.create(path_temp1_sub10, recursive = TRUE) }
  #myDiff_overdispersion_Chisq_g(  methobj2=meth_2two,   path2=path_temp1_sub10  )

  
  path_temp1_sub11 = paste(path_temp1, "11_calculateDiffMethDSS", sep="/")
  if( ! file.exists(path_temp1_sub11) ) { dir.create(path_temp1_sub11, recursive = TRUE) }
  #myDiff_calculateDiffMethDSS_g(  methobj2=meth_2two,   path2=path_temp1_sub11  )
  
} 




#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myMainFunction_2_g  <- function(  myobj_temp1,   path_temp1,   binSize_temp1, binBases_temp1, dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }

  myobj_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  write.table(dataFrame_temp1 , 
              file = paste(path_temp1_sub1,   "0_dataFrame.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(myobj_2two) )
  print( getSampleID(myobj_2two) )
  print( getTreatment(myobj_2two) )
  sink()
  
  tiles_2two = tileMethylCounts( myobj_2two,   win.size=binSize_temp1,   step.size=binSize_temp1,   cov.bases = binBases_temp1  )   
  meth_2two  = unite( tiles_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2A_dimensions-tiles-eachSample.txt", sep="/")  )
  print( tiles_2two )
  sink()

  sink( file=paste(path_temp1_sub1, "2B_only-dimensions-tiles-eachSample.txt", sep="/")  )
  for( i in c(1:length(tiles_2two)) ) {
    print(   dim(tiles_2two[[i]])  )
  }
  sink()
  
  sink( file=paste(path_temp1_sub1 , "2C_dimensions-tiles-merged.txt", sep="/")  )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(getData(meth_2two)[,1:4] , 
              file = paste(path_temp1_sub1,   "3C_regions-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf( file=paste(path_temp1_sub1, "4A_MethylationStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getMethylationStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4B_MethylationStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste(as.vector(dataFrame_temp1$mysampleID)[i],  ":", sep="") )
    print( getMethylationStats( tiles_2two[[i]] )  )
  }
  sink()
  
  pdf( file=paste(path_temp1_sub1, "5A_CoverageStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getCoverageStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  sink( file=paste(path_temp1_sub1, "5B_CoverageStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste( as.vector(dataFrame_temp1$mysampleID)[i],   ":", sep="") )
    print( getCoverageStats( tiles_2two[[i]] )  )
  }
  sink()
  
  #sink( file=paste(path_temp1_sub1, "6A_Correlation-tiles.txt", sep="/")  )
  #pdf( file=paste(path_temp1_sub1, "6A_Correlation-tiles.pdf", sep="/")  )
  #getCorrelation(meth_2two, method = "pearson",   plot=TRUE  )
  #getCorrelation(meth_2two, method = "spearman",  plot=TRUE  )
  #dev.off()
  #sink()
  
  sink( file=paste(path_temp1_sub1, "6B_pearsonCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "pearson",   plot=FALSE  )
  sink()
  
  sink( file=paste(path_temp1_sub1, "6C_spearmanCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "spearman",  plot=FALSE  )
  sink()
  
  
  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  tryCatch(
      MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=30,   height2=10 ),
      error = function(err){"000"}
  )
  
  path_temp1_sub3 = paste(path_temp1, "3_PCA_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub3) ) { dir.create(path_temp1_sub3, recursive = TRUE) }
  MyPCA_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub3,     file2="PCA_byMethylKit_",   width2=10,   height2=10 )
  
  path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=10,   height2=10,  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub5 = paste(path_temp1, "5_PCA_DirectlyByPrcomp", sep="/")
  if( ! file.exists(path_temp1_sub5) ) { dir.create(path_temp1_sub5, recursive = TRUE) }
  PCA_2two_sub5 <- prcomp( t(mat_2two) )
  tryCatch(
      MyPrcompObj_1_g(  prcompObj2=PCA_2two_sub5,   path2=path_temp1_sub5,   file2="PCA_DirectlyByPrcomp",  dataFrame_temp2=dataFrame_temp1   ) ,
      error = function(err){"000"}
  )
  
  path_temp1_sub6 = paste(path_temp1, "6_PCA_byFactoMineR", sep="/")
  if( ! file.exists(path_temp1_sub6) ) { dir.create(path_temp1_sub6, recursive = TRUE) }
  PCA_2two_sub6 <- PCA( t(mat_2two) , graph=FALSE)
  MyPCAobj_FactoMineR_g(  PCAobj2=PCA_2two_sub6,   path2=path_temp1_sub6,   file2="PCA_byFactoMineR",  dataFrame_temp2=dataFrame_temp1   )
                                                                   
  path_temp1_sub7 = paste(path_temp1, "7_HierarchicalClustering", sep="/")
  if( ! file.exists(path_temp1_sub7) ) { dir.create(path_temp1_sub7, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=mat_2two,   path_temp1=path_temp1_sub7,   width1=30, height1=10  )
                                                                    
  path_temp1_sub8 = paste(path_temp1, "8_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )

  
  path_temp1_sub9 = paste(path_temp1, "9_DMR_overdispersion_Ftest", sep="/")
  if( ! file.exists(path_temp1_sub9) ) { dir.create(path_temp1_sub9, recursive = TRUE) }
  myDiff_overdispersion_Ftest_g(  methobj2=meth_2two,   path2=path_temp1_sub9  )

  
  path_temp1_sub10 = paste(path_temp1, "10_DMR_overdispersion_Chisq", sep="/")
  if( ! file.exists(path_temp1_sub10) ) { dir.create(path_temp1_sub10, recursive = TRUE) }
  myDiff_overdispersion_Chisq_g(  methobj2=meth_2two,   path2=path_temp1_sub10  )

  
  path_temp1_sub11 = paste(path_temp1, "11_calculateDiffMethDSS", sep="/")
  if( ! file.exists(path_temp1_sub11) ) { dir.create(path_temp1_sub11, recursive = TRUE) }
  myDiff_calculateDiffMethDSS_g(  methobj2=meth_2two,   path2=path_temp1_sub11  )
  
} 



#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myMainFunction1bp_1_g  <- function(  myobj_temp1,   path_temp1,   dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }
  
  myobj_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  write.table(dataFrame_temp1 , 
              file = paste(path_temp1_sub1,   "0_dataFrame.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(myobj_2two) )
  print( getSampleID(myobj_2two) )
  print( getTreatment(myobj_2two) )
  sink()
  
  meth_2two  = unite( myobj_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2A_dimensions-tiles-eachSample.txt", sep="/")  )
  print( myobj_2two )
  sink()
  
  sink( file=paste(path_temp1_sub1, "2B_only-dimensions-tiles-eachSample.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print(   dim(myobj_2two[[i]])  )
  }
  sink()
  
  sink( file=paste(path_temp1_sub1 , "2C_dimensions-tiles-merged.txt", sep="/")  )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
  
  path_temp1_sub1_1 = paste(path_temp1_sub1, "tile_methylation", sep="/")
  if( ! file.exists(path_temp1_sub1_1) ) { dir.create(path_temp1_sub1_1, recursive = TRUE) }
  for( i in c(1:length(myobj_2two)) ) {
    file_name = Files_All_vector_g[i]
    file_name = gsub(inputDir_g, i, file_name)
    file_name = gsub("/", "__", file_name)
    write.table( myobj_2two[[i]] , 
                 file = paste(path_temp1_sub1_1,   file_name,  sep="/"), 
                 append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
                 row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  }
  
  
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(getData(meth_2two)[,1:4] , 
              file = paste(path_temp1_sub1,   "3C_regions-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf( file=paste(path_temp1_sub1, "4A_MethylationStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getMethylationStats(myobj_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4B_MethylationStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste(as.vector(dataFrame_temp1$mysampleID)[i],  ":", sep="") )
    print( getMethylationStats( myobj_2two[[i]] )  )
  }
  sink()
  
  pdf( file=paste(path_temp1_sub1, "5A_CoverageStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getCoverageStats(myobj_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  sink( file=paste(path_temp1_sub1, "5B_CoverageStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste( as.vector(dataFrame_temp1$mysampleID)[i],   ":", sep="") )
    print( getCoverageStats( myobj_2two[[i]] )  )
  }
  sink()
  
  #sink( file=paste(path_temp1_sub1, "6A_Correlation-tiles.txt", sep="/")  )
  #pdf( file=paste(path_temp1_sub1, "6A_Correlation-tiles.pdf", sep="/")  )
  #getCorrelation(meth_2two, method = "pearson",   plot=TRUE  )
  #getCorrelation(meth_2two, method = "spearman",  plot=TRUE  )
  #dev.off()
  #sink()
  
  sink( file=paste(path_temp1_sub1, "6B_pearsonCorrelation-tiles.txt", sep="/")  )
  #getCorrelation(meth_2two, method = "pearson",   plot=FALSE  )
  sink()
  
  sink( file=paste(path_temp1_sub1, "6C_spearmanCorrelation-tiles.txt", sep="/")  )
  #getCorrelation(meth_2two, method = "spearman",  plot=FALSE  )
  sink()
  
  
  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=30,   height2=10 )
  
  path_temp1_sub3 = paste(path_temp1, "3_PCA_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub3) ) { dir.create(path_temp1_sub3, recursive = TRUE) }
  #MyPCA_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub3,     file2="PCA_byMethylKit_",   width2=10,   height2=10 )
  
  path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  #MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=10,   height2=10,  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub5 = paste(path_temp1, "5_PCA_DirectlyByPrcomp", sep="/")
  if( ! file.exists(path_temp1_sub5) ) { dir.create(path_temp1_sub5, recursive = TRUE) }
  PCA_2two_sub5 <- prcomp( t(mat_2two) )
  tryCatch(
      MyPrcompObj_1_g(  prcompObj2=PCA_2two_sub5,   path2=path_temp1_sub5,   file2="PCA_DirectlyByPrcomp",  dataFrame_temp2=dataFrame_temp1   ) ,
      error = function(err){"000"}
  )
  
  path_temp1_sub6 = paste(path_temp1, "6_PCA_byFactoMineR", sep="/")
  if( ! file.exists(path_temp1_sub6) ) { dir.create(path_temp1_sub6, recursive = TRUE) }
  PCA_2two_sub6 <- PCA( t(mat_2two) , graph=FALSE)
  #MyPCAobj_FactoMineR_g(  PCAobj2=PCA_2two_sub6,   path2=path_temp1_sub6,   file2="PCA_byFactoMineR",  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub7 = paste(path_temp1, "7_HierarchicalClustering", sep="/")
  if( ! file.exists(path_temp1_sub7) ) { dir.create(path_temp1_sub7, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=mat_2two,   path_temp1=path_temp1_sub7,   width1=30, height1=10  )
  
  path_temp1_sub8 = paste(path_temp1, "8_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )
  
  
  path_temp1_sub9 = paste(path_temp1, "9_DMR_overdispersion_Ftest", sep="/")
  if( ! file.exists(path_temp1_sub9) ) { dir.create(path_temp1_sub9, recursive = TRUE) }
  #myDiff_overdispersion_Ftest_g(  methobj2=meth_2two,   path2=path_temp1_sub9  )
  
  
  path_temp1_sub10 = paste(path_temp1, "10_DMR_overdispersion_Chisq", sep="/")
  if( ! file.exists(path_temp1_sub10) ) { dir.create(path_temp1_sub10, recursive = TRUE) }
  #myDiff_overdispersion_Chisq_g(  methobj2=meth_2two,   path2=path_temp1_sub10  )
  
  
  path_temp1_sub11 = paste(path_temp1, "11_calculateDiffMethDSS", sep="/")
  if( ! file.exists(path_temp1_sub11) ) { dir.create(path_temp1_sub11, recursive = TRUE) }
  #myDiff_calculateDiffMethDSS_g(  methobj2=meth_2two,   path2=path_temp1_sub11  )
  
} 





#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myMainFunction1bp_2_g  <- function(  myobj_temp1,   path_temp1,   dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }
  
  myobj_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  write.table(dataFrame_temp1 , 
              file = paste(path_temp1_sub1,   "0_dataFrame.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(myobj_2two) )
  print( getSampleID(myobj_2two) )
  print( getTreatment(myobj_2two) )
  sink()
  
  meth_2two  = unite( myobj_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2A_dimensions-tiles-eachSample.txt", sep="/")  )
  print( myobj_2two )
  sink()
  
  sink( file=paste(path_temp1_sub1, "2B_only-dimensions-tiles-eachSample.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print(   dim(myobj_2two[[i]])  )
  }
  sink()
  
  sink( file=paste(path_temp1_sub1 , "2C_dimensions-tiles-merged.txt", sep="/")  )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
  
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(getData(meth_2two)[,1:4] , 
              file = paste(path_temp1_sub1,   "3C_regions-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf( file=paste(path_temp1_sub1, "4A_MethylationStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getMethylationStats(myobj_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4B_MethylationStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste(as.vector(dataFrame_temp1$mysampleID)[i],  ":", sep="") )
    print( getMethylationStats( myobj_2two[[i]] )  )
  }
  sink()
  
  pdf( file=paste(path_temp1_sub1, "5A_CoverageStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getCoverageStats(myobj_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  sink( file=paste(path_temp1_sub1, "5B_CoverageStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste( as.vector(dataFrame_temp1$mysampleID)[i],   ":", sep="") )
    print( getCoverageStats( myobj_2two[[i]] )  )
  }
  sink()
  
  #sink( file=paste(path_temp1_sub1, "6A_Correlation-tiles.txt", sep="/")  )
  #pdf( file=paste(path_temp1_sub1, "6A_Correlation-tiles.pdf", sep="/")  )
  #getCorrelation(meth_2two, method = "pearson",   plot=TRUE  )
  #getCorrelation(meth_2two, method = "spearman",  plot=TRUE  )
  #dev.off()
  #sink()
  
  sink( file=paste(path_temp1_sub1, "6B_pearsonCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "pearson",   plot=FALSE  )
  sink()
  
  sink( file=paste(path_temp1_sub1, "6C_spearmanCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "spearman",  plot=FALSE  )
  sink()
  
  
  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=30,   height2=10 )
  
  path_temp1_sub3 = paste(path_temp1, "3_PCA_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub3) ) { dir.create(path_temp1_sub3, recursive = TRUE) }
  MyPCA_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub3,     file2="PCA_byMethylKit_",   width2=10,   height2=10 )
  
  path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=10,   height2=10,  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub5 = paste(path_temp1, "5_PCA_DirectlyByPrcomp", sep="/")
  if( ! file.exists(path_temp1_sub5) ) { dir.create(path_temp1_sub5, recursive = TRUE) }
  PCA_2two_sub5 <- prcomp( t(mat_2two) )
  tryCatch(
      MyPrcompObj_1_g(  prcompObj2=PCA_2two_sub5,   path2=path_temp1_sub5,   file2="PCA_DirectlyByPrcomp",  dataFrame_temp2=dataFrame_temp1   ) ,
      error = function(err){"000"}
  )
  
  path_temp1_sub6 = paste(path_temp1, "6_PCA_byFactoMineR", sep="/")
  if( ! file.exists(path_temp1_sub6) ) { dir.create(path_temp1_sub6, recursive = TRUE) }
  PCA_2two_sub6 <- PCA( t(mat_2two) , graph=FALSE)
  MyPCAobj_FactoMineR_g(  PCAobj2=PCA_2two_sub6,   path2=path_temp1_sub6,   file2="PCA_byFactoMineR",  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub7 = paste(path_temp1, "7_HierarchicalClustering", sep="/")
  if( ! file.exists(path_temp1_sub7) ) { dir.create(path_temp1_sub7, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=mat_2two,   path_temp1=path_temp1_sub7,   width1=30, height1=10  )
  
  path_temp1_sub8 = paste(path_temp1, "8_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )
  
  
  path_temp1_sub9 = paste(path_temp1, "9_DMR_overdispersion_Ftest", sep="/")
  if( ! file.exists(path_temp1_sub9) ) { dir.create(path_temp1_sub9, recursive = TRUE) }
  myDiff_overdispersion_Ftest_g(  methobj2=meth_2two,   path2=path_temp1_sub9  )
  
  
  path_temp1_sub10 = paste(path_temp1, "10_DMR_overdispersion_Chisq", sep="/")
  if( ! file.exists(path_temp1_sub10) ) { dir.create(path_temp1_sub10, recursive = TRUE) }
  myDiff_overdispersion_Chisq_g(  methobj2=meth_2two,   path2=path_temp1_sub10  )
  
  
  path_temp1_sub11 = paste(path_temp1, "11_calculateDiffMethDSS", sep="/")
  if( ! file.exists(path_temp1_sub11) ) { dir.create(path_temp1_sub11, recursive = TRUE) }
  myDiff_calculateDiffMethDSS_g(  methobj2=meth_2two,   path2=path_temp1_sub11  )
  
} 

##################################################################################################################





##################################################################################################################
## Annotating differentially methylated bases or regions
myRefSeqGenes_g       = "/public3/yongpeng/AnnotationBED/hg38/otherRegions/RefSeq_Genes.bed"
myRepeats_g           = "/public3/yongpeng/AnnotationBED/hg38/otherRegions/Repeats_rmsk.bed"
myImprintedRegions1_g = "/public3/yongpeng/AnnotationBED/hg38/imprintedRegions/67.Regions.PlosGenetics.ImprintedGenes.hg38.bed"
myImprintedRegions2_g = "/public3/yongpeng/AnnotationBED/hg38/imprintedRegions/75Regions.GR.hg38.bed"
myImprintedRegions3_g = "/public3/yongpeng/AnnotationBED/hg38/imprintedRegions/369Regions.GR.hg38.bed"
myImprintedRegions4_g = "/public3/yongpeng/AnnotationBED/hg38/imprintedRegions/merge1.imprintedRegions.hg38.bed"
myImprintedRegions5_g = "/public3/yongpeng/AnnotationBED/hg38/imprintedRegions/merge2.imprintedRegions.hg38.bed"
region_1_H3K4me1_g        = "/public3/yongpeng/AnnotationBED/hg38/1-H3K4me1.bed"
region_2_H3K4me3_g        = "/public3/yongpeng/AnnotationBED/hg38/2-H3K4me3.bed"
region_3_H3K27ac_g        = "/public3/yongpeng/AnnotationBED/hg38/3-H3K27ac.bed"
region_4_H3K27me3_g       = "/public3/yongpeng/AnnotationBED/hg38/4-H3K27me3.bed"
region_5_380imprint_g     = "/public3/yongpeng/AnnotationBED/hg38/5-380imprint.bed"
region_6_centromeres_g    = "/public3/yongpeng/AnnotationBED/hg38/6-centromeres.bed"
region_7_CpGislands_g     = "/public3/yongpeng/AnnotationBED/hg38/7-CpGislands.bed"
region_8A_activeEn_g      = "/public3/yongpeng/AnnotationBED/hg38/8A-enhancers-active-H3K4me1-H3K27ac.bed"
region_8B_otherEn_g       = "/public3/yongpeng/AnnotationBED/hg38/8B-enhancers-H3K4me1-H3K27ac-H3K27me3.bed"
region_8C_posiedEn_g      = "/public3/yongpeng/AnnotationBED/hg38/8C-enhancers-posied-H3K4me1-H3K27me3.bed"
region_8D_primedEn_g      = "/public3/yongpeng/AnnotationBED/hg38/8D-enhancers-primed-H3K4me1.bed"
region_9_HouseKeeping_g   = "/public3/yongpeng/AnnotationBED/hg38/9-HouseKeepingGenes.bed"
region_10_rRNA_g          = "/public3/yongpeng/AnnotationBED/hg38/10-rRNA-Genes.bed"
region_11A_activeEn_g      = "/public3/yongpeng/AnnotationBED/hg38/11A-active-H3K4me1-H3K27ac.bed"
region_11B_otherEn_g       = "/public3/yongpeng/AnnotationBED/hg38/11B-H3K4me1-H3K27ac-H3K27me3.bed"
region_11C_posiedEn_g      = "/public3/yongpeng/AnnotationBED/hg38/11C-posied-H3K4me1-H3K27me3.bed"
region_11D_primedEn_g      = "/public3/yongpeng/AnnotationBED/hg38/11D-primed-H3K4me1.bed"

gene.obj_g     = readTranscriptFeatures(myRefSeqGenes_g)
myrepeat.obj_g = readFeatureFlank(myRepeats_g,           flank=2000, feature.flank.name=c("Repeats",          "shores"))
imprint1.obj_g = readFeatureFlank(myImprintedRegions1_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj_g = readFeatureFlank(myImprintedRegions2_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj_g = readFeatureFlank(myImprintedRegions3_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj_g = readFeatureFlank(myImprintedRegions4_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj_g = readFeatureFlank(myImprintedRegions5_g, flank=2000, feature.flank.name=c("ImprintedRegions", "shores"))
region_1_H3K4me1.obj_g        = readFeatureFlank(region_1_H3K4me1_g,          flank=2000, feature.flank.name=c("H3K4me1",           "shores"))
region_2_H3K4me3.obj_g        = readFeatureFlank(region_2_H3K4me3_g,          flank=2000, feature.flank.name=c("H3K4me3",           "shores"))
region_3_H3K27ac.obj_g        = readFeatureFlank(region_3_H3K27ac_g,          flank=2000, feature.flank.name=c("H3K27ac",           "shores"))
region_4_H3K27me3.obj_g       = readFeatureFlank(region_4_H3K27me3_g,         flank=2000, feature.flank.name=c("H3K27me3",          "shores"))
region_5_380imprint.obj_g     = readFeatureFlank(region_5_380imprint_g,       flank=2000, feature.flank.name=c("imprint",           "shores"))
region_6_centromeres.obj_g    = readFeatureFlank(region_6_centromeres_g,      flank=2000, feature.flank.name=c("centromere",        "shores"))
region_7_CpGislands.obj_g     = readFeatureFlank(region_7_CpGislands_g,       flank=2000, feature.flank.name=c("CpGislands",        "shores"))
region_8A_activeEn.obj_g      = readFeatureFlank(region_8A_activeEn_g,        flank=2000, feature.flank.name=c("activeEnhancers",   "shores"))
region_8B_otherEn.obj_g       = readFeatureFlank(region_8B_otherEn_g,         flank=2000, feature.flank.name=c("otherEnhancers",    "shores"))
region_8C_posiedEn.obj_g      = readFeatureFlank(region_8C_posiedEn_g,        flank=2000, feature.flank.name=c("posiedEnhancers",   "shores"))
region_8D_primedEn.obj_g      = readFeatureFlank(region_8D_primedEn_g,        flank=2000, feature.flank.name=c("primedEnhancers",   "shores"))
region_9_HouseKeeping.obj_g   = readFeatureFlank(region_9_HouseKeeping_g,     flank=2000, feature.flank.name=c("HouseKeeping",      "shores"))
region_10_rRNA.obj_g          = readFeatureFlank(region_10_rRNA_g,            flank=2000, feature.flank.name=c("rRNA_Genes",        "shores"))
region_11A_activeEn.obj_g      = readFeatureFlank(region_11A_activeEn_g,        flank=2000, feature.flank.name=c("activeEnhancers",   "shores"))
region_11B_otherEn.obj_g       = readFeatureFlank(region_11B_otherEn_g,         flank=2000, feature.flank.name=c("otherEnhancers",    "shores"))
region_11C_posiedEn.obj_g      = readFeatureFlank(region_11C_posiedEn_g,        flank=2000, feature.flank.name=c("posiedEnhancers",   "shores"))
region_11D_primedEn.obj_g      = readFeatureFlank(region_11D_primedEn_g,        flank=2000, feature.flank.name=c("primedEnhancers",   "shores"))
##################################################################################################################





##################################################################################################################
myOutDir_sub1_g = paste(outDir_g, "/0_ReadRawFiles",  sep="") 
if( ! file.exists(myOutDir_sub1_g) ) { dir.create(myOutDir_sub1_g, recursive = TRUE) }

sink( file=paste(myOutDir_sub1_g, "1_length-variables.txt", sep="/") )
length( Files_All_vector_g )
length( Files_All_list_g )
length( mySampleID_All_vector_g )
length( mySampleID_All_list_g )
length( myTreatment_All_vector_g )
length( myTreatment_All_list_g )
length( mySex_All_vector_g )
length( mySex_All_list_g )
length( myTech_All_vector_g )
length( myTech_All_list_g )
print( "#################### Files_All_vector_g ####################" )
print( Files_All_vector_g )
print( "#################### Files_All_list_g ####################" )
print( Files_All_list_g )
print( "#################### mySampleID_All_vector_g ####################" )
print( mySampleID_All_vector_g )
print( "#################### mySampleID_All_list_g ####################" )
print( mySampleID_All_list_g )
print( "#################### myTreatment_All_vector_g ####################" )
print( myTreatment_All_vector_g )
print( "#################### myTreatment_All_list_g ####################" )
print( myTreatment_All_list_g )
print( "#################### mySex_All_vector_g ####################" )
print( mySex_All_vector_g )
print( "#################### mySex_All_list_g ####################" )
print( mySex_All_list_g )
print( "#################### myTech_All_vector_g ####################" )
print( myTech_All_vector_g )
print( "#################### myTech_All_list_g ####################" )
print( myTech_All_list_g )
sink()


# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1_g, "2_theLog-of-read-inputFiles.txt", sep="/") )
myobj_g = methRead( Files_All_list_g,
               sample.id = mySampleID_All_list_g,
               assembly  = "hg38",
               treatment = myTreatment_All_vector_g,
               context   = "CpG",
               pipeline  = "bismarkCoverage",
               mincov    = lowestCoverage_g,       ## >= n
               header    = FALSE
)
sink()


sink( file=paste(myOutDir_sub1_g, "3_all-rawFiles.txt", sep="/") )
    print(Files_All_vector_g)
    print("#########################")
    print(myobj_g)
sink()


sink( file=paste(myOutDir_sub1_g, "4_dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  print( "######################" )
  print(   Files_All_vector_g[i]  )
  print(   dim(myobj_g[[i]])  )
}
sink()


sink( file=paste(myOutDir_sub1_g, "5_dimensions-of-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(myobj_g[[i]])[1], "\t", dim(myobj_g[[i]])[2], "\t", Files_All_vector_g[i], "\n"  )    
}
sink()



continue_on_error_g <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error_g) 




filtered.myobj_g = filterByCoverage( myobj_g,  lo.count = lowestCoverage_g,  lo.perc = NULL,  hi.count = NULL,  hi.perc = 99.99,  chunk.size = 1e+06,  save.db = FALSE )        

sink( file=paste(myOutDir_sub1_g, "6_dimensions-of-eachCov.filtered.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(filtered.myobj_g[[i]])[1], "\t", dim(filtered.myobj_g[[i]])[2], "\t", Files_All_vector_g[i], "\n"  )    
}
sink()


myobj_nor_g <- filtered.myobj_g
##myobj_nor_g <- normalizeCoverage(filtered.myobj_g,  method = "mean")




removed.myobj_g = filterByCoverage( myobj_g,  lo.count = lowestCoverage_g,  lo.perc = 99.99,  hi.count = NULL,  hi.perc = NULL,  chunk.size = 1e+06,  save.db = FALSE )
myOutDir_sub1_g2 = paste(outDir_g, "/0_ReadRawFiles/removed",  sep="") 
if( ! file.exists(myOutDir_sub1_g2) ) { dir.create(myOutDir_sub1_g2, recursive = TRUE) }

for( i in c(1:length(Files_All_vector_g)) ) {
  file_name = Files_All_vector_g[i]
  file_name = gsub(inputDir_g, i, file_name)
  file_name = gsub("/", "__", file_name)
  write.table(getData(removed.myobj_g[[i]])  , 
              file = paste(myOutDir_sub1_g2,   file_name,  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
}



sink( file=paste(myOutDir_sub1_g, "7_dimensions-removed.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(   dim(removed.myobj_g[[i]])[1], "\t", dim(removed.myobj_g[[i]])[2], "\t", Files_All_vector_g[i], "\n"  )    
}
sink()

sink( file=paste(myOutDir_sub1_g, "8_dimensions-all.txt", sep="/")  )
cat("Raw",  "Kept",  "Removed", "Files", "\n")
for( i in c(1:length(Files_All_vector_g)) ) {
  cat(    dim(myobj_g[[i]])[1], dim(filtered.myobj_g[[i]])[1], dim(removed.myobj_g[[i]])[1] ,  Files_All_vector_g[i], "\n"  )    
}
sink()





pdf( file=paste(myOutDir_sub1_g, "9A_MethylationStats-raw.pdf", sep="/")  )
for( i in c(1:length(myobj_g)) ) {
  getMethylationStats(myobj_g[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub1_g, "9B_MethylationStats-raw.txt", sep="/")  )
for( i in c(1:length(myobj_g)) ) {
  print("##############################################")
  print( paste(Files_All_vector_g[i],  ":", sep="") )
  print( getMethylationStats( myobj_g[[i]] )  )
}
sink()
pdf( file=paste(myOutDir_sub1_g, "10A_CoverageStats-raw.pdf", sep="/")  )
for( i in c(1:length(myobj_g)) ) {
  getCoverageStats(myobj_g[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub1_g, "10B_CoverageStats-raw.txt", sep="/")  )
for( i in c(1:length(myobj_g)) ) {
  print("##############################################")
  print( paste( Files_All_vector_g[i],   ":", sep="") )
  print( getCoverageStats( myobj_g[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub1_g, "11A_MethylationStats-kept.pdf", sep="/")  )
for( i in c(1:length(filtered.myobj_g)) ) {
  getMethylationStats(filtered.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub1_g, "11B_MethylationStats-kept.txt", sep="/")  )
for( i in c(1:length(filtered.myobj_g)) ) {
  print("##############################################")
  print( paste(Files_All_vector_g[i],  ":", sep="") )
  print( getMethylationStats( filtered.myobj_g[[i]] )  )
}
sink()
pdf( file=paste(myOutDir_sub1_g, "12A_CoverageStats-kept.pdf", sep="/")  )
for( i in c(1:length(filtered.myobj_g)) ) {
  getCoverageStats(filtered.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub1_g, "12B_CoverageStats-kept.txt", sep="/")  )
for( i in c(1:length(filtered.myobj_g)) ) {
  print("##############################################")
  print( paste( Files_All_vector_g[i],   ":", sep="") )
  print( getCoverageStats( filtered.myobj_g[[i]] )  )
}
sink()





pdf( file=paste(myOutDir_sub1_g, "13A_MethylationStats-removed.pdf", sep="/")  )
for( i in c(1:length(removed.myobj_g)) ) {
  getMethylationStats(removed.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub1_g, "13B_MethylationStats-removed.txt", sep="/")  )
for( i in c(1:length(removed.myobj_g)) ) {
  print("##############################################")
  print( paste(Files_All_vector_g[i],  ":", sep="") )
  print( getMethylationStats( removed.myobj_g[[i]] )  )
}
sink()
pdf( file=paste(myOutDir_sub1_g, "14A_CoverageStats-removed.pdf", sep="/")  )
for( i in c(1:length(removed.myobj_g)) ) {
  getCoverageStats(removed.myobj_g[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub1_g, "14B_CoverageStats-removed.txt", sep="/")  )
for( i in c(1:length(removed.myobj_g)) ) {
  print("##############################################")
  print( paste( Files_All_vector_g[i],   ":", sep="") )
  print( getCoverageStats( removed.myobj_g[[i]] )  )
}
sink()

##################################################################################################################






##################################################################################################################
  mySampleID_temp1_ART  <- c(mySampleID_IVF_fresh_g,   mySampleID_ICSI_fresh_g ,  mySampleID_IVF_frozen_g,   mySampleID_ICSI_frozen_g )
  myTreatment_temp1_ART <- c(myTreatment_IVF_fresh_g,  myTreatment_ICSI_fresh_g,  myTreatment_IVF_frozen_g,  myTreatment_ICSI_frozen_g )
  Sex_temp1_ART         <- c(Sex_IVF_fresh_g,   Sex_ICSI_fresh_g,   Sex_IVF_frozen_g,   Sex_ICSI_frozen_g )
  Tech_temp1_ART        <- c(Tech_IVF_fresh_g,  Tech_ICSI_fresh_g,  Tech_IVF_frozen_g,  Tech_ICSI_frozen_g )
  
  for(i in c(1:length(mySampleID_temp1_ART)) ) {
    myTreatment_temp1_ART[i] = 1
  }
  
  dataFrame_temp1_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_temp1_ART),
    mytreatment = c(myTreatment_NC_g, myTreatment_temp1_ART),
    mysex       = c(Sex_NC_g,         Sex_temp1_ART),
    mytech      = c(Tech_NC_g,        Tech_temp1_ART)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/1A_allSamples_100bpBin_1Cs",  sep="") ,   
                       binSize_temp1 = 100,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp1_A  )
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/1B_allSamples_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp1_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp2_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
    mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
    mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
  )
  
 
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/2_NC-vs-IVFfresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp2_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp3_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_frozen_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_IVF_frozen_g),
    mysex       = c(Sex_NC_g,  Sex_IVF_frozen_g),
    mytech      = c(Tech_NC_g, Tech_IVF_frozen_g)    
  )
  
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/3_NC-vs-IVFfrozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp3_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp4_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_ICSI_fresh_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_ICSI_fresh_g),
    mysex       = c(Sex_NC_g,  Sex_ICSI_fresh_g),
    mytech      = c(Tech_NC_g, Tech_ICSI_fresh_g)    
  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/4_NC-vs-ICSIfresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp4_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp5_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_NC_g,  Sex_ICSI_frozen_g),
    mytech      = c(Tech_NC_g, Tech_ICSI_frozen_g)    
  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/5_NC-vs-ICSIfrozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp5_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp6_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_fresh_g,  mySampleID_IVF_frozen_g),
    mytreatment = c(myTreatment_IVF_fresh_g, myTreatment_IVF_frozen_g),
    mysex       = c(Sex_IVF_fresh_g,         Sex_IVF_frozen_g),
    mytech      = c(Tech_IVF_fresh_g,        Tech_IVF_frozen_g)    
  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/6_IVF_fresh_vs_IVF_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp6_A  ) 
##################################################################################################################





##################################################################################################################
  dataFrame_temp7_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_fresh_g,  mySampleID_ICSI_fresh_g),
    mytreatment = c(myTreatment_IVF_fresh_g, myTreatment_ICSI_fresh_g),
    mysex       = c(Sex_IVF_fresh_g,         Sex_ICSI_fresh_g),
    mytech      = c(Tech_IVF_fresh_g,        Tech_ICSI_fresh_g)    
  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/7_IVF_fresh_vs_ICSI_fresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp7_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp8_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_fresh_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_IVF_fresh_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_IVF_fresh_g,         Sex_ICSI_frozen_g),
    mytech      = c(Tech_IVF_fresh_g,        Tech_ICSI_frozen_g)    
  )
  

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/8_IVF_fresh_vs_ICSI_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp8_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp9_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_frozen_g,  mySampleID_ICSI_fresh_g),
    mytreatment = c(myTreatment_IVF_frozen_g, myTreatment_ICSI_fresh_g),
    mysex       = c(Sex_IVF_frozen_g,         Sex_ICSI_fresh_g),
    mytech      = c(Tech_IVF_frozen_g,        Tech_ICSI_fresh_g)    
  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/9_IVF_frozen_vs_ICSI_fresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp9_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp10_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_frozen_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_IVF_frozen_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_IVF_frozen_g,         Sex_ICSI_frozen_g),
    mytech      = c(Tech_IVF_frozen_g,        Tech_ICSI_frozen_g)    
  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/10_IVF_frozen_vs_ICSI_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp10_A  )
##################################################################################################################





##################################################################################################################
  dataFrame_temp11_A <- data.frame(
    mysampleID  = c(mySampleID_ICSI_fresh_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_ICSI_fresh_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_ICSI_fresh_g,         Sex_ICSI_frozen_g),
    mytech      = c(Tech_ICSI_fresh_g,        Tech_ICSI_frozen_g)    
  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/11_ICSI_fresh_vs_ICSI_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp11_A  )
##################################################################################################################





##################################################################################################################
  mySampleID_temp12_IVF  <- c( mySampleID_IVF_fresh_g,  mySampleID_IVF_frozen_g )
  myTreatment_temp12_IVF <- c( myTreatment_IVF_fresh_g,  myTreatment_IVF_frozen_g )
  Sex_temp12_IVF         <- c( Sex_IVF_fresh_g,  Sex_IVF_frozen_g )
  Tech_temp12_IVF        <- c( Tech_IVF_fresh_g,  Tech_IVF_frozen_g )

  for(i in c(1:length(mySampleID_temp12_IVF)) ) {
    myTreatment_temp12_IVF[i] = 1
  }

  dataFrame_temp12_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_temp12_IVF),
    mytreatment = c(myTreatment_NC_g, myTreatment_temp12_IVF),
    mysex       = c(Sex_NC_g,         Sex_temp12_IVF),
    mytech      = c(Tech_NC_g,        Tech_temp12_IVF)    
  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/12_NC_vs_IVF_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp12_A  )
##################################################################################################################





##################################################################################################################
  mySampleID_temp13_ICSI  <- c( mySampleID_ICSI_fresh_g,  mySampleID_ICSI_frozen_g )
  myTreatment_temp13_ICSI <- c( myTreatment_ICSI_fresh_g,  myTreatment_ICSI_frozen_g )
  Sex_temp13_ICSI         <- c( Sex_ICSI_fresh_g,  Sex_ICSI_frozen_g )
  Tech_temp13_ICSI        <- c( Tech_ICSI_fresh_g,  Tech_ICSI_frozen_g )

  for(i in c(1:length(mySampleID_temp13_ICSI)) ) {
    myTreatment_temp13_ICSI[i] = 1
  }

  dataFrame_temp13_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_temp13_ICSI),
    mytreatment = c(myTreatment_NC_g, myTreatment_temp13_ICSI),
    mysex       = c(Sex_NC_g,         Sex_temp13_ICSI),
    mytech      = c(Tech_NC_g,        Tech_temp13_ICSI)    
  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/13_NC_vs_ICSI_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp13_A  )
##################################################################################################################





##################################################################################################################
  mySampleID_temp14_fresh  <- c( mySampleID_IVF_fresh_g,  mySampleID_ICSI_fresh_g )
  myTreatment_temp14_fresh <- c( myTreatment_IVF_fresh_g,  myTreatment_ICSI_fresh_g )
  Sex_temp14_fresh         <- c( Sex_IVF_fresh_g,  Sex_ICSI_fresh_g )
  Tech_temp14_fresh        <- c( Tech_IVF_fresh_g,  Tech_ICSI_fresh_g )

  for(i in c(1:length(mySampleID_temp14_fresh)) ) {
    myTreatment_temp14_fresh[i] = 1
  }

  dataFrame_temp14_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_temp14_fresh),
    mytreatment = c(myTreatment_NC_g, myTreatment_temp14_fresh),
    mysex       = c(Sex_NC_g,         Sex_temp14_fresh),
    mytech      = c(Tech_NC_g,        Tech_temp14_fresh)    
  )
  

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/14_NC_vs_fresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp14_A  ) 
##################################################################################################################





##################################################################################################################
  mySampleID_temp15_frozen  <- c( mySampleID_IVF_frozen_g,  mySampleID_ICSI_frozen_g )
  myTreatment_temp15_frozen <- c( myTreatment_IVF_frozen_g,  myTreatment_ICSI_frozen_g )
  Sex_temp15_frozen         <- c( Sex_IVF_frozen_g,  Sex_ICSI_frozen_g )
  Tech_temp15_frozen        <- c( Tech_IVF_frozen_g,  Tech_ICSI_frozen_g )

  for(i in c(1:length(mySampleID_temp15_frozen)) ) {
    myTreatment_temp15_frozen[i] = 1
  }

  dataFrame_temp15_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_temp15_frozen),
    mytreatment = c(myTreatment_NC_g, myTreatment_temp15_frozen),
    mysex       = c(Sex_NC_g,         Sex_temp15_frozen),
    mytech      = c(Tech_NC_g,        Tech_temp15_frozen)    
  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/15_NC_vs_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp15_A  )
##################################################################################################################





##################################################################################################################
  mySampleID_temp16_ART  <- c(mySampleID_IVF_fresh_g,  mySampleID_ICSI_fresh_g , mySampleID_IVF_frozen_g,  mySampleID_ICSI_frozen_g )
  myTreatment_temp16_ART <- c(myTreatment_IVF_fresh_g,  myTreatment_ICSI_fresh_g,  myTreatment_IVF_frozen_g,  myTreatment_ICSI_frozen_g )
  Sex_temp16_ART         <- c(Sex_IVF_fresh_g,  Sex_ICSI_fresh_g,  Sex_IVF_frozen_g,  Sex_ICSI_frozen_g )
  Tech_temp16_ART        <- c(Tech_IVF_fresh_g,  Tech_ICSI_fresh_g,  Tech_IVF_frozen_g,  Tech_ICSI_frozen_g )

  for(i in c(1:length(mySampleID_temp16_ART)) ) {
    myTreatment_temp16_ART[i] = 1
  }

  dataFrame_temp16_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_temp16_ART),
    mytreatment = c(myTreatment_NC_g, myTreatment_temp16_ART),
    mysex       = c(Sex_NC_g,         Sex_temp16_ART),
    mytech      = c(Tech_NC_g,        Tech_temp16_ART)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/16_NC_vs_ART_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp16_A  )
##################################################################################################################





##################################################################################################################
  mySampleID_temp17_IVF  <- c(mySampleID_IVF_fresh_g,   mySampleID_IVF_frozen_g )
  myTreatment_temp17_IVF <- c(myTreatment_IVF_fresh_g,   myTreatment_IVF_frozen_g )
  Sex_temp17_IVF         <- c(Sex_IVF_fresh_g,  Sex_IVF_frozen_g )
  Tech_temp17_IVF        <- c(Tech_IVF_fresh_g,   Tech_IVF_frozen_g )

  mySampleID_temp17_ICSI  <- c( mySampleID_ICSI_fresh_g ,  mySampleID_ICSI_frozen_g )
  myTreatment_temp17_ICSI <- c( myTreatment_ICSI_fresh_g,   myTreatment_ICSI_frozen_g )
  Sex_temp17_ICSI         <- c( Sex_ICSI_fresh_g,   Sex_ICSI_frozen_g )
  Tech_temp17_ICSI        <- c( Tech_ICSI_fresh_g,   Tech_ICSI_frozen_g )

  for(i in c(1:length(mySampleID_temp17_IVF)) ) {
    myTreatment_temp17_IVF[i] = 0
  }

  for(i in c(1:length(mySampleID_temp17_ICSI)) ) {
    myTreatment_temp17_ICSI[i] = 1
  }

  dataFrame_temp17_A <- data.frame(
    mysampleID  = c(mySampleID_temp17_IVF,  mySampleID_temp17_ICSI),
    mytreatment = c(myTreatment_temp17_IVF, myTreatment_temp17_ICSI),
    mysex       = c(Sex_temp17_IVF,         Sex_temp17_ICSI),
    mytech      = c(Tech_temp17_IVF,        Tech_temp17_ICSI)    
  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/17_IVF_vs_ICSI_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp17_A  )
##################################################################################################################





##################################################################################################################
  mySampleID_temp18_fresh  <- c(mySampleID_IVF_fresh_g,  mySampleID_ICSI_fresh_g   )
  myTreatment_temp18_fresh <- c(myTreatment_IVF_fresh_g,  myTreatment_ICSI_fresh_g  )
  Sex_temp18_fresh         <- c(Sex_IVF_fresh_g,  Sex_ICSI_fresh_g )
  Tech_temp18_fresh        <- c(Tech_IVF_fresh_g,  Tech_ICSI_fresh_g )

  mySampleID_temp18_frozen  <- c(  mySampleID_IVF_frozen_g,  mySampleID_ICSI_frozen_g )
  myTreatment_temp18_frozen <- c(  myTreatment_IVF_frozen_g,  myTreatment_ICSI_frozen_g )
  Sex_temp18_frozen         <- c(  Sex_IVF_frozen_g,  Sex_ICSI_frozen_g )
  Tech_temp18_frozen        <- c(  Tech_IVF_frozen_g,  Tech_ICSI_frozen_g )

  for(i in c(1:length(mySampleID_temp18_fresh)) ) {
    myTreatment_temp18_fresh[i] = 0
  }

  for(i in c(1:length(mySampleID_temp18_frozen)) ) {
    myTreatment_temp18_frozen[i] = 1
  }

  dataFrame_temp18_A <- data.frame(
    mysampleID  = c(mySampleID_temp18_fresh,  mySampleID_temp18_frozen),
    mytreatment = c(myTreatment_temp18_fresh, myTreatment_temp18_frozen),
    mysex       = c(Sex_temp18_fresh,         Sex_temp18_frozen),
    mytech      = c(Tech_temp18_fresh,        Tech_temp18_frozen)    
  )
  

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/18_fresh_vs_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp18_A  )
##################################################################################################################





##################################################################################################################
  mySampleID_temp100_ART  <- c(mySampleID_IVF_fresh_g,   mySampleID_ICSI_fresh_g ,  mySampleID_IVF_frozen_g,   mySampleID_ICSI_frozen_g )
  myTreatment_temp100_ART <- c(myTreatment_IVF_fresh_g,  myTreatment_ICSI_fresh_g,  myTreatment_IVF_frozen_g,  myTreatment_ICSI_frozen_g )
  Sex_temp100_ART         <- c(Sex_IVF_fresh_g,   Sex_ICSI_fresh_g,   Sex_IVF_frozen_g,   Sex_ICSI_frozen_g )
  Tech_temp100_ART        <- c(Tech_IVF_fresh_g,  Tech_ICSI_fresh_g,  Tech_IVF_frozen_g,  Tech_ICSI_frozen_g )
  
  for(i in c(1:length(mySampleID_temp100_ART)) ) {
    myTreatment_temp100_ART[i] = 1
  }
  
  dataFrame_temp100_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_temp100_ART),
    mytreatment = c(myTreatment_NC_g, myTreatment_temp100_ART),
    mysex       = c(Sex_NC_g,         Sex_temp100_ART),
    mytech      = c(Tech_NC_g,        Tech_temp100_ART)    
  )
  
  
  ###
  myMainFunction1bp_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/AllSamples_1bpBin",  sep="") ,
                          dataFrame_temp1 = dataFrame_temp100_A    )
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/AllSamples_300bpBin_1Cs/1-matrix-all",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp100_A  )
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/AllSamples_300bpBin_2Cs/1-matrix-all",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 2, dataFrame_temp1 = dataFrame_temp100_A  )
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/AllSamples_300bpBin_3Cs/1-matrix-all",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp100_A  )
  

  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/AllSamples_100bpBin_3Cs/1-matrix-all",  sep="") ,   
                       binSize_temp1 = 100,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp100_A  )
##################################################################################################################




  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/1A_allSamples_100bpBin_2Cs",  sep="") ,   
                       binSize_temp1 = 100,   binBases_temp1 = 2, dataFrame_temp1 = dataFrame_temp1_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/1B_allSamples_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp1_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/2_NC-vs-IVFfresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp2_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/3_NC-vs-IVFfrozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp3_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/4_NC-vs-ICSIfresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp4_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/5_NC-vs-ICSIfrozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp5_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/6_IVF_fresh_vs_IVF_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp6_A  ) 
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/7_IVF_fresh_vs_ICSI_fresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp7_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/8_IVF_fresh_vs_ICSI_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp8_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/9_IVF_frozen_vs_ICSI_fresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp9_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/10_IVF_frozen_vs_ICSI_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp10_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/11_ICSI_fresh_vs_ICSI_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp11_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/12_NC_vs_IVF_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp12_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/13_NC_vs_ICSI_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp13_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/14_NC_vs_fresh_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp14_A  ) 
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/15_NC_vs_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp15_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/16_NC_vs_ART_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp16_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/17_IVF_vs_ICSI_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp17_A  )
  myMainFunction_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/18_fresh_vs_frozen_300bpBin_3Cs",  sep="") ,   
                       binSize_temp1 = 300,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp18_A  )
  myMainFunction1bp_2_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g,  "/FullResults",   "/AllSamples_1bpBin",  sep="") ,
                          dataFrame_temp1 = dataFrame_temp100_A    )




