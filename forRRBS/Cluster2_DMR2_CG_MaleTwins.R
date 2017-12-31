##################################################################################################################
## All samples must be overlapped. (100% overlap)
## Suffixes of all self-defined global variables must be "_g".
##
## Example:  
## Rscript  Cluster_DMR_CG_singletons_Pups.R     11B_splitXY/100A-rmXY/5_cov50reads     1011B_splitXY/100A-rmXY/5_cov50reads/Cluster_DMR_CG_singletons_Pups          
         
args_g <- commandArgs(TRUE)
print("args: ")
print(args_g[1])   
print(args_g[2])     
print("#############")

inputDir_g = args_g[1];     ## the path of input files
outDir_g   = args_g[2];     ## the path of output files
# inputDir_g =  "11B_splitXY/100A-rmXY/3_cov30reads"
# outDir_g   =  "1011B_splitXY/100A-rmXY/3_cov30reads/Cluster_DMR_CG_singletons_Pups"
print(inputDir_g)   
print(outDir_g)

if( ! file.exists(inputDir_g) ) { print("##### Error-1 #####") }
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }
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




Files_NC_g <- c(
  paste( inputDir_g,   "67_E24C-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "67_E24D-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56C-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56D-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "69_NC-E123-C-Boy_Rep1.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "69_NC-E123-D-Boy_Rep1.bismark.cov",  sep="/" ) 
)

mySampleID_NC_g <- rep( x="NC", times=length(Files_NC_g) )
Tech_NC_g <- mySampleID_NC_g
  
for(i in c(1:length(mySampleID_NC_g)) ) {
  mySampleID_NC_g[i] = paste(mySampleID_NC_g[i], i, sep="_")
}

myTreatment_NC_g <- rep( x=0,  times=length(Files_NC_g) )

Sex_NC_g = rep( x="boy",  times=length(Files_NC_g) )  
for(i in c(1:length(Sex_NC_g)) ) {
  Sex_NC_g[i] = "boy"
  if(i>=4) { Sex_NC_g[i] = "boy" }
}

  

Files_IVF_fresh_g <- c(
  paste( inputDir_g,   "70_E113C-boy-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "70_E113D-boy-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-C-Boy_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-D-Boy_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "72_ART-W779-C-Boy_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "72_W779D-ART-boy_Rep1.bismark.cov",    sep="/" ) 
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
  if(i>=4) { Sex_IVF_fresh_g[i] = "boy" }
}



Files_IVF_frozen_g <- c(
  paste( inputDir_g,   "9_W1365C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "9_Q5-W1365D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "11_W1398C-boy-IVF-frozen_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "11_W1398D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ) 
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
  if(i>=4) { Sex_IVF_frozen_g[i] = "boy" }
}



Files_ICSI_fresh_g <- c(
  paste( inputDir_g,   "12_W1579C-boy-ICSI-fresh-merge_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "12_Q17-W1579D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_W1647C-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_W1647D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "14_W1719C-boy-ICSI-fresh_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "14_W1719D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ) 
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
  if(i>=4) { Sex_ICSI_fresh_g[i] = "boy" }
}



Files_ICSI_frozen_g <- c(
  paste( inputDir_g,   "15_Q6-W871C-boy-ICSI-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "15_Q4-W871D-boy-ICSI-frozen_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_ICSI_frozen_g <- rep( x="ICSI_frozen", times=length(Files_ICSI_frozen_g) )
Tech_ICSI_frozen_g <- mySampleID_ICSI_frozen_g

for(i in c(1:length(mySampleID_ICSI_frozen_g)) ) {
  mySampleID_ICSI_frozen_g[i] = paste(mySampleID_ICSI_frozen_g[i], i, sep="_")
}

myTreatment_ICSI_frozen_g <- rep( x=4,  times=length(Files_ICSI_frozen_g) )

Sex_ICSI_frozen_g = rep( x="boy",      times=length(Files_ICSI_frozen_g) )  
for(i in c(1:length(Sex_ICSI_frozen_g)) ) {
  Sex_ICSI_frozen_g[i] = "boy"
  if(i>=4) { Sex_ICSI_frozen_g[i] = "boy" }
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
myTech_All_color_g = c( "NC"="black", "IVF_fresh"="blue", "IVF_frozen"="green",  "ICSI_fresh"="red", "ICSI_frozen"="purple" )
mySex_All_shape_g
myTech_All_color_g

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
    if(myTech_color[i] == "NC")          { myTech_color[i] = c("NC"="black") }
    if(myTech_color[i] == "IVF_fresh")   { myTech_color[i] = c("IVF_fresh"="blue") }
    if(myTech_color[i] == "IVF_frozen")  { myTech_color[i] = c("IVF_frozen"="green") }
    if(myTech_color[i] == "ICSI_fresh")  { myTech_color[i] = c("ICSI_fresh"="red") }
    if(myTech_color[i] == "ICSI_frozen") { myTech_color[i] = c("ICSI_frozen"="purple") }
  }
  names(myTech_color) = myTech_vector
  return(myTech_color)
}

MySex_Shape_g( mySex_All_vector_g )
MyTech_color_g( myTech_All_vector_g )
## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple


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


Files_NC_parents_g <- c(
  paste( inputDir_g,   "67_NC-E24-F-Father_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56F-father-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "69_NC-E123-F-Father_Rep1.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "67_NC-E24-M-Mother_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56M-mother-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "69_NC-E123-M-Mother-merge_Rep1.bismark.cov",  sep="/" ) 
)

mySampleID_NC_parents_g <- rep( x="NC", times=length(Files_NC_parents_g) )
Tech_NC_parents_g <- mySampleID_NC_parents_g

for(i in c(1:length(mySampleID_NC_parents_g)) ) {
  mySampleID_NC_parents_g[i] = paste(mySampleID_NC_parents_g[i], i, sep="_")
}

myTreatment_NC_parents_g <- rep( x=0,  times=length(Files_NC_parents_g) )

Sex_NC_parents_g = rep( x="father",  times=length(Files_NC_parents_g) )  
for(i in c(1:length(Sex_NC_parents_g)) ) {
  Sex_NC_parents_g[i] = "father"
  if(i>=4) { Sex_NC_parents_g[i] = "mother" }
}



Files_IVF_fresh_parents_g <- c(
  paste( inputDir_g,   "70_E113F-father-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-F-Father_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "72_ART-W779-F-Father_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "70_E113M-mother-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-M-Mother_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "72_ART-W779-M-Mother_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_IVF_fresh_parents_g <- rep( x="IVF_fresh", times=length(Files_IVF_fresh_parents_g) )
Tech_IVF_fresh_parents_g <- mySampleID_IVF_fresh_parents_g

for(i in c(1:length(mySampleID_IVF_fresh_parents_g)) ) {
  mySampleID_IVF_fresh_parents_g[i] = paste(mySampleID_IVF_fresh_parents_g[i], i, sep="_")
}

myTreatment_IVF_fresh_parents_g <- rep( x=1,  times=length(Files_IVF_fresh_parents_g) )

Sex_IVF_fresh_parents_g = rep( x="father",      times=length(Files_IVF_fresh_parents_g) ) 
for(i in c(1:length(Sex_IVF_fresh_parents_g)) ) {
  Sex_IVF_fresh_parents_g[i] = "father"
  if(i>=4) { Sex_IVF_fresh_parents_g[i] = "mother" }
}



Files_IVF_frozen_parents_g <- c(
  paste( inputDir_g,   "9_W1365F-father-IVF-frozen-merge_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733F-father-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "11_W1398F-father-IVF-frozen_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "9_Q1-W1365M-mother-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733M-Mother-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "12_Q19-W1579M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_IVF_frozen_parents_g <- rep( x="IVF_frozen", times=length(Files_IVF_frozen_parents_g) )
Tech_IVF_frozen_parents_g <- mySampleID_IVF_frozen_parents_g

for(i in c(1:length(mySampleID_IVF_frozen_parents_g)) ) {
  mySampleID_IVF_frozen_parents_g[i] = paste(mySampleID_IVF_frozen_parents_g[i], i, sep="_")
}

myTreatment_IVF_frozen_parents_g <- rep( x=2,  times=length(Files_IVF_frozen_parents_g) )

Sex_IVF_frozen_parents_g = rep( x="father",      times=length(Files_IVF_frozen_parents_g) )
for(i in c(1:length(Sex_IVF_frozen_parents_g)) ) {
  Sex_IVF_frozen_parents_g[i] = "father"
  if(i>=4) { Sex_IVF_frozen_parents_g[i] = "mother" }
}



Files_ICSI_fresh_parents_g <- c(
  paste( inputDir_g,   "12_Q18-W1579F-father-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_W1647F-father-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "14_W1719F-father-ICSI-fresh_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "12_Q19-W1579M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_Q23-W1647M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "14_Q15-W1719M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_ICSI_fresh_parents_g <- rep( x="ICSI_fresh", times=length(Files_ICSI_fresh_parents_g) )
Tech_ICSI_fresh_parents_g <- mySampleID_ICSI_fresh_parents_g

for(i in c(1:length(mySampleID_ICSI_fresh_parents_g)) ) {
  mySampleID_ICSI_fresh_parents_g[i] = paste(mySampleID_ICSI_fresh_parents_g[i], i, sep="_")
}

myTreatment_ICSI_fresh_parents_g <- rep( x=3,  times=length(Files_ICSI_fresh_parents_g) )

Sex_ICSI_fresh_parents_g = rep( x="father",      times=length(Files_ICSI_fresh_parents_g) )  
for(i in c(1:length(Sex_ICSI_fresh_parents_g)) ) {
  Sex_ICSI_fresh_parents_g[i] = "father"
  if(i>=4) { Sex_ICSI_fresh_parents_g[i] = "mother" }
}



Files_ICSI_frozen_parents_g <- c(
  paste( inputDir_g,   "15_Q7-W871F-father-ICSI-frozen_Rep1.bismark.cov",     sep="/" ),
  paste( inputDir_g,   "15_Q10-W871M-mother-ICSI-frozen_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_ICSI_frozen_parents_g <- rep( x="ICSI_frozen", times=length(Files_ICSI_frozen_parents_g) )
Tech_ICSI_frozen_parents_g <- mySampleID_ICSI_frozen_parents_g

for(i in c(1:length(mySampleID_ICSI_frozen_parents_g)) ) {
  mySampleID_ICSI_frozen_parents_g[i] = paste(mySampleID_ICSI_frozen_parents_g[i], i, sep="_")
}

myTreatment_ICSI_frozen_parents_g <- rep( x=4,  times=length(Files_ICSI_frozen_parents_g) )

Sex_ICSI_frozen_parents_g = rep( x="father",      times=length(Files_ICSI_frozen_parents_g) )  
for(i in c(1:length(Sex_ICSI_frozen_parents_g)) ) {
  Sex_ICSI_frozen_parents_g[i] = "father"
  if(i>=2) { Sex_ICSI_frozen_parents_g[i] = "mother" }
}









Files_All_vector_parents_g <- c(
  Files_NC_parents_g,
  Files_IVF_fresh_parents_g,
  Files_IVF_frozen_parents_g,
  Files_ICSI_fresh_parents_g,
  Files_ICSI_frozen_parents_g  
)
Files_All_list_parents_g <- as.list( Files_All_vector_parents_g )


mySampleID_All_vector_parents_g <- c( 
  mySampleID_NC_parents_g,
  mySampleID_IVF_fresh_parents_g,
  mySampleID_IVF_frozen_parents_g,
  mySampleID_ICSI_fresh_parents_g,
  mySampleID_ICSI_frozen_parents_g  
)
mySampleID_All_list_parents_g <- as.list( mySampleID_All_vector_parents_g )


myTreatment_All_vector_parents_g <- c( 
  myTreatment_NC_parents_g,
  myTreatment_IVF_fresh_parents_g,
  myTreatment_IVF_frozen_parents_g,
  myTreatment_ICSI_fresh_parents_g,
  myTreatment_ICSI_frozen_parents_g 
)       
myTreatment_All_list_parents_g <- as.list( myTreatment_All_vector_parents_g )


mySex_All_vector_parents_g <- c( 
  Sex_NC_parents_g,
  Sex_IVF_fresh_parents_g,
  Sex_IVF_frozen_parents_g,
  Sex_ICSI_fresh_parents_g,
  Sex_ICSI_frozen_parents_g 
)       
mySex_All_list_parents_g <- as.list( mySex_All_vector_parents_g )


myTech_All_vector_parents_g <- c( 
  Tech_NC_parents_g,
  Tech_IVF_fresh_parents_g,
  Tech_IVF_frozen_parents_g,
  Tech_ICSI_fresh_parents_g,
  Tech_ICSI_frozen_parents_g 
)       
myTech_All_list_parents_g <- as.list( myTech_All_vector_parents_g )


## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple
mySex_All_shape_parents_g  = c( "boy"=16, "girl"=17, "father"=1, "mother"=2 ) 
myTech_All_color_parents_g = c( "NC"="black", "IVF_fresh"="blue", "IVF_frozen"="green",  "ICSI_fresh"="red", "ICSI_frozen"="purple" )
mySex_All_shape_parents_g
myTech_All_color_parents_g

MySex_Shape_parents_g <- function(  mySex_vector   ) {
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

MyTech_color_parents_g <- function(  myTech_vector  ) {
  myTech_color  = myTech_vector
  for(i in c(1:length(myTech_color)) ) {
    if(myTech_color[i] == "NC")          { myTech_color[i] = c("NC"="black") }
    if(myTech_color[i] == "IVF_fresh")   { myTech_color[i] = c("IVF_fresh"="blue") }
    if(myTech_color[i] == "IVF_frozen")  { myTech_color[i] = c("IVF_frozen"="green") }
    if(myTech_color[i] == "ICSI_fresh")  { myTech_color[i] = c("ICSI_fresh"="red") }
    if(myTech_color[i] == "ICSI_frozen") { myTech_color[i] = c("ICSI_frozen"="purple") }
  }
  names(myTech_color) = myTech_vector
  return(myTech_color)
}

MySex_Shape_parents_g( mySex_All_vector_parents_g )
MyTech_color_parents_g( myTech_All_vector_parents_g )
## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple


length( Files_All_vector_parents_g )
length( Files_All_list_parents_g )
length( mySampleID_All_vector_parents_g )
length( mySampleID_All_list_parents_g )
length( myTreatment_All_vector_parents_g )
length( myTreatment_All_list_parents_g )
length( mySex_All_vector_parents_g )
length( mySex_All_list_parents_g )
length( myTech_All_vector_parents_g )
length( myTech_All_list_parents_g )
##################################################################################################################





##################################################################################################################
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
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	    ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	    ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	            ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
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
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	  ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	  ## facet labels along horizontal direction (element_text; inherits from strip.text)
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


MyCluster_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  dev.off()
}


MyCluster_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  dev.off()
}


MyCluster_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
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
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2L.sd0.35.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.35 )                      
}


MyPCA_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  dev.off()
}


MyPCA_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  dev.off()
}


MyPCA_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
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
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2L.sd0.35.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.35 )                      
}



myHierarchicalClustering_1_g  <- function(  mat_3three,   path_temp1,   dataFrame_temp1  )  {
  res.dist1_3three <- get_dist( t(mat_3three) ,   method = "euclidean")
  res.dist2_3three <- get_dist( t(mat_3three) ,   method = "maximum"  )
  res.dist3_3three <- get_dist( t(mat_3three) ,   method = "manhattan")
  res.dist4_3three <- get_dist( t(mat_3three) ,   method = "canberra" )
  res.dist5_3three <- get_dist( t(mat_3three) ,   method = "binary"   )
  res.dist6_3three <- get_dist( t(mat_3three) ,   method = "minkowski")
  res.dist7_3three <- get_dist( t(mat_3three) ,   method = "pearson"  )
  res.dist8_3three <- get_dist( t(mat_3three) ,   method = "spearman" )
  
  pdf( file = paste(path_temp1, "1A_visualizing-distance-matrix.pdf",  sep="/") )
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
  
  pdf( file = paste(path_temp1, "2A_hierarchical-ward.D-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="ward.D, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="ward.D, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="ward.D, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="ward.D, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="ward.D, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="ward.D, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="ward.D, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="ward.D, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "2B_hierarchical-ward.D-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="ward.D, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="ward.D, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="ward.D, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="ward.D, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="ward.D, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="ward.D, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="ward.D, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="ward.D, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  # Compute hierarchical clustering by "ward.D2"
  res.hc1_3three <- hclust(res.dist1_3three, method = "ward.D2"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "ward.D2"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "ward.D2"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "ward.D2"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "ward.D2"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "ward.D2"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "ward.D2"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "ward.D2"  ) 
  
  pdf( file = paste(path_temp1, "3A_hierarchical-ward.D2-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="ward.D2, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="ward.D2, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="ward.D2, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="ward.D2, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="ward.D2, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="ward.D2, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="ward.D2, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="ward.D2, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "3B_hierarchical-ward.D2-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="ward.D2, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="ward.D2, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="ward.D2, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="ward.D2, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="ward.D2, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="ward.D2, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="ward.D2, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="ward.D2, spearman"   , repel = TRUE)  )          
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
  
  pdf( file = paste(path_temp1, "4A_hierarchical-single-rectangle.pdf",  sep="/") )
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

  # Compute hierarchical clustering by "complete"
  res.hc1_3three <- hclust(res.dist1_3three, method = "complete"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "complete"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "complete"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "complete"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "complete"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "complete"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "complete"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "complete"  ) 
  
  pdf( file = paste(path_temp1, "5A_hierarchical-complete-rectangle.pdf",  sep="/") )
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

  # Compute hierarchical clustering by "average"
  res.hc1_3three <- hclust(res.dist1_3three, method = "average"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "average"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "average"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "average"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "average"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "average"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "average"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "average"  ) 
  
  pdf( file = paste(path_temp1, "6A_hierarchical-average-rectangle.pdf",  sep="/") )
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

  # Compute hierarchical clustering by "mcquitty"
  res.hc1_3three <- hclust(res.dist1_3three, method = "mcquitty"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "mcquitty"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "mcquitty"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "mcquitty"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "mcquitty"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "mcquitty"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "mcquitty"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "mcquitty"  ) 
  
  pdf( file = paste(path_temp1, "7A_hierarchical-mcquitty-rectangle.pdf",  sep="/") )
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

  # Compute hierarchical clustering by "median"
  res.hc1_3three <- hclust(res.dist1_3three, method = "median"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "median"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "median"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "median"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "median"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "median"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "median"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "median"  ) 
  
  pdf( file = paste(path_temp1, "8A_hierarchical-median-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="median, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="median, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="median, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="median, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="median, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="median, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="median, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="median, spearman"   )  )
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
  
  pdf( file = paste(path_temp1, "9A_hierarchical-centroid-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="centroid, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="centroid, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="centroid, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="centroid, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="centroid, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="centroid, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="centroid, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="centroid, spearman"   )  )
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
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, "filterByQuantile_1", sep="/")
  path1_2 = paste(path1, "filterByQuantile_2", sep="/")
  path1_3 = paste(path1, "filterByQuantile_3", sep="/")
  path1_4 = paste(path1, "filterByQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
}


MyPCA_2A_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1,  dataFrame_temp1   ) {  
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, "filterNoQuantile_1", sep="/")
  path1_2 = paste(path1, "filterNoQuantile_2", sep="/")
  path1_3 = paste(path1, "filterNoQuantile_3", sep="/")
  path1_4 = paste(path1, "filterNoQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
}



MyPCA_3A_g <- function(  mymeth2 ,  path2,   file2, width2, height2,  dataFrame_temp2   )  {
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1P.kept0.5percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.995 ,dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1Q.kept0.1percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.999 ,dataFrame_temp1=dataFrame_temp2   )                     
  
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.18.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.22.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.22 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.28.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.28 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.30.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2L.sd0.35.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.35 ,  dataFrame_temp1=dataFrame_temp2   )                       
}


##########################
## Annotating differentially methylated bases or regions
myRefSeqGenes_g = "/home/yongpeng/AnnotationBED/hg38_RefSeq_Genes.bed"
myCpGIslands_g  = "/home/yongpeng/AnnotationBED/hg38_CpG_islands.bed"
myRepeats_g     = "/home/yongpeng/AnnotationBED/hg38_Repeats_rmsk.bed"
myImprintedRegions1_g = "/home/yongpeng/AnnotationBED/67.Regions.PlosGenetics.ImprintedGenes.hg38.bed"
myImprintedRegions2_g = "/home/yongpeng/AnnotationBED/75Regions.GR.hg38.bed"
myImprintedRegions3_g = "/home/yongpeng/AnnotationBED/369Regions.GR.hg38.bed"
myImprintedRegions4_g = "/home/yongpeng/AnnotationBED/merge1.imprintedRegions.hg38.bed"
myImprintedRegions5_g = "/home/yongpeng/AnnotationBED/merge2.imprintedRegions.hg38.bed"

gene.obj_g     = readTranscriptFeatures(myRefSeqGenes_g)
cpg.obj_g      = readFeatureFlank(myCpGIslands_g, flank=10000, feature.flank.name=c("CpGi", "shores"))
myrepeat.obj_g = readFeatureFlank(myRepeats_g,    flank=10000, feature.flank.name=c("Repeats", "shores"))
imprint1.obj_g = readFeatureFlank(myImprintedRegions1_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj_g = readFeatureFlank(myImprintedRegions2_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj_g = readFeatureFlank(myImprintedRegions3_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj_g = readFeatureFlank(myImprintedRegions4_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj_g = readFeatureFlank(myImprintedRegions5_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
##########################


#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myDiff_DMC_DMR_g <- function(  methobj2,   path2  ) {
  
  myDiff_2two_sub2 = calculateDiffMeth(methobj2,  num.cores=16)
  dim(myDiff_2two_sub2)
  names(myDiff_2two_sub2)
  head(myDiff_2two_sub2)
  
  myQvalue_2two_sub2 = myDiff_2two_sub2$qvalue
  myMethDi_2two_sub2 = myDiff_2two_sub2$meth.diff
  
  pdf(paste(path2,  "1A_qvalue_distribution.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1), freq=FALSE) )
  print( hist( myQvalue_2two_sub2[myQvalue_2two_sub2<0.1],  nclass=20, xlim=c(0, 0.1), freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1), freq=TRUE) )
  print( hist( myQvalue_2two_sub2[myQvalue_2two_sub2<0.1],  nclass=20, xlim=c(0, 0.1), freq=TRUE) )
  dev.off() 
  
  pdf(paste(path2,  "1B_methDiff_distribution.pdf",  sep="/"))
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100), freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 50), freq=FALSE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=20], nclass=81, xlim=c(20, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=10], nclass=91, xlim=c(10, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2<=10], nclass=100, xlim=c(0, 10),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100), freq=TRUE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 50), freq=TRUE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=20], nclass=81, xlim=c(20, 100),  freq=TRUE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=10], nclass=91, xlim=c(10, 100),  freq=TRUE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2<=10], nclass=100, xlim=c(0, 10),   freq=TRUE) )
  dev.off() 
  
  myColor1_2two_sub2 <- rep( "no",   times= length(myQvalue_2two_sub2) )
  myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>10) & (myQvalue_2two_sub2<0.001) ]  <- "yes"
  length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  
  DataFrame2_2two_sub2 <- data.frame(myx1 = myMethDi_2two_sub2,   
                                     myy1 =  -log10(myQvalue_2two_sub2),  
                                     mycolor1 = myColor1_2two_sub2 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  
  ggsave( filename = paste(path2,  "1C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  
  ggsave( filename = paste(path2,  "1C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)  
  ggsave( filename = paste(path2,  "1C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + ylim(0, 200)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  + ylim(0, 100)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)   + ylim(0, 30)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )

  myDiff25p.hypo_2two_sub2  = getMethylDiff(myDiff_2two_sub2, difference=10, qvalue=0.001, type="hypo" )  ## less enrich in ART
  myDiff25p.hyper_2two_sub2 = getMethylDiff(myDiff_2two_sub2, difference=10, qvalue=0.001, type="hyper")  ## more enrich in ART
  myDiff25p_2two_sub2       = getMethylDiff(myDiff_2two_sub2, difference=10, qvalue=0.05)
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
  print( diffMethPerChr(myDiff_2two_sub2,  plot=FALSE,  qvalue.cutoff=0.001, meth.cutoff=10) )
  sink()
  
  pdf( file=paste(path2, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=TRUE,  qvalue.cutoff=0.001, meth.cutoff=10) )
  dev.off()
  
  
  ########################## annotation for hypo sites.
  diffGeneAnn_2two_sub2_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub2,  "GRanges"),  gene.obj_g)
  
  sink( file=paste(path2, "4A-distribution-onGenes-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_hypo,  percentage=TRUE)
  print(diffGeneAnn_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "4B-distribution-onGenes-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffGeneAnn_2two_sub2_hypo,precedence=TRUE, main="differential methylation annotation") )
  dev.off()
  
  
  
  ##
  diffCpGann_2two_sub2_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2,"GRanges"),
                                                       cpg.obj_g$CpGi,  cpg.obj_g$shores,
                                                       feature.name="CpGi",flank.name="shores")
  
  sink( file=paste(path2, "4C-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffCpGann_2two_sub2_hypo,  percentage=TRUE)
  print(diffCpGann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "4D-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffCpGann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffrepeatann_2two_sub2_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2,"GRanges"),
                                                          myrepeat.obj_g$Repeats,  myrepeat.obj_g$shores,
                                                          feature.name="Repeats",flank.name="shores")
  
  sink( file=paste(path2, "4E-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffrepeatann_2two_sub2_hypo,  percentage=TRUE)
  print(diffrepeatann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "4F-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffrepeatann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted1Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint1.obj_g$ImprintedRegions,  flank=imprint1.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5A-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted1Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5B-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted1Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted2Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint2.obj_g$ImprintedRegions,  flank=imprint2.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5C-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted2Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5D-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted2Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted3Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint3.obj_g$ImprintedRegions,  flank=imprint3.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5E-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted3Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5F-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted3Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted4Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint4.obj_g$ImprintedRegions,  flank=imprint4.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5G-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted4Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5H-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted4Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted5Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint5.obj_g$ImprintedRegions,  flank=imprint5.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5I-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted5Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5J-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted5Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  ####################
  
  
  
  
  
  
  ########################## annotation for hyper sites.
  diffGeneAnn_2two_sub2_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub2,"GRanges"),  gene.obj_g)
  
  sink( file=paste(path2, "10A-distribution-onGenes-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_hyper,  percentage=TRUE)
  print(diffGeneAnn_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffGeneAnn_2two_sub2_hyper,precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffCpGann_2two_sub2_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2,"GRanges"),
                                                        cpg.obj_g$CpGi,  cpg.obj_g$shores,
                                                        feature.name="CpGi",flank.name="shores")
  
  sink( file=paste(path2, "10C-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffCpGann_2two_sub2_hyper,  percentage=TRUE)
  print(diffCpGann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffCpGann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffrepeatann_2two_sub2_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2,"GRanges"),
                                                           myrepeat.obj_g$Repeats,  myrepeat.obj_g$shores,
                                                           feature.name="Repeats",flank.name="shores")
  
  sink( file=paste(path2, "10E-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffrepeatann_2two_sub2_hyper,  percentage=TRUE)
  print(diffrepeatann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffrepeatann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted1Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint1.obj_g$ImprintedRegions,  flank=imprint1.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11A-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted1Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted1Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted2Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint2.obj_g$ImprintedRegions,  flank=imprint2.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11C-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted2Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted2Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted3Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint3.obj_g$ImprintedRegions,  flank=imprint3.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11E-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted3Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted3Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted4Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint4.obj_g$ImprintedRegions,  flank=imprint4.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11G-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted4Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted4Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted5Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint5.obj_g$ImprintedRegions,  flank=imprint5.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11I-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted5Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted5Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()

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
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(myobj_2two) )
  print( getSampleID(myobj_2two) )
  print( getTreatment(myobj_2two) )
  sink()
  
  tiles_2two = tileMethylCounts( myobj_2two,   win.size=binSize_temp1,   step.size=binSize_temp1,   cov.bases = binBases_temp1  )    
  meth_2two  = unite( tiles_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2_dimensions-tiles.txt", sep="/")  )
  print( tiles_2two )
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
  MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=8,   height2=5 )
  
  path_temp1_sub3 = paste(path_temp1, "3_PCA_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub3) ) { dir.create(path_temp1_sub3, recursive = TRUE) }
  MyPCA_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub3,     file2="PCA_byMethylKit_",   width2=4,   height2=4 )
  
  #path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  #if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  #MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=4,   height2=4,  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub5 = paste(path_temp1, "5_PCA_DirectlyByPrcomp", sep="/")
  if( ! file.exists(path_temp1_sub5) ) { dir.create(path_temp1_sub5, recursive = TRUE) }
  PCA_2two_sub5 <- prcomp( t(mat_2two) )
  MyPrcompObj_1_g(  prcompObj2=PCA_2two_sub5,   path2=path_temp1_sub5,   file2="PCA_DirectlyByPrcomp",  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub6 = paste(path_temp1, "6_PCA_byFactoMineR", sep="/")
  if( ! file.exists(path_temp1_sub6) ) { dir.create(path_temp1_sub6, recursive = TRUE) }
  PCA_2two_sub6 <- PCA( t(mat_2two) , graph=FALSE)
  MyPCAobj_FactoMineR_g(  PCAobj2=PCA_2two_sub6,   path2=path_temp1_sub6,   file2="PCA_byFactoMineR",  dataFrame_temp2=dataFrame_temp1   )
                                                                   
  path_temp1_sub7 = paste(path_temp1, "7_HierarchicalClustering", sep="/")
  if( ! file.exists(path_temp1_sub7) ) { dir.create(path_temp1_sub7, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=mat_2two,   path_temp1=path_temp1_sub7,   dataFrame_temp1=dataFrame_temp1  )
                                                                    
  path_temp1_sub8 = paste(path_temp1, "8_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )
                                                                             
  
} 

##################################################################################################################





##################################################################################################################
myOutDir_sub1_g = paste(outDir_g, "/1_ReadRawFiles",  sep="") 
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
               mincov    = 1,       ## >= n
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
  print(   dim(myobj_g[[i]])  )
}
sink()


continue_on_error_g <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error_g) 


myobj_nor_g <- myobj_g
##myobj_nor_g <- normalizeCoverage(myobj_g)
##################################################################################################################











##################################################################################################################
myOutDir_sub1_parents_g = paste(outDir_g, "/1_ReadRawFiles_parents",  sep="") 
if( ! file.exists(myOutDir_sub1_parents_g) ) { dir.create(myOutDir_sub1_parents_g, recursive = TRUE) }

sink( file=paste(myOutDir_sub1_parents_g, "1_length-variables.txt", sep="/") )
length( Files_All_vector_parents_g )
length( Files_All_list_parents_g )
length( mySampleID_All_vector_parents_g )
length( mySampleID_All_list_parents_g )
length( myTreatment_All_vector_parents_g )
length( myTreatment_All_list_parents_g )
length( mySex_All_vector_parents_g )
length( mySex_All_list_parents_g )
length( myTech_All_vector_parents_g )
length( myTech_All_list_parents_g )
print( "#################### Files_All_vector_parents_g ####################" )
print( Files_All_vector_parents_g )
print( "#################### Files_All_list_parents_g ####################" )
print( Files_All_list_parents_g )
print( "#################### mySampleID_All_vector_parents_g ####################" )
print( mySampleID_All_vector_parents_g )
print( "#################### mySampleID_All_list_parents_g ####################" )
print( mySampleID_All_list_parents_g )
print( "#################### myTreatment_All_vector_parents_g ####################" )
print( myTreatment_All_vector_parents_g )
print( "#################### myTreatment_All_list_parents_g ####################" )
print( myTreatment_All_list_parents_g )
print( "#################### mySex_All_vector_parents_g ####################" )
print( mySex_All_vector_parents_g )
print( "#################### mySex_All_list_parents_g ####################" )
print( mySex_All_list_parents_g )
print( "#################### myTech_All_vector_parents_g ####################" )
print( myTech_All_vector_parents_g )
print( "#################### myTech_All_list_parents_g ####################" )
print( myTech_All_list_parents_g )
sink()


# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1_parents_g, "2_theLog-of-read-inputFiles.txt", sep="/") )
myobj_parents_g = methRead( Files_All_list_parents_g,
                            sample.id = mySampleID_All_list_parents_g,
                            assembly  = "hg38",
                            treatment = myTreatment_All_vector_parents_g,
                            context   = "CpG",
                            pipeline  = "bismarkCoverage",
                            mincov    = 1,       ## >= n
                            header    = FALSE
)
sink()


sink( file=paste(myOutDir_sub1_parents_g, "3_all-rawFiles.txt", sep="/") )
print(Files_All_vector_parents_g)
print("#########################")
print(myobj_parents_g)
sink()


sink( file=paste(myOutDir_sub1_parents_g, "4_dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_parents_g)) ) {
  print( "######################" )
  print(   Files_All_vector_parents_g[i]  )
  print(   dim(myobj_parents_g[[i]])  )
}
sink()


sink( file=paste(myOutDir_sub1_parents_g, "5_dimensions-of-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_parents_g)) ) {
  print(   dim(myobj_parents_g[[i]])  )
}
sink()


continue_on_error_parents_g <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error_parents_g) 


myobj_nor_parents_g <- myobj_parents_g
##myobj_nor_parents_g <- normalizeCoverage(myobj_parents_g)
##################################################################################################################





tiles_2two_parents_g = tileMethylCounts( myobj_nor_parents_g,   win.size=1000,   step.size=1000,   cov.bases = 3  )    
meth_2two_parents_g  = unite( tiles_2two_parents_g, destrand=FALSE, mc.cores=16   )   ## 100% overlap
mat_2two_parents   = percMethylation( meth_2two )



















print("###########sign-2##########")
##################################################################################################################
  dataFrame_temp2_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
    mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
    mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
  )
  
 
  ###
  print("###########sign-2A##########")
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/2A-1_NC-vs-IVFfresh_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp2_A  )

  print("###########sign-2B##########")
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/2B-3_NC-vs-IVFfresh_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp2_A  )
  
  print("###########sign-2C##########")
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/2C-5_NC-vs-IVFfresh_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp2_A  )   
##################################################################################################################





print("###########sign-3##########")
##################################################################################################################
  dataFrame_temp3_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_frozen_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_IVF_frozen_g),
    mysex       = c(Sex_NC_g,  Sex_IVF_frozen_g),
    mytech      = c(Tech_NC_g, Tech_IVF_frozen_g)    
  )
  
  
  ###
  print("###########sign-3A##########")
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/3A-1_NC-vs-IVFfrozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp3_A  )
  print("###########sign-3B##########") 
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/3B-3_NC-vs-IVFfrozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp3_A  )
  print("###########sign-3C##########")
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/3C-5_NC-vs-IVFfrozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp3_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp4_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_ICSI_fresh_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_ICSI_fresh_g),
    mysex       = c(Sex_NC_g,  Sex_ICSI_fresh_g),
    mytech      = c(Tech_NC_g, Tech_ICSI_fresh_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/4A-1_NC-vs-ICSIfresh_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp4_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/4B-3_NC-vs-ICSIfresh_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp4_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/4C-5_NC-vs-ICSIfresh_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp4_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp5_A <- data.frame(
    mysampleID  = c(mySampleID_NC_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_NC_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_NC_g,  Sex_ICSI_frozen_g),
    mytech      = c(Tech_NC_g, Tech_ICSI_frozen_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/5A-1_NC-vs-ICSIfrozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp5_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/5B-3_NC-vs-ICSIfrozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp5_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/5C-5_NC-vs-ICSIfrozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp5_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp6_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_fresh_g,  mySampleID_IVF_frozen_g),
    mytreatment = c(myTreatment_IVF_fresh_g, myTreatment_IVF_frozen_g),
    mysex       = c(Sex_IVF_fresh_g,         Sex_IVF_frozen_g),
    mytech      = c(Tech_IVF_fresh_g,        Tech_IVF_frozen_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/6A-1_IVF_fresh_vs_IVF_frozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp6_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/6B-3_IVF_fresh_vs_IVF_frozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp6_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/6C-5_IVF_fresh_vs_IVF_frozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp6_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp7_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_fresh_g,  mySampleID_ICSI_fresh_g),
    mytreatment = c(myTreatment_IVF_fresh_g, myTreatment_ICSI_fresh_g),
    mysex       = c(Sex_IVF_fresh_g,         Sex_ICSI_fresh_g),
    mytech      = c(Tech_IVF_fresh_g,        Tech_ICSI_fresh_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/7A-1_IVF_fresh_vs_ICSI_fresh_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp7_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/7B-3_IVF_fresh_vs_ICSI_fresh_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp7_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/7C-5_IVF_fresh_vs_ICSI_fresh_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp7_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp8_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_fresh_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_IVF_fresh_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_IVF_fresh_g,         Sex_ICSI_frozen_g),
    mytech      = c(Tech_IVF_fresh_g,        Tech_ICSI_frozen_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/8A-1_IVF_fresh_vs_ICSI_frozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp8_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/8B-3_IVF_fresh_vs_ICSI_frozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp8_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/8C-5_IVF_fresh_vs_ICSI_frozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp8_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp9_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_frozen_g,  mySampleID_ICSI_fresh_g),
    mytreatment = c(myTreatment_IVF_frozen_g, myTreatment_ICSI_fresh_g),
    mysex       = c(Sex_IVF_frozen_g,         Sex_ICSI_fresh_g),
    mytech      = c(Tech_IVF_frozen_g,        Tech_ICSI_fresh_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/9A-1_IVF_frozen_vs_ICSI_fresh_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp9_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/9B-3_IVF_frozen_vs_ICSI_fresh_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp9_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/9C-5_IVF_frozen_vs_ICSI_fresh_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp9_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp10_A <- data.frame(
    mysampleID  = c(mySampleID_IVF_frozen_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_IVF_frozen_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_IVF_frozen_g,         Sex_ICSI_frozen_g),
    mytech      = c(Tech_IVF_frozen_g,        Tech_ICSI_frozen_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/10A-1_IVF_frozen_vs_ICSI_frozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp10_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/10B-3_IVF_frozen_vs_ICSI_frozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp10_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/10C-5_IVF_frozen_vs_ICSI_frozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp10_A  )   
##################################################################################################################





##################################################################################################################
  dataFrame_temp11_A <- data.frame(
    mysampleID  = c(mySampleID_ICSI_fresh_g,  mySampleID_ICSI_frozen_g),
    mytreatment = c(myTreatment_ICSI_fresh_g, myTreatment_ICSI_frozen_g),
    mysex       = c(Sex_ICSI_fresh_g,         Sex_ICSI_frozen_g),
    mytech      = c(Tech_ICSI_fresh_g,        Tech_ICSI_frozen_g)    
  )
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/11A-1_ICSI_fresh_vs_ICSI_frozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp11_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/11B-3_ICSI_fresh_vs_ICSI_frozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp11_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/11C-5_ICSI_fresh_vs_ICSI_frozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp11_A  )   
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
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/12A-1_NC_vs_IVF_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp12_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/12B-3_NC_vs_IVF_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp12_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/12C-5_NC_vs_IVF_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp12_A  )   
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
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/13A-1_NC_vs_ICSI_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp13_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/13B-3_NC_vs_ICSI_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp13_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/13C-5_NC_vs_ICSI_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp13_A  )   
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
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/14A-1_NC_vs_fresh_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp14_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/14B-3_NC_vs_fresh_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp14_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/14C-5_NC_vs_fresh_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp14_A  )   
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
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/15A-1_NC_vs_frozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp15_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/15B-3_NC_vs_frozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp15_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/15C-5_NC_vs_frozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp15_A  )   
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
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/16A-1_NC_vs_ART_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp16_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/16B-3_NC_vs_ART_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp16_A  )

  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/16C-5_NC_vs_ART_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp16_A  )
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
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/17A-1_IVF_vs_ICSI_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp17_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/17B-3_IVF_vs_ICSI_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp17_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/17C-5_IVF_vs_ICSI_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp17_A  )   
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
  
  
  ###
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/18A-1_fresh_vs_frozen_1kb_1bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 1, dataFrame_temp1 = dataFrame_temp18_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/18B-3_fresh_vs_frozen_1kb_3bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 3, dataFrame_temp1 = dataFrame_temp18_A  )
  
  myMainFunction_1_g(  myobj_temp1 = myobj_nor_g,   path_temp1 = paste(outDir_g, "/18C-5_fresh_vs_frozen_1kb_5bases",  sep="") ,   
                       binSize_temp1 = 1000,   binBases_temp1 = 5, dataFrame_temp1 = dataFrame_temp18_A  )   
##################################################################################################################



