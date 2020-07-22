##  http://compbio.ucsd.edu/making-heat-maps-r/

options(jupyter.plot_mimetypes = c("text/plain", "image/png" ))


gLogCpmData = as.matrix(read.table("heatmap_test_matrix.txt"))
gLogCpmData

gAnnotationData = read.table("heatmap_test_annotation.txt")
gAnnotationData


dim(gLogCpmData)
dim(gAnnotationData)


# Make helper function to map metadata category to color
mapDrugToColor<-function(annotations){
  colorsVector = ifelse(annotations["subject_drug"]=="MiracleDrugA", "blue", 
                 ifelse(annotations["subject_drug"]=="MiracleDrugB", "green", "red"))
  return(colorsVector)
}



mapDrugToColor(gAnnotationData)






## heatmap, heatmap is the built-in option for heat maps in R
# Test heatmap with column annotations
testHeatmap<-function(logCPM, annotations) {    
  sampleColors = mapDrugToColor(annotations)
  heatmap(logCPM, margins=c(5,8), ColSideColors=sampleColors)
}
testHeatmap(gLogCpmData, gAnnotationData)







## heatmap.2 is an "enhanced" heat map function from the add-on package gplots (not to be confused with ggplot!)
library(gplots)
# Test heatmap.2 with column annotations and custom legend text
testHeatmap2<-function(logCPM, annotations) {    
  sampleColors = mapDrugToColor(annotations)
  heatmap.2(logCPM, margins=c(5,8), ColSideColors=sampleColors,
            key.xlab="log CPM",
            key=TRUE, symkey=FALSE, density.info="none", trace="none")
}
testHeatmap2(gLogCpmData, gAnnotationData)









##aheatmap, which stands for "annotated heatmap", is a heat map plotting function from the NMF package
library(NMF)
# Test aheatmap with column annotations
testAheatmap<-function(logCPM, annotations) {    
  aheatmap(logCPM, annCol=annotations[
    "subject_drug"])
}
testAheatmap(gLogCpmData, gAnnotationData)









## pheatmap, where the "p" stands for "pretty", is the sole function of the package pheatmap
library(pheatmap)
# Test pheatmap with two annotation options
testPheatmap<-function(logCPM, annotations) {    
  drug_info = data.frame(annotations[,"subject_drug"])
  rownames(drug_info) = annotations[["sample_name"]]
  
  # Assign the column annotation straight from 
  # the input annotation dataframe
  pheatmap(logCPM, annotation_col=drug_info, 
           annotation_names_row=FALSE,
           annotation_names_col=FALSE,
           fontsize_col=5)
  
  # Assign the column annotation to an intermediate
  # variable first in order to change the name 
  # pheatmap uses for its legend
  subject_drug = annotations[["subject_drug"]]
  drug_df = data.frame(subject_drug)
  rownames(drug_df) = annotations[["sample_name"]]
  
  pheatmap(logCPM, annotation_col=drug_df, 
           annotation_names_row=FALSE,
           annotation_names_col=FALSE,
           fontsize_col=5)
}
testPheatmap(gLogCpmData, gAnnotationData)








##heatmap3 is the central function of the heatmap3 package. 
##Beware that this is different from "heatmap.3", of which there are numerous versions

library(heatmap3)
# Test heatmap3 with several annotation options
testHeatmap3<-function(logCPM, annotations) {    
  sampleColors = mapDrugToColor(annotations)
  
  # Assign just column annotations
  heatmap3(logCPM, margins=c(5,8), ColSideColors=sampleColors) 
  
  # Assign column annotations and make a custom legend for them
  heatmap3(logCPM, margins=c(5,8), ColSideColors=sampleColors, 
           legendfun=function()showLegend(legend=c("MiracleDrugA", 
                                                   "MiracleDrugB", "?"), col=c("blue", "green", "red"), cex=1.5))
  
  # Assign column annotations as a mini-graph instead of colors,
  # and use the built-in labeling for them
  ColSideAnn <- data.frame(Drug=annotations[["subject_drug"]])
  heatmap3(logCPM,ColSideAnn=ColSideAnn,
           #ColSideFun=function(x)showAnn(x),
           ColSideWidth=0.8)
}

testHeatmap3(logCPM=gLogCpmData, annotations=gAnnotationData)












