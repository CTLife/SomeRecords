library(ggplot2) 
library(scales)





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





###########################################################################################################################################################################
rawMatrix_1 <- read.table("1.9.1_BA9.permutation.txt",  header=F,   sep=" " )  
rawMatrix_2 <- read.table("1.9.2_BA24.permutation.txt", header=F,   sep=" " ) 
rawMatrix_3 <- read.table("1.3.3_C.permutation.txt",    header=F,   sep=" " ) 
rawMatrix_4 <- read.table("1.10.4_H.permutation.txt",    header=F,   sep=" " ) 
rawMatrix_5 <- read.table("1.3.5_T.permutation.txt",    header=F,   sep=" " ) 

dim(rawMatrix_1 )
dim(rawMatrix_2 )
dim(rawMatrix_3 )
dim(rawMatrix_4 )
dim(rawMatrix_5 )


shared <- read.table("shared.phenotype.txt",    header=F,   sep=" " ) 
dim(shared)
shared = as.vector(shared[,1])
length(shared)



############# compare row order
phe_id_1 = as.vector( rawMatrix_1[,1] )
phe_id_2 = as.vector( rawMatrix_2[,1] )
phe_id_3 = as.vector( rawMatrix_3[,1] )
phe_id_4 = as.vector( rawMatrix_4[,1] )
phe_id_5 = as.vector( rawMatrix_5[,1] )

bool_2 = (phe_id_1 == phe_id_2)
bool_3 = (phe_id_1 == phe_id_3)
bool_4 = (phe_id_1 == phe_id_4)
bool_5 = (phe_id_1 == phe_id_5)

length( bool_2 )
length( bool_3 )
length( bool_4 )
length( bool_5 )

length( bool_2[bool_2] )
length( bool_3[bool_3] )
length( bool_4[bool_4] )
length( bool_5[bool_5] )


phe_id_1 = as.vector( rawMatrix_1[,21] )
phe_id_2 = as.vector( rawMatrix_2[,21] )
phe_id_3 = as.vector( rawMatrix_3[,21] )
phe_id_4 = as.vector( rawMatrix_4[,21] )
phe_id_5 = as.vector( rawMatrix_5[,21] )

phe_id_1[  phe_id_1 == "NA" ]
rawMatrix_1_noNA = rawMatrix_1[(!is.na(phe_id_1)),]
rawMatrix_2_noNA = rawMatrix_2[(!is.na(phe_id_2)),]
rawMatrix_3_noNA = rawMatrix_3[(!is.na(phe_id_3)),]
rawMatrix_4_noNA = rawMatrix_4[(!is.na(phe_id_4)),]
rawMatrix_5_noNA = rawMatrix_5[(!is.na(phe_id_5)),]

dim(rawMatrix_1_noNA )
dim(rawMatrix_2_noNA )
dim(rawMatrix_3_noNA )
dim(rawMatrix_4_noNA )
dim(rawMatrix_5_noNA )






phe_id = as.vector( rawMatrix_1_noNA[,1] )
sharedBOOL  =  phe_id   %in%   shared

length(phe_id)
length(sharedBOOL)

length(sharedBOOL[sharedBOOL])
length(phe_id[sharedBOOL])
length(  shared  )
 




distance_1 = as.numeric( rawMatrix_1_noNA[,7] )
distance_2 = as.numeric( rawMatrix_2_noNA[,7] )
distance_3 = as.numeric( rawMatrix_3_noNA[,7] )
distance_4 = as.numeric( rawMatrix_4_noNA[,7] )
distance_5 = as.numeric( rawMatrix_5_noNA[,7] )

nom_pval_1 = as.numeric( rawMatrix_1_noNA[,16] )
nom_pval_2 = as.numeric( rawMatrix_2_noNA[,16] )
nom_pval_3 = as.numeric( rawMatrix_3_noNA[,16] )
nom_pval_4 = as.numeric( rawMatrix_4_noNA[,16] )
nom_pval_5 = as.numeric( rawMatrix_5_noNA[,16] )

r_squared_1 = as.numeric( rawMatrix_1_noNA[,17] )
r_squared_2 = as.numeric( rawMatrix_2_noNA[,17] )
r_squared_3 = as.numeric( rawMatrix_3_noNA[,17] )
r_squared_4 = as.numeric( rawMatrix_4_noNA[,17] )
r_squared_5 = as.numeric( rawMatrix_5_noNA[,17] )

slope_1 = as.numeric( rawMatrix_1_noNA[,18] )
slope_2 = as.numeric( rawMatrix_2_noNA[,18] )
slope_3 = as.numeric( rawMatrix_3_noNA[,18] )
slope_4 = as.numeric( rawMatrix_4_noNA[,18] )
slope_5 = as.numeric( rawMatrix_5_noNA[,18] )

slope_se_1 = as.numeric( rawMatrix_1_noNA[,19] )
slope_se_2 = as.numeric( rawMatrix_2_noNA[,19] )
slope_se_3 = as.numeric( rawMatrix_3_noNA[,19] )
slope_se_4 = as.numeric( rawMatrix_4_noNA[,19] )
slope_se_5 = as.numeric( rawMatrix_5_noNA[,19] )

adj_emp_pval_1 = as.numeric( rawMatrix_1_noNA[,20] )
adj_emp_pval_2 = as.numeric( rawMatrix_2_noNA[,20] )
adj_emp_pval_3 = as.numeric( rawMatrix_3_noNA[,20] )
adj_emp_pval_4 = as.numeric( rawMatrix_4_noNA[,20] )
adj_emp_pval_5 = as.numeric( rawMatrix_5_noNA[,20] )

adj_beta_pval_1 = as.numeric( rawMatrix_1_noNA[,21] )
adj_beta_pval_2 = as.numeric( rawMatrix_2_noNA[,21] )
adj_beta_pval_3 = as.numeric( rawMatrix_3_noNA[,21] )
adj_beta_pval_4 = as.numeric( rawMatrix_4_noNA[,21] )
adj_beta_pval_5 = as.numeric( rawMatrix_5_noNA[,21] )



length( nom_pval_1[nom_pval_1>0.1] )
length( nom_pval_1[nom_pval_1>0.05] )
nom_pval_1[nom_pval_1>0.05]  = 0.05
nom_pval_2[nom_pval_2>0.05]  = 0.05
nom_pval_3[nom_pval_3>0.05]  = 0.05
nom_pval_4[nom_pval_4>0.05]  = 0.05
nom_pval_5[nom_pval_5>0.05]  = 0.05






AllResults_g <- "Z-FinalFigures"
if( ! file.exists(AllResults_g) ) { dir.create(path=AllResults_g, recursive = TRUE) }





my_line <- function(x,y,  ...){
  points(x,y,  ...)
  abline(a = 0,b = 1, col="red", ...)
  abline(h = 1, col="red", ...)
  abline(v = 1, col="red", ...)
}

matrix1 = cbind(distance_1, nom_pval_1, r_squared_1, slope_1, slope_se_1, adj_emp_pval_1, adj_beta_pval_1)
matrix2 = cbind(distance_2, nom_pval_2, r_squared_2, slope_2, slope_se_2, adj_emp_pval_2, adj_beta_pval_2)
matrix3 = cbind(distance_3, nom_pval_3, r_squared_3, slope_3, slope_se_3, adj_emp_pval_3, adj_beta_pval_3)
matrix4 = cbind(distance_4, nom_pval_4, r_squared_4, slope_4, slope_se_4, adj_emp_pval_4, adj_beta_pval_4)
matrix5 = cbind(distance_5, nom_pval_5, r_squared_5, slope_5, slope_se_5, adj_emp_pval_5, adj_beta_pval_5)

pdf( paste(AllResults_g, "1.same-group.different-parameters.pdf" , sep="/" )  )
pairs( matrix1  ,   col = alpha("red", 0.1),    pch = 19,   cex = 0.01 )  
pairs( matrix2  ,   col = alpha("red", 0.1),    pch = 19,   cex = 0.01 ) 
pairs( matrix3  ,   col = alpha("red", 0.1),    pch = 19,   cex = 0.01 ) 
pairs( matrix4  ,   col = alpha("red", 0.1),    pch = 19,   cex = 0.01 ) 
pairs( matrix5  ,   col = alpha("red", 0.1),    pch = 19,   cex = 0.01 ) 
dev.off()








###############################################################################3
plot(distance_1, nom_pval_1   , type = "p",  xlim = NULL, ylim = c(0, 0.05), col="red", pch = 19, cex=0.2)
plot(distance_1, r_squared_1  , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2) 
plot(distance_1, slope_1      , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(distance_1, slope_se_1      , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(distance_1, adj_emp_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(distance_1, adj_beta_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)

plot(nom_pval_1, r_squared_1  , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2) 
plot(nom_pval_1, r_squared_1  , type = "p",  xlim = c(0, 0.05), ylim = NULL, col="red", pch = 19, cex=0.2) 
plot(nom_pval_1, slope_1      , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(nom_pval_1, slope_1      , type = "p",  xlim = c(0, 0.01), ylim = NULL, col="red", pch = 19, cex=0.2)
plot(nom_pval_1, slope_se_1      , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(nom_pval_1, slope_se_1      , type = "p",  xlim = c(0, 0.01), ylim = NULL, col="red", pch = 19, cex=0.2)
plot(nom_pval_1, adj_emp_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(nom_pval_1, adj_emp_pval_1     , type = "p",  xlim = c(0, 0.01), ylim = NULL, col="red", pch = 19, cex=0.2)
plot(nom_pval_1, adj_beta_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(nom_pval_1, adj_beta_pval_1     , type = "p",  xlim = c(0, 0.01), ylim = NULL, col="red", pch = 19, cex=0.2)

plot(r_squared_1, slope_1      , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(r_squared_1, slope_se_1      , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(r_squared_1, adj_emp_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(r_squared_1, adj_beta_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)

plot(slope_1, slope_se_1      , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(slope_1, adj_emp_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(slope_1, adj_beta_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)

plot(slope_se_1, adj_emp_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
plot(slope_se_1, adj_beta_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)

plot(adj_emp_pval_1, adj_beta_pval_1     , type = "p",  xlim = NULL, ylim = NULL, col="red", pch = 19, cex=0.2)
###############################################################################







library(ggplot2)
xaxis = seq(from=0, to=1, by=0.001)
length(xaxis)



##### number of kept g-m6A peaks
num_1 = xaxis
num_2 = xaxis
num_3 = xaxis
num_4 = xaxis
num_5 = xaxis

for(i in c(1:length(xaxis) ) ) {
    Bool_1 =  adj_beta_pval_1 < xaxis[i]
    Bool_2 =  adj_beta_pval_2 < xaxis[i]
    Bool_3 =  adj_beta_pval_3 < xaxis[i]
    Bool_4 =  adj_beta_pval_4 < xaxis[i]
    Bool_5 =  adj_beta_pval_5 < xaxis[i]
    num_1[i] = length( adj_beta_pval_1[ Bool_1 ] )
    num_2[i] = length( adj_beta_pval_2[ Bool_2 ] )
    num_3[i] = length( adj_beta_pval_3[ Bool_3 ] )
    num_4[i] = length( adj_beta_pval_4[ Bool_4 ] )
    num_5[i] = length( adj_beta_pval_5[ Bool_5 ] )
}

DF_A = data.frame( x=c(xaxis, xaxis, xaxis, xaxis, xaxis), y=c(num_1, num_2, num_3, num_4, num_5), type=c( rep("1_BA9", length(xaxis)) ,   rep("2_BA24", length(xaxis)) ,   rep("3_C", length(xaxis)) ,   rep("4_H", length(xaxis)) ,   rep("5_T", length(xaxis)) )  )
FigureTemp3 = ggplot(DF_A, aes(x=x, y=y,  color=type)) + geom_point(size=0.1)
MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3,  path1=AllResults_g, fileName1= "number.of.kept.g-m6A.peaks",  height1=5, width1=7)

DF_B = data.frame( x=c(xaxis[1:200], xaxis[1:200], xaxis[1:200], xaxis[1:200], xaxis[1:200]), y=c(num_1[1:200], num_2[1:200], num_3[1:200], num_4[1:200], num_5[1:200]), type=c( rep("1_BA9", 200) ,   rep("2_BA24",200) ,   rep("3_C", 200) ,   rep("4_H", 200) ,   rep("5_T", 200) )  )
FigureTemp3 = ggplot(DF_B, aes(x=x, y=y,  color=type)) + geom_point(size=0.1)
MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3,  path1=AllResults_g, fileName1= "number.of.kept.g-m6A.peaks.2",  height1=5, width1=7)





##### number of shared g-m6A peaks in 2 groups at least
num_1 = xaxis
num_2 = xaxis
num_3 = xaxis
num_4 = xaxis
num_5 = xaxis

for(i in c(1:length(xaxis) ) ) {
  Bool_1 =  adj_beta_pval_1 < xaxis[i]
  Bool_2 =  adj_beta_pval_2 < xaxis[i]
  Bool_3 =  adj_beta_pval_3 < xaxis[i]
  Bool_4 =  adj_beta_pval_4 < xaxis[i]
  Bool_5 =  adj_beta_pval_5 < xaxis[i]
  num_1[i] = length( adj_beta_pval_1[ Bool_1  &  ( Bool_2 | Bool_3 | Bool_4 | Bool_5 ) ] ) / length( adj_beta_pval_1[ Bool_1 ] )
  num_2[i] = length( adj_beta_pval_2[ Bool_2  &  ( Bool_1 | Bool_3 | Bool_4 | Bool_5 ) ] ) / length( adj_beta_pval_2[ Bool_2 ] )
  num_3[i] = length( adj_beta_pval_3[ Bool_3  &  ( Bool_2 | Bool_1 | Bool_4 | Bool_5 ) ] ) / length( adj_beta_pval_3[ Bool_3 ] )
  num_4[i] = length( adj_beta_pval_4[ Bool_4  &  ( Bool_2 | Bool_3 | Bool_1 | Bool_5 ) ] ) / length( adj_beta_pval_4[ Bool_4 ] )
  num_5[i] = length( adj_beta_pval_5[ Bool_5  &  ( Bool_2 | Bool_3 | Bool_4 | Bool_1 ) ] ) / length( adj_beta_pval_5[ Bool_5 ] )
}

DF_A = data.frame( x=c(xaxis, xaxis, xaxis, xaxis, xaxis), y=c(num_1, num_2, num_3, num_4, num_5), type=c( rep("1_BA9", length(xaxis)) ,   rep("2_BA24", length(xaxis)) ,   rep("3_C", length(xaxis)) ,   rep("4_H", length(xaxis)) ,   rep("5_T", length(xaxis)) )  )
FigureTemp3 = ggplot(DF_A, aes(x=x, y=y,  color=type)) + geom_point(size=0.1)
MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3,  path1=AllResults_g, fileName1= "shared_in_2_groups",  height1=5, width1=7)

DF_B = data.frame( x=c(xaxis[1:200], xaxis[1:200], xaxis[1:200], xaxis[1:200], xaxis[1:200]), y=c(num_1[1:200], num_2[1:200], num_3[1:200], num_4[1:200], num_5[1:200]), type=c( rep("1_BA9", 200) ,   rep("2_BA24",200) ,   rep("3_C", 200) ,   rep("4_H", 200) ,   rep("5_T", 200) )  )
FigureTemp3 = ggplot(DF_B, aes(x=x, y=y,  color=type)) + geom_point(size=0.1)
MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3,  path1=AllResults_g, fileName1= "shared_in_2_groups.2",  height1=5, width1=7)







##### number of shared g-m6A peaks in 5 groups at least
num_1 = xaxis
num_2 = xaxis
num_3 = xaxis
num_4 = xaxis
num_5 = xaxis

for(i in c(1:length(xaxis) ) ) {
  Bool_1 =  adj_beta_pval_1 < xaxis[i]
  Bool_2 =  adj_beta_pval_2 < xaxis[i]
  Bool_3 =  adj_beta_pval_3 < xaxis[i]
  Bool_4 =  adj_beta_pval_4 < xaxis[i]
  Bool_5 =  adj_beta_pval_5 < xaxis[i]
  num_1[i] = length( adj_beta_pval_1[ Bool_1  &  ( Bool_2 & Bool_3 & Bool_4 & Bool_5 ) ] )  
  num_2[i] = length( adj_beta_pval_2[ Bool_2  &  ( Bool_1 & Bool_3 & Bool_4 & Bool_5 ) ] )  
  num_3[i] = length( adj_beta_pval_3[ Bool_3  &  ( Bool_2 & Bool_1 & Bool_4 & Bool_5 ) ] )  
  num_4[i] = length( adj_beta_pval_4[ Bool_4  &  ( Bool_2 & Bool_3 & Bool_1 & Bool_5 ) ] )  
  num_5[i] = length( adj_beta_pval_5[ Bool_5  &  ( Bool_2 & Bool_3 & Bool_4 & Bool_1 ) ] )  
}

DF_A = data.frame( x=c(xaxis, xaxis, xaxis, xaxis, xaxis), y=c(num_1, num_2, num_3, num_4, num_5), type=c( rep("1_BA9", length(xaxis)) ,   rep("2_BA24", length(xaxis)) ,   rep("3_C", length(xaxis)) ,   rep("4_H", length(xaxis)) ,   rep("5_T", length(xaxis)) )  )
FigureTemp3 = ggplot(DF_A, aes(x=x, y=y,  color=type)) + geom_point(size=0.1)
MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3,  path1=AllResults_g, fileName1= "shared_in_5_groups",  height1=5, width1=7)

DF_B = data.frame( x=c(xaxis[1:200], xaxis[1:200], xaxis[1:200], xaxis[1:200], xaxis[1:200]), y=c(num_1[1:200], num_2[1:200], num_3[1:200], num_4[1:200], num_5[1:200]), type=c( rep("1_BA9", 200) ,   rep("2_BA24",200) ,   rep("3_C", 200) ,   rep("4_H", 200) ,   rep("5_T", 200) )  )
FigureTemp3 = ggplot(DF_B, aes(x=x, y=y,  color=type)) + geom_point(size=0.1)
MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3,  path1=AllResults_g, fileName1= "shared_in_5_groups.2",  height1=5, width1=7)













###################################################################
t1=0.2
length(adj_beta_pval_1[adj_beta_pval_1<t1]) 
length(adj_beta_pval_2[adj_beta_pval_2<t1]) 
length(adj_beta_pval_3[adj_beta_pval_3<t1]) 
length(adj_beta_pval_4[adj_beta_pval_4<t1]) 
length(adj_beta_pval_5[adj_beta_pval_5<t1]) 

length(adj_beta_pval_1[adj_beta_pval_1<t1 & sharedBOOL]) 
length(adj_beta_pval_2[adj_beta_pval_2<t1 & sharedBOOL]) 
length(adj_beta_pval_3[adj_beta_pval_3<t1 & sharedBOOL]) 
length(adj_beta_pval_4[adj_beta_pval_4<t1 & sharedBOOL]) 
length(adj_beta_pval_5[adj_beta_pval_5<t1 & sharedBOOL]) 

length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_4<t1  &   adj_beta_pval_3<t1  &   adj_beta_pval_2<t1  &   adj_beta_pval_1<t1    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_4<t1  &   adj_beta_pval_3<t1  &   adj_beta_pval_2<t1      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_4<t1  &   adj_beta_pval_3<t1      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_4<t1       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_3<t1       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_2<t1       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_1<t1       ]) 


length(adj_beta_pval_1[ adj_beta_pval_1<t1  &   (adj_beta_pval_4<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_5<t1)    ]) / length(adj_beta_pval_1[adj_beta_pval_1<t1]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t1  &   (adj_beta_pval_4<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_5<t1  |   adj_beta_pval_1<t1)    ]) / length(adj_beta_pval_2[adj_beta_pval_2<t1])                         
length(adj_beta_pval_3[ adj_beta_pval_3<t1  &   (adj_beta_pval_4<t1  |   adj_beta_pval_5<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_1<t1)    ]) / length(adj_beta_pval_3[adj_beta_pval_3<t1]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t1  &   (adj_beta_pval_5<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_1<t1)    ]) / length(adj_beta_pval_4[adj_beta_pval_4<t1]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   (adj_beta_pval_4<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_1<t1)    ]) / length(adj_beta_pval_5[adj_beta_pval_5<t1]) 

 
hist( distance_1[  adj_beta_pval_5<t1  &   adj_beta_pval_4<t1  &   adj_beta_pval_3<t1  &   adj_beta_pval_2<t1  &   adj_beta_pval_1<t1  ]    )
hist( distance_1[  adj_beta_pval_1<t1  &   (adj_beta_pval_4<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_5<t1)   ] , breaks=200  )

hist( distance_1[  adj_beta_pval_1 < 2 ]/1000 , breaks=200  )
hist( distance_1[  adj_beta_pval_1 < 0.5 ]/1000 , breaks=200  )
hist( distance_1[  adj_beta_pval_1 < 0.1 ]/1000 , breaks=200  )
hist( distance_1[  adj_beta_pval_1 < 0.01 ]/1000 , breaks=200  )



###################################################################
t2=0.1
length(adj_beta_pval_1[adj_beta_pval_1<t2]) 
length(adj_beta_pval_2[adj_beta_pval_2<t2]) 
length(adj_beta_pval_3[adj_beta_pval_3<t2]) 
length(adj_beta_pval_4[adj_beta_pval_4<t2]) 
length(adj_beta_pval_5[adj_beta_pval_5<t2]) 

length(adj_beta_pval_1[adj_beta_pval_1<t2 & sharedBOOL]) 
length(adj_beta_pval_2[adj_beta_pval_2<t2 & sharedBOOL]) 
length(adj_beta_pval_3[adj_beta_pval_3<t2 & sharedBOOL]) 
length(adj_beta_pval_4[adj_beta_pval_4<t2 & sharedBOOL]) 
length(adj_beta_pval_5[adj_beta_pval_5<t2 & sharedBOOL]) 

length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   adj_beta_pval_4<t2  &   adj_beta_pval_3<t2  &   adj_beta_pval_2<t2  &   adj_beta_pval_1<t2    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   adj_beta_pval_4<t2  &   adj_beta_pval_3<t2  &   adj_beta_pval_2<t2      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   adj_beta_pval_4<t2  &   adj_beta_pval_3<t2      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   adj_beta_pval_4<t2       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   adj_beta_pval_3<t2       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   adj_beta_pval_2<t2       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   adj_beta_pval_1<t2       ]) 

length(adj_beta_pval_1[ adj_beta_pval_1<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_5<t2)    ]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_5<t2  |   adj_beta_pval_1<t2)    ]) 
length(adj_beta_pval_3[ adj_beta_pval_3<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_5<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_1<t2)    ]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t2  &   (adj_beta_pval_5<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_1<t2)    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_1<t2)    ]) 

t1 = 0.1
length(adj_beta_pval_1[ adj_beta_pval_1<t2  &   (adj_beta_pval_4<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_5<t1)    ]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t2  &   (adj_beta_pval_4<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_5<t1  |   adj_beta_pval_1<t1)    ]) 
length(adj_beta_pval_3[ adj_beta_pval_3<t2  &   (adj_beta_pval_4<t1  |   adj_beta_pval_5<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_1<t1)    ]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t2  &   (adj_beta_pval_5<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_1<t1)    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   (adj_beta_pval_4<t1  |   adj_beta_pval_3<t1  |   adj_beta_pval_2<t1  |   adj_beta_pval_1<t1)    ]) 

length(adj_beta_pval_1[ adj_beta_pval_1<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_5<t2)    ]) / length(adj_beta_pval_1[adj_beta_pval_1<t2]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_5<t2  |   adj_beta_pval_1<t2)    ]) / length(adj_beta_pval_2[adj_beta_pval_2<t2])                         
length(adj_beta_pval_3[ adj_beta_pval_3<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_5<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_1<t2)    ]) / length(adj_beta_pval_3[adj_beta_pval_3<t2]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t2  &   (adj_beta_pval_5<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_1<t2)    ]) / length(adj_beta_pval_4[adj_beta_pval_4<t2]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t2  &   (adj_beta_pval_4<t2  |   adj_beta_pval_3<t2  |   adj_beta_pval_2<t2  |   adj_beta_pval_1<t2)    ]) / length(adj_beta_pval_5[adj_beta_pval_5<t2]) 










###################################################################
t3=0.05
length(adj_beta_pval_1[adj_beta_pval_1<t3]) 
length(adj_beta_pval_2[adj_beta_pval_2<t3]) 
length(adj_beta_pval_3[adj_beta_pval_3<t3]) 
length(adj_beta_pval_4[adj_beta_pval_4<t3]) 
length(adj_beta_pval_5[adj_beta_pval_5<t3]) 

length(adj_beta_pval_1[adj_beta_pval_1<t3 & sharedBOOL]) 
length(adj_beta_pval_2[adj_beta_pval_2<t3 & sharedBOOL]) 
length(adj_beta_pval_3[adj_beta_pval_3<t3 & sharedBOOL]) 
length(adj_beta_pval_4[adj_beta_pval_4<t3 & sharedBOOL]) 
length(adj_beta_pval_5[adj_beta_pval_5<t3 & sharedBOOL]) 

length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   adj_beta_pval_4<t3  &   adj_beta_pval_3<t3  &   adj_beta_pval_2<t3  &   adj_beta_pval_1<t3    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   adj_beta_pval_4<t3  &   adj_beta_pval_3<t3  &   adj_beta_pval_2<t3      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   adj_beta_pval_4<t3  &   adj_beta_pval_3<t3      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   adj_beta_pval_4<t3       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   adj_beta_pval_3<t3       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   adj_beta_pval_2<t3       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   adj_beta_pval_1<t3       ]) 

length(adj_beta_pval_1[ adj_beta_pval_1<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_5<t3)    ]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_5<t3  |   adj_beta_pval_1<t3)    ]) 
length(adj_beta_pval_3[ adj_beta_pval_3<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_5<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_1<t3)    ]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t3  &   (adj_beta_pval_5<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_1<t3)    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_1<t3)    ]) 

length(adj_beta_pval_1[ adj_beta_pval_1<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_5<t3)    ]) / length(adj_beta_pval_1[adj_beta_pval_1<t3]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_5<t3  |   adj_beta_pval_1<t3)    ]) / length(adj_beta_pval_2[adj_beta_pval_2<t3])                         
length(adj_beta_pval_3[ adj_beta_pval_3<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_5<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_1<t3)    ]) / length(adj_beta_pval_3[adj_beta_pval_3<t3]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t3  &   (adj_beta_pval_5<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_1<t3)    ]) / length(adj_beta_pval_4[adj_beta_pval_4<t3]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t3  &   (adj_beta_pval_4<t3  |   adj_beta_pval_3<t3  |   adj_beta_pval_2<t3  |   adj_beta_pval_1<t3)    ]) / length(adj_beta_pval_5[adj_beta_pval_5<t3]) 









###################################################################
t4=0.01
length(adj_beta_pval_1[adj_beta_pval_1<t4]) 
length(adj_beta_pval_2[adj_beta_pval_2<t4]) 
length(adj_beta_pval_3[adj_beta_pval_3<t4]) 
length(adj_beta_pval_4[adj_beta_pval_4<t4]) 
length(adj_beta_pval_5[adj_beta_pval_5<t4]) 

length(adj_beta_pval_1[adj_beta_pval_1<t4 & sharedBOOL]) 
length(adj_beta_pval_2[adj_beta_pval_2<t4 & sharedBOOL]) 
length(adj_beta_pval_3[adj_beta_pval_3<t4 & sharedBOOL]) 
length(adj_beta_pval_4[adj_beta_pval_4<t4 & sharedBOOL]) 
length(adj_beta_pval_5[adj_beta_pval_5<t4 & sharedBOOL]) 

length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   adj_beta_pval_4<t4  &   adj_beta_pval_3<t4  &   adj_beta_pval_2<t4  &   adj_beta_pval_1<t4    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   adj_beta_pval_4<t4  &   adj_beta_pval_3<t4  &   adj_beta_pval_2<t4      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   adj_beta_pval_4<t4  &   adj_beta_pval_3<t4      ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   adj_beta_pval_4<t4       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   adj_beta_pval_3<t4       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   adj_beta_pval_2<t4       ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   adj_beta_pval_1<t4       ]) 

length(adj_beta_pval_1[ adj_beta_pval_1<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_5<t4)    ]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_5<t4  |   adj_beta_pval_1<t4)    ]) 
length(adj_beta_pval_3[ adj_beta_pval_3<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_5<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_1<t4)    ]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t4  &   (adj_beta_pval_5<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_1<t4)    ]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_1<t4)    ]) 

length(adj_beta_pval_1[ adj_beta_pval_1<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_5<t4)    ]) / length(adj_beta_pval_1[adj_beta_pval_1<t4]) 
length(adj_beta_pval_2[ adj_beta_pval_2<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_5<t4  |   adj_beta_pval_1<t4)    ]) / length(adj_beta_pval_2[adj_beta_pval_2<t4])                         
length(adj_beta_pval_3[ adj_beta_pval_3<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_5<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_1<t4)    ]) / length(adj_beta_pval_3[adj_beta_pval_3<t4]) 
length(adj_beta_pval_4[ adj_beta_pval_4<t4  &   (adj_beta_pval_5<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_1<t4)    ]) / length(adj_beta_pval_4[adj_beta_pval_4<t4]) 
length(adj_beta_pval_5[ adj_beta_pval_5<t4  &   (adj_beta_pval_4<t4  |   adj_beta_pval_3<t4  |   adj_beta_pval_2<t4  |   adj_beta_pval_1<t4)    ]) / length(adj_beta_pval_5[adj_beta_pval_5<t4]) 










###################################################################
slope_se_1 = abs(slope_se_1) 
slope_se_2 = abs(slope_se_2) 
slope_se_3 = abs(slope_se_3) 
slope_se_4 = abs(slope_se_4) 
slope_se_5 = abs(slope_se_5) 

summary(slope_se_1)

hist(slope_se_1, breaks = 100  )
hist(slope_se_1[adj_beta_pval_1<t1], breaks = 100 , xlim=c(0, 1.5) )
slope_se_1A = slope_se_1[adj_beta_pval_1<t1]
slope_se_1B = slope_se_1[adj_beta_pval_1<0.01]
slope_se_1C = slope_se_1[adj_beta_pval_1<0.001]

length(slope_se_1[slope_se_1>0.5]) / length(slope_se_1)
length(slope_se_1A[slope_se_1A>0.5]) / length(slope_se_1A)
length(slope_se_1B[slope_se_1B>0.5]) / length(slope_se_1B)
length(slope_se_1C[slope_se_1C>0.5]) / length(slope_se_1C)

plot(c(1,2,3,4,5,6), c(2,4,6,8,10,12), xlim=c(0, 12), ylim=c(0,12)  )
plot( c(2,4,6,8,10,12),  c(1,2,3,4,5,6) , xlim=c(0, 12), ylim=c(0,12) )


t11=0.1
length(slope_se_1[slope_se_1<t11]) 
length(slope_se_2[slope_se_2<t11]) 
length(slope_se_3[slope_se_3<t11]) 
length(slope_se_4[slope_se_4<t11]) 
length(slope_se_5[slope_se_5<t11]) 

length(slope_se_1[slope_se_1<t11 & sharedBOOL]) 
length(slope_se_2[slope_se_2<t11 & sharedBOOL]) 
length(slope_se_3[slope_se_3<t11 & sharedBOOL]) 
length(slope_se_4[slope_se_4<t11 & sharedBOOL]) 
length(slope_se_5[slope_se_5<t11 & sharedBOOL]) 

length(slope_se_5[ slope_se_5<t11  &   slope_se_4<t11  &   slope_se_3<t11  &   slope_se_2<t11  &   slope_se_1<t11    ]) 
length(slope_se_5[ slope_se_5<t11  &   slope_se_4<t11  &   slope_se_3<t11  &   slope_se_2<t11      ]) 
length(slope_se_5[ slope_se_5<t11  &   slope_se_4<t11  &   slope_se_3<t11      ]) 
length(slope_se_5[ slope_se_5<t11  &   slope_se_4<t11       ]) 
length(slope_se_5[ slope_se_5<t11  &   slope_se_3<t11       ]) 
length(slope_se_5[ slope_se_5<t11  &   slope_se_2<t11       ]) 
length(slope_se_5[ slope_se_5<t11  &   slope_se_1<t11       ]) 

length(slope_se_1[ slope_se_1<t11  &   (slope_se_4<t11  |   slope_se_3<t11  |   slope_se_2<t11  |   slope_se_5<t11)    ]) 
length(slope_se_2[ slope_se_2<t11  &   (slope_se_4<t11  |   slope_se_3<t11  |   slope_se_5<t11  |   slope_se_1<t11)    ]) 
length(slope_se_3[ slope_se_3<t11  &   (slope_se_4<t11  |   slope_se_5<t11  |   slope_se_2<t11  |   slope_se_1<t11)    ]) 
length(slope_se_4[ slope_se_4<t11  &   (slope_se_5<t11  |   slope_se_3<t11  |   slope_se_2<t11  |   slope_se_1<t11)    ]) 
length(slope_se_5[ slope_se_5<t11  &   (slope_se_4<t11  |   slope_se_3<t11  |   slope_se_2<t11  |   slope_se_1<t11)    ]) 

length(slope_se_1[ slope_se_1<t11  &   (slope_se_4<t11  |   slope_se_3<t11  |   slope_se_2<t11  |   slope_se_5<t11)    ]) / length(slope_se_1[slope_se_1<t11]) 
length(slope_se_2[ slope_se_2<t11  &   (slope_se_4<t11  |   slope_se_3<t11  |   slope_se_5<t11  |   slope_se_1<t11)    ]) / length(slope_se_2[slope_se_2<t11])                         
length(slope_se_3[ slope_se_3<t11  &   (slope_se_4<t11  |   slope_se_5<t11  |   slope_se_2<t11  |   slope_se_1<t11)    ]) / length(slope_se_3[slope_se_3<t11]) 
length(slope_se_4[ slope_se_4<t11  &   (slope_se_5<t11  |   slope_se_3<t11  |   slope_se_2<t11  |   slope_se_1<t11)    ]) / length(slope_se_4[slope_se_4<t11]) 
length(slope_se_5[ slope_se_5<t11  &   (slope_se_4<t11  |   slope_se_3<t11  |   slope_se_2<t11  |   slope_se_1<t11)    ]) / length(slope_se_5[slope_se_5<t11]) 




t1=0.2
t12=0.5
length(adj_beta_pval_5[ adj_beta_pval_5<t1  &   adj_beta_pval_4<t1  &   adj_beta_pval_3<t1  &   adj_beta_pval_2<t1  &   adj_beta_pval_1<t1  & slope_se_5>t12  &   slope_se_4>t12  &   slope_se_3>t12  &   slope_se_2>t12  &   slope_se_1>t12  ]) 

length(adj_beta_pval_1[ (adj_beta_pval_1<t1   &   slope_se_1>t12 ) &    ( (adj_beta_pval_4<t1 &   slope_se_4>t12 )  |   ( adj_beta_pval_3<t1 &   slope_se_3>t12 )  |   ( adj_beta_pval_2<t1 &   slope_se_2>t12 )  |   ( adj_beta_pval_5<t1 &   slope_se_5>t12) )      ]) / length(adj_beta_pval_1[adj_beta_pval_1<t1  &   slope_se_1>t12]) 
length(adj_beta_pval_2[ (adj_beta_pval_2<t1   &   slope_se_2>t12 ) &    ( (adj_beta_pval_4<t1 &   slope_se_4>t12 )  |   ( adj_beta_pval_3<t1 &   slope_se_3>t12 )  |   ( adj_beta_pval_5<t1 &   slope_se_5>t12 )  |   ( adj_beta_pval_1<t1 &   slope_se_1>t12) )      ]) / length(adj_beta_pval_2[adj_beta_pval_2<t1  &   slope_se_2>t12])                         
length(adj_beta_pval_3[ (adj_beta_pval_3<t1   &   slope_se_3>t12 ) &    ( (adj_beta_pval_4<t1 &   slope_se_4>t12 )  |   ( adj_beta_pval_5<t1 &   slope_se_5>t12 )  |   ( adj_beta_pval_2<t1 &   slope_se_2>t12 )  |   ( adj_beta_pval_1<t1 &   slope_se_1>t12) )      ]) / length(adj_beta_pval_3[adj_beta_pval_3<t1  &   slope_se_3>t12]) 
length(adj_beta_pval_4[ (adj_beta_pval_4<t1   &   slope_se_4>t12 ) &    ( (adj_beta_pval_5<t1 &   slope_se_5>t12 )  |   ( adj_beta_pval_3<t1 &   slope_se_3>t12 )  |   ( adj_beta_pval_2<t1 &   slope_se_2>t12 )  |   ( adj_beta_pval_1<t1 &   slope_se_1>t12) )      ]) / length(adj_beta_pval_4[adj_beta_pval_4<t1  &   slope_se_4>t12]) 
length(adj_beta_pval_5[ (adj_beta_pval_5<t1   &   slope_se_5>t12 ) &    ( (adj_beta_pval_4<t1 &   slope_se_4>t12 )  |   ( adj_beta_pval_3<t1 &   slope_se_3>t12 )  |   ( adj_beta_pval_2<t1 &   slope_se_2>t12 )  |   ( adj_beta_pval_1<t1 &   slope_se_1>t12) )      ]) / length(adj_beta_pval_5[adj_beta_pval_5<t1  &   slope_se_5>t12]) 









################################################
t30=0.1
t45=0.5

mybool_1 = adj_beta_pval_1<t30  &  slope_se_1<t45 
mybool_2 = adj_beta_pval_2<t30  &  slope_se_2<t45 
mybool_3 = adj_beta_pval_3<t30  &  slope_se_3<t45 
mybool_4 = adj_beta_pval_4<t30  &  slope_se_4<t45 
mybool_5 = adj_beta_pval_5<t30  &  slope_se_5<t45 

 
mybool_1 = adj_beta_pval_1<t30  &  r_squared_1> t45 
mybool_2 = adj_beta_pval_2<t30  &  r_squared_2> t45 
mybool_3 = adj_beta_pval_3<t30  &  r_squared_3> t45 
mybool_4 = adj_beta_pval_4<t30  &  r_squared_4> t45 
mybool_5 = adj_beta_pval_5<t30  &  r_squared_5> t45 






length(mybool_1[mybool_1])
length(mybool_2[mybool_2])
length(mybool_3[mybool_3])
length(mybool_4[mybool_4])
length(mybool_5[mybool_5])
length(mybool_5[mybool_5  &  mybool_4  &  mybool_3  &  mybool_2  &  mybool_1 ])

length(mybool_1[mybool_1 & sharedBOOL])
length(mybool_2[mybool_2 & sharedBOOL])
length(mybool_3[mybool_3 & sharedBOOL])
length(mybool_4[mybool_4 & sharedBOOL])
length(mybool_5[mybool_5 & sharedBOOL])



rawMatrix_1B = rawMatrix_1_noNA[ mybool_1 , ]
rawMatrix_2B = rawMatrix_2_noNA[ mybool_2 , ]
rawMatrix_3B = rawMatrix_3_noNA[ mybool_3 , ]
rawMatrix_4B = rawMatrix_4_noNA[ mybool_4 , ]
rawMatrix_5B = rawMatrix_5_noNA[ mybool_5 , ]

dim(rawMatrix_1B )
dim(rawMatrix_2B )
dim(rawMatrix_3B )
dim(rawMatrix_4B )
dim(rawMatrix_5B )




write.table(x=rawMatrix_1B, file = "finalPermutation.rawMatrix_1B.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")


write.table(x=rawMatrix_2B, file = "finalPermutation.rawMatrix_2B.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")


write.table(x=rawMatrix_3B, file = "finalPermutation.rawMatrix_3B.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")


write.table(x=rawMatrix_4B, file = "finalPermutation.rawMatrix_4B.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")


write.table(x=rawMatrix_5B, file = "finalPermutation.rawMatrix_5B.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")














