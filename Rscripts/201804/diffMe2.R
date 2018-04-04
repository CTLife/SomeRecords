matrix1 = read.table(file="9_DMR_overdispersion_Ftest/2E_AlldiffMesites_q0.05_diff0.txt", header = TRUE, sep = "\t", quote = "\"'",
                     dec = ".")

dim(matrix1)
# matrix1[1:10,]

myQvalue = matrix1$qvalue
myDiff = matrix1$meth.diff   




print("#####################################################")
nrow( matrix1[( (myQvalue<0.001) & (abs(myDiff)>10) ),]  )
nrow( matrix1[( (myQvalue<0.01) & (abs(myDiff)>10) ),]  )
nrow( matrix1[( (myQvalue<0.05) & (abs(myDiff)>10) ),]  )
nrow( matrix1[( (myQvalue<0.05) & (abs(myDiff)>5) ),]  )
nrow( matrix1[( (myQvalue<0.001) & (abs(myDiff)>5) ),]  )
print("#####################################################")



diff1_all = matrix1[( (myQvalue<0.001) & (abs(myDiff)>20) ),]
diff2_all = matrix1[( (myQvalue<0.01)  & (abs(myDiff)>20) ),]
diff3_all = matrix1[( (myQvalue<0.05)  & (abs(myDiff)>20) ),]
print("#####################################################")
dim(diff1_all) 
dim(diff2_all) 
dim(diff3_all) 
print("#####################################################")



diff1_name = c()
for(i in c(1:nrow(diff1_all)) ) {
  diff1_name[i] = paste("diff1_all_", i, sep=""  )     
}


diff2_name = c()
for(i in c(1:nrow(diff2_all)) ) {
  diff2_name[i] = paste("diff2_all_", i, sep=""  )     
}


diff3_name = c()
for(i in c(1:nrow(diff3_all)) ) {
  diff3_name[i] = paste("diff3_all_", i, sep=""  )     
}




diff1_all = cbind(diff1_all[,1:4],  diff1_name, rep(1, nrow(diff1_all)) , diff1_all[,5:7] )   
diff2_all = cbind(diff2_all[,1:4],  diff2_name, rep(1, nrow(diff2_all)) , diff2_all[,5:7] )    
diff3_all = cbind(diff3_all[,1:4],  diff3_name, rep(1, nrow(diff3_all)) , diff3_all[,5:7] )    

colnames(diff1_all) = c("chr", "start",   "end",  "strand", "name",  "strength",  "pvalue",    "qvalue",  "meth.diff")
colnames(diff2_all) = c("chr", "start",   "end",  "strand", "name",  "strength",  "pvalue",    "qvalue",  "meth.diff")
colnames(diff3_all) = c("chr", "start",   "end",  "strand", "name",  "strength",  "pvalue",    "qvalue",  "meth.diff")



diff1_all[,3] = diff1_all[,3] + 1
diff2_all[,3] = diff2_all[,3] + 1
diff3_all[,3] = diff3_all[,3] + 1



diff1_hyper = diff1_all[ diff1_all$meth.diff > 0, ]
diff1_hypo  = diff1_all[ diff1_all$meth.diff < 0, ]
diff2_hyper = diff2_all[ diff2_all$meth.diff > 0, ]
diff2_hypo  = diff2_all[ diff2_all$meth.diff < 0, ]
diff3_hyper = diff3_all[ diff3_all$meth.diff > 0, ]
diff3_hypo  = diff3_all[ diff3_all$meth.diff < 0, ]

print("#####################################################")
print("diff1:")
dim(diff1_all) 
dim(diff1_hyper) 
dim(diff1_hypo) 
print("diff2:")
dim(diff2_all) 
dim(diff2_hyper) 
dim(diff2_hypo) 
print("diff3:")
dim(diff3_all) 
dim(diff3_hyper) 
dim(diff3_hypo) 
print("#####################################################")


Part1_g = "9_DMR_3conditions/1-rawFiles"
if( ! file.exists(Part1_g) )  { dir.create(Part1_g, recursive = TRUE) }

write.table(x=diff1_all,   file = paste(Part1_g, "diff1_all.txt", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff1_hyper, file = paste(Part1_g, "diff1_hyper.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff1_hypo,  file = paste(Part1_g, "diff1_hypo.txt", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
  
write.table(x=diff2_all,   file = paste(Part1_g, "diff2_all.txt", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff2_hyper, file = paste(Part1_g, "diff2_hyper.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff2_hypo,  file = paste(Part1_g, "diff2_hypo.txt", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )

write.table(x=diff3_all,   file = paste(Part1_g, "diff3_all.txt", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff3_hyper, file = paste(Part1_g, "diff3_hyper.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
write.table(x=diff3_hypo,  file = paste(Part1_g, "diff3_hypo.txt", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )











