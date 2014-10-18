setwd("/media/drozdovapb/big/Studies//bioinformaticsinstitute/bioinfo-lab-14/bi2014-coffea-canephora/annotation/")

#install.packages("hoardeR")
library(hoardeR)

byTomato <- importGFF3("./chr1_by_tomato_Augustus.txt")
byArabidopsis <- importGFF3("./chr1_by_Arabidopsis_augustus.gff")

byTomatoGenes <- subset(byTomato, byTomato$V3=="gene")
byArabidopsisGenes <- subset(byArabidopsis, byArabidopsis$V3 == "gene")

diffStart<-setdiff(byArabidopsisGenes$V4, byTomatoGenes$V4)
diffStop<-setdiff(byArabidopsisGenes$V5, byTomatoGenes$V5)

comStart <- intersect(byArabidopsisGenes$V4, byTomatoGenes$V4) 
comStop <- intersect(byArabidopsisGenes$V5, byTomatoGenes$V5) 

length(comStart)
length(comStop)

merge(byTomatoGenes, byArabidopsisGenes)

byTomatoExon <- subset(byTomato, byTomato$V3=="exon")
byArabidopsisExon <- subset(byArabidopsis, byArabidopsis$V3 == "gene")

intersect(byArabidopsisExon$V4, byTomatoExon$V4) 
intersect(byArabidopsisExon$V5, byTomatoExon$V5) 


inters<- importGFF3("./chr1_intersectBed_f_90.gff")

nrow(subset(inters, inters$V3=="gene"))

inters2 <- importGFF3("./chr1_intersectBed_r_f_90.gff")
nrow(subset(inters2, inters2$V3=="gene"))
