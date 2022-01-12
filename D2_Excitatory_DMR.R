#!/usr/bin/Rscript
library(DSS)
library(bsseq)

         dat_D2_W1 = read.table("wgbs_D2_Rep1.DSS.GT3", header=F)
         colnames(dat_D2_W1) <- c("chr",     "pos",     "N",       "X")

         dat_D2_W2 = read.table("wgbs_D2_Rep2.DSS.GT3", header=F)
         colnames(dat_D2_W2) <- c("chr",     "pos",     "N",       "X")

         dat_D2_E1 = read.table("EM_D2_Rep1.DSS.GT3", header=F)
         colnames(dat_D2_E1) <- c("chr",     "pos",     "N",       "X")

         dat_D2_E2 = read.table("EM_D2_Rep2.DSS.GT3", header=F)
        colnames(dat_D2_E2) <- c("chr",     "pos",     "N",       "X")
        
        dat_E1 = read.table("E1.CG.DSS.GT3", header=F)
        colnames(dat_E1) <- c("chr",     "pos",     "N",       "X")


BSobj = makeBSseqData( list(dat_D2_W1,dat_D2_W2,dat_D2_E1,dat_D2_E2,dat_E1,dat_E2),c("P1","P2","P3","P4","N1","N2"))
    dmlTest = DMLtest(BSobj, group2=c("P1","P2","P3","P4"), group1=c("N1","N2"),smoothing=TRUE)

dmls = callDML(dmlTest, delta=0.2)
write.table(dmls,"DML_D2_vs_Excit.20per_defaultP.txt",sep="\t",row.name=FALSE,quote=F)
dmrs = callDMR(dmlTest,delta=0.1,p.threshold=1e-5,minlen=100,minCG=5,dis.merge=500)
write.table(dmrs,"DMR_D2_vs_Excit_10per_1e5_withsmooth.txt",sep="\t",quote=F,row.name=FALSE)
system("cat DMR_D2_vs_Excit_10per_1e5_withsmooth.txt | sed '1d'| sort -k1,1 -k2,2n | awk 'BEGIN{OFS=\"\t\"}{print $1,$2,$3,\"D2_vs_Excit\"NR,\".\",\"+\"}' > DMR_D2_vs_Excit_10per_1e5_withsmooth.bed")
