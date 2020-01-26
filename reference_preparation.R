library(readr)
library(edgeR)
fantom_c <- read.table("human.srna.cpm.tsv")
fantom_c<-fantom_c[-grep(rownames(fantom_c),pattern = "novelmiR"),]
human_srna_cellontology <- read_delim("human.srna.cellontology.tsv", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
Par<-c(44,53,100,143,55,94,142,39,17,111,59,56)
Par2<-list()
for(i in 1:length(Par)){
  
  Par2[[i]]<-cbind(strsplit(x = as.character(human_srna_cellontology[Par[i],2]),split = ",")[[1]],human_srna_cellontology[Par[i],1])
}       
Par2<-do.call(rbind,Par2)
colnames(Par2)[1]<-"Sample"
Par2<-Par2[order(as.character(Par2$Sample)),]
EXP<-list()
for (i in 1:length(Par)){
  EXP2<-list()
  for (j in 1:length(Par)){
    print(i)
    print(j)
    
    fantom_o<-fantom_c[,strsplit(x = as.character(human_srna_cellontology[Par[i],2]),split = ",")[[1]]]
    fantom_o2<-fantom_c[strsplit(x = as.character(human_srna_cellontology[Par[j],2]),split = ",")[[1]]]
    rownames(fantom_o)<-rownames(fantom_c)
    rownames(fantom_o2)<-rownames(fantom_c)
    group <- factor(c(rep(0,ncol(fantom_o)),rep(1,ncol(fantom_o2))))
    y <- DGEList(counts=cbind(fantom_o,fantom_o2),group=group)
    keep <- rowSums(cpm(y)>10) >= 1
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    design <- model.matrix(~group)
    y <- estimateDisp(y,design)
    
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,coef=2)
    #topTags(qlf)
    
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    #topTags(lrt)
    FDR <- p.adjust(qlf$table$PValue, method="BH")
    
    FDR2 <- p.adjust(lrt$table$PValue, method="BH")
    Dif_table<-cbind(lrt$table,FDR2)
    write.table(Dif_table,file = paste0("cell_reference/",human_srna_cellontology[Par[i],1]," vs ",human_srna_cellontology[Par[j],1],".txt"),sep = "\t",quote = T)
    EXP2[[j]]<-list(qlf,lrt,FDR,FDR2,sum(FDR<0.05),sum(FDR2<0.05),sum(qlf$table$PValue<0.05),sum(lrt$table$PValue<0.05),rownames(lrt$table)[FDR2<0.05])
    
    
  }
  EXP[[i]]<-EXP2
}

FF2<-do.call(rbind,lapply(EXP, function(x){unlist(lapply(x, function(y){y[[6]]}))}))
colnames(FF2)<-human_srna_cellontology$ontologyID[Par]
rownames(FF2)<-human_srna_cellontology$ontologyID[Par]
write.table(FF2,file = "cell_reference/Tables_Differential_Cell_Expression.txt",sep = "\t",quote = T)


