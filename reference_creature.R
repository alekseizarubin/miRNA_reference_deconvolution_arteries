Par<-c(44,53,100,143,94,55,142,39,17,111,59,56)
Par2<-list()
for(i in 1:length(Par)){
  
  Par2[[i]]<-cbind(strsplit(x = as.character(human_srna_cellontology[Par[i],2]),split = ",")[[1]],human_srna_cellontology[Par[i],1])
}       
Par2<-do.call(rbind,Par2)
colnames(Par2)[1]<-"Sample"
Par2$group<-c(rep("SMC_Fib",11),rep("SMC_Fib",3),rep("EC",9),rep("Macr",3),rep("B",4),rep("T",6),rep("NK",3),rep("Neut",3))
Par2<-Par2[!(Par2$Sample %in% c("SRhi10011.TAGCTT.11437","SRhi10010.GTAGAG.11284")),]

EXP<-list()
for (i in 1:length(unique(Par2$group))){
  EXP2<-list()
  for (j in 1:length(unique(Par2$group))){
    print(i)
    print(j)
    
    fantom_o<-fantom_c[,colnames(fantom_c)%in%Par2$Sample[Par2$group %in% unique(Par2$group)[i]]]
    fantom_o2<-fantom_c[,colnames(fantom_c)%in%Par2$Sample[Par2$group %in% unique(Par2$group)[j]]]
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
    write.table(Dif_table,file = paste0("cell_group_reference/",unique(Par2$group)[i]," vs ",unique(Par2$group)[j],".txt"),sep = "\t",quote = T)
    EXP2[[j]]<-list(qlf,lrt,FDR,FDR2,sum(FDR<0.05),sum(FDR2<0.05),sum(qlf$table$PValue<0.05),sum(lrt$table$PValue<0.05),rownames(lrt$table)[FDR2<0.05])
    
    
  }
  EXP[[i]]<-EXP2
}

FF2<-do.call(rbind,lapply(EXP, function(x){unlist(lapply(x, function(y){y[[6]]}))}))
colnames(FF2)<-unique(Par2$group)
rownames(FF2)<-unique(Par2$group)
#DIFF_mir<-unique(unlist(lapply(EXP, function(x){unlist(lapply(x, function(y){y[[9]]}))})))
write.table(FF2,file = "cell_group_reference/Tables_Differential_Cell_Expression.txt",sep = "\t",quote = T)

Par2<-Par2[order(as.character(Par2$Sample)),]
fantom_o<-fantom_c[,colnames(fantom_c)%in%Par2$Sample]
fantom_o<-fantom_o[,order(colnames(fantom_o))]
Par2$name<-paste0(Par2$group,"_",1:nrow(Par2))
colnames(fantom_o)<-Par2$name
rownames(fantom_o)<-rownames(fantom_c)
group <- factor(Par2$group)
y <- DGEList(counts=fantom_o,group=group)
keep <- rowSums(cpm(y)>10) >= length(group)/5
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~group)
y <- estimateDisp(y,design)


fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR<0.05)
FDR2 <- p.adjust(lrt$table$PValue, method="BH")
sum(FDR<0.05)
plotMDS(y)
plotBCV(y)
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
plotMD(qlf)
abline(h=c(-1, 1), col="blue")
heatmap(log2(cpm(y)+1))


DIFF_mir<-row.names(lrt$table[FDR2<0.05,])







Par2<-Par2[order(as.character(Par2$Sample)),]
fantom_o<-fantom_c[DIFF_mir,colnames(fantom_c)%in%Par2$Sample]
fantom_o<-fantom_o[,order(colnames(fantom_o))]
Par2$name<-paste0(Par2$group,"_",1:nrow(Par2))
colnames(fantom_o)<-Par2$name
rownames(fantom_o)<-rownames(fantom_c)
group <- factor(Par2$group)
y <- DGEList(counts=fantom_o,group=group)
keep <- rowSums(cpm(y)>10) >= length(group)/5
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~group)
y <- estimateDisp(y,design)


fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR<0.05)
FDR2 <- p.adjust(lrt$table$PValue, method="BH")
sum(FDR<0.05)
plotMDS(y)
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
heatmap(log2(cpm(y)+1))
