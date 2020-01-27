
library(edgeR)
library(readr)
library(debCAM)

fantom_o<-fantom_c[,colnames(fantom_c)%in%Par2$Sample]
fantom_o<-fantom_o[,order(colnames(fantom_o))]
colnames(fantom_o)<-Par2$name
rownames(fantom_o)<-rownames(fantom_c)
fantom_CAM <- CAM(as.data.frame(fantom_o), K = 2:10, thres.low = 0.10, thres.high = 0.95)

png(filename ="images/MDL_reference.png",width = 6000,height = 3000,res = 400 )
plot(MDL(fantom_CAM), data.term = TRUE)
dev.off()
Aest <- Amat(fantom_CAM, 7)
row.names(Aest)<-colnames(as.data.frame(fantom_o))
write.table(Aest,file = "debCAM_7component.txt",sep = "\t",quote = T)

library(DeconRNASeq)

Ref_fantom<-list()
#for (i in 1:length(unique(Par2$group))){
#  fantom_p<-cpm(y)[,colnames(cpm(y))%in%Par2$name[Par2$group %in% unique(Par2$group)[i]]]
#  Med_fantom[[i]]<-rowMedians(fantom_p)
#}
for (i in 1:length(unique(Par2$group))){
  fantom_p<-cpm(Reference_all)[,colnames(cpm(Reference_all))%in%Par2$name[Par2$group %in% unique(Par2$group)[i]]]
  Ref_fantom[[i]]<-rowMedians(fantom_p)
}
Ref_fantom<-do.call(cbind,Ref_fantom)
colnames(Ref_fantom)<-unique(Par2$group)

rownames(Ref_fantom)<-rownames(cpm(Reference_all))



DEC_fantom<-decon.bootstrap(as.data.frame(fantom_o), as.data.frame(Ref_fantom),round(nrow(Ref_fantom)*0.8), 1000)
write.table(DEC_fantom[,,1],file = "DeconRNASeq.txt",sep = "\t",quote = T)

