library(debCAM)


fantom_o<-fantom_c[,colnames(fantom_c)%in%Par2$Sample]
fantom_o<-fantom_o[,order(colnames(fantom_o))]
colnames(fantom_o)<-Par2$name
rownames(fantom_o)<-rownames(fantom_c)
fantom_CAM <- CAM(as.data.frame(fantom_o), K = 2:10, thres.low = 0.10, thres.high = 0.95)

tiff(filename ="MDL_reference.tiff",width = 6000,height = 3000,compression = "lzw",res = 400 )
plot(MDL(fantom_CAM), data.term = TRUE)
dev.off()
Aest <- Amat(fantom_CAM, 7)
row.names(Aest)<-colnames(as.data.frame(fantom_o))
