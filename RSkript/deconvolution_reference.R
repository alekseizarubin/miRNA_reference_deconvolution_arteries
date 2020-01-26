library(debCAM)

fantom_CAM <- CAM(as.data.frame(fantom_o), K = 2:10, thres.low = 0.25, thres.high = 0.95)
plot(MDL(fantom_CAM), data.term = TRUE)
Aest <- Amat(fantom_CAM, 6)
row.names(Aest)<-colnames(as.data.frame(fantom_o))
