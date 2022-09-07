
pcafunc<-function(boxallnor){
  rownam<-boxallnor[,1]  #accession col
  boxallnor0<-boxallnor[,seq(2,dim(boxallnor)[2])]  #w/o accession col
  boxallnor0<-(apply(boxallnor0,2,as.numeric))  #change to numeric
  if(sum(is.na(boxallnor0))>0){  
    boxallnor0[is.na(boxallnor0)]<-0
    warning("The data contain missing values and are changed to zero")
  }
  
  pcaallnor0<-t(boxallnor0)
  colnames(pcaallnor0)<-rownam
  write.csv(pcaallnor0,file="pca matrix.csv",row.names = F)
  pca<- prcomp(pcaallnor0, scale=T)
  return(pca)
}
