limmafunc<-function(paircomp,colnam1,path0,thresh,bioinfo){
  ord <- c("MHC/W",unique(as.vector(colnam1))[!unique(as.vector(colnam1)) %in% c("MHC/W")])  #WT and the other gene treatment
  colnam2<-colnam1[order(match(colnam1,ord))]
  rownam<-paircomp[,1]  #accession
  paircomp0<-paircomp[,seq(2,dim(paircomp)[2])]  #w/o accession
  paircomp0<-(apply(paircomp0,2,as.numeric))  #change to numeric
  if(sum(is.na(paircomp0))>0){
    paircomp0[is.na(paircomp0)]<-0
    warning("The data contain missing values and are changed to zero")
  }
  colnam2<-relevel(colnam2, ref = "MHC/W")
  design<-model.matrix(~colnam2)
  rownames(paircomp0)<-rownam
  # Ordinary fit
  fit <- lmFit(paircomp0,design)  #linear model in limma
  fit2 <- eBayes(fit)  #get p-value
  table1<-topTable(fit2,coef=2,n=Inf,confint = TRUE)  #table in limma
  table1[,"Accession"]<-rownames(table1)  #add back accession
  setwd(file.path(root0,path0))
  table0<-merge(bioinfo,table1,by="Accession")  #add back bioinfo
  write.csv(table0,str_replace(paste0(paste(as.character(unique(colnam1)),collapse ="_"),".csv"),"/","_"),row.names=F)
  # dim(fit)
  # colnames(fit)
  # rownames(fit)[1:10]
  # names(fit)
  # # Fold-change thresholding
  # fit2 <- treat(fit,lfc=0.1)
  # topTreat(fit2,coef=2)
  # Volcano plot
  table1$updown<-'NO'    #create a col to show up/down regulation
  table1$updown[table1$logFC>thresh & table1$adj.P.Val<0.05]<-'Up'  #set up regulation
  table1$updown[table1$logFC<(-thresh) & table1$adj.P.Val<0.05] <-'Down'  #set down regulation
  
  table0$updown<-'NO'
  table0$updown[table0$logFC>thresh & table0$adj.P.Val<0.05]<-'Up'
  table0$updown[table0$logFC<(-thresh) & table0$adj.P.Val<0.05] <-'Down'
  table0$'Gene Symbol'[table0$updown =='NO']<-NA
  
  mycolors <- c("red", "blue","grey")
  names(mycolors) <- c("Up","Down","NO")
  # table1$lab<-NA
  # table1$lab[table1$updown !='NO']<-rownames(table1)[table1$updown !='NO']
  # colnames(table1)[colnames(table1)=="lab"]<-'Accession'
  table1$Accession[table1$updown =='NO']<-NA
  tablesig<-table1[table1$updown !='NO',]
  tablesig<-as.data.frame(tablesig)
  table2<-tablesig
  table3<-table2[(table2[,"adj.P.Val"]<0.05 & abs(table2[,"logFC"])>thresh),]
  
  table4<-merge(bioinfo,table3,by="Accession",all.y=TRUE)
  write.csv(table4,file=str_replace(paste0(paste(as.character(unique(colnam1)),collapse ="_")," signif.csv"),"/","_"),row.names=F)
  
  tiff(str_replace(paste0(paste(as.character(unique(colnam1)),collapse =" vs "),".tiff"),"/","_"),units="in", width=11, height=6,res=200)
  print(ggplot(data=table0,aes(x=logFC,y=-log10(adj.P.Val),col=updown,label=table0$'Gene Symbol'))+ geom_vline(xintercept=c(-thresh, thresh), col="purple") +
          ylim(0,6)+
          theme_minimal()+             #volcano
          scale_colour_manual(name="Gene Expression",values=mycolors)+
          geom_point()+
          geom_hline(yintercept=-log10(0.05), col="purple")+
          # geom_text(col='black',size=5,position=position_jitter(width=1,height=1))+
          geom_label_repel(label.size=NA,show.legend = FALSE,
                           box.padding   = 0.3, 
                           point.padding = 0.02,)+
          labs(title=paste(as.character(unique(colnam1)),collapse =" vs "),
               x ="Log2(FC)", y = "-Log10(adj.P.value)")+
          theme(plot.title = element_text(hjust = 0.5),
                axis.title = element_text(size=10),
                axis.text=element_text(size=10),
                legend.text=element_text(size=10))+
          coord_cartesian(clip='off')+
          scale_shape_manual(values = c(17, 17)))
          # geom_count(show.legend = FALSE))
  # ggsave(str_replace(paste0(paste(as.character(unique(colnam1)),collapse =" vs "),".tiff"),"/","_"),p)
  dev.off()
  
}
