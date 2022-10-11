library(readr)
root0<-'/Users/rushi/Desktop/FYP/R_code'
if(!dir.exists(root0)){
  dir.create(root0)
}

batch1_renamed <- read_csv("/Users/rushi/Desktop/FYP/Yi's data/batch1_renamed.csv")
batch2_renamed <- read_csv("/Users/rushi/Desktop/FYP/Yi's data/batch2_renamed.csv")
Accession<-batch1_renamed$Accession
batch1<-batch1_renamed[,11:22]  #w/o descriptive info
batch2<-batch2_renamed[,11:23]
bioinfo<-batch1_renamed[,1:10]

b1e<-cbind(Accession,batch1)  #accession + value
b2e<-cbind(Accession,batch2)

group0<-c(rep('batch1',12),rep('batch2',13))
batchname<-c(colnames(b1e)[-1],colnames(b2e)[-1]) #batch names of batch1&2
dict<-cbind(group0,batchname)

###########Get batch number
library(stringr)
batchnumber<-as.numeric(as.vector(sapply(batchname,function(x)str_extract(x,"(?i)(?<=batch)\\d+"))))

boxall<-cbind('Accession'=Accession,sapply(batch1,as.numeric), #accession+batch1+batch2
              sapply(batch2,as.numeric))
boxall<-as.data.frame(boxall)  #convert to df
boxallnumb<-sapply(boxall[,-1],function(x)log2(as.numeric(x)))  #remove accession col and apply log2 on each value
boxallnumb0 = as.data.frame(boxallnumb)
boxallnumb0[is.na(boxallnumb0)] = 0 #Convert missing to 0

######Normalization
boxallnor<-apply(boxallnumb,2,function(x)(x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE))  #normalise each col, margin = 2 means col, na.rm removes missing value
######Row-wise normalization
boxallnorf<-t(apply(boxallnumb,1,function(x)(x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE))) #normalise each row, margin = 1, transpose

######Normalization NA to 0
boxallnor0<-apply(boxallnumb0,2,function(x)(x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE))
boxallnor0 = data.frame("Accession" = boxall[,1], boxallnor0)

######Find missing >=0.5 and missing combination
mispercent<-apply(boxallnor,1, function(x,y=0.5){  #for each row, if less than half is NA, label as 1, or else as 0
  temp<-sum(is.na(x))/length(x)
  a<-c()
  if (temp<y){
    a<-1
  }else{
    a<-0}
  return(a)
})

mispercent2<-apply(boxallnor,1, function(x,y=0.4){  #for each row, if less than half is NA, label as 1, or else as 0
  temp<-sum(is.na(x))/length(x)
  a<-c()
  if (temp<y){
    a<-1
  }else{
    a<-0}
  return(a)
})

combicheck<-c(1,1,1,2,2,2,3,3,3,4,5,6,7,7,7,4,4,4,5,5,5,6,6,6,8)
#seq(1,9)
 #seq(1,8)
mispercentcombi<-apply(boxallnor,1, function(x,y=0.5){
  a<-c()
  if (sum(unique(combicheck[!is.na(x)]))!=sum(seq(1,8))){
    a<-0
  }else{a<-1}
  return(a)
})

combicheck1<-c(1,1,1,2,2,2,3,3,3,9,9,9,7,7,7,4,4,4,5,5,5,6,6,6,8)
mispercentcombi1<-apply(boxallnor,1, function(x,y=0.5){
  a<-c()
  if (sum(unique(combicheck1[!is.na(x)]))!=sum(seq(1,9))){  #at least got a value for each type
    a<-0
  }else{a<-1}
  return(a)
})

combicheck2<-c(1,1,1,2,2,2,3,3,3,8,8,8,7,7,7,4,4,4,5,5,5,6,6,6,0)
mispercentcombi2<-apply(boxallnor,1, function(x,y=0.5){
  a<-c()
  if (sum(unique(combicheck2[!is.na(x)]))!=sum(seq(1,8))){  #at least got a value for each type
    a<-0
  }else{a<-1}
  return(a)
})

######Add accession number
boxallc2<-data.frame('Accession'=boxall[,1],boxallnor)  #add back accession col, boxallc2 = Accession + boxallnor
missclean<-boxallc2[(mispercent & mispercentcombi),]  #filter out some rows
missclean2 = boxallc2[(mispercent & mispercentcombi1),]  #50% NA & 1-9
missclean3 = boxallc2[(mispercent & mispercentcombi2),]  #50% NA 1-8
path1<-"Missing clean Zscore 0.5 missing"
setwd(root0)
if(!dir.exists(path1)){
  dir.create(path1)
}
write.csv(missclean,file=file.path(path1,'Clean data Zscore 0.5 missing.csv'),row.names = FALSE)
diff = subset(missclean2, !(Accession %in% missclean$Accession))
diff$gene = comp_gene_info2[, 'Gene Symbol'][match(diff[,'Accession'], comp_gene_info2[,'Accession'])]
write.csv(diff, file = file.path(root0,"difference.csv"),row.names = FALSE)
#check NA % in diff
naa = c()
for (i in 1:nrow(diff)) {
  nas = length(which(is.na(diff[,-1][i,])))
  perc = nas*100/ncol(diff[,-1]) 
  naa = append(naa, perc)
}
diff$NApercent = naa


######Combat correction design  ##ruv package
library(sva)
combatdesign<-as.factor(c(rep('MHC/W',3),rep('USP10-PA OE',3),rep('USP10 sh',3),
                          'FMR1 RNAi-2','RIN OE','RIN RNAi',rep('FMR1 OE',3),rep('FMR1 RNAi-2',3),
                          rep('RIN OE',3),rep('RIN RNAi',3),'MHC/W'))
t1<-table(combatdesign)
twobatch<-names(t1)[t1==4]
design1<-model.matrix(~combatdesign)
misscleannumb<-missclean[,-1]   #w/o accession
batchcorrect<-ComBat(misscleannumb,batch = group0,mod = design1)
final<-cbind(missclean[,1],batchcorrect)  #add back accession col
path0<-'Batch_removed'
setwd(root0)
if (!file.exists(path0)){
  dir.create(path0)
}
write.csv(final,file=file.path(path0,'Z score 0.5 missing.csv'),row.names = F)


######PCA check
originalcolname<-colnames(final)
#colnames(final)<-originalcolname
pcalabel<-str_replace_all(originalcolname,"_batch\\d*","")
pcalabell<-paste0("    ",pcalabel)
colnames(final)<-pcalabell
groups = c(rep('MHC/W',3),rep('USP10-PA OE',3),rep('USP10 sh',3),
           'FMR1 RNAi-2','RIN OE','RIN RNAi',rep('FMR1 OE',3),rep('FMR1 RNAi-2',3),
           rep('RIN OE',3),rep('RIN RNAi',3),'MHC/W')
# groups<-c(rep('Batch1',12),rep('Batch2',13))
# batchname<-c(colnames(b1e)[-1],colnames(b2e)[-1])
# groups[batchname %in% intersectset]<-"CommenBatch"
# groups[c(2,10,11,12)]<-"CommonBatch1"
# groups[c(18,19,22,25)]<-"CommonBatch2"

source("/Users/rushi/Desktop/FYP/Yi's data/pcafunc.R")
pcacombat<-pcafunc(final)  #use prcomp to do PCA
library(factoextra)
if (!dir.exists(file.path(root0,'PCA_combat'))){
  dir.create(file.path(root0,'PCA_combat'))
}
setwd(file.path(root0,'PCA_combat'))
tiff("Combat Z score normalization 0.5 missing.tiff",units="in", width=11, height=10,res=250)  #safe as high resolution picture
options(ggrepel.max.overlaps=Inf)
fviz_pca_ind(pcacombat,         #plot PCA
             col.ind = groups, # color by groups
             palette = c("#00AFBB","#FC4E07", 'orange','blue','#336600','#FF33FF','purple'), #"grey50","grey50",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             labelsize=5,
             ggtheme = theme(axis.text = element_text(size = 13),
                             title=element_text(size=16),
                             plot.title=element_text(hjust=0.5),
                             legend.text = element_text(size=15)))+
  labs(title ="PCA Scores", x = "PC1", y = "PC2")
dev.off()  #stop plotting PCA

#####50% NA,1-9 combat, PCA
batchcorrect<-ComBat(missclean2[,-1],batch = group0,mod = design1)
final1_9<-cbind(missclean2[,1],batchcorrect)

originalcolname = colnames(final1_9)
pcalabel<-str_replace_all(originalcolname,"_batch\\d*","")
pcalabell<-paste0("    ",pcalabel)
colnames(final1_9)<-pcalabel
groups = c(rep('MHC/W',3),rep('USP10-PA OE',3),rep('USP10 sh',3),
           'FMR1 RNAi-2','RIN OE','RIN RNAi',rep('FMR1 OE',3),rep('FMR1 RNAi-2',3),
           rep('RIN OE',3),rep('RIN RNAi',3),'MHC/W')
source("/Users/rushi/Desktop/FYP/Yi's data/pcafunc.R")
pcacombat1_9<-pcafunc(final1_9)  #use prcomp to do PCA
library(factoextra)
if (!dir.exists(file.path(root0,'PCA_combat'))){
  dir.create(file.path(root0,'PCA_combat'))
}
setwd(file.path(root0,'PCA_combat'))
tiff("Combat PCA 1-9.tiff",units="in", width=11, height=10,res=250)  #safe as high resolution picture
options(ggrepel.max.overlaps=Inf)
fviz_pca_ind(pcacombat1_9,         #plot PCA
             col.ind = groups, # color by groups
             palette = c("#00AFBB","#FC4E07", 'orange','blue','#336600','#FF33FF','purple'), #"grey50","grey50",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             labelsize=5,
             ggtheme = theme(axis.text = element_text(size = 13),
                             title=element_text(size=16),
                             plot.title=element_text(hjust=0.5),
                             legend.text = element_text(size=15)))+
  labs(title ="PCA Analysis", x = "PC1", y = "PC2")
dev.off()

##### 50% NA 1-8 Combat, PCA
batchcorrect2<-ComBat(missclean3[,-1],batch = group0,mod = design1)
final1_8<-cbind(missclean3[,1],batchcorrect2)
originalcolname = colnames(final1_8)
pcalabel<-str_replace_all(originalcolname,"_batch\\d*","")
pcalabell<-paste0("    ",pcalabel)
colnames(final1_8)<-pcalabel
groups = c(rep('MHC/W',3),rep('USP10-PA OE',3),rep('USP10 sh',3),
           'FMR1 RNAi-2','RIN OE','RIN RNAi',rep('FMR1 OE',3),rep('FMR1 RNAi-2',3),
           rep('RIN OE',3),rep('RIN RNAi',3),'MHC/W')
source("/Users/rushi/Desktop/FYP/Yi's data/pcafunc.R")
pcacombat1_8<-pcafunc(final1_8)  #use prcomp to do PCA
library(factoextra)
if (!dir.exists(file.path(root0,'PCA_combat'))){
  dir.create(file.path(root0,'PCA_combat'))
}
setwd(file.path(root0,'PCA_combat'))
tiff("Combat PCA 1-8.tiff",units="in", width=11, height=10,res=250)  #safe as high resolution picture
options(ggrepel.max.overlaps=Inf)
fviz_pca_ind(pcacombat1_8,         #plot PCA
             col.ind = groups, # color by groups
             palette = c("#00AFBB","#FC4E07", 'orange','blue','#336600','#FF33FF','purple'), #"grey50","grey50",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             labelsize=5,
             ggtheme = theme(axis.text = element_text(size = 13),
                             title=element_text(size=16),
                             plot.title=element_text(hjust=0.5),
                             legend.text = element_text(size=15)))+
  labs(title ="PCA Analysis", x = "PC1", y = "PC2")
dev.off()


furtherint<-final
intre<-c('FMR1_RNAi','USP10_sh','RIN_RNAi','MHC_W')
furint<-furtherint[,grepl(paste(intre,collapse="|"),c('Accession',batchname))]
intgrp<-groups[grepl(paste(intre,collapse="|"),batchname)]
furint<-cbind(furtherint[,1],furint)
furintpca<-pcafunc(furint)
if (!dir.exists(file.path(root0,"BatchControl"))){
  dir.create(file.path(root0,"BatchControl"))
}
setwd(file.path(root0,"BatchControl"))
tiff("Interested Z score normalization 0.5 missing.tiff",units="in", width=11, height=6,res=200)
options(ggrepel.max.overlaps=Inf)
fviz_pca_ind(furintpca,
             col.ind = intgrp, # color by groups
             palette = c("grey","grey","#00AFBB","#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
dev.off()

path1<-"missing_values"
if (!dir.exists(file.path(root0,path1))){
  dir.create(file.path(root0,path1))
}
checkmis<-final[!complete.cases(final),]  #compile all rows w missing values
write.csv(checkmis,file=file.path(root0,path1,"Missing after clean and z score trans.csv"),row.names = FALSE)

#####Chromosome position & CG number
id_list = read.csv('/Users/rushi/Desktop/FYP/id_validation_table_115027.csv')
chr_data = read.csv('/Users/rushi/Desktop/FYP/FlyBase_Fields_download.csv')
colnames(id_list) = c('#submitted_item', 'FBID_KEY',	'current_symbol')
comp_gene_info = merge(id_list, chr_data, by='FBID_KEY')
colnames(comp_gene_info)[2] = "Accession"; colnames(comp_gene_info)[4] = "Submitted ID"
comp_gene_info2 = merge(comp_gene_info, bioinfo, by='Accession')
write.csv(comp_gene_info2,file=file.path(root0,"Complete_gene_info_CG_Chr.csv"),row.names = FALSE)

######Remove Controls
# y<-str_match(batchname, ".*[_](\\d{2}).*")
# z<-as.numeric(y[,2])
# pos<-which(z==4)
pos<-c(11,12,13,26)
path1<-"BatchRemovedWithoutControl"
if (!dir.exists(file.path(root0,path1))){
  dir.create(file.path(root0,path1))
}
clean0<-final[,-pos]  #updated final
colnames(clean0)[1]<-'Accession'
clean1 = as.data.frame(clean0)
clean1[,'CG_number'] = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(clean0[,'Accession'], comp_gene_info2[,'Accession'])]  #add CG number
clean1[, 'Description'] = comp_gene_info2[, 'Gene Symbol'][match(clean0[,'Accession'], comp_gene_info2[,'Accession'])]  #add gene symbol
write.csv(clean1,file=file.path(root0,path1,"BatchRemoveWithoutControl.csv"),row.names = FALSE)

#Remove controls: 1-9
clean1_9<-final1_9[,-pos]  #updated final
colnames(clean1_9)[1]<-'Accession'
clean1_9.2 = as.data.frame(clean1_9)
clean1_9.2[,'CG_number'] = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(clean1_9[,'Accession'], comp_gene_info2[,'Accession'])]  #add CG number
clean1_9.2[, 'Description'] = comp_gene_info2[, 'Gene Symbol'][match(clean1_9[,'Accession'], comp_gene_info2[,'Accession'])]  #add gene symbol
write.csv(clean1_9.2,file=file.path(root0,path1,"BatchRemoveWithoutControl 1-9.csv"),row.names = FALSE)

#Remove controls: 1-8
clean1_8<-final1_8[,-pos]  #updated final
colnames(clean1_8)[1]<-'Accession'
clean1_8.2 = as.data.frame(clean1_8)
clean1_8.2[,'CG_number'] = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(clean1_8[,'Accession'], comp_gene_info2[,'Accession'])]  #add CG number
clean1_8.2[, 'Description'] = comp_gene_info2[, 'Gene Symbol'][match(clean1_8[,'Accession'], comp_gene_info2[,'Accession'])]  #add gene symbol
write.csv(clean1_8.2,file=file.path(root0,path1,"BatchRemoveWithoutControl 1-8.csv"),row.names = FALSE)

#####PCA after removing controls 1-9
originalcolname2<-colnames(clean1_9)
pcalabel2<-str_replace_all(colnames(clean1_9),"_batch\\d*","")
#pcalabel2<-paste0("    ",pcalabel2)
colnames(clean1_9)<-pcalabel2
groups2<-c(rep('MHC/W',3),rep('USP10-PA OE',3),rep('USP10 sh',3),
           rep('FMR1 OE',3),rep('FMR1 RNAi-2',3),
           rep('RIN OE',3),rep('RIN RNAi',3))

library(factoextra)
pcacombat2<-pcafunc(clean1_9)  #use prcomp to do PCA
if (!dir.exists(file.path(root0,'PCA_combat'))){
  dir.create(file.path(root0,'PCA_combat'))
}
setwd(file.path(root0,'PCA_combat'))
tiff("Combat PCA after removing control 1-9",units="in", width=11, height=10,res=250)  #safe as high resolution picture
options(ggrepel.max.overlaps=Inf) #???
fviz_pca_ind(pcacombat2,         #plot PCA
             col.ind = groups2, # color by groups
             palette = c("#00AFBB","#FC4E07", 'orange','blue','#336600','#FF33FF','purple'),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE,
             labelsize=5,
             ggtheme = theme(axis.text = element_text(size = 13),
                             title=element_text(size=16),
                             plot.title=element_text(hjust=0.5),
                             legend.text = element_text(size=15)))+
  labs(title ="PCA Analysis after removing controls 1-9", x = "PC1", y = "PC2")
dev.off()  #stop plotting PCA

#####PCA after removing controls 1-8
originalcolname2<-colnames(clean1_8)
pcalabel2<-str_replace_all(colnames(clean1_8),"_batch\\d*","")
#pcalabel2<-paste0("    ",pcalabel2)
colnames(clean1_8)<-pcalabel2
groups2<-c(rep('Control',3),rep('USP10 OE',3),rep('USP10 sh',3),
           rep('FMR1 OE',3),rep('FMR1 RNAi',3),
           rep('RIN OE',3),rep('RIN RNAi',3))
groups3 = factor(
  groups2,
  levels = c(rep('Control',3),rep('USP10 OE',3),rep('USP10 sh',3),
             rep('FMR1 OE',3),rep('FMR1 RNAi',3),
             rep('RIN OE',3),rep('RIN RNAi',3)),
  labels = c(rep('Control',3),rep('USP10 OE',3),rep('USP10 sh',3),
             rep('FMR1 OE',3),rep('FMR1 RNAi',3),
             rep('RIN OE',3),rep('RIN RNAi',3))
)
library(factoextra)
pcacombat3<-pcafunc(clean1_8)  #use prcomp to do PCA
if (!dir.exists(file.path(root0,'PCA_combat'))){
  dir.create(file.path(root0,'PCA_combat'))
}
setwd(file.path(root0,'PCA_combat'))
tiff("Combat PCA after removing control 1-8.tiff",units="in", width=11, height=10,res=250)  #safe as high resolution picture
options(ggrepel.max.overlaps=Inf) #???
fviz_pca_ind(pcacombat3,         #plot PCA
             col.ind = groups3, # color by groups
             palette = c("#00AFBB","#FC4E07", 'orange','blue','#336600','#FF33FF','purple'),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             pointsize = 5,
             label = 'None',
             invisible="quali",
             ggtheme = theme(axis.text = element_text(size = 13),
                             title=element_text(size=16),
                             plot.title=element_text(hjust=0.5),
                             legend.text = element_text(size=15)))+
  labs(title ="PCA Scores", x = "PC1", y = "PC2") +
  scale_shape_manual(values=rep(19,7))
dev.off()  #stop plotting PCA


######Volcano Plot
library(limma)
library(ggrepel)
source("/Users/rushi/Desktop/FYP/Yi's data/limmafunc2.R")
combatdesign0<-combatdesign[-(pos-1)]  #remove10,11,12,25 controls
treatname<-as.vector(unique(combatdesign0)[-1])  #w/o wildtype MHC
##Volcano with adj.p.val
path0<-"Volcano plot"
if(!file.exists(file.path(root0,path0))){
  dir.create(file.path(root0,path0))}
for (i in seq(1,length(treatname))){   #parse treatname
  pos<-combatdesign0 %in% c("MHC/W",treatname[i])
  colnam<-combatdesign0[combatdesign0 %in% c("MHC/W",treatname[i])]  #colnam = combatdesign0p[pos]
  colnam1<-droplevels(colnam)
  paircomp<-cbind(clean1_8[,1],clean1_8[,2:dim(clean1_8)[2]][,pos])  #cbind accession and treatname[i], eg. USP10
  thresh<-0.322
  limmafunc(paircomp, colnam1,path0,thresh,bioinfo)
}

##Dual volcano plot with adj.p.val
# source("/Users/rushi/Desktop/FYP/Yi's data/limmafunc_dual.R")
path6 = "Dual volcano plot"
if(!file.exists(file.path(root0,"Dual volcano plot"))){
  dir.create(file.path(root0,"Dual volcano plot"))}
setwd("/Users/rushi/Desktop/FYP/R_code/Dual volcano plot")
volcano.files <- list.files(pattern="^.*\\.csv")
for (i in 1:length(volcano.files)) {
  file = read.csv(paste0("/Users/rushi/Desktop/FYP/R_code/Dual volcano plot/", volcano.files[i]))
  file$updown<-'No significant change'
  file$updown[file$logFC>0.322 & file$adj.P.Val<0.05]<-'Upregulation (adj. p < 0.05)'
  file$updown[file$logFC<(-0.322) & file$adj.P.Val<0.05] <-'Downregulation (adj. p < 0.05)'
  #file$Gene.Symbol[file$updown =='NO']<-NA
  file$labels = ''
  file$green = ''
  if (i==1) {
    file$labels = ifelse(file$Gene.Symbol == 'Fmr1', 'purple', FALSE)
    file$green = ifelse(file$Gene.Symbol %in% c('CG11876','SdhC','l(1)G0156','SdhB','SdhA','slgA','ND4L','CG9090','CG14209; Mpc1','CG14757','CG9399','ND-B8'), 'Mitochondrial genes (adj. p < 0.05)',FALSE)
  } else if (i==2) {
    file$labels = ifelse(file$Gene.Symbol == 'Fmr1', 'purple', FALSE)
    file$updown[which(file$Gene.Symbol == 'Fmr1')] = 'Downregulation (adj. p < 0.05)'
  } else if (i==3) {
    file$labels = ifelse(file$Gene.Symbol == 'rin', 'purple',FALSE)
    file$green = ifelse(file$Gene.Symbol %in% c('ATPsyn-gamma; ATPsyngamma','CG12079; NDUFS3; ND-30','ATPsyn-b; ATPsynB','SdhB','ND4','CG12400; ND-B14.5B','CG5703; ND-24','l(3)neo18; ND-SGDH','mtacp1; ND-ACP',
                                                'ATP6','CG32230; ND-MLRQ','COX2','CG10320; ND-B12','l(1)G0230; ATPsyndelta','CG3214; ND-B17.2','CG4769; Cyt-c1','l(2)06225; ATPsynG','CG8680; ND-13A','CG6020; ND-39',
                                                'CG3192; ND-ASHI','CG4692; ATPsynF','CG4169; UQCR-C2','blw','ND5', 'CG1746; ATPsynC','CG34439; ND-MWFE','CG12203; ND-18','CG9090','Oscp; ATPsynO','Cyt-c-p',
                                                'CG7712; ND-B14','CG3560; UQCR-14','CG3621; ND-B14.5A','CG14290; Mpc1', 'CG40002; ND-AGGG; CG40472', 'CG3321; ATPsynE','CG9399','Mlc1','NDL4','CG15434; ND-B8'),'Mitochondrial genes (adj. p < 0.05)',FALSE)
  } else if (i==4){
    file$labels = ifelse(file$Gene.Symbol == 'rin', 'purple',FALSE)
  } else if (i==6) {
    file$labels = ifelse(file$Gene.Symbol %in% c('Usp10','Atg9','ref(2)P','Clbn','Fmr1'),'purple',FALSE)
    file$green = ifelse(file$Gene.Symbol %in% c('CG4769; Cyt-c1','mfrn','CG12859; ND-B15','CG7580; UQCR-Q','CYTB','CG3731; UQCR-C1','CG3560; UQCR-14','CG12400; ND-B14.5B',
                                                'CG9306; ND-B22','l(2)35Di; ND-B17','CG4169; UQCR-C2','CG11015; CoVb; COX5B','CG3321; ATPsynE','CG3192; ND-ASHI',
                                                'l(3)neo18; ND-SGDH','CG3214; ND-B17.2','ATPsyn-d; ATPsynD','ND4','SdhC','CG9350; ND-B14.7','CG11876','CG10664; CoIV; COX4',
                                                'COX1','ATPsyn-gamma; ATPsyngamma','CG3446; ND-B16.6','CG14235; CoVIb; COX6B','l(2)06225; ATPsynG','Sirt4','Nc73EF',
                                                'CG3621; ND-B14.5A','l(1)G0230; ATPsyndelta','CG6463; ND-13B','CG32230; ND-MLRQ','ATPsyn-b; ATPsynB','CG11779','ND5',
                                                'ND42; ND-42','ATP6','CG11110','CG9140; ND-51','CG4692; ATPsynF','CG1970; ND-49','CG1746; ATPsynC','ATPsyn-beta; ATPsynbeta','CG8680; ND-13A',
                                                'blw','COX2','ND23; ND-23','ND75; ND-75','CG5703; ND-24','CoVa; COX5A','SdhA','CG34439; ND-MWFE','Oscp; ATPsynO','CG6020; ND-39',
                                                'CG12079; NDUFS3; ND-30','CG12203; ND-18','CG7506','CG7712; ND-B14','Sdh8','ATP8','mtacp1; ND-ACP','l(1)G0136; MagR',
                                                'CG14290; Mpc1','Acon','ND4L','CG18317; Rim2','CG9399','Cyt-c-p','CG9090'),'Mitochondrial genes (adj. p < 0.05)',FALSE)
  } else if (i==5){
    file$labels = ifelse(file$Gene.Symbol == 'Usp10','purple',FALSE)
    file$updown[which(file$Gene.Symbol == 'Usp10')] = 'Downregulation (adj. p < 0.05)'
  }
  write.csv(file, file = file.path(root0, path6, volcano.files[i]))
}

for (i in 1:length(volcano.files)) {
  for (j in 1:length(volcano.files[-i])) {
    volcano1 = read.csv(paste0("/Users/rushi/Desktop/FYP/R_code/Dual volcano plot/", volcano.files[i]))
    volcano1[, "log10(adj.P.Val)"] = -(log10(volcano1$adj.P.Val))
    volcano2 = read.csv(paste0("/Users/rushi/Desktop/FYP/R_code/Dual volcano plot/", volcano.files[-i][j]))
    volcano2[, "log10(adj.P.Val)"] = log10(volcano2$adj.P.Val)
    volcano = rbind(volcano1, volcano2)
    name1 = unlist(strsplit(volcano.files[i], split='.', fixed=TRUE))[1]
    name2 = unlist(strsplit(volcano.files[-i][j], split='.', fixed=TRUE))[1]
    name = paste(name1, 'vs', name2)
    mycolors <- c("red", "blue","grey","green")
    names(mycolors) <- c("Upregulation (adj. p < 0.05)","Downregulation (adj. p < 0.05)","No significant change","Mitochondrial genes (adj. p < 0.05)")
    green.volcano = dplyr::filter(volcano, (volcano$green == 'Mitochondrial genes (adj. p < 0.05)' & volcano$updown != 'No significant change'))
    purple.volcano = dplyr::filter(volcano, (volcano$labels == 'purple'))
    tiff(paste0(name, '.tiff'), units="in", width=11, height=6,res=200)
    print(ggplot(data=volcano,aes(x=logFC,y=`log10(adj.P.Val)`,col=updown)) +
            xlim(-2.5,2.5)+
            theme_minimal()+             
            scale_colour_manual(name="Gene Expression",values=mycolors)+
            geom_point(size = 1)+
            geom_point(data = green.volcano, aes(x=logFC,y=`log10(adj.P.Val)`),color='green',size=1.5)+
            geom_point(data = purple.volcano, aes(x=logFC,y=`log10(adj.P.Val)`),color='purple',size=1.5)+
            geom_hline(yintercept=-log10(0.05), col="purple", size = 0.2)+
            geom_hline(yintercept=log10(0.05), col="purple", size = 0.2)+
            geom_vline(xintercept=c(-0.322,0.322), col="purple", size = 0.2) +
            geom_vline(xintercept = 0, col = "black") +
            geom_hline(yintercept = 0, col = "black") +
            geom_text_repel(label = ifelse(volcano$labels == 'purple', volcano$Gene.Symbol, ""),box.padding = unit(0.45, "lines"),hjust = 1,col = 'purple')+
            labs(title=name, x ="Log2(Fold Change)", y = paste(name2,'                              ', name1, '\n', "-Log10(adjusted P value)"))+
            theme(plot.title = element_text(hjust = 0.5),
                  axis.title = element_text(size=10),
                  axis.text=element_text(size=10),
                  legend.text=element_text(size=10))+
            theme(axis.line.x = element_line(color="black", size = 1),
                  axis.line.y = element_line(color="black", size = 1)) +
            geom_text(aes(x=-2,y=6), label = '')+geom_text(aes(x=-2,y=-6), label = '')+
            #scale_shape_manual(values = c(17, 17)) +
            # scale_y_continuous(labels = abs))
            scale_y_continuous(limits = c(-5,5),breaks = c(5,4,3,2,1,0,-1,-2,-3,-4,-5),
                               labels = c('5','4','3','2','1','0','1','2','3','4','5')))
    dev.off()
  }
}

##Volcano with raw P.Value
source("/Users/rushi/Desktop/FYP/Yi's data/limmafunc_p.R")
path0<-"Volcano plot with p-value"
if(!file.exists(file.path(root0,path0))){
  dir.create(file.path(root0,path0))}
for (i in seq(1,length(treatname))){   #parse treatname
  pos<-combatdesign0 %in% c("MHC/W",treatname[i])
  colnam<-combatdesign0[combatdesign0 %in% c("MHC/W",treatname[i])]  #colnam = combatdesign0p[pos]
  colnam1<-droplevels(colnam)
  paircomp<-cbind(clean0[,1],clean0[,2:dim(clean0)[2]][,pos])  #cbind accession and treatname[i], eg. USP10
  thresh<-0.6
  limmafunc_p(paircomp, colnam1,path0,thresh,bioinfo)
}


# #Remove mitochondrial genome
# comp_gene_info2 = as.data.frame(comp_gene_info2);
# gene_no_mito = comp_gene_info2[which(comp_gene_info2$LOCATION_ARM != "mitochondrion_genome"),]

#####Manhattan plot 
##Manhattan plot with P Value
library(ggplot2)
library(dplyr)
path2 = "Manhattan plot with p-value"
if(!file.exists(file.path(root0,path2))){
  dir.create(file.path(root0,path2))}
setwd("/Users/rushi/Desktop/FYP/R_code/Volcano plot with p-value")
wd<-getwd()
files0 <- list.files(pattern="^.*\\.csv")
setwd("/Users/rushi/Desktop/FYP/R_code/Manhattan plot with p-value")
mycolors <- c("red", "blue","grey")
names(mycolors) <- c("Up","Down","NO")
for (i in 1:length(files0)) {
  manhattan_file = read.csv(paste0("/Users/rushi/Desktop/FYP/R_code/Volcano plot with p-value/",files0[i]))
  manhattan_file[,"Chromosome"] = comp_gene_info2$LOCATION_ARM[match(manhattan_file$Accession, comp_gene_info2$Accession)]  #add chromosome position
  manhattan_file[,"Gene_symbol"] = comp_gene_info2$`Gene Symbol`[match(manhattan_file$Accession, comp_gene_info2$Accession)]
  manhattan_file$updown =  "NO"
  manhattan_file$updown[manhattan_file$logFC>0.6 & manhattan_file$P.Value<0.05] = "Up"  #up or down regulation
  manhattan_file$updown[manhattan_file$logFC< -0.6 & manhattan_file$P.Value<0.05] = "Down"
  manhattan_file = filter(manhattan_file, Chromosome != "-" & !is.na(Chromosome))  ## remove NA and -
  l = unlist(strsplit(files0[i], split='_', fixed=TRUE))
  manhattan_name = paste('MHC_W vs', unlist(strsplit(l[3], split = '.',fixed=TRUE))[1])
  title0 = paste('MHC/W vs', unlist(strsplit(l[3], split = '.',fixed=TRUE))[1])
  # manhattan_name = unlist(strsplit(files0[i], split='.', fixed=TRUE))[1]
  write.csv(manhattan_file, file = paste0(manhattan_name, '.csv'), row.names = F)
  tiff(paste0(manhattan_name,'.tiff'),units="in", width=11, height=6,res=200)
  print(ggplot(manhattan_file, aes(x = Chromosome, y = -log10(P.Value), col = updown)) +
    geom_jitter() + 
    #geom_point() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(name="Gene Expression",values=mycolors)+
    geom_text(aes(label=ifelse(updown=='NO','',manhattan_file$Gene_symbol))) +
    ylim(0.0, 10.0) +
    labs(title = title0)); print('success1')
  dev.off(); print('success2')
}

##Manhattan with adj.P.Value
path2 = "Manhattan plot with adj.p.value"
if(!file.exists(file.path(root0,path2))){
  dir.create(file.path(root0,path2))}
setwd("/Users/rushi/Desktop/FYP/R_code/Volcano plot")
wd<-getwd()
files1 <- list.files(pattern="^.*\\.csv")
setwd("/Users/rushi/Desktop/FYP/R_code/Manhattan plot with adj.p.value")
for (i in 1:length(files1)) {
  manhattan_file = read.csv(paste0("/Users/rushi/Desktop/FYP/R_code/Volcano plot/",files1[i]))
  manhattan_file[,"Chromosome"] = gene_no_mito$LOCATION_ARM[match(manhattan_file$Accession, gene_no_mito$Accession)]  #add chromosome position
  manhattan_file[,"Gene_symbol"] = gene_no_mito$`Gene Symbol`[match(manhattan_file$Accession, gene_no_mito$Accession)]
  manhattan_file$updown =  "NO"
  manhattan_file$updown[manhattan_file$logFC>0.6 & manhattan_file$adj.P.Val<0.05] = "Up"  #up or down regulation
  manhattan_file$updown[manhattan_file$logFC< -0.6 & manhattan_file$adj.P.Val<0.05] = "Down"
  manhattan_file = filter(manhattan_file, Chromosome != "-" & !is.na(Chromosome))  ## remove NA and -
  l = unlist(strsplit(files1[i], split='_', fixed=TRUE))
  manhattan_name = paste('MHC_W vs', unlist(strsplit(l[3], split = '.',fixed=TRUE))[1])
  title0 = paste('MHC/W vs', unlist(strsplit(l[3], split = '.',fixed=TRUE))[1])
  # manhattan_name = unlist(strsplit(files0[i], split='.', fixed=TRUE))[1]
  write.csv(manhattan_file, file = paste0(manhattan_name, '.csv'), row.names = F)
  tiff(paste0(manhattan_name,'.tiff'),units="in", width=11, height=6,res=200)
  print(ggplot(manhattan_file, aes(x = Chromosome, y = -log10(adj.P.Val), col = updown)) +
          geom_jitter() + 
          #geom_point() +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_colour_manual(name="Gene Expression",values=mycolors)+
          geom_text(aes(label=ifelse(updown=='NO','',manhattan_file$Gene_symbol))) +
          ylim(0.0, 10.0) +
          labs(title = title0)); print('success1')
  dev.off(); print('success2')
}


#####Preparation for GSEA Rushi
##Convert BatchRemoveWithoutControl
BatchRemoveWithoutControl <- read.csv("/Users/rushi/Desktop/FYP/R_code/BatchRemovedWithoutControl/BatchRemoveWithoutControl.csv")
treatment = colnames(BatchRemoveWithoutControl)[-c(1:4, 23,24)]
setwd("/Users/rushi/Desktop/FYP/GSEA")
for (i in 1:6) {
  gseafile = cbind(BatchRemoveWithoutControl[, c('CG.number', 'Description')],BatchRemoveWithoutControl[,1:4],BatchRemoveWithoutControl[, treatment[(3*i-2):(3*i)]])
  gseafile = gseafile[,-3]
  colnames(gseafile)[c(1,2)] = c("NAME", "DESCRIPTION")
  write.table(gseafile, paste0('MHC_W vs ', unlist(strsplit(treatment[3*i], split='_0', fixed=TRUE))[1], " GSEA.txt"), append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, quote = FALSE)
}

##Convert BatchRemoveWithoutControl 1-9
path4 = "GSEA 1-9"
if(!file.exists(file.path("/Users/rushi/Desktop/FYP/GSEA",path4))){
  dir.create(file.path("/Users/rushi/Desktop/FYP/GSEA",path4))}
setwd("/Users/rushi/Desktop/FYP/GSEA/GSEA 1-9")
for (i in 1:6) {
  gseafile = cbind(clean1_9.2[, c('CG_number', 'Description')],clean1_9.2[,1:4],clean1_9.2[, treatment[(3*i-2):(3*i)]])
  gseafile = gseafile[,-3]
  colnames(gseafile)[c(1,2)] = c("NAME", "DESCRIPTION")
  write.table(gseafile, paste0('MHC_W vs ', unlist(strsplit(treatment[3*i], split='_0', fixed=TRUE))[1], " GSEA 1-9.txt"), append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, quote = FALSE)
}

##Convert BatchRemoveWithoutControl 1-8
path5 = "GSEA 1-8"
if(!file.exists(file.path("/Users/rushi/Desktop/FYP/GSEA",path5))){
  dir.create(file.path("/Users/rushi/Desktop/FYP/GSEA",path5))}
setwd("/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8")
for (i in 1:6) {
  gseafile = cbind(clean1_8.2[, c('CG_number', 'Description')],clean1_8.2[,1:4],clean1_8.2[, treatment[(3*i-2):(3*i)]])
  gseafile = gseafile[,-3]
  colnames(gseafile)[c(1,2)] = c("NAME", "DESCRIPTION")
  write.table(gseafile, paste0('MHC_W vs ', unlist(strsplit(treatment[3*i], split='_0', fixed=TRUE))[1], " GSEA 1-8.txt"), append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, quote = FALSE)
}

##Check boxallnor against new Reactome mtor signalling pathway
setwd("/Users/rushi/Desktop/FYP/GSEA")
library(GSA)
boxallc3 = boxallc2
boxallc3[,"CG_number"] = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(boxallc3[,'Accession'], comp_gene_info2[,'Accession'])]
reactome <- GSA.read.gmt("~/Desktop/FYP/GSEA/new_Drosophila_Reactome.gmt")
mtor = unlist(reactome$genesets[which(reactome$geneset.names == "MTOR signalling")])[1:30]
mtor_gene_before_remove_na = mtor[sapply(mtor, function(x) x %in% boxallc3$CG_number)]
genesymbol1 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_before_remove_na, comp_gene_info2$ANNOTATION_SYMBOL)]
expression.info = boxallc3[, c(2:26)][match(mtor_gene_before_remove_na, boxallc3$CG_number),]
df_mtor_in_boxallc2 = cbind(data.frame(mtor_gene = mtor_gene_before_remove_na, gene_symbol = genesymbol1), expression.info)
mtor_gene_not_in_boxallc2 = mtor[sapply(mtor, function(x) !(x %in% boxallc3$CG_number))]
write.csv(df_mtor_in_boxallc2, file = 'mtor genes in boxallc2 before NA removal.csv')
# genesymbol2 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_not_in_boxallc2, comp_gene_info2$ANNOTATION_SYMBOL)]
# df_mtor_not_in_boxallc2 = data.frame(mtor_gene = mtor_gene_not_in_boxallc2, gene_symbol = genesymbol2)

##Check final against new KEGG mtor signalling pathway
colnames(final)[1] = 'Accession'
final2 = as.data.frame(missclean2)
final2$CG_number = NA
final2$CG_number = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(missclean2[,'Accession'], comp_gene_info2[,'Accession'])]
mtor_gene_after_remove_na = mtor[sapply(mtor, function(x) x %in% final2$CG_number)]
genesymbol2 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_after_remove_na, comp_gene_info2$ANNOTATION_SYMBOL)]
expression.info2 = boxallc3[, c(2:26)][match(mtor_gene_after_remove_na, boxallc3$CG_number),]
df_mtor_after_remove_na = cbind(data.frame(mtor_gene = mtor_gene_after_remove_na, gene_symbol = genesymbol2), expression.info2)

mtor_gene_removed = mtor_gene_before_remove_na[!(mtor_gene_before_remove_na %in% mtor_gene_after_remove_na)]
genesymbol3 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_removed, comp_gene_info2$ANNOTATION_SYMBOL)]
expression.info3 = boxallc3[, c(2:26)][match(mtor_gene_removed, boxallc3$CG_number),]
df_mtor_removed = cbind(data.frame(mtor_gene = mtor_gene_removed, gene_symbol = genesymbol3), expression.info3)

write.csv(df_mtor_after_remove_na, file = 'mtor genes in final after NA removal.csv')
write.csv(df_mtor_removed, file = 'mtor genes removed.csv')






#####Preparation for GSEA Analysis

setwd("/Users/lily2019e/Dropbox/Wu/Genome_New/Volcano plot/AC_FBgn")
wd<-getwd()
files <- list.files(pattern="^.*\\.csv")
fileorig<-files[!as.vector(sapply(files,function(x)grepl("AC_FB",x,fixed=TRUE)))]
fileorig1<-fileorig[as.vector(sapply(fileorig,function(x)grepl("MHC_W_",x,fixed=TRUE)))]
filedict<-files[as.vector(sapply(files,function(x)grepl("AC_FB",x,fixed=TRUE)))]

for (i in fileorig1){
  datao<-read_csv(i)
  genesymbol<-datao$`Gene Symbol`
  compos<-sapply(datao$`Gene Symbol`,function(x)grepl(";",x))
  genesymbol[compos]<-sapply(genesymbol[compos],function(x)str_extract(x,"[^;]+"))
  datao[,"Genesym"]<-genesymbol
  data1<-data.frame("Genesym"=genesymbol,"Accession"=datao$Accession)
  write.csv(data1,file=file.path(wd,paste0(i," AC_FB.csv")),row.names = FALSE)
}

######Unique
uni<-function(x){
  n_occur <- data.frame(table(x))
  y<-n_occur[n_occur$Freq > 1,1]
  z<-n_occur[n_occur$Freq == 1,1]
  re<-list('uniq1'=z,'rep1'=y)
  return(re)
}
########Unique Full
uniful<-function(x){
  n_occur <- data.frame(table(x))
  y<-n_occur[n_occur$Freq > 1,]
  z<-n_occur[n_occur$Freq == 1,]
  re<-list("uniqc"=z,"repc"=y)
  return(re)
}


######Based on Gene Symbols

# n_occur <- data.frame(table(MHC_W_FMR1_OE_AC_FB$From))
# n_occur[n_occur$Freq > 1,]

standa<-try[,1:4]
head(standa)
standa<-standa[complete.cases(standa$Accession), ]
acces<-try[,5:6]
pos<-which(is.na(standa$FB))
suppFB<-standa$Accession[pos]
length(pos)
length(unique(suppFB))
FBmis<-acces$AC_FB[acces$Accession_F %in% suppFB]
n_occur <- data.frame(table(FBmis))
chec<-n_occur[n_occur$Freq>1,]
chec
makupFB<-n_occur[n_occur$Freq==1,]
repFB<-acces$Accession_F[acces$AC_FB %in% as.vector(chec[,1])]
nuni<-try$OriginGenesym[try$Accession %in% repFB]
makup<-acces[acces$AC_FB %in% makupFB[,1],]

for (j in seq(1,dim(makup)[1])){
  standa$FB[standa$Accession==unlist(makup[j,1])]<-unlist(makup[j,2])
}

standa[is.na(standa$FB),]


#######Based on Accession
standa<-try[,1:4]
head(standa)
standa<-standa[complete.cases(standa$Accession), ]
acces<-try[,5:6]
pos<-uni(acces$Accession_F)
uniqAC_FB<-acces$AC_FB[acces$Accession_F %in% pos$uniq1]
length(uniqAC_FB)==length(unique(uniqAC_FB))

mat<-standa[!(standa$Accession %in% pos$uniq1),]

to1AC<-as.vector(uni(acces$Accession_F)$uniq1)
to2AC<-as.vector(uni(acces$Accession_F)$rep1)
accesuniq<-acces[acces$Accession_F %in% to1AC,]
dim(accesuniq)

##START HERE##
library(readr)
library(dplyr)
BatchRemoveWithoutControl <- read.csv("/Users/rushi/Desktop/FYP/R_code/BatchRemovedWithoutControl/BatchRemoveWithoutControl.csv")
AC_nCG <- read.csv("/Users/rushi/Desktop/FYP/Yi's data/AC_nCG.csv", header = FALSE)
colnames(AC_nCG)<-c("Accession","FB","GeneGrp","Name") #change col names
supAC_CG<-AC_nCG[,c("Accession","Name")]  #only want col accession and name
supAC_CG<-supAC_CG[complete.cases(supAC_CG),]  #remove NA rows
supAC_CG2<-supAC_CG[!duplicated(supAC_CG$Accession), ]  #remove duplicated rows
# supAC_CG %>% group_by(Accession) %>% filter(row_number(Name) == 1)
name1<-colnames(BatchRemoveWithoutControl)
grpinfo<-str_replace_all(name1,"_\\d{2}_batch\\d{2}","")  #remove 'batch'
grpinfo<-grpinfo[!(grpinfo %in% "Accession")]  #remove 'Accession' from grpinfo
unigrp<-unique(grpinfo)
norep<-"sele1"
# norep<-"all"
for (i in seq(1,length(unigrp))){
  if (unigrp[i]!="MHC_W"){
    pos<-sapply(name1,function(x)grepl(paste0("MHC_W|",unigrp[i]),x))
    data0<-as.data.frame(BatchRemoveWithoutControl[,pos])
    data1<-data.frame("Accession"=BatchRemoveWithoutControl$Accession,"Description"=NA,data0)
    data1[match(convertdict$Accession,data1$Accession),"Name"]<-convertdict$CG  # what is converdict? not defined
    
    if (norep=="sele1"){
      AC_CG2<-merge(data1,supAC_CG2,by="Accession",all.y=TRUE)
    }else if(norep=="all"){
      AC_CG2<-merge(data1,supAC_CG,by="Accession",all.y=TRUE)
    }
    
    colnames(AC_CG2)<-c(colnames(AC_CG2)[seq(1,dim(AC_CG2)[2]-2)],"Remove","Name")
    AC_CG2<-AC_CG2[,!(colnames(AC_CG2) %in% "Remove")]
    AC_CG2<-AC_CG2[,c("Name",names(data1)[names(data1) != "Name" & names(data1) !="Accession"])]
    
    AC_CG1<-data1[!is.na(data1$Name),]
    AC_CG1<-AC_CG1[,c("Name",names(data1)[names(data1) != "Name" & names(data1) !="Accession"])]
    
    data2<-rbind(AC_CG1,AC_CG2)
    data3<-aggregate(data2[,2:8], by=list(Name=data2$Name), FUN=sum)
    dirname1<-paste0(norep," MHC_W vs ",unigrp[i])
    if (!dir.exists(file.path("/Users/lily2019e/Dropbox/Wu/Genome_New/GSEA",dirname1))){
      dir.create(file.path("/Users/lily2019e/Dropbox/Wu/Genome_New/GSEA",dirname1))
    }
    path1<-file.path("/Users/lily2019e/Dropbox/Wu/Genome_New/GSEA",dirname1)
    # write.csv(AC_CG2,file=file.path(path1,"check AC_CG2.csv"),row.names = FALSE)
    write.csv(data3,file=file.path(path1,paste0(norep," MHC_W vs ",unigrp[i]," GSEA.csv")),row.names=FALSE)
    write.table(data3, file=file.path(path1,paste0(norep," MHC_W vs ",unigrp[i]," GSEA.txt")), sep = "\t",
                row.names = FALSE) #, col.names = NA)
    
    data4<-data3
    data4[is.na(data4)]<-0
    sum(is.na(data4))
    write.table(data4, file=file.path(path1,paste0(norep," MHC_W vs ",unigrp[i]," GSEA nona.txt")), sep = "\t",
                row.names = FALSE)
    
  }
  
}


no_matching_CG_ying <- read_csv("no_matching_CG_ying.csv")

df1<-AC_FB_CG
df2<-no_matching_CG_ying
df1[match(df2$Accession_F, df1$Accession_F),"CG"]<-df2$CG
write.csv(df1,file="AC_CG1.csv",row.names = FALSE)
sum(is.na(df1$CG))
dim(uniful(df1$CG)$repc)
sum(is.na(df1$CG))
dim(uniful(AC_FB_CG$CG[!is.na(AC_FB_CG$CG)])$repc)
sum(is.na(AC_FB_CG$CG))
convertdict<-data.frame("Accession"=df1$Accession_F,"FB"=df1$FB,"CG"=df1$CG)
write.csv(convertdict,file="convert dictionary.csv",row.names = FALSE)


########Rename GSEA gene library name########

filename<-"Drosophila_Reactome.csv"
# filename<-"Drosophila_KEGG.txt"
gmtf<-read.csv(filename, header = FALSE)
# gmtf<-read.table(filename, sep = '\t',header = FALSE)
View(gmtf)
pos<-2
gmtfnew<-gmtf
gmtfnew[pos:dim(gmtf)[1],"V2"]<-gmtfnew[pos:dim(gmtf)[1],"V1"]
View(gmtfnew)
#write.csv(gmtfnew,file=file.path( "/Users/lily2019e/Dropbox/Wu/Genome_New/GSEA",paste0("new_",sub("\\.txt", "",filename),".csv")),row.names = FALSE)
write.table(gmtfnew,               # Write CSV file without header
            file.path( "/Users/lily2019e/Dropbox/Wu/Genome_New/GSEA",paste0("new_",sub("\\.csv", "",filename),".csv")),
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)

#########Pathway enrichment plot##########
library(forcats)
setwd('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final')
gseafile = list.files(pattern="^.*\\.csv")
fmr1oe.file = read.csv(paste0('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final/', gseafile[1]))
fmr1oe = data.frame(pathways=fmr1oe.file$NAME[which(fmr1oe.file$plot==TRUE)], NES=fmr1oe.file$NES[which(fmr1oe.file$plot==TRUE)], 
                    FDR.q.value=fmr1oe.file$FDR.q.val[which(fmr1oe.file$plot==TRUE)], Size=fmr1oe.file$SIZE[which(fmr1oe.file$plot==TRUE)])
fmr1oe$pathways[7] = "ACTIVATION OF THE MRNA UPON BINDING OF THE CAP-BINDING \n COMPLEX AND EIFS, AND SUBSEQUENT BINDING TO 43S"
fmr1rnai.file = read.csv(paste0('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final/', gseafile[2]))
fmr1rnai = data.frame(pathways=fmr1rnai.file$NAME[which(fmr1rnai.file$plot==TRUE)], NES=fmr1rnai.file$NES[which(fmr1rnai.file$plot==TRUE)], 
                      FDR.q.value=fmr1rnai.file$FDR.q.val[which(fmr1rnai.file$plot==TRUE)], Size=fmr1rnai.file$SIZE[which(fmr1rnai.file$plot==TRUE)])
rinoe.file = read.csv(paste0('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final/', gseafile[3]))
rinoe = data.frame(pathways=rinoe.file$NAME[which(rinoe.file$plot==TRUE)], NES=rinoe.file$NES[which(rinoe.file$plot==TRUE)], 
                   FDR.q.value=rinoe.file$FDR.q.val[which(rinoe.file$plot==TRUE)], Size=rinoe.file$SIZE[which(rinoe.file$plot==TRUE)])
rinrnai.file = read.csv(paste0('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final/', gseafile[4]))
rinrnai = data.frame(pathways=rinrnai.file$NAME[which(rinrnai.file$plot==TRUE)], NES=rinrnai.file$NES[which(rinrnai.file$plot==TRUE)], 
                     FDR.q.value=rinrnai.file$FDR.q.val[which(rinrnai.file$plot==TRUE)], Size=rinrnai.file$SIZE[which(rinrnai.file$plot==TRUE)])
usp10sh.file = read.csv(paste0('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final/', gseafile[6]))
usp10sh = data.frame(pathways=usp10sh.file$NAME[which(usp10sh.file$plot==TRUE)], NES=usp10sh.file$NES[which(usp10sh.file$plot==TRUE)], 
                     FDR.q.value=usp10sh.file$FDR.q.val[which(usp10sh.file$plot==TRUE)], Size=usp10sh.file$SIZE[which(usp10sh.file$plot==TRUE)])
usp10oe.file = read.csv(paste0('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final/', gseafile[5]))
usp10oe = data.frame(pathways=usp10oe.file$NAME[which(usp10oe.file$plot==TRUE)], NES=usp10oe.file$NES[which(usp10oe.file$plot==TRUE)], 
                     FDR.q.value=usp10oe.file$FDR.q.val[which(usp10oe.file$plot==TRUE)], Size=usp10oe.file$SIZE[which(usp10oe.file$plot==TRUE)])
usp10oe$pathways[4] = "ACTIVATION OF THE MRNA UPON BINDING OF THE CAP-BINDING \n COMPLEX AND EIFS, AND SUBSEQUENT BINDING TO 43S"
plots = list(fmr1oe,fmr1rnai,rinoe,rinrnai,usp10oe,usp10sh)
for (i in 1:length(plots)){
  plotname = unlist(strsplit(gseafile[i], split='.', fixed=TRUE))[1]
  tiff(paste0(plotname, '.tiff'), units="in", width=11, height=6,res=200)
  print(ggplot(plots[[i]], aes(x = NES, y = fct_reorder(pathways, negate(FDR.q.value)))) + 
          xlim(1,2.5)+
          geom_point(aes(size = Size, color = FDR.q.value)) +
          theme_bw(base_size = 14) +
          scale_colour_gradient(limits=c(0, 0.05), low="red", high = 'blue') +
          ylab(NULL) +
          labs(title=plotname, x ="Normalized Enrichment Score (NES)", y='')+
          ggtitle(plotname))
  dev.off()
}

#####heatmaps for USP10 OE
library(pheatmap)
library(hash)
negate = function(x) {return(-x)}
heat.file = read.delim("~/Desktop/FYP/GSEA/GSEA 1-8/MHC_W vs USP10_PA_OE GSEA 1-8.txt")
tp53.genes = c('lic','DOR','CkIIbeta','caz','CG17471; PIP4K','Akt1','SNF4Agamma','alc','CG6707','Caf1; Caf1','dod','CG9246',
               'SNF1A; AMPKalpha','CkIIalpha','Mi-2','Usp7','RpS27A','CG10673; Prpk','Cdk5')
mtor.genes = c('alph','Akt1','RagC; RagC-D','SNF4Agamma','eIF-4E','alc','eIF4G','Tor','SNF1A; AMPKalpha','S6k','RpS6','Rheb','FK506-bp2','14-3-3zeta','RagA; RagA-B')
auto.genes = c('ref(2)P','CkIIbeta','sw','shrb','RagC; RagC-D','HDAC6','SNF4Agamma','Dlic','alc','Atg18; Atg18a','Atg4; Atg4a','CG11975','Aut1; Atg3','Vps20','Tor',
               'Atg8a','SNF1A; AMPKalpha','vps2','CkIIalpha','Mrp4','Atg7','vps24','RpS27A','CG7627','ben','Uev1A','Rheb','CG31793','CG5498','Dhc64C','Marf','porin',
               'RagA; RagA-B','dj-1beta','ctp','CG5676','CG5789')
macro.genes = auto.genes
colnames(heat.file)[2] ='Gene Name'
heat = hash()
heat[['Regulation of TP53 Activity']] = tp53.genes; heat[['mTOR signalling']]=mtor.genes;heat[['Autophagy']]=auto.genes;heat[['Macroautophagy']]=macro.genes
setwd('/Users/rushi/Desktop/FYP/GSEA/GSEA 1-8/GSEA results final')
for (i in keys(heat)) {
  tp53 = heat.file[which(heat.file$`Gene Name`%in% heat[[i]]),][,-1]
  rownames(tp53)=tp53$`Gene Name`
  tp53 = tp53[,-1]
  #tp53 = cbind(negate(tp53[,c(1,2,3)]), tp53[,c(4,5,6)])
  if (i=='Autophagy' | i=='Macroautophagy') {
    tiff(paste0(i," heatmap.tiff"),units="in", width=6, height=10,res=250)
    print(pheatmap(tp53, color =  colorRampPalette(c("dark blue", "white","red"))(20), scale = 'row', cluster_cols = F,
                   angle_col = 45, main = i, cellwidth = 30,cellheight = 15))
    dev.off()
  } else {
    tiff(paste0(i," heatmap.tiff"),units="in", width=6, height=6,res=250)
    print(pheatmap(tp53, color =  colorRampPalette(c("dark blue", "white","red"))(20), scale = 'row', cluster_cols = F,
                   angle_col = 45, main = i, cellwidth = 30,cellheight = 15))
    dev.off()
  }
}




