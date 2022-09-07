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

######Normalization
boxallnor<-apply(boxallnumb,2,function(x)(x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE))  #normalise each col, margin = 2 means col, na.rm removes missing value
######Row-wise normalization
boxallnorf<-t(apply(boxallnumb,1,function(x)(x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE))) #normalise each row, margin = 1, transpose

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

combicheck<-c(1,1,1,2,2,2,3,3,3,4,5,6,7,7,7,8,8,8,9,9,9,10,10,10,11)  #classify each col into diff types
mispercentcombi<-apply(boxallnor,1, function(x,y=0.5){
  a<-c()
  if (sum(unique(combicheck[!is.na(x)]))!=sum(seq(1,11))){  #at least got a value for each type
    a<-0
  }else{a<-1}
  return(a)
})

######Add accession number
boxallc2<-data.frame('Accession'=boxall[,1],boxallnor)  #add back accession col
missclean<-boxallc2[(mispercent & mispercentcombi),]  #filter out some rows
path1<-"Missing clean Zscore 0.5 missing"
setwd(root0)
if(!dir.exists(path1)){
  dir.create(path1)
}
write.csv(missclean,file=file.path(path1,'Clean data Zscore 0.5 missing.csv'),row.names = FALSE)

######Combat correction design
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
#originalcolname<-colnames(final)
colnames(final)<-originalcolname
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
  labs(title ="PCA Analysis", x = "PC1", y = "PC2")
dev.off()  #stop plotting PCA

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
clean1[,'CG number'] = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(clean0[,'Accession'], comp_gene_info2[,'Accession'])]  #add CG number
clean1[, 'Description'] = comp_gene_info2[, 'Gene Symbol'][match(clean0[,'Accession'], comp_gene_info2[,'Accession'])]  #add gene symbol
write.csv(clean1,file=file.path(root0,path1,"BatchRemoveWithoutControl.csv"),row.names = FALSE)

#####PCA after removing controls
originalcolname2<-colnames(clean0)
pcalabel2<-str_replace_all(colnames(clean0),"_batch\\d*","")
pcalabel2<-paste0("    ",pcalabel2)
colnames(clean0)<-pcalabel2
groups2<-c(rep('MHC/W',3),rep('USP10-PA OE',3),rep('USP10 sh',3),
           rep('FMR1 OE',3),rep('FMR1 RNAi-2',3),
           rep('RIN OE',3),rep('RIN RNAi',3))

library(factoextra)
pcacombat2<-pcafunc(clean0)  #use prcomp to do PCA
if (!dir.exists(file.path(root0,'PCA_combat'))){
  dir.create(file.path(root0,'PCA_combat'))
}
setwd(file.path(root0,'PCA_combat'))
tiff("Combat after removing control.tiff",units="in", width=11, height=10,res=250)  #safe as high resolution picture
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
  labs(title ="PCA Analysis after Removing controls", x = "PC1", y = "PC2")
dev.off()  #stop plotting PCA

######Volcano Plot
library(limma)
library(ggrepel)
source("/Users/rushi/Desktop/FYP/Yi's data/limmafunc2.R")
combatdesign0<-combatdesign[-(pos-1)]  #remove10,11,12,25 why?
treatname<-as.vector(unique(combatdesign0)[-1])  #w/o wildtype MHC
##Volcano with adj.p.val
path0<-"Volcano plot"
if(!file.exists(file.path(root0,path0))){
  dir.create(file.path(root0,path0))}
for (i in seq(1,length(treatname))){   #parse treatname
  pos<-combatdesign0 %in% c("MHC/W",treatname[i])
  colnam<-combatdesign0[combatdesign0 %in% c("MHC/W",treatname[i])]  #colnam = combatdesign0p[pos]
  colnam1<-droplevels(colnam)
  paircomp<-cbind(clean0[,1],clean0[,2:dim(clean0)[2]][,pos])  #cbind accession and treatname[i], eg. USP10
  thresh<-0.6
  limmafunc(paircomp, colnam1,path0,thresh,bioinfo)
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

##Check boxallnor against new KEGG mtor signalling pathway
library(GSA)
boxallc2[,"CG_number"] = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(boxallc2[,'Accession'], comp_gene_info2[,'Accession'])]
kegg <- GSA.read.gmt("~/Desktop/FYP/GSEA/new_Drosophila_KEGG.gmt")
mtor = unlist(k$genesets[which(k$geneset.names == "mTOR signaling pathway")])[1:100]
mtor_gene_before_remove_na = mtor[sapply(mtor, function(x) x %in% boxallc2$CG_number)]
genesymbol1 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_before_remove_na, comp_gene_info2$ANNOTATION_SYMBOL)]
df_mtor_in_boxallc2 = data.frame(mtor_gene = mtor_gene_before_remove_na, gene_symbol = genesymbol1)
mtor_gene_not_in_boxallc2 = mtor[sapply(mtor, function(x) !(x %in% boxallc2$CG_number))]
write.csv(df_mtor_in_boxallc2, file = 'mtor genes in boxallc2 before NA removal.csv')
# genesymbol2 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_not_in_boxallc2, comp_gene_info2$ANNOTATION_SYMBOL)]
# df_mtor_not_in_boxallc2 = data.frame(mtor_gene = mtor_gene_not_in_boxallc2, gene_symbol = genesymbol2)

##Check final against new KEGG mtor signalling pathway
colnames(final)[1] = 'Accession'
final = as.data.frame(final)
final$CG_number = NA
final$CG_number = comp_gene_info2[,'ANNOTATION_SYMBOL'][match(final[,'Accession'], comp_gene_info2[,'Accession'])]
mtor_gene_after_remove_na = mtor[sapply(mtor, function(x) x %in% final$CG_number)]
genesymbol2 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_after_remove_na, comp_gene_info2$ANNOTATION_SYMBOL)]
df_mtor_after_remove_na = data.frame(mtor_gene = mtor_gene_after_remove_na, gene_symbol = genesymbol2)
mtor_gene_removed = mtor_gene_before_remove_na[!(mtor_gene_before_remove_na %in% mtor_gene_after_remove_na)]
genesymbol3 = comp_gene_info2$`Gene Symbol`[match(mtor_gene_removed, comp_gene_info2$ANNOTATION_SYMBOL)]
df_mtor_removed = data.frame(mtor_gene = mtor_gene_removed, gene_symbol = genesymbol3)
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
