library(ggplot2)
library(pheatmap)
library(dplyr)

height <- function(df){
  if(nrow(df) >= 1 & nrow(df) <= 20 ){
    height=7
  } else if(nrow(df) > 20 & nrow(df) <= 60 ) {
    height=10
  } else if(nrow(df) > 60 & nrow(df) <= 200 ) {
    height=15
  } else if (nrow(df) > 200 ) {
    height=140
  }
  return(height)
} 

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


EXP <- read.csv('/home/disk/ZF/mTOR_0720/res/FPKM/combined_fpkm.csv')
#op="/home/disk/ZF/mTOR_0720/res/Heatmap"
op="/home/disk/ZF/mTOR_0720/res/Heatmap_noIg"

FC <- 1
Qvalue <- 0.1
CUT=paste('FC',FC,'_Q',Qvalue,sep='')
#DEG_file='DEG_new'
DEG_file='DEG_new_noIg'

## 1. 两两比
cond=c('WT-DMSO','WT-200nMAZD','WT-500nMAZD','KI-DMSO','KI-200nMAZD','KI-500nMAZD')
combn_lst <- t(combn(cond,2))

if(!file.exists(file.path(op,CUT))){
  dir.create(file.path(op,CUT),recursive = TRUE)
}

## offtarget
for(i in 13:nrow(combn_lst)){
  cond1 <- combn_lst[i,1]
  cond2 <- combn_lst[i,2]
  TAG <- paste(cond1,cond2,sep='_')
  
  sigUP_lst <- read.csv(paste('/home/disk/ZF/mTOR_0720/res/',DEG_file,'/',CUT,'/',cond1,'_',cond2,'/sigUP_gene.csv',sep=''))[,1]
  sigDN_lst <- read.csv(paste('/home/disk/ZF/mTOR_0720/res/',DEG_file,'/',CUT,'/',cond1,'_',cond2,'/sigDN_gene.csv',sep=''))[,1]
  
  df <- EXP[EXP$gene_short_name %in% c(sigUP_lst,sigDN_lst) ,
            c(colnames(EXP)=='gene_short_name',grep(sub('-','_',cond1),colnames(EXP)),grep(sub('-','_',cond2),colnames(EXP)))]
  print(dim(df))
  rownames(df) <- df$gene_short_name
  df <- df[,-1]
  
  h=height(df)
  tmp=df
  stmp=apply((tmp+1),2,log,2)
  stmp=t(apply(tmp,1,scale))
  rownames(stmp)=rownames(df)
  colnames(stmp)=colnames(df)
  this_diff=apply(stmp[,1:3],1,mean)-apply(stmp[,4:6],1,mean)
  stmp=stmp[order(this_diff),]
  
  bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))
  p <- pheatmap(stmp,scale = "none",cexRow=0.5,fontsize=8,fontsize_row = 8,fontsize_col=8,cluster_col = F, 
                cluster_row = F,treeheight_row=0, treeheight_col=0,border_color=NA,
                color = colorRampPalette(colors = c("yellow2","blue2"))(400),
                breaks=seq(-2,2,0.01),
                legend_labels = c( "-2", "-1", "0", "1", "2"))
  save_pheatmap_pdf(p, paste(op,'/',CUT, '/',TAG,'.pdf',sep=''),4,height = h)
}

### 2. 3列比
# FC <- 2
# Qvalue <- 0.05
# CUT=paste('FC',FC,'_Q',Qvalue,sep='')

sigUP_lst <- read.csv(paste('/home/disk/ZF/mTOR_0720/res/',DEG_file,'/',CUT,'/KI200_500_overlap/sigUP_gene.csv',sep=''))[,1]
sigDN_lst <- read.csv(paste('/home/disk/ZF/mTOR_0720/res/',DEG_file,'/',CUT,'/KI200_500_overlap/sigDN_gene.csv',sep=''))[,1]

df <- EXP[EXP$gene_short_name %in% c(sigUP_lst,sigDN_lst) ,
          c(colnames(EXP)=='gene_short_name',
            grep('KI_DMSO',colnames(EXP)),grep('KI_200nMAZD',colnames(EXP)),grep('KI_500nMAZD',colnames(EXP)))]

rownames(df) <- df$gene_short_name
df <- df[,-1]
df <- as.data.frame(t(df))
genes <- c(sigUP_lst,sigDN_lst)
df <- as.data.frame(t(select(df,genes)))
dim(df)

h=10
tmp=df
stmp=apply((tmp+1),2,log,2)
stmp=t(apply(tmp,1,scale))
rownames(stmp)=rownames(df)
colnames(stmp)=colnames(df)

p <- pheatmap(stmp,scale = "none",cexRow=0.5,fontsize=8,fontsize_row = 8,fontsize_col=8,cluster_col = F,
              cluster_row = F,treeheight_row=0, treeheight_col=0,border_color=NA,
              color = colorRampPalette(colors = c("yellow2","blue2"))(400),
              breaks=seq(-2,2,0.01),
              legend_labels = c( "-2", "-1", "0", "1", "2"))
save_pheatmap_pdf(p, paste(op,'/',CUT, '/KI200_500_overlap.pdf',sep=''),6,height = h)



### 4. on-target gene
op_DEG=paste('/home/disk/ZF/mTOR_0720/res/',DEG_file,'/',sep='')
#op_heatmap='/home/disk/ZF/mTOR_0720/res/Heatmap/'
op_heatmap='/home/disk/ZF/mTOR_0720/res/Heatmap_noIg/'


AZD=500

# FC=1.3
# Qvalue=0.05 
# CUT=paste('FC',FC,'_Q',Qvalue,sep='')


if(AZD==200){
  TAG='WT_KI_200AZD_overlap'
}else if(AZD==500){
  TAG='WT_KI_500AZD_overlap'
}



TYPE <- c('ALL','UP','DN')
for(i in 1:length(TYPE)){
  type=TYPE[i] #DN/UP
  sigG_lst <- read.csv(paste(op_DEG,'FC',FC,'_Q',Qvalue,'/',TAG,'/sig',type,'_gene.csv',sep=''))[,1]
  
  df <- EXP[EXP$gene_short_name %in% sigG_lst ,
            c(colnames(EXP)=='gene_short_name',
              grep('KI_DMSO',colnames(EXP)),grep('KI_200nMAZD',colnames(EXP)),grep('KI_500nMAZD',colnames(EXP)))]
  
  rownames(df) <- df$gene_short_name
  df <- df[,-1]
  df <- as.data.frame(t(df))
  df <- as.data.frame(t(select(df,sigG_lst)))
  dim(df)
  
  h=height(df)
  tmp=df
  stmp=apply((tmp+1),2,log,2)
  stmp=t(apply(tmp,1,scale))
  rownames(stmp)=rownames(df)
  colnames(stmp)=colnames(df)
  
  p <- pheatmap(stmp,scale = "none",cexRow=0.5,fontsize=8,fontsize_row = 8,fontsize_col=8,cluster_col = F,
                cluster_row = F,treeheight_row=0, treeheight_col=0,border_color=NA,
                color = colorRampPalette(colors = c("yellow2","blue2"))(400),
                breaks=seq(-2,2,0.01),
                legend_labels = c("-2", "-1", "0", "1", "2"))
  save_pheatmap_pdf(p, paste(op_heatmap,'/',CUT, '/',TAG,'_',type,'.pdf',sep=''),6,height = h)
}

