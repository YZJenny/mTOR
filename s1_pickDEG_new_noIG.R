##删除Ig开头的基因

##1.先处理一下offtarget 相关的就好，应该是每个cutoff的前4组，KI-DMSO vs 200nM, KI-DMSO vs 500nM, overlap的， KI-200nM vs 500nM的。
##然后用剩下的差异基因画heatmap，做一下David与EnrichR富集分析，富集上次说也可以只用差异基因GSEA做一下看看,就是用三种方法看看有么有pathway

##2.再处理一下ontarget 相关,WT_KI_200AZD_overlap,WT_KI_200AZD_overlap
## mkdir
wp='/home/disk/ZF/mTOR_0720/res/DEG_new/'
op='/home/disk/ZF/mTOR_0720/res/DEG_new_noIg/'

FC <- 1
Q <- 0.1
CUT=paste('FC',FC,'_Q',Q,sep='')

condi=list.files(path=paste(wp,CUT,sep=''))
for(i in 1:length(condi)){
  if(!file.exists(file.path(op, paste(CUT,'/',condi[i],sep='')))){
    dir.create(file.path(op, paste(CUT,'/',condi[i],sep='')),recursive = TRUE)
  }
}

path=paste('/home/disk/ZF/mTOR_0720/res/DEG_new_noIG/FC',FC,'_Q',Q,sep='')
for(i in 1:length(condi)){
  print(condi[i])
  for(f in c('sigALL_gene.csv','sigDN_gene.csv','sigUP_gene.csv')){
    data <- read.csv(paste(wp,CUT,'/',condi[i],'/',f,sep=''))
    index <- grep('^Ig',data[,1])
    data <- data[-index,]
    write.csv(data,paste(op,CUT,'/',condi[i],'/',f,sep=''),row.names = F)
  }
}

