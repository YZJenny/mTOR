### 画 barplot

library(RColorBrewer)
color <- brewer.pal(9,"Set1");
mycol <- colorRampPalette(color)(9)

get_barplot <- function(up_file,dn_file,path,type){
  if(file.exists(up_file)){
    up_res <- read.table(up_file,sep='\t',stringsAsFactors =F, header =T,fill = T,quote = "")[,c(2,13)]
    up_res <- up_res[up_res$FDR < 0.05,]
    if(nrow(up_res)==0){
      up_res <- c()
    }else if(nrow(up_res) > 0 & nrow(up_res) <= 20){
      up_res <- up_res
      rownames(up_res) <- up_res$Term
    }else{
      up_res <- up_res[1:20,]
      rownames(up_res) <- up_res$Term
    }
  }else{
    up_res <- c()
  }
  
  if(file.exists(dn_file)){
    dn_res <- read.table(dn_file,sep='\t',stringsAsFactors =F, header =T,fill = T,quote = "")[,c(2,13)]
    dn_res <- dn_res[dn_res$FDR < 0.05,]
    if(nrow(dn_res)==0){
      dn_res <- c()
    }else if(nrow(dn_res) > 0 & nrow(dn_res) <= 20){
      dn_res <- dn_res
      rownames(dn_res) <- dn_res$Term
    }else{
      dn_res <- dn_res[1:20,]
      rownames(dn_res) <- dn_res$Term
    }
  }else{
    dn_res <- c()
  }
  
  
  if(length(dn_res) > 0 & length(up_res) > 0){
    mycolor<-c(rep(mycol[2],nrow(dn_res)),rep(mycol[1],nrow(up_res)))
    end_X <- round(max(-log(dn_res[,2],10),-log(up_res[,2],10)))+2
    
    pdf(paste0(path,'/',condi[i],'/sig',type,'.pdf'),width=10,height=7)
    bb<-barplot(c(-log(dn_res[,2],10),-log(up_res[,2],10)),
                horiz=TRUE,axisnames=F,main=condi[i],cex.axis=0.9,col=mycolor,border=NA,
                xlim=c(0,end_X),space=0.3,xlab="-log10(qvalue)",ylab=type)
    legend("bottomright",legend=c("up-regulated","down-regulated"),fill=c(mycol[1],mycol[2]))
    text(0.001,bb,c(rownames(dn_res),rownames(up_res)),cex=0.8,xpd=T,srt=0,pos=4)
    dev.off()  
  }else if(length(dn_res) > 0 & length(up_res) == 0){
    mycolor<-rep(mycol[2],nrow(dn_res))
    end_X <- round(max(-log(dn_res[,2],10)))+2
    
    pdf(paste0(path,'/',condi[i],'/sig',type,'.pdf'),width=10,height=7)
    bb<-barplot(-log(dn_res[,2],10),
                horiz=TRUE,axisnames=F,main=condi[i],cex.axis=0.9,col=mycolor,border=NA,xlim=c(0,end_X),space=0.3,xlab="-log10(qvalue)",ylab=type)
    legend("bottomright",legend=c("up-regulated","down-regulated"),fill=c(mycol[1],mycol[2]))
    text(0.001,bb,rownames(dn_res),cex=0.8,xpd=T,srt=0,pos=4)
    dev.off() 
  }else if(length(dn_res) == 0 & length(up_res) > 0){
    mycolor<-rep(mycol[1],nrow(up_res))
    end_X <- round(max(-log(up_res[,2],10)))+2
    
    pdf(paste0(path,'/',condi[i],'/sig',type,'.pdf'),width=10,height=7)
    bb<-barplot(-log(up_res[,2],10),
                horiz=TRUE,axisnames=F,main=condi[i],cex.axis=0.9,col=mycolor,border=NA,xlim=c(0,end_X),space=0.3,xlab="-log10(qvalue)",ylab=type)
    legend("bottomright",legend=c("up-regulated","down-regulated"),fill=c(mycol[1],mycol[2]))
    text(0.001,bb,rownames(up_res),cex=0.8,xpd=T,srt=0,pos=4)
    dev.off() 
  }
}

get_barplot_all <- function(david_file,path,type){
  if(file.exists(david_file)){
    david_res <- read.table(david_file,sep='\t',stringsAsFactors =F, header =T,fill = T,quote = "")[,c(2,13)]
    david_res <- david_res[david_res$FDR < 0.05,]
    if(nrow(david_res)==0){
      david_res <- c()
    }else if(nrow(david_res) > 0 & nrow(david_res) <= 20){
      david_res <- david_res
      rownames(david_res) <- david_res$Term
    }else{
      david_res <- david_res[1:20,]
      rownames(david_res) <- david_res$Term
    }
  }else{
    david_res <- c()
  }
  
  if(length(david_res) > 0){
    mycolor<-c(rep(mycol[2],nrow(david_res)))
    end_X <- round(max(-log(david_res[,2],10)))+2
    
    pdf(paste0(path,'/',condi[i],'/sig',type,'_All.pdf'),width=10,height=7)
    bb<-barplot(c(-log(david_res[,2],10)),
                horiz=TRUE,axisnames=F,main=condi[i],cex.axis=0.9,col=mycolor,border=NA,xlim=c(0,end_X),space=0.3,xlab="-log10(qvalue)",ylab=type)
    legend("bottomright",legend=c("sig-regulated"),fill=c(mycol[2]))
    text(0.001,bb,c(rownames(david_res)),cex=0.8,xpd=T,srt=0,pos=4)
    dev.off()  
  }
}


## mkdir
# wp='/home/disk/ZF/mTOR_0720/res/DEG/'
# op='/home/disk/ZF/mTOR_0720/res/DAVID'

wp='/home/disk/ZF/mTOR_0720/res/DEG_new_noIg/'
op='/home/disk/ZF/mTOR_0720/res/DAVID_new_noIg'


FC <- 1
Q <- 0.05
CUT=paste('FC',FC,'_Q',Q,sep='')

condi=list.files(path=paste(wp,CUT,sep=''))
for(i in 1:length(condi)){
  if(!file.exists(file.path(op, paste(CUT,'/',condi[i],sep='')))){
    dir.create(file.path(op, paste(CUT,'/',condi[i],sep='')),recursive = TRUE)
  }
}


path=paste('/home/disk/ZF/mTOR_0720/res/DAVID_new_noIg/FC',FC,'_Q',Q,sep='')
condi=list.files(path=path)
for(i in 1:length(condi)){
  print(condi[i])
  
  ## GO BP
  upBP_file <- paste(path,'/',condi[i],'/sigUP_GO.txt',sep='')
  dnBP_file <- paste(path,'/',condi[i],'/sigDN_GO.txt',sep='')
  print(file.exists(upBP_file))
  print(file.exists(dnBP_file))
  
  get_barplot(upBP_file,dnBP_file,path,'GO')
  
  BP_file <- paste(path,'/',condi[i],'/sigALL_GO.txt',sep='')
  print(file.exists(BP_file))
  get_barplot_all(BP_file,path,'GO')
  
  ## KEGG  
  upKEGG_file <- paste(path,'/',condi[i],'/sigUP_KEGG.txt',sep='')
  dnKEGG_file <- paste(path,'/',condi[i],'/sigDN_KEGG.txt',sep='')
  print(file.exists(upKEGG_file))
  print(file.exists(dnKEGG_file))
  get_barplot(upKEGG_file,dnKEGG_file,path,'KEGG')
  
  KEGG_file <- paste(path,'/',condi[i],'/sigALL_KEGG.txt',sep='')
  print(file.exists(KEGG_file))
  get_barplot_all(KEGG_file,path,'KEGG')
  
}


### 2. on-target gene
## mkdir
op='/home/disk/ZF/mTOR_0720/res/DAVID'

FC <- 1
Q <- 0.05
CUT=paste('FC',FC,'_Q',Q,sep='')

condi=c('WT_KI_200AZD_overlap','WT_KI_500AZD_overlap')
for(i in 1:length(condi)){
  if(!file.exists(file.path(op, paste(CUT,'/',condi[i],sep='')))){
    dir.create(file.path(op, paste(CUT,'/',condi[i],sep='')),recursive = TRUE)
  }
}
path=paste('/home/disk/ZF/mTOR_0720/res/DAVID/FC',FC,'_Q',Q,sep='')

for(i in 1:length(condi)){
  print(condi[i])
  ## GO BP
  upBP_file <- paste(path,'/',condi[i],'/sigUP_GO.txt',sep='')
  dnBP_file <- paste(path,'/',condi[i],'/sigDN_GO.txt',sep='')
  print(file.exists(upBP_file))
  print(file.exists(dnBP_file))
  
  get_barplot(upBP_file,dnBP_file,path,'GO')
  
  BP_file <- paste(path,'/',condi[i],'/sigALL_GO.txt',sep='')
  print(file.exists(BP_file))
  get_barplot_all(BP_file,path,'GO')
  
  ## KEGG  
  upKEGG_file <- paste(path,'/',condi[i],'/sigUP_KEGG.txt',sep='')
  dnKEGG_file <- paste(path,'/',condi[i],'/sigDN_KEGG.txt',sep='')
  print(file.exists(upKEGG_file))
  print(file.exists(dnKEGG_file))
  get_barplot(upKEGG_file,dnKEGG_file,path,'KEGG')
  
  KEGG_file <- paste(path,'/',condi[i],'/sigALL_KEGG.txt',sep='')
  print(file.exists(KEGG_file))
  get_barplot_all(KEGG_file,path,'KEGG')
  
}
