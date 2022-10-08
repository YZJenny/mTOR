.libPaths("/home/yzj/R/x86_64-pc-linux-gnu-library/4.0")
library(enrichR)
library(RColorBrewer)
library("VennDiagram")
color <- brewer.pal(9,"Set1");
mycol <- colorRampPalette(color)(9)
dbs <- listEnrichrDbs() 

height <- function(df){
  if(nrow(df) >= 1 & nrow(df) <= 20 ){
    height=7
  } else if(nrow(df) > 20 & nrow(df) <= 60 ) {
    height=15
  } else if(nrow(df) > 60 & nrow(df) <= 200 ) {
    height=18
  } else if (nrow(df) > 200 ) {
    height=140
  }
  return(height)
} 


#wp='/home/disk/ZF/mTOR_0720/res/DEG_new/'
#op='/home/disk/ZF/mTOR_0720/res/EnrichR_new/'

wp='/home/disk/ZF/mTOR_0720/res/DEG_new_noIg/'
op='/home/disk/ZF/mTOR_0720/res/EnrichR_new_noIg/'

FC <- 2
Qvalue <- 0.1
CUT=paste('FC',FC,'_Q',Qvalue,sep='')

## 1. 两两比
if(!file.exists(file.path(op,CUT))){
  dir.create(file.path(op,CUT),recursive = TRUE)
}


condi=list.files(path=paste(wp,CUT,sep=''))
#for(i in 1:length(condi)){
for(i in 1:6){
  TAG <- condi[i]
  print(TAG)
  
  TYPE <- c('ALL','UP','DN')
  for(i in 1:length(TYPE)){
    type=TYPE[i] 
    sigG_lst <- read.csv(paste(wp,'FC',FC,'_Q',Qvalue,'/',TAG,'/sig',type,'_gene.csv',sep=''))[,1]
    
    if(length(sigG_lst)>0){
      ## GO BP
      enrichr_res=enrichr(sigG_lst,databases=c('GO_Biological_Process_2018'))
      tmp=enrichr_res[[1]][,c(1,4)]
      print(dim(tmp))
      if(nrow(tmp) > 0){
        tmp=tmp[which(tmp$Adjusted.P.value < 0.05),];rownames(tmp)=tmp[,1]
        tmp=tmp[order(tmp[,2],decreasing=F),]
        print(dim(tmp))
        
        if(nrow(tmp) > 0){
          write.csv(tmp,paste0(op,CUT,'/',TAG,"_BP_",type,".csv"),quote=F, row.names=F )
          if(nrow(tmp) >50){
            tmp <-tmp[1:50,]
          }else{
            tmp <- tmp
          }
          h=height(tmp)
          mycolor<-rep(mycol[2],nrow(tmp))
          pdf(paste0(op,CUT,'/',TAG,"_BP_",type,".pdf"),width=13,height=h)
          bb<-barplot(c(-log(tmp[,2],10)),horiz=TRUE,axisnames=F,main=TAG,cex.axis=0.9,col=mycolor,border=NA,xlim=c(0,25),space=0.3,xlab="-log10(qvalue)",ylab='GO BP')
          text(0.001,bb,rownames(tmp),cex=0.8,xpd=T,srt=0,pos=4)
          dev.off() 
        }else{
          print('no sigBP!')
        }
      }else{
        print('no BP!')
      }
      
      ## KEGG
      enrichr_res=enrichr(sigG_lst,databases=c('KEGG_2019_Human'))
      tmp=enrichr_res[[1]][,c(1,4)]
      print(dim(tmp))
      if(nrow(tmp) > 0){
        tmp=tmp[which(tmp$Adjusted.P.value < 0.05),];rownames(tmp)=tmp[,1]
        tmp=tmp[order(tmp[,2],decreasing=F),]
        print(dim(tmp))
        
        if(nrow(tmp) > 0){
          write.csv(tmp,paste0(op,CUT,'/',TAG,"_KEGG_",type,".csv"),quote=F, row.names=F )
          if(nrow(tmp) >50){
            tmp <-tmp[1:50,]
          }else{
            tmp <- tmp
          }
          h=height(tmp)
          mycolor<-rep(mycol[2],nrow(tmp))
          pdf(paste0(op,CUT,'/',TAG,"_KEGG_",type,".pdf"),width=13,height=h)
          bb<-barplot(c(-log(tmp[,2],10)),horiz=TRUE,axisnames=F,main=TAG,cex.axis=0.9,col=mycolor,border=NA,xlim=c(0,25),space=0.3,xlab="-log10(qvalue)",ylab='GO BP')
          text(0.001,bb,rownames(tmp),cex=0.8,xpd=T,srt=0,pos=4)
          dev.off() 
        }else{
          print('no sigKEGG!')
        }
      }else{
        print('no KEGG!')
      }
      
    }else{
      print('no genes')
    }
    
  }
}


