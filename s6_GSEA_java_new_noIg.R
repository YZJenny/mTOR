library(ggplot2)
library(forcats)

wp='/home/disk/ZF/mTOR_0720/res/DIFF/'
op="/home/disk/ZF/mTOR_0720/res/GSEA_new_noIg/"
EXP_all <- read.csv('/home/disk/ZF/mTOR_0720/res/FPKM/combined_fpkm.csv') 

go <- 'GO_mouse.gmt'
kegg <- 'KEGG_mouse.gmt'


FC <- 2
Qvalue <- 0.1
CUT=paste('FC',FC,'_Q',Qvalue,sep='')

if(!file.exists(file.path(op,CUT))){
  dir.create(file.path(op,CUT),recursive = TRUE)
}

condi=list.files(path=wp)
#for(i in 1:length(condi)){
for(i in 1:3){
  print(condi[i])
  cond1 = unlist(strsplit(condi[i],'_'))[1]
  cond2 = unlist(strsplit(condi[i],'_'))[2]
  ## mkdir
  if(!file.exists(file.path(paste(op,CUT,sep='/'), condi[i]))){
    dir.create(file.path(paste(op,CUT,sep='/'), condi[i]),recursive = TRUE)
  }
  
  OUT_PATH=paste(op,CUT,'/',condi[i],'/',sep='')
  print(OUT_PATH)
  
  ## Construct EXP file
  # if(!file.exists(paste(OUT_PATH,'EXP.txt',sep=''))){
  #   data<-read.table(paste(wp,condi[i],"/gene_exp.diff",sep=""),header = T)
  #   EXP <- data[,c(3,1,8,9)]
  #   colnames(EXP) <- c('NAME','DESCRIPTION',Cond1,Cond2)
  #   write.table(EXP,paste(OUT_PATH,'EXP.txt',sep=''),col.names = T,row.names = F,quote = F,sep='\t')
  # }
  
  if(!file.exists(paste(OUT_PATH,'EXP.txt',sep=''))){
    sigUP_lst <- read.csv(paste('/home/disk/ZF/mTOR_0720/res/DEG_new_noIg/',CUT,'/',cond1,'_',cond2,'/sigUP_gene.csv',sep=''))[,1]
    sigDN_lst <- read.csv(paste('/home/disk/ZF/mTOR_0720/res/DEG_new_noIg/',CUT,'/',cond1,'_',cond2,'/sigDN_gene.csv',sep=''))[,1]
    
    EXP <- EXP_all[]
    EXP <- EXP_all[EXP_all$gene_short_name %in% c(sigUP_lst,sigDN_lst) ,
                   c(which(colnames(EXP_all) %in% c('gene_short_name','gene_id')),
                     grep(sub('-','_',cond1),colnames(EXP_all)),
                     grep(sub('-','_',cond2),colnames(EXP_all)))]
    colnames(EXP)[1:2] <- c('NAME','DESCRIPTION')
    write.table(EXP,paste(OUT_PATH,'EXP.txt',sep=''),col.names = T,row.names = F,quote = F,sep='\t')
  }
  
  ## Construct CLS fileu
  if(!file.exists(paste(OUT_PATH,'CLS.cls',sep=''))){
    L1 <- matrix(c(6,2,1),nrow = 1)
    L2 <- matrix(c('#',cond1,cond2),nrow = 1)
    L3 <-  matrix(c(rep(cond1,3),rep(cond2,3)),nrow = 1)
    CLS <- list(L1,L2,L3)
    for (i in 1:length(CLS)) {
      write.table(CLS[[i]], paste(OUT_PATH,'CLS.cls',sep=''),row.names = F,col.names = F,quote = F,sep=' ', append = T)
    }
  }
  
  ###################
  ## Run GSEA
  ###################
  for(GMT in c('GO_mouse.gmt','KEGG_mouse.gmt')){
    system(
      paste('/home/disk/fyh/tools/GSEA_Linux_4.0.2/gsea-cli GSEA -set_min 3 -set_max 2000 -out ',OUT_PATH,' -rnd_seed 13579 -res ',OUT_PATH,'EXP.txt -cls ',OUT_PATH,'CLS.cls -gmx /home/yzj/publicData/GMT/',GMT,' -plot_top_x 2000 -permute gene_set -collapse false',sep='')
    )
    system(
      paste('mv ',OUT_PATH,'my_analysis* ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],sep=''))
    system(
      paste('mv ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',cond1,'_*.xls ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',cond1,'.xls',sep='')
    )
    system(
      paste('mv ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',cond2,'_*.xls ',OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',cond2,'.xls',sep='')
    )
    
    
    ###################
    ## DotPlot
    ###################
    if(file.exists(paste(OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],sep=''))){
      A<-read.table(file=paste(OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',cond1,'.xls',sep=''),sep="\t",head=T,row.names=1,fill=T)
      A<- A[A$FDR.q.val < 0.05,]
      B<-read.table(file=paste(OUT_PATH,'Gsea.',unlist(strsplit(GMT,'\\.'))[1],'/gsea_report_for_',cond2,'.xls',sep=''),sep="\t",head=T,row.names=1,fill=T)
      B<- B[B$FDR.q.val < 0.05,]
      
      gsea_dot <- rbind(A,B)
      print(dim(gsea_dot))
      print(summary(gsea_dot$NES))
      
      gsea_dot=gsea_dot[order(gsea_dot$FDR.q.val,decreasing=F),]
      if(nrow(gsea_dot) > 0){
        ## mkdir
        
        if(!file.exists(paste(OUT_PATH,'Figure.',unlist(strsplit(GMT,'\\.'))[1],sep=''))){
          dir.create(file.path(OUT_PATH, paste('Figure.',unlist(strsplit(GMT,'\\.'))[1],sep='')),recursive = TRUE)
        }
        
        
        FIG_PATH=paste(OUT_PATH,'Figure.',unlist(strsplit(GMT,'\\.'))[1],'/',sep='')
        if(nrow(gsea_dot) > 20){
          gsea_dot=gsea_dot[1:20,]
        }
        
        p<- ggplot(gsea_dot, aes(x = NES, y = fct_reorder(GS.br..follow.link.to.MSigDB, NES))) +
          geom_point(aes(color = FDR.q.val,size = SIZE)) +
          coord_cartesian(xlim=c(floor(min(gsea_dot$NES)),ceiling(max(gsea_dot$NES))))+
          scale_x_continuous(breaks=floor(min(gsea_dot$NES)):ceiling(max(gsea_dot$NES)))+
          scale_colour_gradientn(limits=c(0, 1), colours=rainbow(6)) +
          theme_bw(base_size = 14) +
          ylab(NULL) +scale_fill_brewer(palette = 'Accent')+theme(panel.background=element_rect(fill='grey95'),panel.border = element_blank())+
          geom_vline(xintercept = 0,size=1)+
          theme(axis.text=element_text(size=9),axis.text.x = element_text(size = 9))+
          geom_segment(aes(x=-2,y=0,xend=2,yend=0))
        
        pdf(paste(FIG_PATH,"DotPlot.pdf",sep=''),8,10)
        print(p)
        dev.off()
      }else{
        print('no GSEA res!')
      } 
    }else{
      print('no Files!')
    } 
  }
}


