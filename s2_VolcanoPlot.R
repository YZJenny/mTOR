## FC的阈值是1/1.3，q值是0.05/0.1

library(ggrepel)
library(ggplot2)

FC=1
Qvalue=0.1

wp='/home/disk/ZF/mTOR_0720/res/DIFF/'
#op="/home/disk/ZF/mTOR_0720/res/Volcano/"
op="/home/disk/ZF/mTOR_0720/res/Volcano_noIg/"

condi=list.files(path=wp)
for(i in 1:length(condi)){
  print(condi[i])
  ## mkdir
  if(!file.exists(file.path(op, condi[i]))){
    dir.create(file.path(op, condi[i]),recursive = TRUE)
  }
  ##
  data<-read.table(paste(wp,condi[i],"/gene_exp.diff",sep=""),header = T)
  
  ## 删除Ig开头基因
  index <- grep('^Ig',data[,'gene'])
  data <- data[-index,]
  ##
  
  data$significant <- as.factor(ifelse(data$q_value < Qvalue & abs(data$log2.fold_change.) > log(FC,2),
                                       ifelse(data$log2.fold_change. > log(FC,2),'UP','DOWN'),'NOT'))
  gene<-data$gene
  gene_log2FC<-data$log2.fold_change.
  gene_logQ<--log10(data$q_value)
  boolean_isInf<-is.infinite(gene_log2FC)
  filter_gene<-as.vector(gene[!boolean_isInf])
  filter_gene_log2FC<-gene_log2FC[!boolean_isInf]
  filter_gene_logQ<-gene_logQ[!boolean_isInf]
  significant<-data$significant[!boolean_isInf]
  filter_data<-data.frame(filter_gene,filter_gene_log2FC,filter_gene_logQ,significant)
  boolean_tag=(filter_gene_log2FC> log(FC,2) | filter_gene_log2FC< -log(FC,2)) & filter_gene_logQ > -log10(Qvalue)
  filter_gene[which(boolean_tag==F)]=""
  pdf(paste(op,condi[i],'/FC',FC,'_Q',Qvalue,'.pdf',sep=""))
  figure<-ggplot(filter_data,aes(x=filter_gene_log2FC,y=filter_gene_logQ))
  xstt <- floor(min(filter_data$filter_gene_log2FC))
  xend <- ceiling(max(filter_data$filter_gene_log2FC))
  x <- max(abs(xstt),abs(xend))
  xstt <- -x
  xend <- x

  ystt <- floor(min(filter_data$filter_gene_logQ))
  yend <- ceiling(max(filter_data$filter_gene_logQ))

  figure1 <- figure+geom_point(aes(color=significant))+xlim(xstt,xend) + ylim(ystt,yend)+
    scale_color_manual(values = c("#377EB8","#999999","#E41A1C"))+
    labs(title=condi[i],x="log2FC",y="-log10(Qvalue)")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept = -log10(Qvalue),linetype=3)+
    geom_vline(xintercept=c(-log(FC,2),log(FC,2)),linetype=3)
  print(figure1)
  dev.off()
}
