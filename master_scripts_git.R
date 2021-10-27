library(DESeq2)
library(Rsubread)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(caret)
library(DescTools)
library(pROC)
library(lmtest)
library(ROCR)
library(glmnet)
library(ggsignif)
library(ggpubr)
library(tidyverse)
library(gridExtra)

theme_set(theme_bw())
set.seed(2021)

############################################
######## Data Ingest ######################
############################################

source("master_functions_Oct2021.R")

dir_path<-"dir"
output_path<-"dir/output/"
dir.create(output_path)

fCounts <- read.csv(paste0(dir_path,"supplemental_table2.csv"),stringsAsFactors = F)
row.names(fCounts)<-fCounts[,1]
fCounts<-fCounts[,-1]
key_df<-read.csv(paste0(dir_path,"supplemental_table1.csv"),stringsAsFactors = F)
coldata <- colnames(fCounts)
coldata<-data.frame(DEIDENT=coldata)
coldata<-merge(coldata,key_df,"DEIDENT")
coldata<-coldata %>% arrange(as.numeric(column_order))
counts<-fCounts


############################################
####### DESeq ##############################
############################################

coldata<-coldata %>% mutate(RaceDescription=ifelse(RaceDescription!="Caucasian","Other",RaceDescription))
coldata$scaled_age<-scale(coldata$Age)
coldata$scaled_bmi<-scale(coldata$BMI)
coldata$smoking_status<-factor(coldata$smoking_status)
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=coldata, 
                              design=~class+Age+factor(smoking_status)+BMI+GenderDescription+RaceDescription)
dds<- DESeq(dds)

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
resLFC_clean <- lfcShrink(ddsClean, coef="class_HCC_vs_Cirrhosis", type="apeglm")
cleanresOrdered <- resLFC_clean[order(resLFC_clean$pvalue),]
res1<-as.data.frame(cleanresOrdered)
write.csv(res1,
          paste0(output_path,"DESeq_HCCvCirrhosis_cleanres.csv"))
pdf(paste0(output_path,"DESeqPlots.pdf"))
plotMA(resLFC_clean)
graphics.off()

p1 <- EnhancedVolcano(resLFC_clean,
                      lab = rownames(resLFC_clean),
                      x = "log2FoldChange",
                      y = "padj",
                      pCutoff = 0.05,
                      FCcutoff = 5,
                      xlim = c(-5.5, 5.5),
                      ylim = c(0, -log10(5E-6)),
                      pointSize = 2,
                      labSize = 3.0,
                      title="",
                      subtitle = "",
                      caption="",
                      legendPosition = "right",
                      legendLabSize = 8,
                      colAlpha = 0.85,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      ylab =bquote("~-Log[10]~ FDR P"))
p1
pdf(paste0(output_path,"fig1a.pdf"),width=8,height=6)
p1
graphics.off()

### Compare with tissue
geo2r<-fread("geo2r.txt",stringsAsFactors = F,data.table = F)
geo2r<-geo2r %>% filter(grepl("hsa-",miRNA_ID_LIST))
dim(geo2r)
sum(geo2r$adj.P.Val<0.05)
res1$ID<-row.names(res1)
clean<-res1

sum(geo2r$miRNA_ID_LIST %in% clean$ID)
sum(geo2r$adj.P.Val<0.05)
sum(geo2r$miRNA_ID_LIST %in% clean$ID&geo2r$adj.P.Val<0.05)

sum(clean$ID %in% geo2r$miRNA_ID_LIST&clean$padj<0.05,na.rm=T)

gc<-merge(geo2r,res1,by.y="ID",by.x="miRNA_ID_LIST")
dim(gc)

gc[which(gc$adj.P.Val<0.05&gc$padj<0.05),] %>% arrange(padj,adj.P.Val)

sum(gc$adj.P.Val<0.05&gc$padj<0.05&gc$logFC<0&gc$log2FoldChange<0,na.rm = T)
sum(gc$adj.P.Val<0.05&gc$padj<0.05&gc$logFC<0&gc$log2FoldChange>0,na.rm = T)
sum(gc$adj.P.Val<0.05&gc$padj<0.05&gc$logFC>0&gc$log2FoldChange>0,na.rm = T)
sum(gc$adj.P.Val<0.05&gc$padj<0.05&gc$logFC>0&gc$log2FoldChange<0,na.rm = T)

gc$miRNA_ID_LIST[which(gc$adj.P.Val<0.05&gc$padj<0.05&gc$logFC>0&gc$log2FoldChange<0)]
gc$miRNA_ID_LIST[which(gc$adj.P.Val<0.05&gc$padj<0.05&gc$logFC<0&gc$log2FoldChange<0)]

gc_20<-gc %>% select(miRNA_ID_LIST,ID,logFC,adj.P.Val,log2FoldChange,padj) %>%
  arrange(padj,adj.P.Val) %>% head(20)

write.csv(gc_20,file.path(output_path,"tissue_saliva_20.csv"),row.names = F)

res1[grepl("hsa-mir-106",res1$ID),]
clean[grepl("hsa-mir-106",clean$ID),]
clean[grepl("hsa-mir-125",clean$ID),]

### CLD Specific Analysis
counts_cld<-counts[,colnames(counts) %in% coldata$DEIDENT[!is.na(coldata$cirrhosis.vs.fibrosis)]]
dds_cld <- DESeqDataSetFromMatrix(countData=counts_cld, 
                              colData=coldata[!is.na(coldata$cirrhosis.vs.fibrosis),], 
                              design=~class+Age+factor(smoking_status)+BMI+GenderDescription+RaceDescription)
dds_cld<- DESeq(dds_cld)

ddsClean_cld <- replaceOutliersWithTrimmedMean(dds_cld)
ddsClean_cld <- DESeq(ddsClean_cld)
resLFC_clean_cld <- lfcShrink(ddsClean_cld, coef="class_HCC_vs_Cirrhosis", type="apeglm")
cleanresOrdered_cld <- resLFC_clean_cld[order(resLFC_clean_cld$pvalue),]
res2<-as.data.frame(cleanresOrdered_cld)
write.csv(res2,
          paste0(output_path,"DESeq_HCCvCirrhosis_cleanres_cld.csv"))

res1<-res1 %>% filter(padj<0.05)
res2<-res2 %>% filter(padj<0.05)

dim(res1)
dim(res2)

sum(row.names(res2) %in% row.names(res1))
sum(res2$log2FoldChange<0)/nrow(res2)
sum(res2$log2FoldChange< -2)/nrow(res2)

res2 %>% arrange(padj) %>% head()
write.csv(res2,
          paste0(output_path,"DESeq_HCCvCirrhosis_cleanres_cld_sig.csv"))


p1 <- EnhancedVolcano(resLFC_clean_cld,
                      lab = rownames(resLFC_clean_cld),
                      x = "log2FoldChange",
                      y = "padj",
                      pCutoff = 0.05,
                      FCcutoff = 5,
                      xlim = c(-5.5, 5.5),
                      ylim = c(0, -log10(5E-6)),
                      pointSize = 2,
                      labSize = 3.0,
                      title="",
                      subtitle = "",
                      caption="",
                      legendPosition = "right",
                      legendLabSize = 8,
                      colAlpha = 0.85,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      ylab =bquote("~-Log[10]~ FDR P"))
p1
pdf(paste0(output_path,"fig1b.pdf"),width=8,height=6)
p1
graphics.off()


######################################################
########## Heatmaps and Clustering ###################
######################################################

### Setup annotations

ramp <- colorRampPalette(colors=c("green", "darkgreen","grey","darkred","red"))(20)
df<-coldata %>% select(DEIDENT,class,cirrhosis.vs.fibrosis) %>%
  mutate(cirrhosis.vs.fibrosis=ifelse(!is.na(cirrhosis.vs.fibrosis),
                                      cirrhosis.vs.fibrosis,"none documented"),
         class=ifelse(class=="HCC","HCC","Control"))
row.names(df)<-df$DEIDENT
df<-df[,-1]
colnames(df)<-c("class","Liver Disease")
var1<-c("#2E86AB","#E88873")
names(var1)<-c("Control","HCC")
var2<-c("#3EC300","#FFCB47","#796293")
names(var2)<-c("cirrhosis","fibrosis","none documented")
ann_colors<-list(class=var1,`Liver Disease`=var2)

### Heatmap of all miRNA

dist_mat<-assay(ddsClean)
out<-pheatmap(dist_mat, 
              color = ramp,
              show_rownames=F,
              show_colnames = F,
              cluster_rows=F,
              cluster_columns=T,
              cutree_cols = 8,
              scale = "row",
              treeheight_col = 50,
              cluster_cols=T, annotation_col=df,main=paste(""),
              annotation_colors = ann_colors)
save_pheatmap_pdf(out,paste(output_path,"all_mirna.pdf",sep=""))

### Clustering

assay_dds_clean<-assay(ddsClean)
clust_method<-"manhattan"
p_vals<-0.05
fold_change<-1.5
clust_num<-4
for(p_vals in c(0.001,0.01,0.05)){
  for(fold_change in c(1,1.5,2)){
    for(clust_num in 2:5){
top.genes <- which(abs(resLFC_clean$log2FoldChange)>fold_change & resLFC_clean$padj<p_vals)
topgenes <- rownames(unique(resLFC_clean[top.genes,]))
out<-pheatmap(assay_dds_clean[top.genes,], 
              color = ramp,
              show_rownames=F,
              show_colnames = T,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = clust_method,
              clustering_method = "ward.D",
              cutree_cols = clust_num,
              cluster_rows=F,
              scale = "row",
              treeheight_col = 50,
              cluster_cols=T, annotation_col=df,main=paste("P-val<",p_vals,"fold_changes>",fold_change,"k=",clust_num),
              annotation_colors = ann_colors)
save_pheatmap_pdf(out,paste(output_path,clust_method,"_",p_vals,"_",fold_change,"_",clust_num,"_ann.pdf",sep=""),
                  height = 5)
    }
  }
}
############################################
######### Modeling #########################
############################################

select.genes <- rownames(resLFC_clean[which(resLFC_clean$padj<0.05),])
df <- ddsClean[select.genes,]
df <- counts(df, normalized=T) %>% t() %>% data.frame()
df$class<-dds$class

df_plot<-df
df_plot<-df_plot %>% mutate(class=ifelse(class=="HCC","HCC","CLD"))
pdf(file.path(output_path,"pca_plot.pdf"),width =8.5,height = 6)
x<-prcomp(df[,-ncol(df)],scale. = T,center = T)
prop_x<-summary(x)
x<-x$x %>% data.frame()
x$ID<-row.names(x)
df_plot$ID<-row.names(df_plot)
x<-merge(x,df_plot,"ID")
comp1<-colnames(x)[2]
comp2<-colnames(x)[2]
plot_pca<-function(x,comp1,comp2,prop_x){
  p<-NULL
  if(comp1==comp2){
    loc_x<-min(x[,comp1])+(max(x[,comp1])-min(x[,comp1]))/2
    p<-ggplot(data=x,aes_string(comp1,comp2))+geom_point(col="white")+
      geom_text(aes(loc_x,loc_x,label=
                                  paste0(round(prop_x$importance[2,comp1]*100,2),"%")),size=5)+
      coord_cartesian(clip = "off")+theme(axis.title = element_blank())
  }
  else{p<-ggplot(data=x,aes_string(comp1,comp2,color="class"))+geom_text(aes(label=ID),size=3)+
    theme(legend.position = "none",axis.title = element_blank())+coord_cartesian(clip = "off")
  }
  return(p)
  }

list_of_plots<-list()
for(comp1 in colnames(x)[2:5]){
  for(comp2 in colnames(x)[2:5]){
    list_of_plots[[paste0(comp1,comp2)]]<-plot_pca(x,comp1,comp2,prop_x)
  }
}
do.call("grid.arrange",c(list_of_plots,ncol=4))
graphics.off()



### Model 
track_models_1<-mirna.model(df,"All",coldata)

## Modeling liver disease samples

samples_with_cirrhosis<-coldata$DEIDENT[!is.na(coldata$cirrhosis.vs.fibrosis)]
train_cirr<-df[row.names(df) %in% samples_with_cirrhosis,]

track_models_2<-mirna.model(train_cirr,"CLD",coldata)

acc_metrics<-rbind.data.frame(track_models_1,track_models_2)
acc_metrics %>% filter(Balanced.Accuracy==`Max Balanced Accuracy`)

## Visualize 


roc_df<-acc_metrics
roc_df<-roc_df %>% mutate(color=ifelse(label=="All","goldenrod","royalblue"))%>%
  arrange(Sensitivities)
p<-ggplot(data=roc_df,aes(x=1-Specificities,y=Sensitivities,color=label))+
      geom_segment(x=0,xend=1,y=0,yend=1,lty=2,color="black")+
      geom_line(size=0.6)+
      geom_point(size=1)+
    scale_color_manual(values=unique(roc_df$color))+
    labs(x="False Positive Rate",y="Sensitivity")+
    facet_wrap(vars(`Demographic Variables Included`))+
  theme(legend.position = "top")
jpeg(paste0(output_path,"roc_curvers_final",".jpg"),width = 5,height = 3,
         units = "in",res=350)
print(p)
graphics.off()

## Table 4

acc_metrics<-acc_metrics %>% filter(Balanced.Accuracy==`Max Balanced Accuracy`)
acc_metrics<-acc_metrics %>% select(mirna,`Demographic Variables Included`,sigma,C,Sensitivities,Specificities,Balanced.Accuracy,`Pos Pred Value`,
                       `Neg Pred Value`,AUC)

acc_metrics<-acc_metrics %>% mutate_if(is.numeric,function(x) round(x,2))
acc_metrics$mirna<-gsub("[+]",", ",acc_metrics$mirna)
acc_metrics$mirna<-gsub("[.]","-",acc_metrics$mirna)
write.csv(acc_metrics,paste0(output_path,"model_metrics.csv"),row.names = F)

## Visualize miRNAs

mirs<-data.frame(miRNA=unlist(strsplit(acc_metrics$mirna,"[+]")))
mirs<-mirs %>% group_by(miRNA) %>% count()

pdf(paste0(output_path,"violin_plots.pdf"),
    width=5,height=5)
for(g in mirs$miRNA){
  sub<-df[,c(g,"class")]
  colnames(sub)<-c("rna","condition")
  print(ggviolin(sub, x="condition", y="rna",
                 xlab = "",
                 ylab="Normalized Expression",
                 title=gsub("[.]","-",g),
                 add = "dotplot", shape="condition",
                 draw_quantiles = 0.5,
                 color="black",
                 fill="condition",
                 add.params = list(fill = "black")) +
          yscale("log10", .format=T) +
          theme_bw() +
          theme(axis.title.y = element_text(size = 25), 
                axis.text = element_text(size = 25),
                plot.title = element_text(hjust = 0.5, size=30),
                axis.title.x = element_text(size = 30)) +
          rremove("legend"))
}
graphics.off()

