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

theme_set(theme_bw())
set.seed(2021)

############################################
######## Data Ingest ######################
############################################

source("scripts_peerJ/master_functions.R")

dir_path<-"data_peerJ/"
output_path<-"data_peerJv5/output/"
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

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=coldata, 
                              design=~class)
dds<- DESeq(dds)

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
resLFC_clean <- lfcShrink(ddsClean, coef="class_HCC_vs_Cirrhosis", type="apeglm")
cleanresOrdered <- resLFC_clean[order(resLFC_clean$pvalue),]
write.csv(as.data.frame(cleanresOrdered),
          paste0(output_path,"DESeq_HCCvCirrhosis_cleanres.csv"))
pdf(paste0(output_path,"DESeqPlots.pdf"))
plotMA(resLFC_clean)
graphics.off()

p1 <- EnhancedVolcano(resLFC_clean,
                      lab = rownames(resLFC_clean),
                      x = "log2FoldChange",
                      y = "pvalue",
                      pCutoff = 5e-6,
                      FCcutoff = 5,
                      xlim = c(-5.5, 5.5),
                      ylim = c(0, -log10(10e-12)),
                      pointSize = 2,
                      labSize = 3.0,
                      title="",
                      subtitle = "",
                      caption="",
                      legendPosition = "right",
                      legendLabSize = 8,
                      colAlpha = 0.85,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5)
p1
pdf(paste0(output_path,"fig1.pdf"),width=8,height=6)
p1
graphics.off()

######################################################
########## Heatmaps and Clustering ###################
######################################################

### Setup annotations

ramp <- colorRampPalette(colors=c("green", "darkgreen","grey","darkred","red"))(20)
df<-coldata %>% select(DEIDENT,class,cirrhosis.vs.fibrosis) %>%
  mutate(cirrhosis.vs.fibrosis=ifelse(!is.na(cirrhosis.vs.fibrosis),
                                      cirrhosis.vs.fibrosis,"none documented"))
row.names(df)<-df$DEIDENT
df<-df[,-1]
colnames(df)<-c("class","Liver Disease")
var1<-c("#2E86AB","#E88873")
names(var1)<-c("Control","HCC")
var2<-c("#3EC300","#FFCB47","#796293")
names(var2)<-c("cirrhosis","fibrosis","none documented")
ann_colors<-list(Condition=var1,`Liver Disease`=var2)

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
p_vals<-0.001
fold_change<-2
clust_num<-4
top.genes <- which(abs(resLFC_clean$log2FoldChange)>fold_change & resLFC_clean$padj<p_vals)
topgenes <- rownames(unique(resLFC_clean[top.genes,]))
out<-pheatmap(assay_dds_clean[top.genes,], cluster_rows=T, 
              color = ramp,
              show_rownames=F,
              show_colnames = F,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = clust_method,
              clustering_method = "ward.D",
              cutree_cols = clust_num,
              cutree_rows = 3,
              scale = "row",
              treeheight_col = 50,
              cluster_cols=T, annotation_col=df,main=paste("P-val<",p_vals,"fold_changes>",fold_change,"k=",clust_num),
              annotation_colors = ann_colors)
save_pheatmap_pdf(out,paste(output_path,clust_method,"_ann.pdf",sep=""),
                  height = 5)

############################################
######### Modeling #########################
############################################

select.genes <- rownames(resLFC_clean[which(resLFC_clean$padj<0.05),])
df <- ddsClean[select.genes,]
df <- counts(df, normalized=T) %>% t() %>% data.frame()
df$class<-dds$class
trainIndex<-which(key_df$Test==F)
train_df<-df[trainIndex,]
test_df<-df[-trainIndex,]


train_df %>% select(class) %>% table()
test_df %>% select(class) %>% table()


### Model 
track_models<-mirna.model(train_df,test_df,coldata,model.type="rf")

## Modeling liver disease samples

samples_with_cirrhosis<-coldata$DEIDENT[!is.na(coldata$cirrhosis.vs.fibrosis)]
train_cirr<-train_df[row.names(train_df) %in% samples_with_cirrhosis,]
test_cirr<-test_df[row.names(test_df) %in% samples_with_cirrhosis,]

track_models_2<-mirna.model(train_cirr,test_cirr,coldata,model.type="rf",sample_type = "cld")

acc_metrics<-rbind.data.frame(track_models[[2]],track_models_2[[2]])
track_models<-rbind.data.frame(track_models[[1]],track_models_2[[1]])

## Visualize 

mirs<-track_models %>% filter(!grepl("T",subset_label)) %>% group_by(sample_type) %>%
  mutate(max_auc=max(AUC)) %>%filter(AUC==max_auc) %>% select(subset_label) %>% unique()


roc_df<-track_models %>% filter(grepl("Test",subset_label) | subset_label %in% mirs$subset_label)

roc_df<-roc_df %>% mutate(color=ifelse(sample_type=="","goldenrod","royalblue"))
for(i in unique(roc_df$sample_type)){
  for(j in unique(roc_df$subset_label[roc_df$sample_type==i])){
    temp_df<-roc_df[roc_df$subset_label==j&roc_df$sample_type==i,] %>%
      arrange(Sensitivities)
    p<-ggplot(data=temp_df,aes(x=1-Specificities,y=Sensitivities))+
      geom_segment(x=0,xend=1,y=0,yend=1,lty=2)+
      geom_line(size=1.2,color=unique(temp_df$color))+
      geom_point(size=2,color=unique(temp_df$color))+labs(x="False Positive Rate",y="Sensitivity",
                                                          title = j)
    jpeg(paste0(output_path,i,"_",j,".jpg"),width = 3,height = 2,
         units = "in",res=350)
    print(p)
    graphics.off()
  }
}

## Table 4

acc_metrics %>% filter(subset_label %in% sapply(strsplit(mirs$subset_label," "),"[",1)) %>%
  filter(Partition=="Test")

acc_metrics$metric<-gsub("[0-9]","",row.names(acc_metrics))

spec_metrics<-c("Sensitivity","Specificity","Balanced Accuracy","Pos Pred Value","Neg Pred Value")

acc_metrics$subset_label<-
  sapply(strsplit(acc_metrics$subset_label,":"),"[",1)

unique_mirnas<-acc_metrics %>% select(miRNAs=subset_label,sample_type) %>% unique()

for(i in spec_metrics){
  t1<-acc_metrics[acc_metrics$metric==i&acc_metrics$Partition=="Test",] %>%
    select("subset_label",".","sample_type")
  colnames(t1)[2]<-i
  unique_mirnas<-merge(unique_mirnas,t1,by.x=c("miRNAs","sample_type"),
                       by.y=c("subset_label","sample_type"))
}

unique_mirnas %>% arrange(sample_type)

write.csv(unique_mirnas,paste0(output_path,"model_metrics.csv"),row.names = F)

unique_mirnas %>% filter(!grepl("miRNAs",miRNAs)) %>%
  group_by(sample_type) %>% summarise_if(is.numeric,function(x) paste(round(min(x),2),
                                                                      round(max(x),2),sep="-"))


## Visualize miRNAs

mirs<-track_models$subset_label[grep("hsa",track_models$subset_label)]
mirs<-sapply(strsplit(mirs," AUC:"), "[",1) %>% unique()

pdf(paste0(output_path,"violin_plots.pdf"),
    width=5,height=5)
for(g in mirs){
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

