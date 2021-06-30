output_path<-"data_peerJv5/output/"
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

plot_gini<-function(training_fit,label_m,sample_type){
  imp<-data.frame(training_fit$finalModel$importance)
  imp_sd<-data.frame(training_fit$finalModel$importanceSD)
  imp<-merge(imp,imp_sd,by=0)
  imp$sd1<-imp$MeanDecreaseAccuracy.x-imp$MeanDecreaseAccuracy.y
  imp$sd2<-imp$MeanDecreaseAccuracy.x+imp$MeanDecreaseAccuracy.y
  imp$mirna<-gsub("[.]","-",imp$Row.names)
  sel_color<-"#395C6B"
  write.csv(imp,paste0(output_path,"gini_scores_",label_m,"_",sample_type,".csv"),row.names = F)
  if(nrow(imp)>30){
    imp<-imp[1:25,]
  }
  p<-ggplot(data=imp,aes(y=reorder(mirna,-MeanDecreaseAccuracy.x)))+
    geom_segment(aes(x=sd1,xend=sd2,yend=mirna),col=sel_color,size=.8)+
    geom_point(aes(x=MeanDecreaseAccuracy.x),col=sel_color,shape=15,size=2.2)+theme_bw()+
    labs(y="miRNA",x="Mean Decrease Accuracy",title="A)")
  
  p2<-ggplot(data=imp,aes(x=MeanDecreaseGini,y=reorder(mirna,-MeanDecreaseGini)
  ))+
    geom_col(fill=sel_color)+theme_bw()+labs(y="miRNA",x="Mean Decrease Gini",
                                             title="B)")+
    scale_fill_viridis_d()
  
  jpeg(paste0(output_path,"rf_",label_m,"_",sample_type,".jpg"),width=6,height=ifelse(nrow(imp)>10,6,2),units="in",res=300)
  gridExtra::grid.arrange(p,p2,ncol=2)
  graphics.off()
}

auc_calc<-function(temp_df,training_fit,subset_label=NA){
  pred <- predict(training_fit,temp_df, type = "prob")
  roc.results <- roc(response = temp_df$class, predictor = pred$HCC)
  auc <- auc(roc.results)
  #print(auc)
  roc_df<-data.frame(Thresholds=roc.results$thresholds,
                     Sensitivities=roc.results$sensitivities,
                     Specificities=roc.results$specificities,
                     `Balanced Accuracy`=
                       (roc.results$sensitivities+roc.results$specificities)/2)
  roc_df$AUC<-auc
  roc_df$`Max Balanced Accuracy`<-max(roc_df$Balanced.Accuracy)
  roc_df$subset_label<-subset_label
  return(roc_df)
}



calculate_acc_metrics<-function(temp_df,training_fit,subset_label=NA,part=""){
  pred <- predict(training_fit,temp_df, type = "prob")
  roc.results <- roc(response = temp_df$class, predictor = pred$HCC)
  roc_df<-data.frame(Thresholds=roc.results$thresholds,
                     Sensitivities=roc.results$sensitivities,
                     Specificities=roc.results$specificities,
                     `Balanced Accuracy`=
                       (roc.results$sensitivities+roc.results$specificities)/2)
  roc_df$`Max Balanced Accuracy`<-max(roc_df$Balanced.Accuracy)
  thr<-roc_df$Thresholds[roc_df$`Max Balanced Accuracy`==roc_df$Balanced.Accuracy]
  train_thr<-ifelse(pred$HCC<thr,"HCC","Cirrhosis") %>% factor()
  temp_cm<-confusionMatrix(train_thr,temp_df$class,"HCC") 
  if(temp_cm$byClass["Balanced Accuracy"]<max(roc_df$Balanced.Accuracy)){
    train_thr<-ifelse(pred$HCC>thr,"HCC","Cirrhosis") %>% factor()
    temp_cm<-confusionMatrix(train_thr,temp_df$class,"HCC")  
  }
  comb<-temp_cm$byClass %>% data.frame()
  comb$Partition<-part
  comb$subset_label<-subset_label
  comb$Threshold<-thr
  return(comb)
}


mirna.model <- function(train_df,test_df,coldata,model.type="rf",
                        sample_type=""){
  cat("Training Model...")
  ctrl <- trainControl( savePredictions = TRUE, 
                        classProbs = TRUE, 
                        summaryFunction = twoClassSummary,
                        search = "grid")
  subset_list<-colnames(train_df)[-ncol(train_df)]
  track_roc<-NULL
  while(length(subset_list)!=0){
    formula=as.formula(paste0("class~",paste0(subset_list, 
                                              collapse = "+"))) 
    training_fit <- train(formula, data=train_df, method=model.type, trControl=ctrl, na.action =     na.pass,importance = TRUE,metric="ROC",resamples="all",tuneLength=10)
    imp<-training_fit$finalModel$importance
    x<-data.frame(num=length(subset_list),mean=mean(training_fit$resample$ROC,
                                                    na.rm = T),
                  sd=sd(training_fit$resample$ROC,na.rm = T))
    x$vars<-list(subset_list)
    train_pred<-auc_calc(train_df,training_fit,"Training")
    x$Train<-unique(train_pred$AUC)
    test_pred<-auc_calc(test_df,training_fit,"Test")
    x$Test<-unique(test_pred$AUC)
    track_roc<-rbind.data.frame(track_roc,x)
    subset_list<-row.names(imp)[imp[,4]>quantile(imp[,4],.05)]
  }
  sel_mod<-track_roc[15:1,] %>% select(mean) %>% unlist() %>% which.max()
  sel_mod<-nrow(track_roc)-sel_mod
  sel_mod<-track_roc$mean[sel_mod]
  track_roc<-track_roc %>% mutate(max_test=ifelse(Test==max(Test,na.rm = T),
                                                  paste0(num,"miRNAs: ",Test),""))
  track_sub<-track_roc %>% filter(max_test!="") %>%
    filter(num==min(num))
  print(track_sub)
  track_models<-NULL
  track_acc<-NULL
  for(m1 in seq(1,nrow(track_sub))){
    subset_list<-unlist(track_sub[m1,"vars"])
    formula=as.formula(paste0("class~",paste0(subset_list, 
                                              collapse = "+"))) 
    training_fit <- train(formula, data=train_df, method=model.type, trControl=ctrl, na.action =na.pass,
                          importance = TRUE,metric="ROC",resamples="all",tuneLength=10)
    cat("Done.\n")
    test_pred<-auc_calc(test_df,training_fit,track_sub$max_test[m1])
    track_models<-rbind.data.frame(track_models,test_pred)
    temp_acc<-calculate_acc_metrics(train_df,training_fit,track_sub$max_test[m1],"Train")
    track_acc<-rbind.data.frame(track_acc,temp_acc)
    temp_acc<-calculate_acc_metrics(test_df,training_fit,track_sub$max_test[m1],"Test")
    track_acc<-rbind.data.frame(track_acc,temp_acc)
    plot_gini(training_fit,track_sub$num[m1],sample_type)
    for(m2 in subset_list){
    formula=as.formula(paste0("class~",paste0(m2, 
                                                collapse = "+"))) 
    training_fit2 <- train(formula, data=train_df, method=model.type, trControl=ctrl, na.action =na.pass,
                            importance = TRUE,metric="ROC",resamples="all",tuneLength=10)
    test_pred<-auc_calc(test_df,training_fit2,m2)
    track_models<-rbind.data.frame(track_models,test_pred)
    temp_acc<-calculate_acc_metrics(train_df,training_fit2,m2,"Train")
    track_acc<-rbind.data.frame(track_acc,temp_acc)
    temp_acc<-calculate_acc_metrics(test_df,training_fit2,m2,"Test")
    track_acc<-rbind.data.frame(track_acc,temp_acc)
    }
  }
  track_acc$sample_type<-sample_type
  track_roc<-track_roc %>% 
    mutate(model_select=ifelse(num==track_sub$num,"Selected",""))
  p<-ggplot(data=track_roc %>%filter(num<70),aes(x=-num,y=mean))+
    geom_segment(aes(xend=-num,yend=mean+sd,y=mean-sd))+
    geom_point(aes(col=model_select))+
    theme(legend.title = element_blank())+
    scale_color_manual(values=c("black","goldenrod"))+
    labs(y="AUC",x="miRNA N")+scale_x_continuous(labels = function(x) abs(x))+
    theme(legend.position = "top")
  jpeg(paste0(output_path,"auc_by_features_",sample_type,".jpg"),width = 3.5,
       height = 2.5,units = "in",res=350)  
  print(p)
  graphics.off()
  
  track_roc$Diff<-track_roc$Train-track_roc$Test
  
  p<-ggplot(data=track_roc %>% filter(num<70),aes(y=Diff,x=-num))+
    geom_point(aes(col=model_select),size=3)+
    scale_color_manual(values=c("black","goldenrod"))+
    geom_smooth(se=F,col="darkgrey")+
    theme(legend.position = "right")+
    labs(y="Overfitting",x="miRNA N",
         color="Models")+
    scale_x_continuous(labels = function(x) abs(x))

  jpeg(paste0(output_path,"overfitting_by_features_",sample_type,".jpg"),width = 4.5,
       height = 2.5,units = "in",res=350)
  print(p)
  graphics.off()
    
  p<-ggplot(data=track_roc %>% filter(num<70),aes(y=Diff,x=-num))+
    geom_line(col="darkgrey")+
    geom_point(aes(col=ifelse(model_select=="Selected","Selected Model","")),size=3)+
    scale_color_manual(values=c("black","goldenrod"))+
    theme(legend.position = "top",legend.title = element_blank())+
    labs(y="Model Fit",x="miRNA N")+
    scale_x_continuous(labels = function(x) abs(x))
  
  jpeg(paste0(output_path,"overfitting_v2_by_features_",sample_type,".jpg"),width = 5,
       height = 3.5,units = "in",res=350)  
  print(p)
  graphics.off()
  
  write_rds(track_roc,paste0(output_path,"overfitting_by_features_",sample_type,".rds"))
  
  track_models$subset_label[grep("RNA",track_models$subset_label)]<-"Test"
  track_models<-rbind.data.frame(auc_calc(train_df,training_fit,"Training"),track_models)
  track_models<-track_models %>% mutate(fpr=1-Specificities) %>%
    arrange(Sensitivities)
  track_models$subset_label<-paste(track_models$subset_label," AUC: ",round(track_models$AUC,2),sep="")
  
  p<-ggplot(data=track_models %>% filter(grepl("T",subset_label)),
            aes(y=Sensitivities,x=fpr,group=subset_label,
                                  col=subset_label))+
    geom_point(size=2)+
    geom_line(size=1.2)+
    theme(legend.title = element_blank(),legend.position = "none")+
    labs(x="False Positive Rate",y="Sensitivity")+
    scale_color_manual(values = c("royalblue","goldenrod"))+
    geom_segment(aes(x=0,xend=1,y=0,yend=1),lty=2,color="black")+
    facet_wrap(~subset_label)
  jpeg(paste0(output_path,"model_roc_curves_",sample_type,".jpg"),width = 3.5,
       height = 2.5,units = "in",res=350)
  print(p)
  graphics.off()
  
  p<-ggplot(data=track_models %>% filter(!grepl("T",subset_label)),
            aes(y=Sensitivities,x=fpr,group=subset_label))+
    geom_point(size=2)+
    geom_line(size=1.2)+
    theme(legend.title = element_blank(),legend.position = "none")+
    labs(x="False Positive Rate",y="Sensitivity")+
    geom_segment(aes(x=0,xend=1,y=0,yend=1),lty=2,color="black")+
    facet_wrap(~subset_label)+
    geom_label(x=.75,y=.25,aes(label=paste("AUC:",round(AUC,2))))
  jpeg(paste0(output_path,"model_single_curves_",sample_type,".jpg"),width = 5.5,
       height = 2.5,units = "in",res=350)
  print(p)
  graphics.off()
  track_models$sample_type<-sample_type
  write.csv(track_models,paste0(output_path,"model_roc_",sample_type,".csv"),row.names = F)
  return( list(track_models,track_acc))
}
