output_path<-"dir/output/"
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

auc_calc<-function(pred){
  roc.results <- roc(response = pred$obs, predictor = pred$HCC)
  auc <- auc(roc.results)
  roc_df<-data.frame(Thresholds=roc.results$thresholds,
                     Sensitivities=roc.results$sensitivities,
                     Specificities=roc.results$specificities,
                     `Balanced Accuracy`=
                       (roc.results$sensitivities+roc.results$specificities)/2)
  roc_df$AUC<-auc
  roc_df$`Max Balanced Accuracy`<-max(roc_df$Balanced.Accuracy)
  return(roc_df)
}

mirna.model <- function(train_df,subset_label="",coldata){
  subset_list<-colnames(train_df)[-ncol(train_df)]
  formula=as.formula(paste0("class~",paste0(subset_list, 
                                              collapse = "+"))) 
  x<-train_df %>% select(-class) %>% as.matrix()
  ctrl <- rfeControl(functions = rfFuncs,
                       method = "LOOCV",
                       verbose = FALSE)
  ctrl$returnResamp<-"all"
  svmProfile <- rfe(x, train_df$class,
                    sizes = c(2,5,7,10),
                    rfeControl =ctrl,preProcess=c("scale","center","nzv"))
  fitControl<-trainControl(method = "LOOCV",
                           classProbs=T,savePredictions = T,
                           summaryFunction=twoClassSummary,
                           seeds = )
  fitControl$returnResamp<-"all"
  auc_df_fin<-NULL
  for(add_demo in c(F,T)){
  formula<-as.formula(paste0("class~",paste0(predictors(svmProfile)[1:10], 
                                               collapse = "+")))
  training_fit<-train(formula, data=train_df, method="svmRadial",
                        preProcess=c("scale","center","nzv"),metric="ROC",trControl=fitControl)
  
  if(add_demo==T){
    temp_df<-train_df
    temp_df$DEIDENT<-row.names(temp_df)
    temp_df<-merge(coldata %>% select(DEIDENT,RaceDescription,GenderDescription,smoking_status,
                                      BMI,Age),temp_df,by="DEIDENT")
    temp_df<-temp_df %>% select(-DEIDENT)
    formula<-as.formula(paste0("class~",paste0(c(predictors(svmProfile)[1:10],
                                                 "RaceDescription","GenderDescription","smoking_status",
                                                 "BMI","Age"), 
                                               collapse = "+")))
    training_fit<-train(formula, data=temp_df, method="svmRadial",
                        preProcess=c("scale","center","nzv"),metric="ROC",trControl=fitControl)
  }
  
  predictions<-training_fit$pred %>% filter(sigma==training_fit$bestTune$sigma&C==training_fit$bestTune$C)  
  cm<-confusionMatrix(predictions$pred,predictions$obs,"HCC")
  print(cm$byClass)
  auc_df<-auc_calc(predictions)
  auc_df$label=subset_label
  auc_df$mirna<-paste0(predictors(svmProfile)[1:10], 
                       collapse = "+")
  auc_df$`Pos Pred Value`<-cm$byClass["Pos Pred Value"]
  auc_df$`Neg Pred Value`<-cm$byClass["Neg Pred Value"]
  auc_df$sigma<-training_fit$bestTune$sigma
  auc_df$C<-training_fit$bestTune$C
  auc_df$`Demographic Variables Included`<-add_demo
  auc_df_fin<-rbind.data.frame(auc_df_fin,auc_df)
  saveRDS(training_fit,file.path(output_path,
                                 paste0("model_",subset_label,"_",add_demo,".RDS")))
  }
  return(auc_df_fin)
  }
