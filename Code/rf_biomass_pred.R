
library(dplyr)
library(caret)
##----------------------------------------------------##
##------ Try with all nonNA metabs first  ------------##
##----------------------------------------------------##
#full_df <- read.csv("Data/Combined_noNA_pheno_geno_data.csv")
full_df <- read.csv("Data/Combined_noNA_pheno_geno_data.csv")
noise_data <- FALSE
roots_only <- FALSE
leaves_only <- FALSE
accession_only <- TRUE
rm_accession <- FALSE
#Options: AG, BG, Total
response <- "Total"

#Define the response chosen ("AG","BG","Total")
if(response=="AG"){
  cols_to_rm <- c("ID","BG_weight","Total_weight")
  full_df <- full_df[,-which(colnames(full_df)%in%cols_to_rm)]
}else if(response =="BG"){
  cols_to_rm <- c("ID","AG_weight","Total_weight")
  full_df <- full_df[,-which(colnames(full_df)%in%cols_to_rm)]
}else{
  #Total
  cols_to_rm <- c("ID","BG_weight","AG_weight")
  full_df <- full_df[,-which(colnames(full_df)%in%cols_to_rm)]
}
colnames(full_df)[3] <- "Biomass"

#Filter down to roots/leaves only as appropriate
if(roots_only){
  full_df <- full_df[,c(1:3,grep("Roots",colnames(full_df)))]
} 

if(leaves_only){
  full_df <- full_df[,c(1:3,grep("Above",colnames(full_df)))]
}

if(accession_only){
  full_df <- full_df[,c(1:3)]
}else if(rm_accession){
  #metab_mat <- data.matrix(full_df[,-c(1:3)])
  lmres <- lm(rnorm(nrow(full_df))~Biomass-1,data=full_df)
  full_df$Biomass <- lmres$residuals
}

#Separate control and drought data

control_df <- filter(full_df,Treatment=="Control")%>%select(-Treatment)
drought_df <- filter(full_df,Treatment=="Drought")%>%select(-Treatment)

#Try the "noise data" idea (scramble the responses so the metabolites and biomass don't map)
if(noise_data){
  set.seed(10)
  control_df$Biomass <- control_df$Biomass[sample(nrow(control_df),replace = FALSE)]
  drought_df$Biomass <- drought_df$Biomass[sample(nrow(drought_df),replace = FALSE)]
}

##----------------------------------------------------##
##------ Create control/drought data sets  -----------##
##----------------------------------------------------##

#RF set up used by both control and testing models
control <- trainControl(method="repeatedcv", number=5, repeats=5)
metric <- "RMSE"
set.seed(9832)
if(accession_only){
  mtry <- seq(5,10,length=3)
}else{
  mtry <- seq(100,300,length=5)
}
tunegrid <- expand.grid(.mtry=mtry)

##-----------------------##
#--- Control analysis (AG) ----#

#Split to testing/training
train_split <- createDataPartition(y=control_df$Biomass, p=.75,list=FALSE)
training_con <- control_df#[train_split,]
testing_con <- control_df#[-train_split,]

# Create model with default paramters
#crf_default_con <- train(Biomass~., data=training_con, method="cforest", metric=metric,trControl=control,  preProcess=c("center"))
rf_default_con <- train(Biomass~., data=training_con, method="ranger", metric=metric, tuneGrid=tunegrid, 
                        trControl=control,preProcess=c("center"),importance='impurity')

#Add importance=TRUE to get importance measures for random forests

plot(rf_default_con$resample$Rsquared,main=mean(rf_default_con$resample$Rsquared))
plot(rf_default_con$resample$RMSE,main=mean(rf_default_con$resample$RMSE))
mean(rf_default_con$resample$Rsquared)
mean(rf_default_con$resample$RMSE)

rf_preds_con <- predict(rf_default_con, newdata = testing_con)
plot(rf_preds_con,testing_con$Biomass);abline(0,1)
SSReg_con <- sum((rf_preds_con-testing_con$Biomass)^2)
1-SSReg_con/(sum((testing_con$Biomass-mean(testing_con$Biomass))^2))
sqrt(mean(unname(rf_preds_con-testing_con$Biomass)^2)) #RMSE

metabs_con <- rf_default_con$finalModel$variable.importance
plot(metabs_con);abline(v=c(29.5,735.5),col='gray50')


##-----------------------##
#--- Drought analysis (AG) ----#

#Split to testing/training
train_split <- createDataPartition(y=drought_df$Biomass, p=.75,list=FALSE)
training_dr <- drought_df#[train_split,]
testing_dr <- drought_df#[-train_split,]

# Create model with default paramters
rf_default_dr <- train(Biomass~., data=training_dr, method="ranger", metric=metric, tuneGrid=tunegrid, 
                        trControl=control,preProcess=c("center"),importance='impurity')


plot(rf_default_dr$resample$Rsquared,main=mean(rf_default_dr$resample$Rsquared))
plot(rf_default_dr$resample$RMSE,main=mean(rf_default_dr$resample$RMSE))
mean(rf_default_dr$resample$Rsquared)
mean(rf_default_dr$resample$RMSE)


rf_preds_dr <- predict(rf_default_dr, newdata = testing_dr)
plot(rf_preds_dr,testing_dr$Biomass);abline(0,1)
SSReg_dr <- sum((rf_preds_dr-testing_dr$Biomass)^2)
1-SSReg_dr/(sum((testing_dr$Biomass-mean(testing_dr$Biomass))^2))
sqrt(mean(unname(rf_preds_dr-testing_dr$Biomass)^2)) #RMSE


#metabs_dr <- rf_default_dr$finalModel$importance[,1]
#plot(metabs_dr);abline(v=c(29.5,735.5),col='gray50')
metabs_dr <- rf_default_dr$finalModel$variable.importance
plot(metabs_dr);abline(v=c(29.5,735.5),col='gray50')

##------------------------------------------------##
##--- Compare control and drought "signatures ----##

plot(metabs_con,metabs_dr,xlab="Control",ylab="Drought");abline(0,1)

#Better plots
all_imports <- data.frame(Metab=rf_default_con$coefnames,Idx=1:length(metabs_con),Control=metabs_con, Drought=metabs_dr)
all_imports$Tissue <- "Accession"
all_imports$Tissue[grep("Above", all_imports$Metab)] <- "Leaf Metabolites"
all_imports$Tissue[grep("Roots", all_imports$Metab)] <- "Root Metabolites"


p1 <- ggplot(data=all_imports)+geom_vline(xintercept=c(29.5,735.5),colour='gray70')+geom_point(aes(Idx,Control,colour=Tissue))+xlab("")+
  ylab("Mean Drop in RMSE")+annotate("text",x=600,y=0.02,label="Leaf\nMetabolites")+annotate("text",x=1000,y=0.02,label="Root\nMetabolites")+
  scale_colour_manual(values=c('black','orange','blue'))+guides(colour=FALSE)+theme_bw()+ggtitle("A - Control")
p1

p2 <- ggplot(data=all_imports)+geom_vline(xintercept=c(29.5,735.5),colour='gray70')+geom_point(aes(Idx,Drought,colour=Tissue))+xlab("")+
  ylab("Mean Drop in RMSE")+annotate("text",x=500,y=0.02,label="Leaf\nMetabolites")+annotate("text",x=1000,y=0.02,label="Root\nMetabolites")+
  scale_colour_manual(values=c('black','orange','blue'))+guides(colour=FALSE)+theme_bw()+ggtitle("B - Drought")
p2

p3 <- ggplot(data=all_imports)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue),alpha=.6)+
  coord_fixed()+xlim(c(0,.026))+ylim(c(0,0.026))+scale_colour_manual(values=c('black','orange','blue'))+guides(colour=FALSE)+theme_bw()+
  ggtitle("C")
p3

layout_mat <- rbind(c(1,1,1,1,3,3),
                    c(2,2,2,2,3,3))

gridExtra::grid.arrange(p1,p2,p3,layout_matrix=layout_mat)
#ggsave(filename = paste0("RF_results_",response,".pdf"),plot = gridExtra::grid.arrange(p1,p2,p3,layout_matrix=layout_mat))

#Save "all_imports" object
all_imports$Tissue <- response
save(all_imports, file = paste0("Results/",response,"_RF_Results.RData"))


all_imports$ConRank <- rank(-all_imports$Control)
all_imports$DrRank <- rank(-all_imports$Drought)

#Filter to top 10 control and rank metabolites
most_imports <- all_imports[unique(c(which(all_imports$ConRank<=20), which(all_imports$DrRank<=20))),]
most_imports <- most_imports%>%select(-Control,-Drought)%>%arrange(ConRank)
bottom_rank <- c(1:20,20+order(most_imports$DrRank[-c(1:20)]))

#write.csv(most_imports[bottom_rank,],"Results/Ranked_Metabs.csv",row.names = FALSE)








