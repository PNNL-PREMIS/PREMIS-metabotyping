library(ggplot2)
library(dplyr)
theme_hm <- function(col=NULL){
  return(theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    axis.line=element_blank(),
    axis.ticks=element_blank(),
    #axis.text.y=element_blank(),
    axis.text.x=element_text(colour='black',face='bold',angle=45,hjust=.8,siz=9),
    panel.background=element_rect(fill=col),
    #legend.position = "none",
    axis.title=element_blank(),
    panel.border=element_blank(),
    plot.background=element_blank()
  ))
}

#Read in metab and phenotype data
geno_data <- data.table::fread("Data/GoBrachy_LCMS_GCMS_DATASET_July_27_2017.csv",data.table=FALSE)
#pheno_data <- data.table::fread("Data/Cleand_GoBrachy_Data.csv",data.table=FALSE)

#Reexpress the phenotype data so it's in a similar format to geno_data
#pheno_data_control <- pheno_data[,c(1:5)]
#pheno_data_drought <- pheno_data[,c(1,2,6:8)]
#colnames(pheno_data_drought) <- colnames(pheno_data_control) <- c("Plant_Line","Replicate","AG_weight","BG_weight","Total_weight")

#pheno_data_control$Treatment <- "Control"
#pheno_data_drought$Treatment <- "Drought"
#all_pheno <- rbind(pheno_data_control,pheno_data_drought)
#all_pheno <- all_pheno[,c(1,2,6,3:5)]

#sum_pheno <- all_pheno%>%group_by(Plant_Line,Treatment)%>%summarize(AG=mean(AG_weight),BG=mean(BG_weight),Total=mean(Total_weight))
#plot(sum_pheno[sum_pheno$Treatment=="Control",]$AG,sum_pheno[sum_pheno$Treatment!="Control",]$AG)
#plot(sum_pheno[sum_pheno$Treatment=="Control",]$BG,sum_pheno[sum_pheno$Treatment!="Control",]$BG,pch=19);abline(0,1)


#Give each "Unknown" a unique columne name
colnames(geno_data)[grep("Unknown",colnames(geno_data))] <- paste("Unknown",1:length(grep("Unknown",colnames(geno_data))),sep="_")

#Select "Roots" or "Above" (aka leaves)
org <- "Above"
geno_data <- dplyr::filter(geno_data,Organ==org)
if(org=="Above"){
  org <- "Leaves"
}

#Need to replace "-" with "_" because dashes cause problems
geno_data$Line <- gsub(pattern = "-",replacement = "_", geno_data$Line)

##------------------------------------------##
#-- Do pMart type analysis of geno_data ----##

library(MSomicsQC)
#The data need to be transposed for MSomics family
metabs <- colnames(geno_data)[-c(1:5)]
geno_data$ID2 <- 1:nrow(geno_data)
samp_names <- paste("ID",geno_data[,"ID2"],geno_data[,"Line"],sep="_")
raw_metab_mat <- unname(data.matrix(geno_data[,-c(1:5)]))
raw_metab_mat <- data.frame(t(raw_metab_mat))
raw_metab_mat <- cbind(Metabolites = metabs,raw_metab_mat)
colnames(raw_metab_mat)[-1] <- samp_names

#Create the f-data object
gob_fdata <- geno_data[,c(1:5)]
gob_fdata[,1] <- samp_names
colnames(gob_fdata)[1] <- "SampleID"
gob_fdata$gp <- paste(gob_fdata$Treatment,gob_fdata$Line,sep="_")

#Put into metab
gob_data <- as.metabData(e_data=raw_metab_mat, f_data = gob_fdata,edata_cname = "Metabolites",fdata_cname = "SampleID", data_scale='abundance')

#------ Log base 2 transform the data ----#
#For some reason, as.metabData automatically days edata is on log2 scale, but it's not
gob_data_trans <- edata_transform(omicsData =  gob_data, data_scale = "log2")

#------ Group the data ----#
gob_data_trans <- group_designation(omicsData = gob_data_trans, main_effects = c("gp"), time_course=NULL)

#------- Filter out biomolecules using IMD-ANOVA filter -----#
imd_afilt <- imdanova_filter(omicsData = gob_data_trans)
gob_data_trans <- applyFilt(filter_object = imd_afilt, omicsData = gob_data_trans, min_nonmiss_anova = 2,min_nonmiss_gtest = 3)

##---------------------##
#- Stats analysis -----##
library(MSomicsSTAT)
library(dplyr)
library(scales)
comp_df <- data.frame(Control=paste("W",unique(geno_data$Line),sep="_"),Test=paste("D",unique(geno_data$Line),sep="_"))
attr(gob_data_trans,"group_DF") <- attr(gob_data_trans,"group_DF")[,c(1:2)]

all(comp_df$Control%in%attr(gob_data_trans,"group_DF")$Group)

gob_res <- imd_anova(omicsData = gob_data_trans, comparisons = comp_df, test_method = 'combined', pval_thresh = 0.01, pval_adjust = "holm")

##---------------------##
#- Plot of counts -----##

comp_df_melt <- reshape2::melt(attr(gob_res,"number_significant"),id.vars="Comparison",value.name="Count",variable.name="Direction")
levels(comp_df_melt$Comparison) <- unique(gob_data_trans$f_data$Line)
levels(comp_df_melt$Comparison) <- gsub(pattern = "_", replacement = " ",levels(comp_df_melt$Comparison))
total <- comp_df_melt%>%group_by(Comparison)%>%summarize(Total=sum(abs(Count)))%>%arrange(desc(Total))
comp_df_melt$Comparison <- factor(comp_df_melt$Comparison,levels=total$Comparison)

#To turn the up direction to green, and down direction to red
#pal <- RColorBrewer::brewer.pal(9, "Set1")
#scale_fill_manual("",values=pal[c(3,1)],labels=c("More in Drought","Less in Drought"))+ #Added to plot below

##Up direction is positive, down direction is negative
comp_df_melt[comp_df_melt$Direction=="Down",]$Count <- (-comp_df_melt[comp_df_melt$Direction=="Down",]$Count)
ggplot(data=comp_df_melt,aes(Comparison,Count,fill=Direction))+geom_bar(stat='identity')+
  geom_hline(aes(yintercept=0),colour='gray50')+theme_bw()+
  theme(axis.text.x=element_text(angle=40,vjust=1,hjust=.9,size=7),axis.text=element_text(colour='black'))+
  scale_fill_manual("",values=c("#4DBBD5FF","#F39B7FFF"),labels=c("More in Drought","Less in Drought"))+
  xlab("")+ylab("Number of Statistically\nSignificant Metabolites")+ylim(c(-200,200))
fig_id <- ifelse(org=="Roots","B","A")
#ggsave(paste0("Extended_Data_Fig4",fig_id,".eps"),width=7,height=4.375)


##----------------------------------------------------------------##
##------------ Full heat map (control means)  --------------------##
##-----------  These are the heamaps published on the shiny app --##
##----------------------------------------------------------------##

line_names <- read.csv("Data/GoBrachy_30_lines.csv",header=FALSE)
metab_means <- gob_res$Full_results
fc_cols <- c(1,grep("Mean_W",colnames(metab_means)))
metab_means <- metab_means[,fc_cols]
colnames(metab_means)[-1] <- paste(line_names$V2,org)
#write.csv(metab_means,"~/Desktop/Leaf_Metabolites.csv")

#Filter to three interesting lines
library(reshape2)
library(dplyr)
colnames(metab_means)[-1] <- as.character(line_names$V2)
#metab_means <- metab_means[-grep("Unknown",metab_means$Metabolites),]
#metab_means$Metabolites <- factor(metab_means$Metabolites)

#Remove meatbolites with missing data
na_rows <- unname(which(apply(apply(metab_means[,-1],1,is.na),2,any)))
metab_means <- metab_means[-na_rows,]

#Cluster to organize rows
cl_metab_clusts <- hclust(dist(data.matrix(metab_means[,-c(1)])),method='ward.D')
matab_order <- metab_means$Metabolites[rev(cl_metab_clusts$order)]

#Cluster to organize the columns
cl_access_clusts <- hclust(dist(t(data.matrix(metab_means[,-c(1)]))),method='ward.D')
access_order <- colnames(metab_means[,-1])[cl_access_clusts$order]

metab_means_melt <- melt(metab_means,id.vars=c("Metabolites"),variable.name='Accession')
metab_means_melt$Accession <- factor(metab_means_melt$Accession,levels=access_order)
metab_means_melt$Metabolites <- factor(metab_means_melt$Metabolites,levels=matab_order)

#Try heatmaply
library(heatmaply)
rownames(metab_means) <- metab_means$Metabolites
tmp <- heatmaply(metab_means[,-1],scale_fill_gradient_fun=
            ggplot2::scale_fill_gradient2(midpoint=15,low="#003366",high="#990000",mid='white',limits=c(0,29)),
          fontsize_col=20,label_names=rev(c("Abundance","Accession","Metabolite")),hclust_method="ward.D",
          margins=c(200,20,20,20),showticklabels = c(TRUE,FALSE),file=paste0("control_heatmap_",org,".png"),
          row_dend_left=TRUE,height=3000,width=1200,subplot_heights=c(.1,.9))%>%layout(showlegend=FALSE)
rm(tmp)

#Standardized-version
metab_means_std <- t(apply(metab_means[,-1],1,function(x)return((x-min(x))/(max(x)-min(x)))))
tmp <- heatmaply(metab_means_std,scale_fill_gradient_fun=
                   ggplot2::scale_fill_gradient2(midpoint=.5,low="#003366",high="#990000",mid='white',limits=c(0,1)),
                 fontsize_col=20,label_names=rev(c("Abundance","Accession","Metabolite")),hclust_method="ward.D",
                 margins=c(200,20,20,20),showticklabels = c(TRUE,FALSE),file=paste0("control_heatmap_std_",org,".png"),
                 row_dend_left=TRUE,height=3000,width=1200,subplot_heights=c(.1,.9))%>%layout(showlegend=FALSE)
rm(tmp)



##------------------------------------------------##
##------------ interesting lines only
##------------------------------------------------##
metab_res <- gob_res$Full_results
fc_cols <- c(1,grep("Fold_change",colnames(metab_res)))
metab_res <- metab_res[,fc_cols]
colnames(metab_res)[-1] <- paste(line_names$V2,org)

#Filter to  interesting lines
library(reshape2)
library(dplyr)
colnames(metab_res)[-1] <- as.character(line_names$V2)
#int_lines <- c("Adi-2","Bd21","Arn1","BdTR1i","Mon3","Per-1")
int_lines <- colnames(metab_res)[-1]
metab_res <- metab_res[,c("Metabolites",int_lines)]
metab_res <- metab_res[-grep("Unknown",metab_res$Metabolites),]
metab_res$Metabolites <- factor(metab_res$Metabolites)

#Cluster to organize the rows
na_rows <- unname(which(apply(apply(metab_res[,-1],1,is.na),2,any)))
if(length(na_rows)>0){
  metab_res <- metab_res[-na_rows,]
}
cl_metab_clusts <- hclust(dist(data.matrix(metab_res[,-c(1)])),method='ward.D')
metab_res$Metabolites <- factor(metab_res$Metabolites,levels=metab_res$Metabolites[cl_metab_clusts$order])

#Denote stat significant ones?
metab_res_melt <- melt(metab_res,id.vars=c("Metabolites"),variable.name='Accession')
xlines_df <- data.frame(x=.5+0:length(unique(metab_res_melt$Accession)))
ylines_df <- data.frame(y=.5+0:length(unique(metab_res_melt$Metabolites)))
#Sig identifiers
sig_mat <- gob_res$Flags
sig_mat <- filter(sig_mat,Metabolites%in%metab_res_melt$Metabolites)
sig_mat <- melt(sig_mat,id.vars=c("Metabolites"),variable.name='Accession',value.name="Sig")
rows <- NULL
sig_mat$Accession <- factor(sig_mat$Accession,labels=int_lines)
metab_res_melt <- merge(metab_res_melt,sig_mat)

#Keep all metabs that are significant for atleast one line
metab_keeps <- filter(sig_mat,Sig==1)$Metabolites
metab_keeps <- factor(metab_keeps)
metab_res <- metab_res[metab_keeps,]


#Abbreviate the metabolite names (only take those before ";")
to_abbrev <- grep(";",levels(metab_res_melt$Metabolites))
levels(metab_res_melt$Metabolites)[to_abbrev] <- gsub(";.*$","",levels(metab_res_melt$Metabolites)[to_abbrev])
levels(metab_res_melt$Metabolites) <- tolower(levels(metab_res_melt$Metabolites))

#Try Albert's standardization idea
metab_res_std <- metab_res[,-1]
metab_res_std <- t(apply(metab_res_std,1,function(x)return(x/max(abs(x)))))
cl_metab_clusts <- hclust(dist(data.matrix(metab_res_std)),method='ward.D')

metab_res_std <- data.frame(Metabolites=metab_res$Metabolites,metab_res_std)
#metab_res_std$Metabolites <- factor(metab_res_std$Metabolites,levels=metab_res_std$Metabolites[cl_metab_clusts$order])

metab_res_std_melt <- melt(metab_res_std,id.vars=c("Metabolites"),variable.name='Accession')
metab_res_std_melt$Accession <- factor(metab_res_std_melt$Accession, labels=int_lines)
metab_res_std_melt <- merge(metab_res_std_melt,sig_mat)

#Abbreviate the metabolite names (only take those before ";")
to_abbrev <- grep(";",levels(metab_res_std_melt$Metabolites))
levels(metab_res_std_melt$Metabolites)[to_abbrev] <- gsub(";.*$","",levels(metab_res_std_melt$Metabolites)[to_abbrev])
levels(metab_res_std_melt$Metabolites) <- tolower(levels(metab_res_std_melt$Metabolites))

metab_res_std_melt$Metabolites <- gsub("1-aminocyclopropane-1-carboxylic acid","1-aminocyclopropane",metab_res_std_melt$Metabolites)

if(org=="Roots"){
  org2 <- "Belowground"
  ht <- 8
}else{
  org2 <- "Aboveground"
  ht <- 6
}

#make the plot
sig_metabs <- filter(metab_res_std_melt,abs(Sig)>0)
sig_metabs <- sig_metabs$Metabolites
metab_res_std_melt_sub <- metab_res_std_melt[metab_res_std_melt$Metabolites%in%sig_metabs,]
qplot(Accession,Metabolites,fill=value,data=metab_res_std_melt_sub,geom='raster')+scale_fill_gradient2()+
  theme(axis.text.y = element_text(color="Black",face='bold',size=10))+labs(fill="Scaled log\nfold change")+theme_hm(col='white')+
  geom_hline(aes(yintercept=y),data=ylines_df,colour='gray75',size=.2)+
  geom_vline(aes(xintercept=x),data=xlines_df,colour='gray75',size=.2)+coord_fixed(ratio=1)+
  geom_tile(aes(x=Accession,y=Metabolites),data=filter(metab_res_std_melt,abs(Sig)>0),colour='yellow',size=1.5)+
  ggtitle(org2)+theme(plot.title=element_text(face='bold',size=15))
#ggsave(paste0("Figure4",fig_id,".eps"),width=7,height=4.375)



