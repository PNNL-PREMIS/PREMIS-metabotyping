library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
source('Code/pheno_funs.R')
line_names <- read.csv("Data/GoBrachy_30_lines.csv",header=FALSE)

#Use results from 'Code/bayes_allometric.R' to order the accessions
line_order <- c("Mon3","Adi-2","Bd1-1","Tek-1","Arn1","BdTR2B","Kah-5", "Adi-10","Bd21-3","BdTR13a", "Bd2-3", "Koz-3", "Kah-1", "BdTR12C", "Bis-1",
                "BdTR9K","Bd3-1", "BdTR5i","Kah-6", "BdTR3C","BdTR2G","Koz-1", "BdTR11A", "BdTR10C", "Gaz-8", "Adi-12","Bd30-1","Per-1", "BdTR1i",'Bd21')

#---- Get data in a reasonable form -----#
gobrach_data <- read.csv("Data/Cleand_GoBrachy_Biomass_Data.csv")
colnames(gobrach_data)[-c(1:2)] <- c(paste("Control",c("AG","BG","Total"),sep="_"),paste("Drought",c("AG","BG","Total"),sep="_"))

control_data <- gobrach_data[c(1:5)]
drought_data <- gobrach_data[-c(3:5)]
colnames(control_data)[c(3:5)] <- colnames(drought_data)[c(3:5)] <- c("AG","BG","Total")
control_data$Treatment <- "Control"
drought_data$Treatment <- "Drought"

all_data <- rbind(control_data,drought_data)
all_data$NumLine <- all_data$Line
all_data$Line <- as.factor(all_data$Line)
all_data$AG_Ratio <- all_data$BG/(all_data$Total)

all_data$Both <- paste0(all_data$Line,all_data$Treatment)
all_data <- all_data%>%arrange(Line,Treatment)

all_data$Treatment <- factor(all_data$Treatment,levels=c("Control","Drought"),labels=c("Control","Drought"))
all_data$Line_Name <- line_names$V2[as.numeric(as.character(all_data$Line))]
all_data$Line_Name <- factor(all_data$Line_Name,levels=line_order)


##--------------------------------------------------##
#-----  Cross hair plot for both treatments      ---##
##--------------------------------------------------##

both_summary <- all_data%>%group_by(Line_Name,Treatment)%>%
  summarize(AG_Mean=mean(AG,na.rm=T), BG_Mean=mean(BG,na.rm=T), Total_Mean=mean(Total,na.rm=T),
            AG_SD=sd(AG,na.rm=T), BG_SD=sd(BG,na.rm=T), Total_SD=sd(Total,na.rm=T))%>%
  mutate(AG_Low=AG_Mean-AG_SD, AG_High=AG_Mean+AG_SD,BG_Low=BG_Mean-BG_SD, BG_High=BG_Mean+BG_SD)

#Data.frame specifying defining which points need to be labeled explicitly
label_df <- both_summary[c(2,59),]
label_df$Line_Name <- as.character(label_df$Line_Name)
label_df$Line_Name[1] <- "Mon 3"
label_df$Line_Name[2] <- "Bd-21"


con_plot <- ggplot(filter(both_summary,Treatment=="Control"))+geom_abline(aes(intercept=0,slope=1),colour='gray50',linetype=2)+
  theme_bw()+geom_errorbarh(aes(x=AG_Mean,xmin=AG_Low,xmax=AG_High,y=BG_Mean),colour='gray70')+
  geom_errorbar(aes(x=AG_Mean,ymin=BG_Low,ymax=BG_High),colour='gray70')+
  geom_point(aes(AG_Mean,BG_Mean,colour=Line_Name),size=2)+guides(colour=FALSE)+coord_fixed()+
  xlab("Aboveground Biomass (g)")+ylab("Belowground Biomass (g)")+
  geom_text(aes(x=AG_Mean,y=BG_Mean,label=Line_Name),filter(label_df,Treatment=="Control"),nudge_x=.03,nudge_y=.015,size=2.5)+
  xlim(c(.05,.55))+ylim(c(0.025,.55))

drt_plot <- ggplot(filter(both_summary,Treatment!="Control"))+geom_abline(aes(intercept=0,slope=1),colour='gray50',linetype=2)+
  theme_bw()+geom_errorbarh(aes(x=AG_Mean,xmin=AG_Low,xmax=AG_High,y=BG_Mean),colour='gray70')+
  geom_errorbar(aes(x=AG_Mean,ymin=BG_Low,ymax=BG_High),colour='gray70')+
  geom_point(aes(AG_Mean,BG_Mean,colour=Line_Name),size=2)+guides(colour=FALSE)+coord_fixed()+
  xlab("Aboveground Biomass (g)")+ylab("Belowground Biomass (g)")+
  geom_text(aes(x=AG_Mean,y=BG_Mean,label=Line_Name),filter(label_df,Treatment!="Control"),nudge_x=.03,nudge_y=.015,size=2.5)+
  xlim(c(.05,.55))+ylim(c(0.025,.55))

gridExtra::grid.arrange(con_plot,drt_plot,nrow=1)
#ggsave(plot=gridExtra::grid.arrange(con_plot,drt_plot,nrow=1),"Figure2_AB.eps",width=7,height=4)


#Extract legend and create own object out of it
point_plot <- ggplot(filter(both_summary,Treatment=="Control"))+geom_abline(aes(intercept=0,slope=1),colour='gray50',linetype=2)+
  theme_bw()+geom_errorbarh(aes(x=AG_Mean,xmin=AG_Low,xmax=AG_High,y=BG_Mean),colour='gray70')+
  geom_errorbar(aes(x=AG_Mean,ymin=BG_Low,ymax=BG_High),colour='gray70')+
  geom_point(aes(AG_Mean,BG_Mean,colour=Line_Name),size=4)+
  theme(legend.position='top',legend.title=element_blank())+
  guides(colour=guide_legend(nrow=3))

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)}
#ggsave(plot=grid::grid.draw(g_legend(point_plot)),"Accession_Legends.png",width=8,height=1.25)

##------------------------------------------------------------------------------##
#-----  Create cross-hair plot of means with sds (extended data figure 2)    ---##
##------------------------------------------------------------------------------##
library(tidyr)
bivar_all_data <- all_data%>%group_by(Treatment,Line_Name)%>%
  summarize(AG = mean(AG, na.rm=T), BG = mean(BG, na.rm=T), Total = mean(Total, na.rm=T))%>%
  gather(Type,Biomass,AG:Total)%>%spread(Treatment,Biomass)

#Add +/- sd bars by computing sd creating new vars in bivar_all_data
bivar_all_data_sd <- all_data%>%group_by(Treatment,Line_Name)%>%
  summarize(AG_sd = sd(AG, na.rm=T), BG_sd = sd(BG, na.rm=T),Total_sd = sd(Total, na.rm=T))%>%
  gather(Type,Biomass,AG_sd:Total_sd)%>%spread(Treatment,Biomass)

bivar_all_data$Control_Low <-  bivar_all_data$Control - bivar_all_data_sd$Control
bivar_all_data$Control_High <-  bivar_all_data$Control + bivar_all_data_sd$Control
bivar_all_data$Drought_Low <-  bivar_all_data$Drought - bivar_all_data_sd$Drought
bivar_all_data$Drought_High <-  bivar_all_data$Drought + bivar_all_data_sd$Drought

bivar_all_data$Type <- factor(bivar_all_data$Type, levels=c("AG","BG","Total"),labels=c("Above Ground","Below Ground","Total"))

#Data.frame specifying defining which points need to be labeled explicitly
label_df <- bivar_all_data[c(2,3,88),]
label_df$Line_Name <- as.character(label_df$Line_Name)
label_df$Line_Name[3] <- "Bd-21"

ggplot(bivar_all_data)+geom_abline(aes(intercept=0,slope=1),colour='gray50',linetype=2)+
  facet_wrap(~Type,scales='free')+theme_bw()+geom_errorbarh(aes(x=Control,xmin=Control_Low,xmax=Control_High,y=Drought),colour='gray70')+
  geom_errorbar(aes(x=Control,ymin=Drought_Low,ymax=Drought_High),colour='gray70')+
  geom_point(aes(Control,Drought,colour=Line_Name),size=2)+guides(colour=FALSE)+
  geom_text(aes(x=Control,y=Drought,label=Line_Name),label_df,nudge_x=c(0.045,0.04,0.05),nudge_y=c(-.012,.015,.02),size=2.5)
#ggsave("Extended_Data_Fig2.eps",width=9,height=9/3)






