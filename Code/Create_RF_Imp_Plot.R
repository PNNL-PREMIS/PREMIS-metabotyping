g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

library(dplyr)
library(tidyr)
library(ggplot2)
load("Results/BG_RF_Results.RData")
bg_res <- all_imports
load("Results/AG_RF_Results.RData")
ag_res <- all_imports
load("Results/Total_RF_Results.RData")
total_res <- all_imports

all_res <- rbind(bg_res,ag_res,total_res)
all_res$Biomass_Type <- all_res$Tissue
all_res$Tissue <- "Genotype"
all_res$Tissue[grep("Above", all_res$Metab)] <- "Aboveground"
all_res$Tissue[grep("Roots", all_res$Metab)] <- "Belowground"

all_res$Biomass_Type <- factor(all_res$Biomass_Type,levels=c("AG","BG","Total"),labels=c("Aboveground biomass (g)","Belowground biomass (g)","Total biomass (g)"))
rownames(all_res) <- NULL

#Three scatter plots of control versus drought variable importance

ggplot(data=all_res)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue),alpha=.8)+
  scale_colour_manual(values=c('#56B4E9','#FF9999','#666666'))+theme_bw()+xlim(c(0,.1))+ylim(c(0,0.1))+
  facet_grid(~Biomass_Type)+labs(colour="")+theme(panel.spacing=unit(2,"lines"))+coord_fixed(1)
#ggsave("~/Documents/iPASS/Manuscripts/GoBrachy MS1/Figures/RF_analysis_results.pdf",width=10,height=3)

#AG only
ago <- filter(all_res,Biomass_Type=="Aboveground biomass (g)",Tissue!="Genotype")
ilm <- max(max(ago$Control),max(ago$Drought))
agoplot <- ggplot(data=ago)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue))+#,alpha=.8)+
  scale_colour_manual(values=c('#56B4E9','#FF9999','#666666'))+theme_bw()+labs(colour="")+#facet_grid(~Biomass_Type)+
  theme(panel.spacing=unit(2,"lines"))+coord_fixed()+xlim(c(0,ilm))+ylim(c(0,ilm))+theme(legend.position='none')+
  xlab("")
#BG only
bgo <- filter(all_res,Biomass_Type=="Belowground biomass (g)",Tissue!="Genotype")
ilmbg <- max(max(bgo$Control),max(bgo$Drought))
bglplt <- ggplot(data=bgo)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue))+#,alpha=.8)+
  scale_colour_manual(values=c('#56B4E9','#FF9999','#666666'))+theme_bw()+labs(colour="")+#facet_grid(~Biomass_Type)+
  theme(panel.spacing=unit(2,"lines"))+coord_fixed()+xlim(c(0,ilmbg))+ylim(c(0,ilmbg))+theme(legend.position='none')+
  ylab("")

#Total only
tgo <- filter(all_res,Biomass_Type=="Total biomass (g)",Tissue!="Genotype")
ilmtot <- max(max(tgo$Control),max(tgo$Drought))
totplt <- ggplot(data=tgo)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue))+#,alpha=.8)+
  scale_colour_manual(values=c('#56B4E9','#FF9999','#666666'))+theme_bw()+labs(colour="")+#facet_grid(~Biomass_Type)+
  theme(panel.spacing=unit(2,"lines"))+coord_fixed()+xlim(c(0,ilmtot))+ylim(c(0,ilmtot))+
  ylab("")+xlab("")

legend <- g_legend(totplt)

gridExtra::grid.arrange(agoplot,bglplt,totplt+theme(legend.position='none'),legend,nrow=1,
                        widths=c(rep(3,3),1))

#ggsave("Figure6.eps",width=10,height=3,plot=gridExtra::grid.arrange(agoplot,bglplt,totplt+theme(legend.position='none'),legend,nrow=1,widths=c(rep(3,3),1.5)))


all_res_melt <- all_res%>%filter(Tissue!="Genotype")%>%gather(Treatment,Value,-Metab,-Idx,-Tissue,-Biomass_Type)
qplot(Idx,Value,data=all_res_melt,colour=Tissue,facets=Biomass_Type~Treatment)+
  scale_colour_manual(values=c('#56B4E9','#FF9999','#666666'))+theme_bw()+labs(colour="")+
  xlab("")+ylab("Mean Drop in RMSE (g)")+theme(legend.text=element_text(size=12))+
  guides(colour=guide_legend(override.aes=list(size=4)))+theme(strip.background = element_blank(),strip.text.x = element_blank(),strip.text.y = element_blank())
#ggsave("~/Documents/iPASS/Manuscripts/GoBrachy MS1/Figures/RF_analysis_results_comparison.pdf",width=8,height=10)
#ggsave("/Volumes/Genotrait/GoBrachy MS1/Manu_Figures/Extended_Data_Fig6.eps",width=8,height=10)


#Make csv for Pubudu
rownames(bg_res) <- NULL
bg_res <- bg_res[,-c(2,5)]
bg_res <- bg_res[-grep("Plant_Line",bg_res$Metab),]
colnames(bg_res)[c(2,3)] <- paste0("BG_Biomass_",colnames(bg_res)[c(2,3)])

rownames(ag_res) <- NULL
ag_res <- ag_res[,-c(2,5)]
ag_res <- ag_res[-grep("Plant_Line",ag_res$Metab),]
colnames(ag_res)[c(2,3)] <- paste0("AG_Biomass_",colnames(ag_res)[c(2,3)])

rownames(total_res) <- NULL
total_res <- total_res[,-c(2,5)]
total_res <- total_res[-grep("Plant_Line",total_res$Metab),]
colnames(total_res)[c(2,3)] <- paste0("Total_Biomass_",colnames(total_res)[c(2,3)])

comb_csv <- merge(ag_res,bg_res)
comb_csv <- merge(comb_csv,total_res)

#write.csv(x=comb_csv,file = "~/Documents/iPASS/Manuscripts/GoBrachy MS1/Tables/RF_Metabolite_importance.csv")


##-------------------------------------##
#Version for poster with more obvious colours
#AG only
ago <- filter(all_res,Biomass_Type=="Aboveground biomass (g)",Tissue!="Genotype")
ilm <- max(max(ago$Control),max(ago$Drought))
agoplot <- ggplot(data=ago)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue),alpha=.8)+
  scale_colour_manual(values=c('black','red','#666666'))+theme_bw()+labs(colour="")+#facet_grid(~Biomass_Type)+
  theme(panel.spacing=unit(2,"lines"))+coord_fixed()+xlim(c(0,ilm))+ylim(c(0,ilm))+theme(legend.position='none')+
  xlab("")
#BG only
bgo <- filter(all_res,Biomass_Type=="Belowground biomass (g)",Tissue!="Genotype")
ilmbg <- max(max(bgo$Control),max(bgo$Drought))
bglplt <- ggplot(data=bgo)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue),alpha=.8)+
  scale_colour_manual(values=c('black','red','#666666'))+theme_bw()+labs(colour="")+#facet_grid(~Biomass_Type)+
  theme(panel.spacing=unit(2,"lines"))+coord_fixed()+xlim(c(0,ilmbg))+ylim(c(0,ilmbg))+theme(legend.position='none')+
  ylab("")

#Total only
tgo <- filter(all_res,Biomass_Type=="Total biomass (g)",Tissue!="Genotype")
ilmtot <- max(max(tgo$Control),max(tgo$Drought))
totplt <- ggplot(data=tgo)+geom_abline(aes(intercept=0,slope=1),colour='gray70')+geom_point(aes(Control,Drought,colour=Tissue),alpha=.8)+
  scale_colour_manual(values=c('black','red','#666666'))+theme_bw()+labs(colour="")+#facet_grid(~Biomass_Type)+
  theme(panel.spacing=unit(2,"lines"))+coord_fixed()+xlim(c(0,ilmtot))+ylim(c(0,ilmtot))+
  ylab("")+xlab("")

legend <- g_legend(totplt)

gridExtra::grid.arrange(agoplot,bglplt,totplt+theme(legend.position='none'),legend,nrow=1,
                        widths=c(rep(3,3),1))

#ggsave("~/Documents/iPASS/Manuscripts/GoBrachy MS1/Figures/RF_analysis_results_updated_blackred.pdf",width=10,height=3,
#       plot=gridExtra::grid.arrange(agoplot,bglplt,totplt+theme(legend.position='none'),legend,nrow=1,widths=c(rep(3,3),1.5)))

