mult_fun_ex <- function(x,b,k){
  return(b*(x^k))
}

add_fun_ex <- function(x,b,k){
  return(log(b)+k*x)
}

library(rstan)
library(dplyr)
library(ggplot2)
line_names <- read.csv("Data/GoBrachy_30_lines.csv",header=FALSE)

#---- Get data in a reasonable form -----#
bg_data <- read.csv("Data/Cleand_GoBrachy_Biomass_Data.csv")
colnames(bg_data)[-c(1:2)] <- c(paste("Control",c("AG","BG","Total"),sep="_"),paste("Drought",c("AG","BG","Total"),sep="_"))
control_data <- bg_data[c(1:5)]
drought_data <- bg_data[-c(3:5)]
colnames(control_data)[c(3:5)] <- colnames(drought_data)[c(3:5)] <- c("AG","BG","Total")
control_data$Treatment <- "Control"
drought_data$Treatment <- "Drought"

bg_data <- rbind(control_data,drought_data)
bg_data$Accession <- as.factor(bg_data$Line)
bg_data$LogAG <- log(bg_data$AG)
bg_data$LogBG <- log(bg_data$BG)
bg_data$LogTotal <- log(bg_data$Total)

#Remove NAs
to_rm <- unique(c(which(is.na(bg_data$AG)),which(is.na(bg_data$BG))))
bg_data <- bg_data[-to_rm,]

#Create drought indicator
bg_data$DIndicator <- 0
bg_data$DIndicator[bg_data$Treatment=="Drought"] <- 1

control_dat <- filter(bg_data,Treatment=="Control")
drought_dat <- filter(bg_data,Treatment=="Drought")

#------ Models -----------#


##----------------------------------------------------------------------------------##
#------ BHM model will accession specific slopes/intercepts on log scale  -----------#
##----------------------------------------------------------------------------------##

#Tau is the measurement error on the shoot biomass on the log scale
bhm_lscale_data <- list(N=nrow(bg_data), log_x = log(bg_data$AG), log_y = log(bg_data$BG), tau = .01,
                       drought=bg_data$DIndicator, accession = as.numeric(as.character(bg_data$Accession)),
                       M=length(unique(bg_data$Accession)))

#Had divergent transitions with default adapt_delta (0.8) so I increased it to 0.95
#Several transitions after warmup exceeded the default maximum treedepth (10) so I increased it to 15

bhm_lscale_res <- stan(file = "Code/BHM_allom_lscale.stan", data= bhm_lscale_data,
                       iter=5000, chains=3,thin = 5, verbose = FALSE, control=list(adapt_delta=0.95, max_treedepth = 15),
                       pars=c("k_drought","k_control","b_drought","b_control","beta_drought","beta_interaction","sigma",
                              "x_mean","x_sd","b_control_pop","b_drought_pop","k_control_pop","k_drought_pop"))

#Extract parameters and look at traceplots, histograms...
bhm_lscale_pars <- extract(bhm_lscale_res)

#Plot all the data along with their allometric lines
kd_post <- colMeans(bhm_lscale_pars$k_drought)
kc_post <- colMeans(bhm_lscale_pars$k_control)
bd_post <- colMeans(bhm_lscale_pars$b_drought)
bc_post <- colMeans(bhm_lscale_pars$b_control)

kd_pop_post <- mean(bhm_lscale_pars$k_drought_pop)
kc_pop_post <- mean(bhm_lscale_pars$k_control_pop)
bd_pop_post <- mean(bhm_lscale_pars$b_drought_pop)
bc_pop_post <- mean(bhm_lscale_pars$b_control_pop)

ests_bhm <- data.frame(Line=1:30,ContronB = bc_post, ControlK = kc_post, DroughtB = bd_post, DroughtK = kd_post, Accession=line_names$V2)
#Original scale
control_plot <- qplot(AG,BG,colour=Accession,data=control_dat)+xlab("Above Ground Biomass (g)")+ylab("Below Ground Biomass (g)")+
  theme_bw()+theme(legend.position="none")+ggtitle("Control Group")
drought_plot <- qplot(AG,BG,colour=Accession,data=drought_dat)+xlab("Above Ground Biomass (g)")+ylab("Below Ground Biomass (g)")+
  theme_bw()+theme(legend.position="none")+ggtitle("Drought Group")

for(i in 1:30){
  control_plot <- control_plot+stat_function(fun=mult_fun_ex,args=list(b=ests_bhm$ContronB[i],k=ests_bhm$ControlK[i]),colour="gray75",alpha=.4)
  drought_plot <- drought_plot+stat_function(fun=mult_fun_ex,args=list(b=ests_bhm$DroughtB[i],k=ests_bhm$DroughtK[i]),colour="gray75",alpha=.4)
}

#Add pop level line
control_plot <- control_plot+stat_function(fun=mult_fun_ex,args=list(b=bc_pop_post,k=kc_pop_post),colour="black")
drought_plot <- drought_plot+stat_function(fun=mult_fun_ex,args=list(b=bd_pop_post,k=kd_pop_post),colour="black")

#Label the obviously different accessions
drought_plot <- drought_plot+annotate('text',label="Mon3",x=.32,y=.47)
drought_plot <- drought_plot+annotate("text",label="Bd30-1",x=.29,y=.36)
control_plot <- control_plot+annotate("text",label="Bd30-1",x=.4,y=.37)
gridExtra::grid.arrange(control_plot,drought_plot,nrow=1)
#ggsave("Allometric_BHM_Accession_Model.png",plot = gridExtra::grid.arrange(control_plot,drought_plot,nrow=1),width=10,height=5)


#Log scale
control_plot <- ggplot(data=control_dat)+ylab("Belowground Biomass (log g)")+xlab("Aboveground Biomass (log g)")+
  theme_bw()+theme(legend.position="none")+xlim(c(-3,0))

drought_plot <- ggplot(data=drought_dat)+ylab("Belowground Biomass (log g)")+xlab("Aboveground Biomass (log g)")+
  theme_bw()+theme(legend.position="none")+xlim(c(-3,0))

#EPS doesn't support 'alpha' arguments so I've commented them out
for(i in 1:30){
  control_plot <- control_plot+stat_function(fun=add_fun_ex,args=list(b=ests_bhm$ContronB[i],k=ests_bhm$ControlK[i]),colour="gray75")#,alpha=0.4)
  drought_plot <- drought_plot+stat_function(fun=add_fun_ex,args=list(b=ests_bhm$DroughtB[i],k=ests_bhm$DroughtK[i]),colour="gray75")#,alpha=0.4)
}

#Add pop level line
control_plot <- control_plot+stat_function(fun=add_fun_ex,args=list(b=bc_pop_post,k=kc_pop_post),colour="black")
drought_plot <- drought_plot+stat_function(fun=add_fun_ex,args=list(b=bd_pop_post,k=kd_pop_post),colour="black")

#Add points on top
control_plot <- control_plot+geom_point(aes(log(AG),log(BG),colour=Accession))
drought_plot <- drought_plot+geom_point(aes(log(AG),log(BG),colour=Accession))


#Label the obviously different accessions
drought_plot <- drought_plot+annotate('text',label="Mon3",x=-2.5,y=-.55)
drought_plot <- drought_plot+annotate("text",label="Bd30-1",x=-.5,y=-.5)
control_plot <- control_plot+annotate("text",label="Bd30-1",x=-.5,y=-.6)
gridExtra::grid.arrange(control_plot,drought_plot,nrow=1)
#ggsave("Figure2_CD.eps",plot = gridExtra::grid.arrange(control_plot,drought_plot,nrow=1),width=10,height=5)

#------- Difference between treats within each accession--------#

#For which lines is there a significant change in b (beta_drought) or k (beta_interaction) in the drought/control treatments
beta_d <- data.frame(bhm_lscale_pars$beta_drought)
beta_i <- data.frame(bhm_lscale_pars$beta_interaction)
colnames(beta_d) <- colnames(beta_i) <- line_names$V2

beta_d_melt <- reshape2::melt(beta_d,variable.name='Accession')
beta_d_melt$Parameter <- "b"
beta_i_melt <- reshape2::melt(beta_i,variable.name='Accession')
beta_i_melt$Parameter <- "k"
beta_melt <- rbind(beta_d_melt,beta_i_melt)

#What percentile is 0 in posterior distribution for each accession?
b_pvals <- apply(beta_d,2,function(x){length(which(x<0))/length(x)})
b_pvals[b_pvals>.5] <- 1-b_pvals[b_pvals>.5]
k_pvals <- apply(beta_i,2,function(x){length(which(x>0))/length(x)})
k_pvals[k_pvals>.5] <- 1-k_pvals[k_pvals>.5]

pval_df <- data.frame(Accession=line_names$V2, b=b_pvals, b_sig = b_pvals<0.05, k=k_pvals, k_sig=k_pvals<0.05)
pval_df <- pval_df[order(pval_df$b),]
pval_df$Significant <- as.logical(pval_df$b_sig+pval_df$k_sig)

#Reorder beta_melt$Accession based on p-value
beta_melt$Accession <- factor(beta_melt$Accession,levels=pval_df$Accession)

#Plot parameter distributions and color based on p-value significance
qplot(value,geom='density',data=beta_melt,fill=Parameter,alpha=I(.5))+facet_wrap(~Accession)+
  geom_vline(aes(xintercept=0,colour=Significant),data=pval_df)+scale_colour_manual(values=c("black","red"),guide=FALSE)+
  xlab("")+ylab("")+theme(legend.text=element_text(face="italic"))
#ggsave("Extended_Data_Fig1.eps",width=8,height=8)


#------- Difference between accessions for each parameter --------#

k_drought <- data.frame(bhm_lscale_pars$k_drought, Treatment = "Drought",Parameter = "k")
k_control <- data.frame(bhm_lscale_pars$k_control, Treatment = "Control",Parameter = "k")
b_drought <- data.frame(bhm_lscale_pars$b_drought, Treatment = "Drought",Parameter = "b")
b_control <- data.frame(bhm_lscale_pars$b_control, Treatment = "Control",Parameter = "b")
colnames(k_drought)[1:30] <- colnames(k_control)[1:30] <- as.character(line_names$V2)
colnames(b_drought)[1:30] <- colnames(b_control)[1:30] <- as.character(line_names$V2)

allom_pars <- rbind(k_drought,k_control,b_drought,b_control)

allom_pars_melt <- reshape2::melt(allom_pars,id.vars=c("Treatment","Parameter"),variable.name="Accession")

#Plot each accession distribution and facet by parameter
qplot(value,geom='density',data=filter(allom_pars_melt,Parameter=="k"),fill=Accession,alpha=I(.5))+facet_grid(~Treatment)+
  ggtitle("Posterior Distributions for k")+theme_bw()+xlab("")+ylab("")+geom_vline(xintercept=.5,colour='gray50')
#ggsave("Allometric_posterior_k.png",width=11,height=5)
qplot(value,geom='density',data=filter(allom_pars_melt,Parameter=="b"),fill=Accession,alpha=I(.5))+facet_grid(~Treatment)+
  ggtitle("Posterior Distributions for b")+theme_bw()+xlab("")+ylab("")
#ggsave("Allometric_posterior_b.png",width=11,height=5)


