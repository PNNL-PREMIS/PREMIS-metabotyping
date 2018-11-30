diff_est_byg <- function(y,des_mat,C_mat,gps,lines){
  ngps <- length(gps)
  res_df <- resids <- hats <- NULL
  for(i in 1:ngps){
    
    if(i<=(ngps-1)){
      if(length(gps[[i+1]])==1){
        gpi <- c(gps[[i]],gps[[i+1]])
        i <- i+1
      }else{
        gpi <- gps[[i]]
      }
    }else{
      gpi <- gps[[i]] 
    }
    
    rowsi <- which(lines%in%gpi)
    all_resi <- diff_est(yvec=y[rowsi],xmat=des_mat[rowsi,],Cmat=C_mat[gpi,])
    resi <- all_resi$resdf
    resids <- c(resids,list(all_resi$resids))
    hats <- c(hats,list(all_resi$fits))
    resi$Line <- gpi
    res_df <- rbind(res_df,resi)
    if(nrow(res_df)==length(unique(lines))){
      break
    }
  }
  res_df <- res_df[order(res_df$Line),]
  return(list(res_df=res_df,es=resids,yhats=hats))
}


#A function that takes the respone of interest (AG, BG, Total Biomass) and returnes estimated diferences for each Line
diff_est <- function(yvec,xmat,Cmat,sig_only=FALSE){
  ##------- Remove NAs from yvec and xmat ----------#
  if(any(is.na(yvec))){
    xmat <- xmat[-which(is.na(yvec)),]
    yvec <- yvec[-which(is.na(yvec))]
  }
  #Projection matrix to get residuals and esimate of sigma^2
  Px_Total <- xmat%*%MASS::ginv(t(xmat)%*%xmat)%*%t(xmat)
  yhat <- Px_Total%*%yvec
  sigma2 <- t(yvec)%*%(diag(length(yvec))-Px_Total)%*%yvec/(length(yvec)-qr(xmat)$rank)
  if(sig_only){
    return(sigma2)
  }
  
  par_vec <- MASS::ginv(t(xmat)%*%xmat)%*%t(xmat)%*%yvec 
  Cbeta <- Cmat%*%par_vec
  CbetaSE <- sqrt(diag(sigma2[1,1]*Cmat%*%MASS::ginv(t(xmat)%*%xmat)%*%t(Cmat)))
  tests <- Cbeta/CbetaSE
  pvals <- pt(tests,length(yvec)-qr(xmat)$rank,lower.tail = FALSE)
  return(list(resdf=data.frame(Line=1:length(Cbeta),Diff_hat=Cbeta,Diff_SE=CbetaSE,Pvals=pvals),
              resids=yvec-yhat,
              fits=yhat))
}

#Group the accessions such that the min and max accession SD is within a factor of 2
group_by_sd <- function(x, max_ratio=2){
  glist <- NULL
  gids <- 1:length(x)
  orig_x <- x
  while(length(x)>0){
    minx <- min(x)
    gi <- which(x<max_ratio*minx)
    glist <- c(glist,list(which(orig_x%in%x[gi])))
    x <- x[-gi]
  }
  return(glist)
}

#Function to create an X matrix based on a covariate data frame
build_x_mat <- function(cov_df){
  
  if(is.null(ncol(cov_df))){
    cov_df <- matrix(cov_df,ncol=1)
  }
  
  #If all covariates are numeric, simply return the same matrix back
  if(is.numeric(cov_df)){
    return(data.matrix(cov_df))
  }
  
  #If the covariates are a mix of numeric, factors and characters, return matrix of group identifiers
  Xmatrix <- NULL
  for(i in 1:ncol(cov_df)){
    
    if(is.numeric(cov_df[,i])){
      #If column i is numeric, append it to X
      Xmatrix <- cbind(Xmatrix,cov_df[,i])
    }else{
      coli_levels <- unique(cov_df[,i])
      n_levels <- length(coli_levels)
      if(n_levels!=length(cov_df[,i])){
        Xcoli <- matrix(0,nrow(cov_df),n_levels)
        
        for(j in 1:n_levels){
          Xcoli[cov_df[,i]==coli_levels[j],j] <- 1
        }
        Xmatrix <- cbind(Xmatrix,Xcoli)
      }
    }
  }
  return(Xmatrix)
}


#FUnction to take root/leaf data and return imd_anova results
run_imd_anova <- function(rl_data, id.cols=1:5){
  
  library(MSomicsQC)
  #The data need to be transposed for MSomics family
  metabs <- colnames(rl_data)[-id.cols]
  samp_names <- rep("ID",nrow(rl_data))
  for(i in 1:length(id.cols)){
    samp_names <- paste(samp_names,rl_data[,i],sep="_")
  }
  samp_names <- gsub(pattern = "-",replacement = "_",x = samp_names)
  root_metab_mat <- unname(data.matrix(rl_data[,-id.cols]))
  root_metab_mat <- data.frame(t(root_metab_mat))
  root_metab_mat <- cbind(Metabolites = metabs,root_metab_mat)
  colnames(root_metab_mat)[-1] <- samp_names
  
  #Create the f-data object
  root_fdata <- rl_data[,id.cols]
  root_fdata[,1] <- samp_names
  colnames(root_fdata)[1] <- "SampleID"
  root_fdata$gp <- paste(root_fdata$Treatment,root_fdata$Plant_Line,sep="_")
  
  #Put into metab
  root_gob_data <- as.metabData(e_data=root_metab_mat, f_data = root_fdata,edata_cname = "Metabolites",fdata_cname = "SampleID", data_scale='abundance')
  
  #------ Log base 2 transform the data ----#
  #For some reason, as.metabData automatically days edata is on log2 scale, but it's not
  rl_data_trans <- edata_transform(omicsData =  root_gob_data, data_scale = "log2")
  
  #------ Group the data ----#
  rl_data_trans <- group_designation(omicsData = rl_data_trans, main_effects = c("gp"), time_course=NULL)
  
  #------- Filter out biomolecules using IMD-ANOVA filter -----#
  imd_afilt <- imdanova_filter(omicsData = rl_data_trans)
  rl_data_trans <- applyFilt(filter_object = imd_afilt, omicsData = rl_data_trans, min_nonmiss_anova = 2,min_nonmiss_gtest = 3)
  #rl_data_trans <- MSomics_filter(filter_object = imd_afilt, omicsData = rl_data_trans, min_nonmiss_anova = 2,min_nonmiss_gtest = 3)
  
  
  ##---------------------##
  #- Stats analysis -----##
  library(MSomicsSTAT)
  
  comp_df <- data.frame(Control=paste("Control",unique(rl_data$Plant_Line),sep="_"),Test=paste("Drought",unique(rl_data$Plant_Line),sep="_"))
  attr(rl_data_trans,"group_DF") <- attr(rl_data_trans,"group_DF")[,c(1:2)]
  
  all(comp_df$Control%in%attr(rl_data_trans,"group_DF")$Group)
  
  root_res <- imd_anova(omicsData = rl_data_trans, comparisons = comp_df, test_method = 'comb', pval_thresh = 0.01, pval_adjust = "holm")
  return(root_res)
  
}



