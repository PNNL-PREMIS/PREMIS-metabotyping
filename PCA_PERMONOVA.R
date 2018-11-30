########################################
## PRINCIPAL COMPONENT ANALYSIS (PCA) ##
########################################

# Load the packages
require(missMDA)
require(ggplot2)
require(FactoMineR)

# Load the dataset.
DATASET <- read.csv("File_path.csv", sep=",", header=T)

DATA.SUBSET <- subset(DATASET)

# REPLACE 0 for NA
DATA.SUBSET[DATA.SUBSET==0] <- NA

# Delete Entire NA columns (if ANY after the DATA.SUBSET)
DATA.SUBSET <- DATA.SUBSET[colSums(!is.na(DATA.SUBSET)) > 0]

FACTOR.1.Shape <- as.factor(DATA.SUBSET$CELL.FACTOR) # Factor.1.Shape will determine the shape of cases
FACTOR.2.Color <- as.factor(DATA.SUBSET$CELL.FACTOR) # Factor.2.Color will determine the color of cases

DATA <- DATA.SUBSET[,-c(1:N)] #where N is the number of categorical factors

#### PCA ANALYSIS ####

## Estimate number of dimensions necessary to impute the dataset.
Estim.ncp <- estim_ncpPCA(DATA, ncp.min=0, ncp.max=15)
ncp.min <- min(Estim.ncp$criterion)

## Impute the values of your dataset with NAS
DATA_imputed_PCA <- imputePCA(DATA, ncp=ncp.min, scale = TRUE, method=c("Regularized"))

## Perform the PCA
PCA.RESULTS <- PCA(DATA_imputed_PCA) 

## Variability for each of the 6 first PCs
PC1.var <- format(round(PCA.RESULTS$eig[,2][1],2))
PC2.var <- format(round(PCA.RESULTS$eig[,2][2],2))
PC3.var <- format(round(PCA.RESULTS$eig[,2][3],2))
PC4.var <- format(round(PCA.RESULTS$eig[,2][4],2))
PC5.var <- format(round(PCA.RESULTS$eig[,2][5],2))
PC6.var <- format(round(PCA.RESULTS$eig[,2][6],2))

## We create the labels for each of the 6 first axes of the plot
PC1.axis <- paste("PC1 ", "(", PC1.var, "%", ")", sep="")
PC2.axis <- paste("PC2 ", "(", PC2.var, "%", ")", sep="")
PC3.axis <- paste("PC3 ", "(", PC3.var, "%", ")", sep="")
PC4.axis <- paste("PC4 ", "(", PC4.var, "%", ")", sep="")
PC5.axis <- paste("PC5 ", "(", PC5.var, "%", ")", sep="")
PC6.axis <- paste("PC6 ", "(", PC6.var, "%", ")", sep="")

#### CASE PLOT PCA (with ggplot2) ####
# Extract the CASE coordinates from the PCA.RESULTS.
CASE.Coordinate.table <- as.data.frame(PCA.RESULTS$ind$coord)

## Create the CASE PLOT through GGPLOT2
ggplot(data = CASE.Coordinate.table, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(colour= FACTOR.2.Color, shape = FACTOR.1.Shape), size=8) +
                 theme_bw() +
                 ggtitle("Case Plot") +
                 theme(legend.position = "right") +
                 labs(x=PC1.axis, y=PC2.axis, size=10)
                 
                 
                 
                 
########################################
### PERMUTATIONAL MANOVA (PERMANOVA) ###
########################################

# Load the packages
require(vegan)

# Load the dataset
DATASET <- read.csv("File_patwh.csv", sep=",", header=T)


#### PREPARATION OF THE DATASET ####

# MODIFY EACH TIME IF NECESSARY ## We can create a subset of the data if necessary (run always anyway)
DATA.SUBSET <- subset(DATASET[,])

# REPLACE 0 for NA
DATA.SUBSET[DATA.SUBSET==0] <- NA
# Delete Entire NA columns (if any after DATA.SUBSET)
DATA.SUBSET <- DATA.SUBSET[colSums(!is.na(DATA.SUBSET)) > 0]

# MODIFY EACH TIME ## Exclude the Categorical Factors from the dataframe to create the matrix DATA
DATA <- DATA.SUBSET[,-c(1:N)] # N is the number of categorical factors


#### PERMANOVA MODEL (MATRIX ~ Facor1*Factor2*Factor3...) ####

adonis(DATA ~ Treatment*Genotype*Tissue, perm=10000, na.rm = T, method="euclidian")
                 