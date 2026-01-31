library(readxl)
library("Hmisc") #function corrr
library("corrplot") #function corrplot
library(tidyverse)  
library("corrr") #functions correlate and network_plot, needed for PCA 
library(ggcorrplot) #functions needed for PCA
library("FactoMineR") #functions needed for PCA
library(factoextra ) #functions needed for PCA
library(ggfortify) #functions needed for PCA

#Import and modify data for the analyses
#Important data from excel sheet with assembly parameters for different species
Species_Data_Input <- read_excel("~/000_changed_documents_20230810/Analyses/InvertOmics/Species_Data_Input.xlsx")
#change tibble format to data frame format
Species_Data_Input <- as.data.frame(Species_Data_Input)
#make first column to names of row
rownames(Species_Data_Input) <- Species_Data_Input$Species
#remove first column
Species_Data_Input <- Species_Data_Input[,-1]
#replace NA values with 0
Data_NA <- Species_Data_Input %>% replace(is.na(.), 0)

## Conduct a correlation analysis to determine highly correlated variables, which effectively measure same across the different steps of the analysis (e.g., BUSCO scores)
#remove columns 1-3 with  character values
Data_Reduced <- Data_NA[,-(1:3)]

#calculate correlation coefficients and plot the distribution of correlation coefficients
Data_Corr <- rcorr(as.matrix(Data_Reduced))
Data_Corr_NA <- Data_Corr$r
diag(Data_Corr_NA) <- NA
plot(density(na.omit(as.numeric(Data_Corr_NA))))

#calculate hierarchical cluster of correlation coefficients
hc <- hclust(dist(Data_Corr_NA), "ave")
#plot the original hierarchical cluster
plot(hc,  hang = -1, cex = 0.3)
#Generate a corellogarm for the same
corrplot(as.matrix(Data_Corr_NA), type="upper", order="hclust", hclust.method = "average", addgrid.col = "NA", tl.cex = 0.3)
#generate a correlation network
Data_Corr_NA %>% correlate() %>% network_plot(min_cor = 0.5, repel = FALSE, curved = FALSE)

##### NOTE OF HIGHLY CORRELATED VARIABLES OF SAME PROPERTY WITHIN A CLUSTER AND HEIGHT BELOW 0.5 #####
##### Busco_comp_A = Busco_comp_F = Busco_comp_P 
##### Error_A = Error_F = Error_P
##### N50_A = N50_F = N50_P
##### Larg_con_A = Larg_con_P
##### Per_corr_phyl_A = Per_corr_phyl_P
##### BUSCO_frag_A = BUSCO_frag_F = BUSCO_frag_P
##### L50_F = L50_P
##### Num_con_F = Num_con_P
##### G_size_A = G_size_F = G_size_P
##### Comp_Pri_F = Comp_Pri_P
##### BUSCO_dup_A = BUSCO_dup_F
##### GC_A = GC_F = GC_P
##### Mean_readlen = Median_readlen

#Exclude redundant variables
Data_Reduced2 <- subset(Data_Reduced, select = -c(Busco_comp_F, Busco_comp_P, Error_F, Error_P, N50_F, N50_P, Larg_con_P, Per_corr_phyl_P, BUSCO_frag_F, BUSCO_frag_P, L50_P, Num_con_P, G_size_F, G_size_P, Comp_Pri_P, BUSCO_dup_F, GC_F, GC_P, Median_readlen))

#Rerun analyses with reduced dataset
#calculate correlation coefficients and plot the distribution of correlation coefficients
Data_Corr2 <- rcorr(as.matrix(Data_Reduced2))
Data_Corr2_NA <- Data_Corr2$r
diag(Data_Corr2_NA) <- NA
plot(density(na.omit(as.numeric(Data_Corr2_NA))))

#calculate hierarchical cluster of correlation coefficients
hc2 <- hclust(dist(Data_Corr2_NA), "ave")
#plot the original hierarchical cluster
plot(hc2,  hang = -1, cex = 0.3)
#Generate a corellogarm for the same
corrplot(as.matrix(Data_Corr2_NA), type="upper", order="hclust", hclust.method = "average", addgrid.col = "NA", tl.cex = 0.3)
#generate a correlation network
Data_Corr2_NA %>% correlate() %>% network_plot(min_cor = 0.5, repel = FALSE, curved = FALSE)

#PCA analysis of the data
#Normalization of the data
Data_Reduced2_Norm <- scale(Data_Reduced2)
#Conduct PCA on the normalized data
Data_Reduced2_PCA <- princomp(na.omit(Data_Reduced2_Norm))
summary(Data_Reduced2_PCA)
#Determine the contributions of each variable with up 95% cumulative proportion
Data_Reduced2_Loadings95 <- Data_Reduced2_PCA$loadings[, 1:18]
#View(Data_Reduced2_Loadings95)
#Generate a scree plot to visualize the contribution of each component
fviz_eig(Data_Reduced2_PCA, addlabels = TRUE)
#Generate a biplot of the contribution of the variables (attributes) to the two first components
fviz_pca_var(Data_Reduced2_PCA, col.var = "black")
#Determine the contribution of each variable to the two first components
fviz_cos2(Data_Reduced2_PCA, choice = "var", axes = 1:2)
#Combine the biplot and the contribution into one visual display
fviz_pca_var(Data_Reduced2_PCA, col.var = "cos2",
             gradient.cols = c("black", "blue", "green"),
             repel = TRUE)
#Generate a biplot with the rows dispayed colored if WGA has been used or not
autoplot(Data_Reduced2_PCA, data = Species_Data_Input, colour = 'WGA',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3, shape = FALSE, label.size = 3, frame = TRUE)
#autoplot(Data_Reduced2_PCA, data = Species_Data_Input, colour = 'WGA', shape = FALSE, label.size = 3)
#autoplot(Data_Reduced2_PCA, data = Species_Data_Input, colour = 'Phylum', label = TRUE, label.size = 3, frame = TRUE)
autoplot(Data_Reduced2_PCA, data = Species_Data_Input, colour = 'Phylum', frame = TRUE)
autoplot(Data_Reduced2_PCA, data = Species_Data_Input, colour = 'Busco_comp_A') +
  scale_colour_gradientn(colours = c("black","blue","green"))
