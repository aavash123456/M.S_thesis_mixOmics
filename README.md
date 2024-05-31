# M.S_thesis_mixOmics
another chapter of my masters thesis, integrating transcriptomics and proteomics datasets. 

library(tidyverse)
library(readxl)
library(dplyr)
library("data.table")
library(edgeR)
library(zoo)
library(pls)
library(igraph)
library(mixOmics)
library(caTools)


protein_data="/Users/aavashadhikari/Desktop/data/mergedResultsMaster.xlsx"
rna_data="/Users/aavashadhikari/Desktop/data/mergedcounts.csv"
meta_file="/Users/aavashadhikari/Desktop/data/design.csv"
ResultDir<-"/Users/aavashadhikari/Desktop/results/"

# rnaseq data
rna_df=fread(rna_data)
colnames(rna_df) <- sub("-.*", "", colnames(rna_df))
rna_df <- rna_df[, -"FC2"]
#rownames(rna_df) <- rna_df$"geneid"
new_rna_df = t(rna_df)
colnames(new_rna_df) <- rna_df$geneid
new_rna_df <- new_rna_df[-1, ]

# Protein data
protein_df<-read_excel(protein_data, sheet = "mergedResultsMaster")
protein_df_1<- protein_df[,c(37:60)]
rownames(protein_df_1) <- protein_df$UniprotID

protein_df2 <- t(protein_df_1)
rownames(protein_df2) <- sub(".y", "", rownames(protein_df2))

protein_df2 <- as.data.frame(protein_df2)

#Identify common rows present in both samples:

common_rows <- intersect(colnames(rna_df), rownames(protein_df2))

# selecting only common samples in both dataset
new_protein_df <- protein_df2[common_rows, ]
new_protein_df<-type.convert(new_protein_df, as.is = TRUE)

new_protein_df=na.aggregate(new_protein_df, FUN = mean)


# meta_data
meta_df <- read.csv("/Users/aavashadhikari/Desktop/data/design.csv")
meta_df$sample <- trimws(meta_df$sample)

# normalize the RNAseq data

# Change the name of the first column to "geneid"

 x <- rna_df %>% column_to_rownames(., var = "geneid")

  group<-factor(meta_df$condition)

 y <- DGEList(counts=x, group=group,remove.zeros = TRUE)

 keep<-filterByExpr(y)
 y<-y[keep,,keep.lib.sizes=FALSE]
 y<-normLibSizes(y)
 design<-model.matrix(~group)
 y<-estimateDisp(y,design)
 fit<-glmQLFit(y,design)
 qlf<-glmQLFTest(fit,coef=2)

# get the normalized counts:
 rnaseq_cpms <- cpm(y, log=FALSE)
 rnaseq_cpms2 <- t(rnaseq_cpms)

 dim(new_protein_df)
 
 
 # Initail observation and exploration of dataset

 # PCA
 result_PCA_transcriptomcis2 <- pca(rnaseq_cpms2)
 plotIndiv(result_PCA_transcriptomcis2, group = meta_df$condition)
 

 result_pca_proteomics <- pca(new_protein_df)
 plotIndiv(result_pca_proteomics, group = meta_df$condition)

 tune.pca.rna <- tune.pca(rnaseq_cpms2, ncomp = 10, scale = TRUE)
 plot(tune.pca.rna)
 
 tune.pca.rna$cum.var
 
 
 tune.pca.protein <- tune.pca(new_protein_df, ncomp = 10, scale = TRUE)
 plot(tune.pca.protein)
 tune.pca.protein$cum.var
 
 # FINAL PCA WITH 3 COMPONENT

final.pca.rna <- pca(rnaseq_cpms2, ncomp = 3, center = TRUE, scale = TRUE)

plotIndiv(final.pca.rna)
#plotVar(final.pca.rna)
selectVar(final.pca.rna, comp =1)
final.pca.rna$var.tot
final.pca.rna$prop_expl_var
final.pca.rna$cum.var

head(selectVar(final.pca.rna, comp = 3)$value)

# sample plots
install.packages("rgl")
library(rgl)

image <- plotIndiv(final.pca.rna, style = '3d', group = meta_df$condition)


plotVar(final.pca.rna , comp = c(1, 2),
        var.names = TRUE,
        cex = 3,         # To change the font size
        # cutoff = 0.5,  # For further cutoff
        title = 'pca rna, PCA comp 1 - 2',
        cutoff = 0.95)

# BIPLOT
biplot(final.pca.rna, group = meta_df$condition, cutoff = c(0.975))

# boxplot

FBgn0261674 <- scale(rnaseq_cpms2[, 'FBgn0261674'], center = TRUE, scale = TRUE)

boxplot(FBgn0261674 ~ meta_df$condition, col = color.mixo(1:9),
        xlab = 'condition', ylab = 'Expression levels, scaled',
        par(cex.axis = 0.5), # Font size
        main = 'FBgn0261674 gene')


# sparse PCA

grid.keepX = c(seq(5, 20, 5))

set.seed(30) # For reproducibility with this handbook, remove otherwise
tune.spca.rna <- tune.spca(rnaseq_cpms2, ncomp = 3, 
                              folds = 3, 
                              test.keepX = grid.keepX, nrepeat = 5)

tune.spca.rna$choice.keepX

plot(tune.spca.rna)

keepX.select <- tune.spca.rna$choice.keepX[1:2]

final.spca.rna <- spca(rnaseq_cpms2, ncomp = 2, keepX = keepX.select)

final.spca.rna$prop_expl_var

plotIndiv(final.spca.rna, comp = c(1, 2), ind.names = TRUE, 
          group = meta_df$condition, legend = TRUE)

 biplot(final.spca.rna, group = meta_df$condition, cutoff = 0.95,
       legend =FALSE)


plotVar(final.spca.rna, comp = c(1, 2), var.names = TRUE, cex = 4, cutoff = 0.9)

# We can extract the variable names and their positive or negative contribution to a given component (here 2), using the selectVar() function:


selectVar(final.spca.rna, comp = 2)$value

plotLoadings(final.spca.rna, comp = 1)

plotLoadings(final.spca.rna, comp = 2)

# PlS2

tune.pls2 <- pls(X= rnaseq_cpms2, Y = new_protein_df, ncomp = 5, mode = 'regression')

set.seed(33)

# long step Q2.pls2 <- perf(tune.pls2, validaton = 'loo', folds = 10, nrepeat = 5)


plot(Q2.pls2, criterion = 'Q2.total')

?tune.spls

list.keepX <- c(seq(5, 20, 5))
list.keepy <- C(3:5)

set.seed(33)


library(BiocParallel)
register(MulticoreParam(workers = 1))

tune.spls1 <- tune.spls(X = rnaseq_cpms2, Y = new_protein_df, test.keepX = list.keepX,
                       test.keepY = list.keepy, ncomp = 2, validation = 'Mfold', folds = 3, 
                        mode = 'regression')
#got stuck here:)

 # PLS-DA

plsda.transcriptomics <- plsda(rnaseq_cpms2, meta_df$condition, ncomp =10)
cim(plsda.transcriptomics)

perf.plsda.transcriptomics <- perf(plsda.transcriptomics, validation = 'Mfold',
                                   nrepeat = 3, folds =10, 
                                   progressBar = TRUE)
plot(perf.plsda.transcriptomics, sd = TRUE, legend.position = 'horizontal')


final.plsa.transcriptomics <- plsda(rnaseq_cpms2, meta_df$condition, ncomp =2)
                                



plotIndiv(final.plsa.transcriptomics, ellipse = TRUE,
          title = 'PLS-DA-Transcriptomics')
        


plsda.proteomics <- plsda(new_protein_df, meta_df$condition, ncomp =10)
perf.plsda.proteomics <- perf(plsda.proteomics, validation = 'Mfold',
                                   nrepeat = 3, folds =10, 
                                   progressBar = TRUE)
plot(perf.plsda.proteomics, sd = TRUE, legend.position = 'horizontal')


final.plsa.proteomics <- plsda(new_protein_df, meta_df$condition, ncomp =2)


plotIndiv(final.plsa.proteomics, ellipse = TRUE,
          title = 'PLS-DA-proteomics')






# sPLS-DA

list.keepX1 <- c(seq(5, 20, 5))

tune.splsda.proteomics <- tune.splsda(new_protein_df, meta_df$condition, ncomp =4,
                                      validation = 'Mfold', 
                                      folds = 5, dist = 'max.dist', 
                                      test.keepX = list.keepX, nrepeat = 10)

head(tune.splsda.proteomics$error.rate)

plot(tune.splsda.proteomics, sd = TRUE)

#Finding the optimal no of components

tune.splsda.proteomics$choice.ncomp$ncomp

tune.splsda.proteomics$choice.keepX


select.keepX <- tune.splsda.proteomics$choice.keepX[1:2]
splsda.proteomics <- splsda(new_protein_df, meta_df$condition, ncomp = 2, 
                            keepX = select.keepX)




perf.splsda.proteomics <- perf(splsda.proteomics, folds =5, validation = 'Mfold',
                               dist = "max.dist", progressBar = TRUE, nrepeat = 10)

perf.splsda.proteomics$error.rate.class


#checking reproducibility for cross validation

par(mfrow=c(1,2))
# For component 1
stable.comp1 <- perf.splsda.proteomics$features$stable$comp1


barplot(stable.comp1, xlab = 'variables selected across CV folds', 
        ylab = 'Stability frequency',
        main = 'Feature stability for comp = 1')

# For component 2
stable.comp2 <- perf.splsda.proteomics$features$stable$comp2

barplot(stable.comp2, xlab = 'variables selected across CV folds', 
        ylab = 'Stability frequency',
        main = 'Feature stability for comp = 2')
par(mfrow=c(1,1))

# First extract the name of selected var:
select.name <- selectVar(splsda.proteomics, comp = 1)$name

#then extract the stability values from perf:
stability <- perf.splsda.proteomics$features$stable$comp1[select.name]

# Just the head of the stability of the selected var:
head(cbind(selectVar(splsda.proteomics, comp = 1)$value, stability))


 # sample visualization
plotIndiv(splsda.proteomics, comp = c(1,2),
          ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          star = TRUE,
          title = 'proteomics, sPLS-DA comp 1 - 2')

  # variable visualization
var.name.short <- substr(new_protein_df$UniprotID[, 2], 1, 10)
plotVar(splsda.proteomics, comp = c(1,2), 
        var.names = list(var.name.short), cex = 3)


# PLOTLOADONGS
plotLoadings(splsda.proteomics, comp = 1, method = 'mean', contrib = 'max')
 
 
plotLoadings(splsda.proteomics, comp = 2, method = 'mean', contrib = 'max')


cim(splsda.proteomics)

# comp = 1 gives clear seperation of results

# N INTEGRATION





#  train and test
#set.seed(123)
samples <- colnames(rna_df)
samples <- samples[2:19]

#samplet <- sample.split(meta_df$sample, 0.75)
#samples_train <- subset(meta_df$sample, samplet==TRUE)
#samples_test <- subset(meta_df$sample, samplet==FALSE)

# Define the rows you want to keep
selected_rows <- c("FA1", "FA2", "FA3", 
                   "FB3", "FB1",
                   "FC1", 
                   "SA1", "SA2",
                   "SB2", "SB3",
                   "SC1", "SC3") 

unselected_rows <- c( "FB2", "FC3", "SA3", "SC2", "SB1")


# Subset the data frame into test and train based on selected rows

#train
rna_df_train <- rnaseq_cpms2[selected_rows, ]

protein_df_train <- new_protein_df[selected_rows, ]


#test
rna_df_test <- rnaseq_cpms2[unselected_rows, ]

protein_df_test <- new_protein_df[unselected_rows, ]

#metadata
metadf_train <- meta_df[meta_df$sample %in% selected_rows, ]

metadf_test <- meta_df[meta_df$sample %in% unselected_rows, ]

# Check if the row names of rna_df_train and protein_df_train are identical
identical(row.names(rna_df_train), row.names(protein_df_train))



#------------------Mixomics--------------------------------------


X <- list(mRNA = rna_df_train,
          protein = protein_df_train)

Y<-metadf_train$condition


# Create a design matrix with two columns
design <- matrix(0.8, nrow = length(X), ncol = 2,
                 dimnames = list(names(X), c("X", "Y")))

# Set diagonal elements to 0
diag(design) <- 0

#perform regression analyses with PLS to further understand the correlation between data sets.

#we run PLS with one component and calculate the cross-correlations between components associated to each data set:

res1.pls <- pls(X$mRNA, X$protein, ncomp = 1)

cor(res1.pls$variates$X, res1.pls$variates$Y)

#YIELDS VALUE OF 90.32%


### 2) Number of components
diablo <- block.plsda(X, Y, ncomp = 5, design = design)

                    
set.seed(123) # For reproducibility
perf.diablo = perf(diablo, validation = 'Mfold', folds = 5, nrepeat =10, progressBar = TRUE)


# Plot of the error rates based on weighted vote
tiff(file=file.path(ResultDir,paste0("model_performence.tiff")), units="in", width=12, height=12, res=800)
plot(perf.diablo)
dev.off()

ncomp <- perf.diablo$choice.ncomp $WeightedVote["Overall.BER", "centroids.dist"]

#output for no of comp = 2


perf.diablo$choice.ncomp$WeightedVote



### 3) Number of variables to select
set.seed(123) # for reproducibility
test.keepX <- list(mRNA = c(seq(10, 100, 10)),
                   proteomics = c(seq(5, 50, 10)))


tune.diablo <- tune.block.splsda(X, Y, ncomp = 2,
                                      test.keepX = test.keepX, design = design,
                                      validation = 'loo', 
                                      BPPARAM = BiocParallel::SnowParam(workers = 2),
                                      dist = "centroids.dist", progressBar = TRUE)

list.keepX <- tune.diablo$choice.keepX

list.keepX



# ########--------------------- Final model -----------------------------------
diablo <- block.splsda(X, Y, ncomp = 2,
                            keepX = list.keepX, design = design)

select_protein_comp1=selectVar(diablo, block = 'protein', comp = 1)$protein$value
select_protein_comp2=selectVar(diablo, block = 'protein', comp = 2)$protein$value

select_mRNA_comp1=selectVar(diablo, block = 'mRNA', comp = 1)$mRNA$value
select_mRNA_comp2=selectVar(diablo, block = 'mRNA', comp = 2)$mRNA$value


# #######---------------------Sample plots---------------------------------------

# plotDiablo
tiff(file=file.path(ResultDir,paste0("correlations_components_comp1.tiff")), units="in", width=12, height=12, res=800)
plotDiablo(diablo, ncomp = 1)
dev.off()

tiff(file=file.path(ResultDir,paste0("correlations_components_comp2.tiff")), units="in", width=12, height=12, res=800)
plotDiablo(diablo, ncomp = 2)
dev.off()


# plotIndiv
tiff(file=file.path(ResultDir,paste0("scatter_plots_for_individuals.tiff")), units="in", width=12, height=12, res=800)
plotIndiv(diablo, ind.names = FALSE, legend = TRUE, title = 'Control, treatment comp 1 - 2')
dev.off()


# plotArrow
#tiff(file=file.path(ResultDir,paste0("scatter_plots_for_individuals_arrow.tiff")), units="in", width=12, height=12, res=800)
plotArrow(diablo, ind.names = FALSE, legend = TRUE, title = 'Control, treatment comp 1 - 2')
dev.off()


# #######--------------------Variable plots-------------------------------------------
# plotVar

tiff(file=file.path(ResultDir,paste0("variables_representation.tiff")), units="in", width=12, height=12, res=800)
plotVar(diablo, var.names = FALSE, style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2,2),
        col = c('darkorchid', 'brown1'),
        title = 'Control, treatment comp 1 - 2')
dev.off()


# circosPlot
tiff(file=file.path(ResultDir,paste0("circos_plot.tiff")), units="in", width=12, height=12, res=800)
circosPlot(diablo, cutoff = 0.7, line = TRUE,
           color.blocks = c('darkorchid', 'brown1'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)
dev.off()

# network
network(diablo, blocks = c(1,2),
        cutoff = 0.4,
        color.node = c('darkorchid', 'brown1'),save = 'png')

# To save the plot, uncomment below line
#save = 'png', name.save = 'diablo-network'

myNetwork <- network(diablo, blocks = c(1,2), cutoff = 0.4)

write_graph(myNetwork$gR, file = "myNetwork.gml", format = "gml")

# #plotLoadings
tiff(file=file.path(ResultDir,paste0("plot_loading_comp1.tiff")), units="in", width=12, height=12, res=800)
plotLoadings(diablo, comp = 1, contrib = 'max', method = 'median')
dev.off()


tiff(file=file.path(ResultDir,paste0("plot_loading_comp2.tiff")), units="in", width=12, height=12, res=800)
plotLoadings(diablo, comp = 2, contrib = 'max', method = 'median')
dev.off()

# cimDiablo #HEATMAP
tiff(file=file.path(ResultDir,paste0("Clustered_Image_Map.tiff")), units="in", width=12, height=12, res=800)
cimDiablo(diablo, color.blocks = c('darkorchid', 'brown1'),
          comp = 1, margin=c(8,20), legend.position = "right")

dev.off()

# taking 3 components
tune.diablo3 <- tune.block.splsda(X, Y, ncomp = 3,
                                  test.keepX = test.keepX, design = design,
                                  validation = 'loo', 
                                  BPPARAM = BiocParallel::SnowParam(workers = 2),
                                  dist = "centroids.dist", progressBar = TRUE)

list.keepX.3 <- tune.diablo3$choice.keepX

list.keepX.3

# final model

diablo3 <- block.splsda(X, Y, ncomp = 3,
                        keepX = list.keepX.3, design = design)


proteins_comp3=selectVar(diablo3, block = 'protein', comp = 3)$protein$value

mRNA_comp3=selectVar(diablo3, block = 'mRNA', comp = 3)$mRNA$value

# plots

plotDiablo(diablo3, ncomp = 3)

# loading plot

plotLoadings(diablo3, comp = 3, contrib = 'max', method = 'median')

heatmap

cimDiablo(diablo3, color.blocks = c('darkorchid', 'brown1'),
          comp = 3, margin=c(8,20), legend.position = "right")





plotMarkers(diablo3, block = "mRNA", comp = 1)dev.off()

# Model performance and prediction
set.seed(123) # For reproducibility with this handbook, remove otherwise
perf.diablo <- perf(diablo,  validation = 'Mfold', folds = 5,
                         nrepeat = 10, dist = 'centroids.dist')

# Performance with Majority vote
perf.diablo$MajorityVote.error.rate

# Performance with Weighted vote
perf.diablo$WeightedVote.error.rate

tiff(file=file.path(ResultDir,paste0("Model_performance_prediction_mRNA.tiff")), units="in", width=12, height=12, res=800)
auc.diablo <- auroc(diablo, roc.block = "mRNA", roc.comp = 2,print = FALSE)

# ROC CURVE
auc.diablo

dev.off()

tiff(file=file.path(ResultDir,paste0("Model_performance_prediction_protein.tiff")), units="in", width=12, height=12, res=800)
auc.diablo <- auroc(diablo, roc.block = "protein", roc.comp = 5,print = FALSE)
dev.off()


# Prepare test set data: here one block (proteins) is missing
data.test <- list(mRNA = rna_df_test,
                       protein = protein_df_test)

predict.diablo <- predict(diablo, newdata = data.test)


confusion.mat <- get.confusion_matrix(truth = metadf_test$condition,
                                           predicted = predict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat

get.BER(confusion.mat)

# CREATING roc CURVE AFTER THIS


# Performance with Weighted vote on the test set
auc.diablo1 <- auroc(diablo, newdata = data.test, outcome.test = as.factor(metadf_test$condition), 
                          roc.block = c("rna", "protein"), roc.comp = 2,
                          print = FALSE)




# To highlight keyrelationships between more than two datasets without any  supervised purpose, we can apply regularized generalized canonical correlation analysis(RGCCA)-related methods by following this workflow:


Data <- list(Transcriptomics = rnaseq_cpms2, Proteomics = new_protein_df)

Result<- wrapper.rgcca(X = Data, ncomp=2)

plotIndiv(Result, group = meta_df$condition, legend = TRUE)


#  RGCCA
data <- list(Transcriptomics = rnaseq_cpms2, Proteomics = new_protein_df)

Result_wrapper <- wrapper.rgcca(X = data, ncomp = 3)

Result_wrapper$loadings$Transcriptomics

dim(Result_wrapper$loadings$Proteomics)

dev.off()
plotIndiv(Result_wrapper, group = meta_df$condition, legend = TRUE)

plotLoadings(Result_wrapper, comp =1, ndisplay = 30, contrib = 'max')


plotVar(Result_wrapper, var.names = TRUE, cutoff = 0.97)

# Specify the file path with the filename and the ".csv" extension
file_path <- "/Users/aavashadhikari/Desktop/result.proteomics.csv"

# Save the data frame as a CSV file
write.csv(Result_wrapper$loadings$Proteomics, file = file_path, row.names = TRUE)





