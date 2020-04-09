# Script for marker selection through multistudy experiments with multiple batches.

####################################################################
################################### PREPARATION OF ENVIRONMENT
####################################################################

# required packages
library(data.table)
library(parallel)
library(doParallel)
library(caret)
library(tuple)
library(vscc)
library(factoextra)
library(FactoMineR)
library(NbClust)
library(cluster)
library(rafalib)
library(ggsci)
library(cowplot)
library(reshape2)
library(ggplot2)
library(affy)
library(ProteoMM)
library(vsn)

# stop parallel
stopCluster(cl)
stopImplicitCluster()
registerDoSEQ()

# start parallel processing
fc <- as.numeric(detectCores(logical = F))
cl <- makePSOCKcluster(fc+1)
registerDoParallel(cl)

############################# LOAD INPUT DATA
setwd("C:/...")

#by data.table package
dsr <-as.data.frame(fread(input = ".csv", header=T))
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
colnames(dsr)[1] <- c("Label")

# Half minimum missing value imputation
# i in 2:ncol to exclude "Label" column
for (i in 2:ncol(dsr)) {
  dsr[,i][is.na(dsr[,i])] <- (min(dsr[,i], na.rm = T))/2 }

# Univariate filtering
# if some error in Shapiro normality test:
# use shapiro.wilk.test function from cwhmisc instead shapiro.test from stats
# library(cwhmisc)
# norm.test <- apply(xx, 2, function(t) shapiro.wilk.test(t)$p)

uvf <- function(x){
  
  norm_homog_tests <- function(x) {
    xx <- x[,-1]
    # normality test
    norm.test <- apply(xx, 2, function(t) shapiro.test(t)$p.value)
    # homogeneity test
    homog.test <- apply(xx, 2, function(t) bartlett.test(t,g = x[,1])$p.value)
    return(as.data.frame(cbind(norm.test, homog.test)))}
  
  res_tests <- norm_homog_tests(x)
  
  
  wilcox_test <- function(x,y) {
    xx <- x[,-1]
    wx.t <- as.vector(which(y[,1] < 0.05))
    wilcox_test <- list()
    ifelse(identical(wx.t, integer(0)), return (wilcox_test <- 1), wx.t)
    wilcox_test <- lapply(as.data.frame(xx[,wx.t]),  function(t) as.vector(pairwise.wilcox.test(x = t, g =  x[,1], p.adjust.method ="BH", paired=F)$p.value))
    names(wilcox_test) <- (colnames(x)[-1])[wx.t]
    return(as.list(wilcox_test))}
  
  wx.t.res <- wilcox_test(x, res_tests)
  
  
  welch_test <- function(x,y) {
    xx <- x[,-1]
    wl.t <- as.vector(which(y[,1] > 0.05 & y[,2] < 0.05))
    welch_test <- list()
    ifelse(identical(wl.t, integer(0)), return (welch_test <- 1), wl.t)
    welch_test <- lapply(as.data.frame(xx[,wl.t]), function(t) as.vector(pairwise.t.test(x = t, g = x[,1], p.adjust.method = "BH", pool.sd = F)$p.value))
    names(welch_test) <- (colnames(x)[-1])[wl.t]
    return(as.list(welch_test))}
  
  wl.t.res <- welch_test(x, res_tests)
  
  
  student_test <- function(x,y) {
    xx <- x[,-1]
    st.t <- as.vector(which(y[,1] > 0.05 & y[,2] > 0.05))
    student_test <- list()
    ifelse(identical(st.t, integer(0)), return (student_test <- 1), st.t)
    student_test <- lapply(as.data.frame(xx[,st.t]), function(t) as.vector(pairwise.t.test(x = t, g = x[,1], p.adjust.method = "BH", pool.sd = T)$p.value))
    names(student_test) <- (colnames(x)[-1])[st.t]
    return(as.list(student_test))}
  
  st.t.res <- student_test(x, res_tests)
  
  filt_p_val <- function(x, y, z, w){
    
    #x = ds
    #y = wx.t.res
    #z = wl.t.res
    #w = st.t.res
    
    wx.t.n <- names(y)
    wx.t.res2 <-as.data.frame(y)
    colnames(wx.t.res2) <- wx.t.n
    wx.t.res2[is.na(wx.t.res2)] <- max(wx.t.res2, na.rm = T)
    
    wl.t.n <- names(z)
    wl.t.res2 <- as.data.frame(z)
    colnames(wl.t.res2) <- wl.t.n
    wl.t.res2[is.na(wl.t.res2)] <- max(wl.t.res2, na.rm = T)
    
    st.t.n <- names(w)
    st.t.res2 <- as.data.frame(w)
    colnames(st.t.res2) <- st.t.n
    st.t.res2[is.na(st.t.res2)] <- max(st.t.res2, na.rm = T)
    
    wxx <- apply(wx.t.res2, 2, min)
    wll <- apply(wl.t.res2, 2, min)
    stt <- apply(st.t.res2, 2, min)
    wxx2 <- as.data.frame(wxx)
    wll2 <- as.data.frame(wll)
    stt2 <- as.data.frame(stt)
    
    wxx3 <- rownames(wxx2)[which(wxx2 <= 0.05)]
    wll3 <- rownames(wll2)[which(wll2 <= 0.05)]
    stt3 <- rownames(stt2)[which(stt2 <= 0.05)]
    aff <- c(wxx3, wll3, stt3)
    
    ds_fil <- cbind(x[,1], x[, aff])
    return(ds_fil)
  }
  return(filt_p_val(x, wx.t.res, wl.t.res, st.t.res))
}

ds <- uvf(dsr)
colnames(ds)[1] <- "Label"

########################################################
################################### DATA MINING PIPELINE
########################################################

#repeated cross validation

trainControl <- trainControl(method="repeatedcv", number=5, repeats=3)
metric <- "Accuracy"

# KNN
set.seed(1234)
fit.knn <- train(Label~ ., data=ds, method="knn", metric=metric, trControl=trainControl, tuneLength = 5)

# SVM
set.seed(1234)
fit.svm <- train(Label~ ., data=ds, method="svmRadial", metric=metric, trControl=trainControl, tuneLength = 5)

# PLS
set.seed(1234)
fit.pls <- train(Label~ ., data=ds, method="pls", metric=metric, trControl=trainControl, tuneLength = 5)

# RF
set.seed(1234)
fit.rf <- train(Label~ ., data=ds, method="rf", metric=metric, trControl=trainControl, tuneLength = 5)

# only accuracy for all models
results <- resamples(list(knn=fit.knn, svm=fit.svm, rf=fit.rf, pls=fit.pls), trControl = trainControl, metric=metric)
results_df <- as.data.frame(results)
results_df_fin <- apply(results_df[,-5], 2, mean)
results_df_fin

# by model specific value
set.seed(1234)
Imp.rf <- varImp(fit.rf, scale = FALSE)
Imp.rf <- Imp.rf$importance
rownames(Imp.rf) <- gsub("`", '', rownames(Imp.rf))

set.seed(1234)
Imp.pls <- varImp(fit.pls, scale = FALSE)
Imp.pls <- Imp.pls$importance
rownames(Imp.pls) <- gsub("`", '', rownames(Imp.pls))

set.seed(1234)
Imp.svm <- varImp(fit.svm, scale = FALSE)
Imp.svm <- Imp.svm$importance
rownames(Imp.svm) <- gsub("`", '', rownames(Imp.svm))

set.seed(1234)
Imp.knn <- varImp(fit.knn, scale = FALSE)
Imp.knn <- Imp.knn$importance
rownames(Imp.knn) <- gsub("`", '', rownames(Imp.knn))

# creating of list with top = n important features by model
n_model = 50
set.seed(1234)

Imp.rf[,c("sum")] <- apply(X = Imp.rf, 1, FUN = sum)
Imp.rf <- Imp.rf[order(Imp.rf$sum, decreasing = T),]
Imp.rf <- rownames(Imp.rf)[c(1:n_model)]

Imp.pls[,c("sum")] <- apply(X = Imp.pls, 1, FUN = sum)
Imp.pls <- Imp.pls[order(Imp.pls$sum, decreasing = T),]
Imp.pls <- rownames(Imp.pls)[c(1:n_model)]

Imp.svm[,c("sum")] <- apply(X = Imp.svm, 1, FUN = sum)
Imp.svm <- Imp.svm[order(Imp.svm$sum, decreasing = T),]
Imp.svm <- rownames(Imp.svm)[c(1:n_model)]

Imp.knn[,c("sum")] <- apply(X = Imp.knn, 1, FUN = sum)
Imp.knn <- Imp.knn[order(Imp.knn$sum, decreasing = T),]
Imp.knn <- rownames(Imp.knn)[c(1:n_model)]

# minimum n time of duplicated 
all <- c(Imp.rf, Imp.svm, Imp.knn, Imp.pls)
n_tuple <- 2
all1 <- all[which(tuplicated(all, n_tuple), T)]
all1 <- unique(all1)
ds_d <- cbind(ds[,1], ds[,all1])

# Feature Selection
set.seed(1234)
vscc <- vscc(ds_d[,-1], G=1:9, automate = "mclust", initial = NULL, train = NULL, forcereduction = FALSE)
plot(vscc)
summary(vscc)
vscc_sel_df <- colnames(as.data.frame(vscc$topselected))[-1]
ds_vscc <- cbind(ds[,1], ds[ ,vscc_sel_df])

######################################
################################### AUROC calculation
######################################

# raw data (ds_raw)
ds_raw <-as.data.frame(fread(input = ".csv", header=T))
rownames(ds_raw) <- ds_raw[,1]
ds_raw <- ds_raw[,-1]
colnames(ds_raw)[1] <- c("Label")

# ds_raw after EigenMS (dsr)
dsr <-as.data.frame(fread(input = ".csv", header=T))
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
colnames(dsr)[1] <- c("Label")

# dsr after UVF+MVI (ds)
ds <-as.data.frame(fread(input = ".csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <- c("Label")

# ds after ML+SFE (ds_d)
ds_d <-as.data.frame(fread(input = ".csv", header=T))
rownames(ds_d) <- ds_d[,1]
ds_d <- ds_d[,-1]
colnames(ds_d)[1] <- c("Label")

# ds_d after RFE (ds_rfe)
ds_rfe <-as.data.frame(fread(input = ".csv", header=T))
rownames(ds_rfe) <- ds_rfe[,1]
ds_rfe <- ds_rfe[,-1]
colnames(ds_rfe)[1] <- c("Label")

# mean AUROC for ds_raw
set.seed(1234)
auroc_ds_raw <- filterVarImp(ds_raw[,-1], as.factor(ds_raw[,1]))
auroc_ds_raw[,c("sum")] <- apply(X = auroc_ds_raw, 1, FUN = sum)
auroc_ds_raw_mean <- round(mean(auroc_ds_raw$sum/(ncol(auroc_ds_raw)-1)),3)

# mean AUROC for dsr
set.seed(1234)
auroc_dsr <- filterVarImp(dsr[,-1], as.factor(dsr[,1]))
auroc_dsr[,c("sum")] <- apply(X = auroc_dsr, 1, FUN = sum)
auroc_dsr_mean <- round(mean(auroc_dsr$sum/(ncol(auroc_dsr)-1)),3)

# mean AUROC for ds
set.seed(1234)
auroc_ds <- filterVarImp(ds[,-1], as.factor(ds[,1]))
auroc_ds[,c("sum")] <- apply(X = auroc_ds, 1, FUN = sum)
auroc_ds_mean <- round(mean(auroc_ds$sum/(ncol(auroc_ds)-1)),3)

# mean AUROC for ds_d
set.seed(1234)
auroc_ds_d <- filterVarImp(ds_d[,-1], as.factor(ds_d[,1]))
auroc_ds_d[,c("sum")] <- apply(X = auroc_ds_d, 1, FUN = sum)
auroc_ds_d_mean <- round(mean(auroc_ds_d$sum/(ncol(auroc_ds_d)-1)),3)

# mean AUROC for ds_rfe
set.seed(1234)
auroc_ds_rfe <- filterVarImp(ds_rfe[,-1], as.factor(ds_rfe[,1]))
auroc_ds_rfe[,c("sum")] <- apply(X = auroc_ds_rfe, 1, FUN = sum)
auroc_ds_rfe_mean <- round(mean(auroc_ds_rfe$sum/(ncol(auroc_ds_rfe)-1)),3)

# final table with AUROC calculations & # of features
auroc <- data.frame(rbind(auroc_ds_raw_mean, auroc_dsr_mean, auroc_ds_mean, auroc_ds_d_mean, auroc_ds_rfe_mean))
vars_roc <- data.frame(cbind(c((ncol(ds_raw)-1), (ncol(dsr)-1),(ncol(ds)-1),(ncol(ds_d)-1),(ncol(ds_rfe)-1)), auroc))
row.names(vars_roc) <- c("ds_raw", "dsr", "ds", "ds_d", "ds_rfe")
colnames(vars_roc) <- c("# of features","mean AUROC")
vars_roc

########################################################################
# UNSUPERVISED PROJECTION METHODS
########################################################################

#input
ds_base_ul <- ds_rfe
ds_ul <- ds_base_ul[,-1]
ds_group <- ds_base_ul[,1]
label_ul <- unique(ds_base_ul[,1])

# PCA
palette_pca <- "jco" # color

pca.ds <- PCA(ds_ul, scale.unit = T, graph = FALSE)
a <- fviz_pca_ind(pca.ds,
                  title = "",
                  geom.ind = "point", # show points only (nbut not "text")
                  col.ind = ds_group, # color by groups
                  palette = palette_pca,
                  addEllipses = T, # Concentration ellipses
                  legend.title = "")


#HCA
Cols = function(vec, ord){
  cols = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(vec)))
  return(cols[as.fumeric(as.character(vec))[ord]])}

k <- length(label_ul)
ds_ul1 <- ds_ul
rownames(ds_ul1) = make.names(ds_group, unique=TRUE)
res.dist1 <- dist(ds_ul1, method = "manhattan")
res.hc1 <- hclust(d = res.dist1, method = "ward.D2")
b <- fviz_dend(res.hc1, k = k, # Cut in k groups
               cex = 0.3, # label size
               k_colors = unique(Cols(ds_group,res.hc1$order)), # "lancet" color "jco" gray JGRAY(k_hc)
               color_labels_by_k = F, # color labels by groups
               label_cols = Cols(ds_group,res.hc1$order),#Cols(ds_rfe[,1])[res.hc1$order], #as.fumeric(ds_rfe[,1])[res.hc1$order],
               rect = T, # Add rectangle around groups
               rect_fill = T,
               rect_border = unique(Cols(ds_group,res.hc1$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
               main = "",
               ylab = "",
               horiz = F,
               lwd = 0.5, # lines size
               show_labels = T)

# Silhouette method
# pam or hcut
c <- fviz_nbclust(ds_ul, clara, linecolor = "red", method = "silhouette")+
  labs(subtitle = "Silhouette method")

# "NbClust" function
nb <- NbClust(ds_ul, distance = "manhattan", min.nc = 2,max.nc = 20, method = "ward.D2")
d <-fviz_nbclust(nb)

############################
# boxplots in one graph
############################

df.m <- melt(ds_rfe, id.var = "Label") # reshape data frame

p <- ggplot(data = df.m, aes(x=variable, y=value)) + xlab("") + ylab("") +
  geom_boxplot(aes(fill=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) 
                                                                                                      
pp <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="bottom") 

############################
# creating & saving multiplot
############################

# Unsupervised Projection
ul_plot <- plot_grid(a,b,  labels = c("A", "B"), ncol = 2, nrow = 1)

# Validation Clustering 
val_clust_plot <- plot_grid(c,d,  labels = c("A", "B"), ncol = 2, nrow = 1)

# saving plots
ggsave("ul plot.tiff", ul_plot, dpi = 350, height = 86, width = 180, limitsize = F, units = "mm") 

ggsave("val clust.tiff", val_clust_plot, dpi = 350, height = 86, width = 180, limitsize = F, units = "mm") 

ggsave("boxplots2.tiff", pp, dpi = 350, height = 86, width = 180, limitsize = F, units = "mm") 

############################
# save outputs
############################

fwrite(dsr, "ds_raw.csv", row.names = T)

fwrite(ds, "ds_uvf.csv", row.names = T)

fwrite(ds_d, "ds_ranking.csv", row.names = T)

fwrite(ds_vscc, "ds_fs_by_vscc.csv", row.names = T)

############################
# Signal drift correction
############################

# EigenMS 

setwd("C:/...")
data <-as.data.frame(fread(input = ".csv", header=T))
row.names(data) <- data[,1]
data <- data[,-1]
colnames(data)[1] <- c("Label")

m_logInts = as.data.frame(t(data[,-1]))
m_logInts = convert_log2(m_logInts)
grps = as.factor(c(data[,1]))
m_prot.info = cbind(colnames(data)[-1], colnames(data)[-1])
m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
m_ints_norm1 = eig_norm2(rv=m_ints_eig1)

data_em <- data
data_em[,-1] <- as.data.frame(t(m_ints_norm1$norm_m))
data_em <- cbind(rownames(data), data_em)
fwrite(data_em, ".csv")

# Quantile

setwd("C:/...")
data <-as.data.frame(fread(input = ".csv", header=T))
row.names(data) <- data[,1]
data <- data[,-1]
colnames(data)[1] <- c("Label")

data_quantile <- as.matrix(t(data[,-1]))

QUANTILE <- function(data) {
  normalize.quantile <- get("normalize.quantiles", en = asNamespace("affy"))
  quantile.data <- normalize.quantile(data)
  rownames(quantile.data) <- rownames(data)
  return(quantile.data)
}

quantile_data <- as.data.frame(QUANTILE(data_quantile))

data_q <- data
data_q[,-1] <- as.data.frame(t(quantile_data))
data_q <- cbind(rownames(data), data_q)
fwrite(data_q, ".csv")

# Cybic Spline

setwd("C:/...")
data <-as.data.frame(fread(input = ".csv", header=T))
row.names(data) <- data[,1]
data <- data[,-1]
colnames(data)[1] <- c("Label")

data_cubic_spline <- as.matrix(t(data[,-1]))

CUBIC <- function(data) {
  spline.data <- normalize.qspline(data, samples = 0.02, 
                                   target = apply(data, 1, mean))
  rownames(spline.data) <- rownames(data)
  return(spline.data)
}

cubic_spline_data <- as.data.frame(QUANTILE(data_cubic_spline))

data_cs <- data
data_cs[,-1] <- as.data.frame(t(cubic_spline_data))
data_cs <- cbind(rownames(data), data_cs)
fwrite(data_cs, ".csv")

# VSN normalization

setwd("C:/...")
data <-as.data.frame(fread(input = ".csv", header=T))
row.names(data) <- data[,1]
data <- data[,-1]
colnames(data)[1] <- c("Label")

data_vsn <- as.matrix(t(data[,-1]))

VSN <- function(data) {
  vsn.model <- vsn2(data)
  vsn.data <- predict(vsn.model, data)
  return(vsn.data)
}

vsn_data <- as.data.frame(VSN(data_vsn))

data_v <- data
data_v[,-1] <- as.data.frame(t(vsn_data))
data_v <- cbind(rownames(data), data_v)
fwrite(data_v, ".csv")