###This is a code to reproduce statistical analysis of proteomics data from Proteomic profiling of carotid atherosclerotic plaques and adjacent intact arterial segments

##Openinig_the_data
#We reccomed to set the woring directory to make easy to reproduce the code
dat <- data.frame(read.table("report.pg_matrix.tsv", sep = '\t', header = TRUE))

head(dat)
str(dat)
dat1 <- dat
rownames(dat1) <- dat[,1]
dat1 <- dat1[,c(6:45)]
head(dat1)

#Opening the sample_info
library(readxl)
fact <- data.frame(read_excel("sample_info.xlsx"))

rownames(fact) <- fact[,1]
fact <- fact[,-1]
fact$Group <- as.factor(fact$Group)
fact$Group

fact$Donor <- as.factor(as.character(fact$Donor))
fact$Donor

fact$Replicate <- as.factor(as.character(fact$Replicate))
fact$Replicate

colnames(dat1) <- rownames(fact)

colSums(is.na(dat1))

dat2 <- dat1[ , -which(names(dat1) %in% c("62_1", "62_2", "65_1", "65_2", "74_1", "74_2", "81_1", "81_2", "63_1", "63_2", "80_1", "80_2"))]

colSums(is.na(dat2))

fact <- fact[ -which(rownames(fact) %in% c("62_1", "62_2", "65_1", "65_2", "74_1", "74_2", "81_1", "81_2", "63_1", "63_2", "80_1", "80_2")), ]

## Qualitative analysis
datq <- dat2
datq[is.na(datq)] <- 0
head(datq)
datq1 <- data.frame(datq[,"64_1"]+datq[,"64_2"], datq[,"66_1"]+datq[,"66_2"], datq[,"67_1"]+datq[,"67_2"], datq[,"68_1"]+datq[,"68_2"], datq[,"69_1"]+datq[,"69_2"], datq[,"70_1"]+datq[,"70_2"], datq[,"71_1"]+datq[,"71_2"], datq[,"72_1"]+datq[,"72_2"], datq[,"73_1"]+datq[,"73_2"], datq[,"75_1"]+datq[,"75_2"], datq[,"76_1"]+datq[,"76_2"], datq[,"77_1"]+datq[,"77_2"], datq[,"78_1"]+datq[,"78_2"], datq[,"79_1"]+datq[,"79_2"])
colnames(datq1) <- c("64_1", "66_1","67_1","68_1","69_1","70_1","71_1","72_1","73_1","75_1","76_1","77_1","78_1","79_1")
rownames(datq1) <- rownames (dat2)

factq <- fact[colnames(datq1),]

datq1[datq1 == 0] <- NA

venn_Plaq <- datq1[which(rowMeans(!is.na(datq1[, rownames(subset(factq,Group=="Intact"))])) >= 6/7), ]
venn_Intact <- datq1[which(rowMeans(!is.na(datq1[, rownames(subset(factq,Group=="Plaque"))])) >= 6/7), ]

#Venn diagramm
library(gplots)
v.table1 <- venn(list(rownames(venn_Plaq), rownames(venn_Intact)))
print(v.table1)



## Quantitative analysis
#Removing rows with a lot of missing values
dat3 <- dat2[which(rowMeans(!is.na(dat2)) >= 24/28), ]
mean(complete.cases(dat3))
colSums(is.na(dat3))

#Raw data
library(RColorBrewer)

#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$Group]
boxplot(dat3, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$Group), fill = pal, bty = "n", xpd = T)

#knn imputation of missng values
library(impute)
dat_knn <- impute.knn(t(dat3), k = 5)
dat_knn <- t(dat_knn$data)
mean(complete.cases(dat_knn))

#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
boxplot(dat_knn, outline = FALSE, col = cols, main = "Data after missed values impuration")
legend("topright", levels(fact$Group), fill = pal, bty = "n", xpd = T)

#dev.off()
colSums(dat_knn)
head(dat_knn)

dat_log <- log2(dat_knn+1)

#Quantile normalizatiomn
library(limma)
data_norm <- normalizeQuantiles(dat_log)
boxplot(data_norm, outline = FALSE, col = cols, main = "Normalized data")
legend("topright", levels(fact$Group), fill = pal, bty = "n", xpd = T)



#MA-plot
maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  scatter.smooth(x = X, y = Y,
                 main = main, pch = pch,
                 xlab = xlab, ylab = ylab,
                 lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

maplot(dat_log[, rownames(subset(fact,Group=="Intact"))], dat_log[, rownames(subset(fact,Group=="Plaque"))], main = "Log data")
maplot(data_norm[, rownames(subset(fact,Group=="Intact"))], data_norm[, rownames(subset(fact,Group=="Plaque"))], main = "Normalized data")


library(mixOmics)
dat_pca <- pca(t(data_norm), ncomp = 8, center = TRUE)

#tiff('PCA_group.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')
plotIndiv(dat_pca, comp = c(1, 2), ind.names = T, 
          group = fact$Donor, legend = TRUE, ellipse = F,
          title = 'PCA')
#dev.off()

library(mixOmics)
#PLS-DA 
ordination.optimum.splsda <- splsda(t(data_norm), fact$Group, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('PLSDA_all.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
#dev.off()
layout(1,1)



#Limma
library(limma)

X <- model.matrix(~ fact$Group)
X

fit <- lmFit(data_norm, design = X, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit <- eBayes(fit)

# Dif_expr_table
topTable(efit, coef = 2)
full_list <- topTable(efit, coef = 2, number = length(data_norm[,2]))
#write.csv(full_list,'Dif_expr_int_vs_plaq.csv')
head(full_list)


#Vulcano
library(EnhancedVolcano)

#tiff('Dif_expr.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(full_list,
                lab = rownames(full_list),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-2, 2.5),
                ylim = c(0, 11),
                title ="Intact vs plaque",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


p_aboveS <- full_list$adj.P.Val <= 0.05
sum(p_aboveS)

head(full_list)
#proteins
boxplot(dat_knn[c("P27169"),] ~ Group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knn[c("P04217"),] ~ Group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knn[c("P05164"),] ~ Group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knn[c("P02763"),] ~ Group, data = fact,
        varwidth = TRUE, log = "y", las = 1)
