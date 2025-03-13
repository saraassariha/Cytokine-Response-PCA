df=read.table("data/cytokine_responses.csv",header=TRUE, sep=";",dec=",")
df[1:4,1:4]
FC=df


# Calculate fold change
for (i in 2:11){
  FC[,i]=df[,i]/df$Control
}


# Compute lo2 of the FC of the drugs
log2_FC=cbind(df$cytokines,log2(FC[,3:11]))
colnames(log2_FC)[1]="cytokine"

#box plot cyclosporineA
boxplot(cyclosporineA  ~ cytokine,data=log2_FC,ylab="Log2FC",main="cyclosporineA",border = 1:3,col=rainbow(9),las = 2)

#box plot tacrolimus
boxplot(tacrolimus   ~ cytokine,data=log2_FC,ylab="Log2FC",main="Tacrolimus",border = 1:3,col=rainbow(9),las = 2)
#box plot tofacitinib
boxplot(tofacitinib   ~ cytokine,data=log2_FC,ylab="Log2FC",main="tofacitinib",border = 1:3,col=rainbow(9),las = 2)
#box plot ruxolitinib
boxplot(ruxolitinib   ~ cytokine,data=log2_FC,ylab="Log2FC",main="Ruxolitinib",border = 1:3,col=rainbow(9),las = 2)

#Q3
#correlation 
cor(log2_FC[ ,2:10],method="pearson")

#Bonferroni corrected p- values
library(Hmisc)
pmin(rcorr(as.matrix(log2_FC[ ,2:10]), type="pearson")$P*36, 1)

#rounding 
round(pmin(rcorr(as.matrix(log2_FC[ ,2:10]), type="pearson")$P*36, 1),2)

#selecting NA in the second half
library(corrplot)
cor_mat=round(pmin(rcorr(as.matrix(log2_FC[ ,2:10]), type="pearson")$P*36, 1),2)
corrplot(cor_mat,type = "lower")
cor_mat=round(pmin(rcorr(as.matrix(log2_FC[ ,2:10]), type="pearson")$P*36, 1),2)
diag(cor_mat) = NA
cor_mat[upper.tri (cor_mat)] = NA
cor_mat

#comparing-not-needed
library(reshape)
cor_melt=melt(cor_mat)
na.omit(cor_melt)

#seperating the p-value lower than 0.05
index=which(cor_melt$value<0.050 )
most_cor=cor_melt[index,] # Extract correlations that are smaller than 0.05 
most_cor[order(most_cor$value,decreasing = TRUE),]# Print sorted correlations

#Q4
#scatter plot 
library(corrplot)
cor_mat=cor(log2_FC[ ,2:10],method="pearson")
corrplot(cor_mat,type = "lower")

plot(log2_FC[,10],log2_FC[,7] ,xlab="tofacitinib",ylab="ruxolitinib")
plot(log2_FC$tofacitinib,log2_FC$ruxolitinib)
head(log2_FC)
#Q5
#PCA
pca=prcomp(log2_FC[,2:10], scale = TRUE)
summary(pca)
var=as.numeric(pca$sdev^2)
var[1]/sum(var)
var[2]/sum(var)

#draw scatter plot 
library(ggplot2)
library(scales)
prop.pca = pca$sdev^2/sum(pca$sdev^2) # % Variance
dataset = data.frame(pca$x) # Create a data frame of the PCs
p1 = ggplot(dataset) + geom_point(aes(PC1,PC2), size = 3.5) +
  geom_text(aes(PC1, PC2,label=1:length(PC1)),size=2.2,col="white") +
  labs(x = paste("PC1 (", percent(prop.pca[1]), ")", sep=""),
       y = paste("PC2 (", percent(prop.pca[2]), ")", sep=""))
p1


#generating eigenvalue
eigenvalues=pca$sdev^2
print(eigenvalues)
mean(pca$sdev^2)

#selecting the values >1
eigenvalues_greater_than_1=eigenvalues[eigenvalues > 1]
number_of_component=length(eigenvalues_greater_than_1)
number_of_component

#calculating the ptoportion by hand 
explained_variance=eigenvalues_greater_than_1 / sum(eigenvalues)
cumulative_variance=cumsum(explained_variance)
cumulative_variance

#plot the values
plot(pca,type="l")

#Q6
library(psych)
varmax=principal(log2_FC[,2:10], nfactors = 2, rotate = "varimax")
varmax

#Q7
ggplot(log2_FC, aes(x = varmax$scores[,2], y = tacrolimus)) +
  geom_point(alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Scatter Plot: Rotated Component 2 Scores vs. Log2 FC for Tacrolimus",
       x = "Rotated Component 2 Scores",
       y = "Log2 Fold Change (Tacrolimus)") +
  theme_minimal()

# Scatter plot between rotated component 2 scores and log2 FC for tacrolimus
plot(varmax$scores[,2], log2_FC$tacrolimus, 
     xlab = "Rotated Component 2 Scores", 
     ylab = "log2 FC for Tacrolimus", 
     main = "Scatter Plot: PC2 vs log2 FC Tacrolimus",
     col = "black", pch = 16)

# Compute and print correlation
correlation <- cor(varmax$scores[,2], log2_FC$tacrolimus)
print(correlation)



# Independent simulated MANOVA example (not the cytokine dataset) 
set.seed(890604)  # Example: 890123
Procalcitonin=rnorm(21,20,2)
CRP=c(rnorm(14,30,4),rnorm(7,37,4))
Group=factor(c(rep("Virus",7),rep("Bacteria",7),rep("Malaria",7)),levels=c("Virus","Bacteria","Malaria"))
df=data.frame(Procalcitonin,CRP,Group)

print(df[1:6,]) 

#boxplot CRP
boxplot(CRP ~ Group,data=df,ylab="CRP",border = 1:3)

#jittering CRP
points(x = jitter(as.numeric(df$Group)), y = df$CRP,col=as.numeric(df$Group))


#boxplot Procalcitonin 
boxplot(Procalcitonin  ~ Group,data=df,ylab="Procalcitonin",border = 1:3)

#jittering Procalcitonin 
points(x = jitter(as.numeric(df$Group)), y = df$Procalcitonin ,col=as.numeric(df$Group))

#MANOVA
fit = manova(cbind(CRP,Procalcitonin) ~ Group,data=df)
summary(fit)


#post hoc
TukeyHSD(aov(CRP ~ Group, data = df))
TukeyHSD(aov(Procalcitonin ~ Group, data = df))



write.csv(pca$x, "results/pca_scores.csv")
write.csv(data.frame(component=seq_along(pca$sdev), variance=pca$sdev^2/sum(pca$sdev^2)), "results/explained_variance.csv", row.names=FALSE)
ggsave("results/pca.png", p1, width=7, height=5, dpi=150)
