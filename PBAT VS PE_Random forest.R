rm(list=ls(all=TRUE))
library(vegan)
library(dplyr)
library(ggplot2)
setwd("D:/Desktop/randomforest_unoise_re8/AE_figure")
otu <- read.table("Galaxy98-[resample_Remove8_16s_ALL_UNOISE.txt].txt", sep="\t", header=T, row.names=1)
tre <- read.table("16S_PBAT_PE.txt",sep="\t", header=T, row.names=1)
classifier <- read.table("Galaxy97-[Classifier_of_16srrna.txt].txt",sep="\t", header=T, row.names=1)


sub_otu = otu[, rownames(tre)]
rowS = rowSums(sub_otu>0)
valid.row = which(rowS>0)
sub_otu = sub_otu[valid.row, ]
sub_otu

library(randomForest)
#tre$group <- factor(tre$group, levels=tre$group)
set.seed(315)
rf = randomForest(t(sub_otu), tre$group, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)
names(rf)
plot(rf)

varImpPlot(rf, main = "Feature OTU importance",n.var = 21, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray" )
## 交叉验证选择Features
set.seed(315) # 随机数据保证结果可重复，必须
# rfcv是随机森林交叉验证函数：Random Forest Cross Validation
result = rfcv(t(sub_otu), tre$group, cv.fold=100)
error.cv = data.frame(num = result$n.var, error.1 =  result$error.cv)
for (i in 316:(319)){
  print(i)
  set.seed(i)
  result= rfcv(t(sub_otu), tre$group, cv.fold=5)
  error.cv = cbind(error.cv, result$error.cv)
}

n.var = error.cv$num
error.cv = error.cv[,2:6]
colnames(error.cv) = paste('err',1:5,sep='.')
err.mean = apply(error.cv,1,mean)
allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
# number of features selected
optimal = 21

write.table(allerr, file = "Order_rfcv.txt", sep = "\t", quote = F, row.names = T, col.names = T)

p = ggplot() + 
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') + 
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + 
  #  geom_hline(yintercept = min(allerr$err.mean), colour='black', lwd=0.36, linetype="dashed") + 
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 2, 5, 20, 50, 200, 600)) + 
  labs(title=paste('Training set(n= ', dim(sub_otu)[1],')', sep = ''), 
       x='number of Orders', y=' error rates of Cross-validations') + 
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep="")) + 
  theme_classic()
ggsave(p, file = "Order_rfcv.pdf", width = 99, height = 70, unit = 'mm')
p

imp <- importance(rf)
imp <- data.frame(predictors = rownames(imp), imp)
write.table(imp,file = "randomforest_feature.txt",quote = F,sep = '\t', row.names = F, col.names = T)

imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.21 <- imp.sort[1:21, ]
imp.21$predictors <- factor(imp.21$predictors, order=T, levels = rev(imp.21$predictors))
sub_classifier = classifier[rownames(imp.21),]
imp.21 = cbind(imp.21,sub_classifier)
write.table(imp.21, file =  "order_注释_feature.txt",quote = F,sep = '\t', row.names = F, col.names = T)

p2 = ggplot(imp.21, aes(x = predictors, y = MeanDecreaseGini, fill = Order)) +
  geom_bar(stat = "identity") + 
  coord_flip() + theme_bw() +
  ggtitle("Most important OTUs for classifying samples") +
  ylab("Mean of decreased accuracy")+xlab(NULL) + labs(fill = "Order")

p2 = p2 +  scale_x_discrete(label = rev(imp.21$Genus))

ggsave(paste("rf_imp_feature",".pdf", sep=""), p2, width = 6, height =4)


# 图4.2. 绘制时间序列热图

# 加载热图绘制包
library(pheatmap)

# 数据筛选23个feature展示
sub_abu = sub_otu[rownames(imp.21),]

# 简化名字
#rownames(sub_abu)=imp[rownames(sub_abu),"Famliy"]

# 直接自动聚类出图
pheatmap(sub_abu, scale = "row")
# 保存结果
pheatmap(sub_abu, scale = "row",angle_col = 0, filename = "heatmap_samples.pdf", width = 6, height = 8)


# 按膜为组合并均值
tre$names= row.names(tre)
#sampFile = data.frame(tre,row.names = row.names(tre))
#colnames(sampFile)[1] = "group"
mat_t =as.data.frame( t(sub_abu))
mat_t$names = row.names(mat_t)
mat_t2 = merge(tre, mat_t, by= "names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = mat_mean$group
pheatmap(df,fontsize = 20,fontface="italic",fontfamily= "Arial",otu_norm_group,scale="row",cluster_cols = F, cluster_rows = F)
pheatmap(otu_norm_group, scale="row",cluster_cols = F,angle_col = 0, cluster_rows = F, filename = "heatmap_films.pdf", width = 3, height = 6)


# 按时间为组合并PBAT均值
group <- read.table("16S_PBAT_group5.txt",sep="\t", header=T, row.names=1)
group$names= row.names(group)
sampFile = data.frame(group,row.names = row.names(group))
#colnames(sampFile)[1] = "group"
mat_t =as.data.frame( t(sub_abu))
mat_t$names = row.names(mat_t)
mat_t2 = merge(group, mat_t, by= "names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = c(1:5)
rownames(otu_norm_group) = sub_classifier$Genus
pheatmap(df,fontsize = 20,fontface="italic",fontfamily= "Arial", otu_norm_group,scale="row",cluster_cols = F, cluster_rows = F)
pheatmap(df,fontsize = 20,fontface="italic",fontfamily= "Arial", otu_norm_group, scale="row",cluster_cols = F, cluster_rows = F, angle_col = 0, filename = "heatmap_PBAT.pdf", width = 5, height = 8)

# 按时间为组合并PE均值
group <- read.table("16S_PE_group5.txt",sep="\t", header=T, row.names=1)
group$names= row.names(group)
sampFile = data.frame(group,row.names = row.names(group))
#colnames(sampFile)[1] = "group"
mat_t =as.data.frame( t(sub_abu))
mat_t$names = row.names(mat_t)
mat_t2 = merge(group, mat_t, by= "names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = c(1:5)
rownames(otu_norm_group) = sub_classifier$Genus
pheatmap(df,fontsize = 20,fontface="italic",fontfamily= "Arial", otu_norm_group,scale="row",cluster_cols = F, cluster_rows = F)
pheatmap(df,fontsize = 20,fontface="italic",fontfamily= "Arial", otu_norm_group, scale="row",cluster_cols = F, cluster_rows = F, angle_col = 0, filename = "heatmap_PE.pdf", width = 5, height = 8)

