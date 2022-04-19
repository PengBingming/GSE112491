
rm(list = ls())
options(stringsAsFactors = F,scipen = 200)

library(limma)
# 数据初步整理 ------------------------------------------------------------------

df1<-read.csv('GSE112491_A1A3_RNAseq_log2normalized_geneCounts.CSV',header = T)


lis<-strsplit(df1$GeneName,'[.]')
for (i in 1:length(lis)) {
  df1$symbol[i]<- lis[[i]][1]
};rm(i)
df1<- avereps(df1[,c(-1,-14)],ID=df1$symbol)

df1<-df1[,colSums(df1)>0]
df1<-df1[rowSums(df1)>0,]
df1<-as.matrix(df1)
group<-c(rep('wt_con',3),rep('wt_dex',3),rep('ko_con',3),rep('ko_dex',3))

# df1<-df1[,1:6];group<-group[1:6] # wt_dex-wt_con
# df1<-df1[,c(1:3,7:9)];group<-group[c(1:3,7:9)] # ko_con-wt_con
# df1<-df1[,c(4:6,10:12)];group<-group[c(4:6,10:12)] # ko_dex-wt_dex
# df1<-df1[,7:12];group<-group[7:12] # ko_dex-ko_con
# df1<-df1[,c(1:3,10:12)];group<-group[c(1:3,10:12)] # ko_dex-wt_con
# df1<-df1[,c(4:9)];group<-group[c(4:9)] # ko_con-wt_dex

save(df1,group,file='dex_brg1.Rdata')

# PCA主成分分析 ----------------------------------------------------------------
rm(list = ls())

load(file = 'dex_brg1.Rdata')
df1[1:4,1:4]#每次都要检测数据
df1<-t(df1)#画PCA图时要求是行名是样本名，列名是探针名，因此此时需要转换 t()
df1<-as.data.frame(df1)#将matrix转换为data.frame
df1<-cbind(df1,group) #cbind横向追加，即将分组信息追加到最后一列

library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
library(ggplot2)
dat.pca <- PCA(df1[,-ncol(df1)], graph = FALSE)
head(dat.pca$ind$coord)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = df1$group, # color by groups
             
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
# ggsave('Expreset1_PCA.png') 


# 层次聚类分析 ------------------------------------------------------------------

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)

load(file = 'dex_brg1.Rdata')

df1[1:4,1:4]#每次都要检测数据

df1<-t(df1)#层次聚类分析要求行名是样本名，列名是探针名，因此此时需要转换
df1[1:4,1:4]

d<-dist(df1) # 欧几里得距离
fit.complete<-hclust(d,method="complete")
plot(fit.complete,hang = -1,cex=0.8)

# 基因的差异性分析 ----------------------------------------------------------------

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)

load(file = 'dex_brg1.Rdata')

df1[1:4,1:4]#每次都要检测数据
boxplot(df1["SMARCA4",]~group,ylab = 'Smarca4') #samrca4以group分组画箱线图
boxplot(df1["ACTB",]~group,ylab = 'ACTB')

# 无意义
boxplot((df1["SMARCA4",]/df1["ACTB",])~group,ylab = 'SMARCA4 / ACTB')

#install.packages("ggpubr")  
bp=function(gene,name,tissue){         #自定义一个函数bp，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=gene,stage=group)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter",
                 ylab = name,xlab = tissue)
  #  Add p-value
  p + stat_compare_means()
}
## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(df1['SMARCA4',],'Smarca4','cell') 


# 差异分析走标准的limma流程 ---------------------------------------------------------
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'dex_brg1.Rdata')

library(limma)
#创建一个分组的矩阵

design=model.matrix(~0+factor(group))#创建一个分组的矩阵
colnames(design)=levels(factor(group))
rownames(design)=colnames(df1)
head(design)
dematrix<-as.data.frame(design)
table(dematrix$ko_con)

# 创建差异比较矩阵 ----------------------------------------------------------------

contrast.matrix<-makeContrasts(ko_con-wt_con,levels = design)
contrast.matrix ##这个矩阵声明，我们要Tumor组和Non_Tumor组进行差异分析比较

#第一步lmFit，#lmFit为每个基因给定一系列的阵列来拟合线性模型
fit<-lmFit(df1,design)
#第二步eBayes，#eBayes给出了一个微阵列线性模型拟合，通过经验贝叶斯调整标准误差到一个共同的值来计算修正后的t统计量、修正后的f统计量和微分表达式的对数概率。
fit1<-contrasts.fit(fit, contrast.matrix)
fit1<-eBayes(fit1)
#第三步topTable,#topTable从线性模型拟合中提取出排名靠前的基因表。
options(digits = 4) #设置全局的数字有效位数为4
#topTable(fit1,coef=2,adjust='BH') 
tempOutput<-topTable(fit1,coef=1,n=Inf) 
tempOutput<-na.omit(tempOutput)#移除NA值
head(tempOutput)
tempOutput["SMARCA4",]



# 上调基因和下调基因 ---------------------------------------------------------------

tempOutput$g=ifelse(tempOutput$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
                    ifelse( tempOutput$logFC>1.2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                            ifelse( tempOutput$logFC < -1.2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)
table(tempOutput$g)

tem<-tempOutput[!tempOutput$g=='stable',]
save(df1,group,tem,tempOutput,file = '差异基因.Rdata')

# 差异基因火山图 -----------------------------------------------------------------

rm(list = ls())
library(ggplot2)
load(file = '差异基因.Rdata')
DEG=tempOutput
DEG<-na.omit(DEG)
colnames(DEG)[7]<-'change'

# 取前两位小数
logFC_cutoff <- 1.2
# 确定上下调表达基因

table(DEG$change)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,2),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='up',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='down',]))

g2= ggplot(data=DEG,
           aes(x=logFC, y=-log10(P.Value),colour=change)) + 
  geom_point(shape = 16, alpha=0.8, size=1)+ 
  ggtitle(this_tile)+ theme(plot.title = element_text(size=9, hjust = 0.5))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black" ),
        plot.margin = unit(rep(4,4),'lines'), axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  scale_colour_manual(values = c('blue3', 'black', 'forestgreen', 'red3')) + #corresponding to the level)
  geom_hline(yintercept = log10(19.95262), linetype = "dashed", color ='grey', size = 0.5)+
  geom_vline(xintercept = logFC_cutoff, linetype = "dashed", color ='grey',size = 0.5)+
  geom_vline(xintercept = -logFC_cutoff, linetype = "dashed", color ='grey',size = 0.5)

print(g2) 



# for heatmap 热图 ----------------------------------------------------------

rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = '差异基因.Rdata')
#每次都要检测数据
df1[1:4,1:4]
table(group)
#挑选出差异最大100个基因和差异最小的100个基因
x=tem$logFC #deg取logFC这列并将其重新赋值给x
names(x)=row.names(tem)#将基因名作为名字，命名给给x
is.vector(x)
x[1:4]
#对x进行从小到大排列，取前100及后100，并取其对应的基因名，作为向量赋值给cg
cg=c(names(head(sort(x),50)),names(tail(sort(x),50)))

library(pheatmap)

#对dat按照cg取行，所得到的矩阵来画热图
head(df1[cg,])
annotation_col<-as.data.frame(matrix(0,ncol = 2,nrow = 3))

annotation_col = data.frame(
  group= factor(group), 
  rep= 1:3)
rownames(annotation_col)<-colnames(df1)

pheatmap(df1[cg,],annotation_col = annotation_col,scale = 'row',cluster_cols = F,show_colnames =F,show_rownames = T) 


# 输出数据 --------------------------------------------------------------------

write.csv(tem,'wt_dex-ko_con差异基因.csv')
write.csv(tempOutput,'wt_dex-ko_con全部基因.csv')

