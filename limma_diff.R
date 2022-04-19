

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

# 差异分析走标准的limma流程 ---------------------------------------------------------
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'DEG.Rdata')

#创建一个分组的矩阵
design=model.matrix(~0+factor(group))#创建一个分组的矩阵
colnames(design)=levels(factor(group))
rownames(design)=colnames(df1)
head(design)
dematrix<-as.data.frame(design)
table(dematrix$ko_dex)

# 创建差异比较矩阵 ----------------------------------------------------------------

contrast.matrix<-makeContrasts(ko_dex-wt_con,levels = design)
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

tempOutput$gene=ifelse(tempOutput$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
                    ifelse( tempOutput$logFC>1.2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                            ifelse( tempOutput$logFC < -1.2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)
table(tempOutput$gene)

tem<-tempOutput[!tempOutput$gene=='stable',]
save(rdat,group_list,tem,tempOutput,file = '差异基因.Rdata')

