# 差异分析走标准的limma流程 ---------------------------------------------------------
rm(list = ls())  ## 魔幻操作，一键清空~
library(limma)
options(stringsAsFactors = F)
load(file = 'dex_brg1.Rdata')

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

tempOutput$g=ifelse(tempOutput$P.Value>0.01,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
                    ifelse( tempOutput$logFC>1.2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                            ifelse( tempOutput$logFC < -1.2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)
table(tempOutput$g)

tem<-tempOutput[!tempOutput$g=='stable',]
save(df1,group,tem,tempOutput,file = '差异基因.Rdata')



# 富集分析 --------------------------------------------------------------------

library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
library(enrichplot)

rm(list = ls()) 
options(stringsAsFactors = F)
load(file = '差异基因.Rdata')

DEG<-tempOutput
DEG<-DEG[,1:6]
# colnames(DEG)[c(2,5)]<-c('logFC','P.Value')
head(DEG)

## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=1
DEG$gene=ifelse(DEG$P.Value>0.01,'stable',
                ifelse( DEG$logFC> logFC_t,'UP',
                          ifelse( DEG$logFC< -logFC_t,'DOWN','stable') )
)#P>0.05输出stable，其中设定当logFC大于1.5为上调输出'UP'，大于-1.5为下调输出'DOWN',，如果都不是则输出'stable'，从而增加了一列g，筛选出了上调和下调的基因
table(DEG$gene)

head(DEG)



# id转换
DEG$symbol=rownames(DEG)
df <- bitr(unique(DEG$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
#bitr功能为ID转换，
#bitr(geneID, fromType, toType, OrgDb, drop = TRUE)；
#geneid ：基因ID输入 ； fromtype ： 输入ID型；toType：输出ID型；orgdb ：注释数据库）
head(df)
deg<-DEG
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
table(DEG$gene)
#把数据DEG,df通过，DEG的'symbol'列，df的'SYMBOL'列连接在一起，转化ID
head(DEG)

gene_up= DEG[DEG$gene== 'UP','ENTREZID'] #选出上调基因ID
gene_down=DEG[DEG$gene== 'DOWN','ENTREZID'] #选出下调基因ID
gene_diff=c(gene_up,gene_down)#得出上下调基因ID
gene_all=as.character(DEG[ ,'ENTREZID'] )#得出所有基因ID

geneList=DEG$logFC#把 DEG 数据logFC列值赋值给数据geneList
names(geneList)=DEG$ENTREZID#把ID赋值给geneList数据的名字
geneList=sort(geneList,decreasing = T)#把数据进行排序
head(geneList)



# 富集画图 --------------------------------------------------------------------
# KEGG、GSEA 分析
kegg_plot <- function(diff_kegg){
  dat<-diff_kegg
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  dat=dat[order(dat$pvalue,decreasing =F),]
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing =T)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+
    theme(text = element_text(size=8),plot.title = element_text(hjust = 0.1))+
    ggtitle("Pathway Enrichment") 
} 
if(T){
  ###   over-representation test
  # KEGG
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff =0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  dotplot(kk.up )
  # ggsave('kk.up.dotplot.png')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = "hsa",
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  dotplot(kk.down)
  # ggsave('kk.down.dotplot.png')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.diff)[,1:6]
  dotplot(kk.diff )
  # ggsave('kk.diff.dotplot.png')
  
  
  kegg_down_dt <- as.data.frame(kk.down)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  print(kegg_plot(down_kegg))
  
  kegg_up_dt <- as.data.frame(kk.up)
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  print(kegg_plot(up_kegg))
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  diff_kegg<-kegg_diff_dt[kegg_diff_dt$pvalue<0.05,];diff_kegg$group=1
  print(kegg_plot(diff_kegg))
  # ggsave(g_kegg,filename = 'diff_kegg.png')
} 
# GSEA
kk_gse <- gseKEGG(geneList          = geneList,
                  organism          = 'hsa',
                  use_internal_data = FALSE,
                  nPerm             = 1000,
                  minGSSize         = 120,
                  pvalueCutoff      = 1,
                  verbose           = FALSE)
kk_gse[,2]

up<-which(kk_gse$pvalue<0.05 & kk_gse$enrichmentScore> 0)
down<-which(kk_gse$pvalue<0.05 & kk_gse$enrichmentScore< 0)
diff_gsea_kegg<-kk_gse[c(up,down),];diff_gsea_kegg$group=c(rep(1,length(up)),rep(-1,length(down)))
print(kegg_plot(diff_gsea_kegg)) # 富集

# gseaplot(kk_gse, geneSetID = rownames(kk_gse[16,]))
gseaplot2(kk_gse, geneSetID = rownames(kk_gse[c(6,34,44,87,74),]),pvalue_table=T,title = 'BRG1_CON', )
  

# go 富集

g_list=list(gene_up=gene_up,gene_down=gene_down,gene_diff=gene_diff)

if(T){
  go_enrich_results <- lapply( g_list , function(gene) {
    lapply( c('BP','MF','CC') , function(ont) {
      cat(paste('Now process ',ont ))
      ego <- enrichGO(gene          = gene,
                      universe      = gene_all,
                      OrgDb         = org.Hs.eg.db,
                      ont           = ont ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.99,
                      qvalueCutoff  = 0.99,
                      readable      = TRUE)
      
      # print( head(ego) )
      return(ego)
    })
  })
  save(go_enrich_results,file = 'go_enrich_results.Rdata')
}

load(file = 'go_enrich_results.Rdata')
n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
# 画图 
for (i in 1:3){
  for (j in 1:3){
    p<- dotplot(go_enrich_results[[i]][[j]])+
      ggtitle(paste0('dotplot_',n1[i],'_',n2[j]))
    print(p)
  } 
}

