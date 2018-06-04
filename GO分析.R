setwd('E:\\Rdata')
CNA <- read.table('Gene2.txt')
GeneID <- as.factor(CNA[,1])


source("https://bioconductor.org/biocLite.R") 
biocLite('org.Hs.eg.db')
biocLite("clusterProfiler")
geneIDannotation <- function(geneLists=c(1,2,9),name=T,map=T,ensemble=F,accnum=F){
  ## input ID type : So far I just accept entrezID or symbol
  ## default, we will annotate the entrezID and symbol, chromosone location and gene name 
  
  suppressMessages(library("org.Hs.eg.db"))
  all_EG=mappedkeys(org.Hs.egSYMBOL) 
  EG2Symbol=toTable(org.Hs.egSYMBOL)
  if( all(! geneLists %in% all_EG) ){
    inputType='symbol'
    geneLists=data.frame(symbol=geneLists)
    results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
  }else{
    inputType='entrezID'
    geneLists=data.frame(gene_id=geneLists)
    results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
  }
  
  if ( name ){
    EG2name=toTable(org.Hs.egGENENAME)
    results=merge(results,EG2name,by='gene_id',all.x=T)
  }
  if(map){
    EG2MAP=toTable(org.Hs.egMAP)
    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  if(ensemble){
    EG2ENSEMBL=toTable(org.Hs.egENSEMBL)
    results=merge(results,EG2ENSEMBL,by='gene_id',all.x=T)
  }
  if(accnum){
    EG2accnum=toTable(org.Hs.egREFSEQ) 
    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  return(results)
}
geneIDannotation()
gene <- geneIDannotation(GeneID)






library(clusterProfiler)
#没有organism="human"，改为OrgDb=org.Hs.eg.db,就不会有问题啦~~~
ego_cc <- enrichGO(gene = gene$gene_id,
                   OrgDb=org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   readable = TRUE)

write.table(as.data.frame(ego_cc@result), file="test_CC.txt")
#KEGG又不一�?,它就可以�?
kk <- enrichKEGG(gene = gene$gene_id, 
                 organism = "human",
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)
write.table(as.data.frame(kk@result), file="test_kk.txt")
#具体的参数在R中？或者help()，讲的巨清楚�?
#结果就是[url=http://www.bio-info-trainee.com/370.html]http://www.bio-info-trainee.com/370.html[/url]那样的�?
#但是当excel 打开有点错行，（我的错行）。改一改第一行不麻烦�?
#画图�?
barplot(ego_cc, showCategory=15,title="EnrichmentGO_cc")#条状图，按p从小到大排的
dotplot(ego_BP,title="EnrichmentGO_CC_dot")#点图，按富集的数从大到小�?
