rm(list = ls())
options(stringsAsFactors = F)
setwd(dir = 'TB/')

library(tidyverse)
library(Seurat)
library(CellChat)
exprset <- read.table('matrix_X.csv.gz',sep = ',',header = T,row.names = 1,check.names = F)
metadata <- read.table('adata.obs.csv',sep = ',',header = T,row.names = 1,check.names = F)
sce <- CreateSeuratObject(exprset,project = 'dis',assay = 'RNA',meta.data = metadata)
sce <- NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 1e4)
saveRDS(sce,'normalize.sce.cellchat.RDS')

data.input <- sce@assays$RNA@data
meta <- sce@meta.data
sce$type = Idents(sce)
cellchat <- createCellChat(object = sce, meta = meta, group.by = "Celltype")

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "cell-cell_communications.all.csv")

df.net1 <- subsetCommunication(cellchat,slot.name = "netP")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

cellchat@netP$pathways
pathways.show <- c("MHC-I")
vertex.receiver = c(1,2,4,5) 
netVisual_aggregate(cellchat, signaling = pathways.show,  
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

saveRDS(cellchat,'cellchatResult.RDS')

sceTB <- subset(x = sce,subset = Group == "TB")
data.input <- sceTB@assays$RNA@data
meta <- sceTB@meta.data
cellchatTB <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
cellchat <- cellchatTB
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

write.csv(df.net, "TB_cell-cell_communications.all.csv")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cellchatTB <- cellchat
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

saveRDS(cellchatTB,'TB_cellchat.RDS')

sceLTB <- subset(x = sce,subset = Group == "LTB")
data.input <- sceLTB@assays$RNA@data
meta <- sceLTB@meta.data
cellchatLTB <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
cellchat <- cellchatLTB
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

write.csv(df.net, "LTB_cell-cell_communications.all.csv")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cellchatLTB <- cellchat
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

saveRDS(cellchatLTB,'LTB_cellchat.RDS')

sceHealthy <- subset(x = sce,subset = Group == "Healthy")
data.input <- sceHealthy@assays$RNA@data
meta <- sceHealthy@meta.data
cellchatHealthy <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")
cellchat <- cellchatHealthy
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

write.csv(df.net, "Healthy_cell-cell_communications.all.csv")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cellchatHealthy <- cellchat
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

saveRDS(cellchatHealthy,'Healthy_cellchat.RDS')

cellchat.list <- list(ATB = cellchatTB,LTB = cellchatLTB)
cellchat <- mergeCellChat(cellchat.list,add.names = names(cellchat.list),cell.prefix = TRUE)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
p <- gg1 + gg2
ggsave("Overview number strength.pdf",p,width = 6, height = 4)

par(mfrow = c(1,3))
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1 + h2

library(sscVis)
library(data.table)
library(grid)
library(cowplot)
library(ggrepel)
library(readr)
library(plyr)
library(ggpubr)
library(tidyverse)
library(viridis)
library(Seurat)
library(pheatmap)

if(T){
  do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                            meta.cluster = cellInfo.tb$meta.cluster,
                            colname.patient = "patient",
                            loc = cellInfo.tb$loc,
                            out.prefix,
                            pdf.width=3,
                            pdf.height=5,
                            verbose=0){
    ##input data 
    library(data.table)
    dir.create(dirname(out.prefix),F,T)
    
    cellInfo.tb = data.table(cellInfo.tb)
    cellInfo.tb$meta.cluster = as.character(meta.cluster)
    
    if(is.factor(loc)){
      cellInfo.tb$loc = loc
    }else{cellInfo.tb$loc = as.factor(loc)}
    
    loc.avai.vec <- levels(cellInfo.tb[["loc"]])
    count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
    freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
    freq.dist.bin <- floor(freq.dist * 100 / 10)
    print(freq.dist.bin)
    
    {
      count.dist.melt.ext.tb <- test.dist.table(count.dist)
      p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
      OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
      OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
      rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
    }
    
    sscVis::plotMatrix.simple(OR.dist.mtx,
                              out.prefix=sprintf("%s.OR.dist",out.prefix),
                              show.number=F,
                              waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                              exp.name=expression(italic(OR)),
                              z.hi=4,
                              palatte=viridis::viridis(7),
                              pdf.width = 4, pdf.height = pdf.height)
    if(verbose==1){
      return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                  "p.dist.tb"=p.dist.tb,
                  "OR.dist.tb"=OR.dist.tb,
                  "OR.dist.mtx"=OR.dist.mtx))
    }else{
      return(OR.dist.mtx)
    }
  }
  
  test.dist.table <- function(count.dist,min.rowSum=0)
  {
    count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid","cid","count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
      this.row <- count.dist.melt.tb$rid[i]
      this.col <- count.dist.melt.tb$cid[i]
      this.c <- count.dist.melt.tb$count[i]
      other.col.c <- sum.col[this.col]-this.c
      this.m <- matrix(c(this.c,
                         sum.row[this.row]-this.c,
                         other.col.c,
                         sum(sum.col)-sum.row[this.row]-other.col.c),
                       ncol=2)
      res.test <- fisher.test(this.m)
      data.frame(rid=this.row,
                 cid=this.col,
                 p.value=res.test$p.value,
                 OR=res.test$estimate)
    }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                    by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
    return(count.dist.melt.ext.tb)
  }
}

meta <- read.table('adata.obs.csv',sep = ',',header = T,check.names = F)
meta$loc <- meta$Group
meta$meta.cluster <- meta$Celltype

out.prefix <- "./"
OR.immune.list <- do.tissueDist(cellInfo.tb=meta,
                                out.prefix=sprintf("%s.Immune_cell",out.prefix),
                                pdf.width=4,pdf.height=8,verbose=1
)
OR.immune.list

a=OR.immune.list[["OR.dist.tb"]]
a <- as.data.frame(a)
rownames(a) <- a$rid
a <- a[,-1]
a <- na.omit(a)
a

b <- OR.immune.list$count.dist.melt.ext.tb[,c(1,2,6)]
b <- spread(b,key = "cid", value = "adj.p.value")
b <- data.frame(b[,-1],row.names = b$rid)
b <- b[rownames(a),]
b

col <- viridis(11,option = "D")
b = ifelse(b >= 0.05&(a>1.5|a<0.5), "",
           ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                  ifelse(b<0.001&(a>1.5|a<0.5),"***",
                         ifelse(b<0.01&(a>1.5|a<0.5),"**",
                                ifelse(b < 0.05&(a>1.5|a<0.5),"*","")))))

bk=c(seq(0,0.99,by=0.01),seq(1,2,by=0.01))

pheatmap(a[,], border_color = "NA", fontsize = 9,cellheight = 12,cellwidth = 20,clustering_distance_rows="correlation",
         display_numbers = b,number_color="black",fontsize_number=10,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, treeheight_row = 20,treeheight_col = 20,
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))

df <- read.table('CelltypeCount.txt',sep = '\t',header = T,check.names = F)
library(ggplot2)
library(ggalluvial)
pp <- ggplot(df,aes(x = 3,y = Per,fill= Celltype)) + geom_col(width=1.5,color ='white') + facet_grid(.~Group)
pp
colorlist <- c("#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18","#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f","#279772","#add387","#d9b71a")
p1 <- pp + coord_polar(theta = "y") + xlim(c(0.2,3.8)) + 
  scale_fill_manual(values = colorlist) + theme_void() + 
  theme(
    strip.text.x = element_text(size = 14),
    legend.title= element_text(size=15),
    legend.text = element_text(size = 14))

p1

df$Celltype <- factor(df$Celltype,levels = unique(df$Celltype))
p <- ggplot(df,aes(x= Group,y=Per,fill=Celltype,stratum=Celltype,alluvium=Celltype)) + 
  scale_fill_manual(values=colorlist) + 
  scale_y_continuous(expand=c(0,0)) + 
  theme_classic()
p5 <- p + geom_col(width = 0.6, color = NA,size = 0.5) + 
  geom_flow(width = 0.6,alpha = 0.22,knot.pos = 0,color = 'white',size = 0.5) + 
  geom_alluvium(width = 0.6,alpha = 1,knot.pos = 0,fill = NA,color = 'white',size = 0.5)
p5

library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)
#devtools::install_github("junjunlab/ClusterGVis")
sce <- readRDS('Cellchat/normalize.sce.cellchat.RDS')
sce@active.ident <- as.factor(sce@meta.data$Celltype)
Idents(sce) <- sce@meta.data[["Celltype"]]
pbmc.markers.all <- Seurat::FindAllMarkers(sce, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
pbmc.markers.all <- read.table('pbmc.markers.all.txt',sep = '\t',header = T,row.names = 1,check.names = F)
pbmc.markers <- pbmc.markers.all %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
st.data <- prepareDataFromscRNA(object = sce,diffData = pbmc.markers,showAverage = TRUE)
enrich <- enrichCluster(object = st.data,OrgDb = org.Hs.eg.db,type = "BP",organism = "hsa",pvalueCutoff = 0.5,topn = 5, seed = 5201314)
markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,replace = F)]
visCluster(object = st.data, plot.type = "line")

pdf('sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 70,
           markGenes = markGenes,
           cluster.order = c(1:13))
dev.off()

pdf('sc2.pdf',height = 16,width = 18,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 70,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:13),
           go.col = rep(jjAnno::useMyCol("stallion",n = 13),each = 5),
           add.bar = T)
dev.off()
saveRDS(sce,"heatmap.rds")
write.table(pbmc.markers.all,'pbmc.markers.all.txt',sep = '\t')