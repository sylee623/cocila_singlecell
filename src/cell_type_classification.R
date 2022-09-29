library("Seurat")
library("dplyr")

setwd("/home/sylee/work/single_cell")
R_SN.data <- Read10X(data.dir = "./raw_data/R_SN")
R_SN <- CreateSeuratObject(counts = R_SN.data, project ="R_SN")
C_SN.data <- Read10X(data.dir = "./raw_data/C_SN")
C_SN <- CreateSeuratObject(counts = C_SN.data, project ="C_SN")

CDh.data <- Read10X(data.dir = "./raw_data/CD_h")
CDh <- CreateSeuratObject(counts = CDh.data, project ="CDh")
CDt.data <- Read10X(data.dir = "./raw_data/CD_t")
CDt <- CreateSeuratObject(counts = CDt.data, project ="CDt")
CDh
CDt
combined_CD <- merge(x=CDh, y=CDt, add.cell.ids = c("CDh", "CDt"))
combined_CD <- PercentageFeatureSet(object = combined_CD, pattern = "^MT-", col.name = "percent.mt")
combined_CD<-subset(x = combined_CD, subset = nCount_RNA< 10000 & nFeature_RNA >300& percent.mt<5)
combined_CD <- NormalizeData(object = combined_CD)
combined_CD<-FindVariableFeatures(object=combined_CD, selection.method="vst", nfeatures=2000)

all.genes<-rownames(x=combined_CD)
combined_CD<-ScaleData(object=combined_CD, features=all.genes)
combined_CD<-RunPCA(object=combined_CD, features=VariableFeatures(object=combined_CD))
combined_CD<-JackStraw(object=combined_CD, num.replicate=100)
combined_CD<-ScoreJackStraw(object=combined_CD, dims=1:20)
JackStrawPlot(combined_CD, dims = 1:15)
ElbowPlot(combined_CD, ndims=20)

CD_umap<-RunUMAP(object=combined_CD, dims=1:10)
CD_umap<-FindNeighbors(object=CD_umap, dims=1:10)
CD_umap<-FindClusters(object=CD_umap, resolution=0.5)
DimPlot(object=CD_umap, reduction="umap")
DimPlot(object=CD_umap, reduction="umap", group.by="orig.ident")

CD_tsne<-RunTSNE(object=combined_CD, dims=1:8)
CD_tsne<-FindNeighbors(object=CD_tsne, dims=1:8)
CD_tsne<-FindClusters(object=CD_tsne, resolution=0.4)
DimPlot(object=CD_tsne, reduction="tsne")
DimPlot(object=CD_tsne, reduction="tsne", group.by="orig.ident")

## Marker gene expression
## astrocytes (AQP4)
## oligo- dendrocytes (MOBP)
## oligodendrocyte precursor cells (OPCs) (VCAN)
## macrophages/microglia (APBB1IP)
## endothelial cells (FLT1)
FeaturePlot(CD_umap, features = c("AQP4", "MOBP", "VCAN", "APBB1IP"))
FeaturePlot(CD_tsne, features = c("AQP4", "MOBP", "VCAN", "APBB1IP"))
VlnPlot(CD_tsne, features = c("AQP4", "MOBP", "VCAN", "APBB1IP"),pt.size = 0)
VlnPlot(CD_umap, features = c("AQP4", "MOBP", "VCAN", "APBB1IP"),pt.size = 0)


## Rename cluster
new.cluster.ids <- c("Olig.1",
                     "Micro.1",
                     "Micro.2",
                     "Micro.3",
                     "Astro.1",
                     "Olig.2",
                     "Micro.4",
                     "Astro.2",
                     "Olig.3",
                     "Micro.5",
                     "Micro.6",
                     "OPC",
                     "Micro.7")
names(new.cluster.ids) <- levels(CD_tsne)
CD_tsne <- RenameIdents(CD_tsne, new.cluster.ids)
DimPlot(CD_tsne, reduction = "tsne", pt.size = 0.4)
