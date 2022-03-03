library(Seurat)
library(filesstrings)
##Rscript scRNA_Violin_plot.R raw_expression.csv cluster.csv gene_list.csv name
##TODO:read in 10X data directly ("filtered_feature_bc_matrix")


args=commandArgs(TRUE)
raw_expr_matrix<-args[1]
cluster<-args[2]
gene_list<-args[3]
name<-args[4]

folder_name=paste(name,"_Violin_plots",sep="")
dir.create(folder_name)

data<-read.csv(raw_expr_matrix,row.name=1,check.names=F,header=T)
cluster<-read.csv(cluster,row.names=1,check.names=F,header=T)
gene_list<-read.csv(gene_list,header=T)
gene_list<-as.vector(gene_list[,1])
print(head(cluster))
colnames(data)<-sub("\\.","-",colnames(data))
data<-data[,rownames(cluster)]
dim(data)
print(data[1:3,1:3])
object1 <- CreateSeuratObject(counts = data, project = "test", min.cells = 0, min.features = 0)
object1@meta.data<-cbind(object1@meta.data[intersect(rownames(cluster),rownames(object1@meta.data)),],cluster[intersect(rownames(cluster),rownames(object1@meta.data)),])
print(head(object1@meta.data))
Idents(object1)<-"Cluster_name"
print(head(object1@meta.data))
object1 <- NormalizeData(object1, normalization.method = "LogNormalize", scale.factor = 10000)
object1 <- FindVariableFeatures(object1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(object1)
object1 <- ScaleData(object1, features = all.genes)


Violin_plot<-function(gene){
  for (i in 1:length(gene)){
    print(gene[i])
    name = paste(gene[i],"Violin_plot.pdf",sep = "_")
    pdf(name,10,10)
    plot2<-VlnPlot(object1, features = gene[i],log = TRUE,group.by = "Cluster_name")
    print(plot2)
    dev.off()
    file.move(name,folder_name,overwrite = TRUE)
  }
}

Violin_plot(gene_list)

all_violin_name=paste(name,"_Violin_all.pdf",sep="")

pdf(all_violin_name,30,30)
VlnPlot(object1, features = gene_list,log = TRUE,group.by = "Cluster_name")
dev.off()

file.move(all_violin_name,folder_name,overwrite = TRUE)
