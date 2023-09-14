library(Seurat)
tpmdata<-read.delim("TPM_matrix_rename.txt",sep="\t", header=T)
rownames(tpmdata)<-tpmdata$gene
tpmdata$gene<-NULL

logtpm<-log2(tpmdata+1)
seu<-CreateSeuratObject(counts=logtpm)
seu$sex<-seu$orig.ident
seu$orig.ident<-colnames(seu)
seu$nCount_RNA<-NULL

# add life-stage, extract from orig.ident
seu$stage<-gsub("(.*_){1}(.*)_.+", "\\2", seu$orig.ident)
colnames(seu@meta.data)<-c("sample_name", "genes", "sex", "stage")

# rename stages
seu<-SetIdent(seu, value="stage")
new_stage<-c("1_Eggs", "2_Miracidia", "3a_Sporocysts.1d", "3b_Sporocysts.5d", "3c_Sporocysts.32d", "4_Cercariae", "5_Somules.2d", "6_Juveniles")
names(new_stage)<-levels(seu)
seu<-RenameIdents(seu, new_stage)
seu$stage<-Idents(seu)

# rename sex
seu<-SetIdent(seu, value="sex")
new_sex<-c("Female", "Male", "Mixed")
names(new_sex)<-levels(seu)
seu<-RenameIdents(seu, new_sex)
seu$sex<-Idents(seu)

seu<-SetIdent(seu, value="stage")

# add pca embedding
seu<-FindVariableFeatures(seu)
seu<-ScaleData(seu)
seu<-RunPCA(seu)
DimPlot(seu, reduction = "pca", group.by="stage")+coord_fixed()
seu$pca<-CreateDimReducObject(embeddings=seu$pca@cell.embeddings[, 1:2], key='PC_', assay='RNA')

saveRDS(seu, file="Sm_lifeCycle_update.rds")