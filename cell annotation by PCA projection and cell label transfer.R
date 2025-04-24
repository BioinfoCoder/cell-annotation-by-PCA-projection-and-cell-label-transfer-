library(Seurat)
library(Matrix)
library(irlba)  # For PCA
library(RcppAnnoy)  # For fast nearest neighbor search
library(dplyr)
# Assuming the PBMC datasets (3k and 10k) are already normalized
# and represented as sparse matrices
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
#AvailableData()
InstallData("pbmc3k")

pbmc3k<-UpdateSeuratObject(pbmc3k)
pbmc3k@meta.data %>% head()
# routine processing pipeline
pbmc3k<- pbmc3k %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE) %>%
  RunUMAP(dims = 1:10, verbose = FALSE)
p1<- DimPlot(pbmc3k, reduction = "umap", label = TRUE, group.by = 
               "RNA_snn_res.0.5")

p2<- DimPlot(pbmc3k, reduction = "umap", label = TRUE, group.by = "seurat_annotations", label.size = 3)

p1 + p2
# download it here curl -Lo pbmc_10k_v3.rds https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1 

pbmc10k<- readRDS("~/pbmc_10k_v3.rds")
pbmc10k<-UpdateSeuratObject(pbmc10k)
pbmc10k@meta.data %>% head()
DimPlot(pbmc10k, label = TRUE, repel = TRUE) + NoLegend()
pbmc3k_genes <- rownames(pbmc3k)
pbmc10k_genes <- rownames(pbmc10k)

# Find common genes
common_genes <- intersect(pbmc3k_genes, pbmc10k_genes)


pbmc3k <- subset(pbmc3k, features = common_genes)
pbmc10k <- subset(pbmc10k, features = common_genes)

all.equal(rownames(pbmc3k), rownames(pbmc10k))
install.packages("irlba")
library(irlba)
# use the scaled matrix 
pbmc10k_scaled <- pbmc10k@assays$RNA@scale.data

dim(pbmc10k_scaled)
# Perform PCA using irlba (for large matrices). We transpose it first to gene x sample
pca_10k <- irlba(t(pbmc10k_scaled), nv = 100)  # Keep 100 PCs. The orginal seurat object kept 100 PCs
# get the gene loadings (V matrix). 
gene_loadings_10k <- pca_10k$v  # Gene loadings (features/genes in rows, PCs in columns)
dim(gene_loadings_10k)
rownames(gene_loadings_10k) <- rownames(pbmc10k_scaled)
colnames(gene_loadings_10k) <- paste0("PC", 1:100)
# 2068 most variable genes (after subsetting the common genes with the pbmc3k data)
# ideally, we should re-run FindVariableFeatures, but I am skipping it
VariableFeatures(pbmc10k) %>% length()

# Get PCA embeddings/cell embeddings (U matrix * D matrix) 
cell_embeddings_10k <- pca_10k$u %*% diag(pca_10k$d)  # Cell embeddings (10k cells in rows)
dim(cell_embeddings_10k)
rownames(cell_embeddings_10k) <- colnames(pbmc10k_scaled)
colnames(cell_embeddings_10k) <- colnames(gene_loadings_10k)

cell_embeddings_10k[1:5, 1:10]
pbmc3k_normalized <- pbmc3k@assays$RNA$data

# Center the 3k PBMC dataset based on 10k dataset's gene means
pbmc3k_scaled <- scale(t(pbmc3k_normalized), 
                       center = rowMeans(pbmc10k@assays$RNA$data), 
                       scale = TRUE)
dim(pbmc3k_scaled)
# subset the same genes for the scaled 
pbmc3k_scaled<- pbmc3k_scaled[, rownames(pbmc10k_scaled)]

dim(pbmc3k_scaled)
library(ggplot2)
cell_embeddings_3k <- as.matrix(pbmc3k_scaled) %*% gene_loadings_10k

cell_embeddings_3k[1:5, 1:5]
all.equal(rownames(cell_embeddings_3k), rownames(pbmc3k@meta.data))
cbind(cell_embeddings_3k, pbmc3k@meta.data) %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color = seurat_annotations)) +
  theme_classic(base_size = 14)
# PCA space based on pbmc3k its own 
DimPlot(pbmc3k, reduction = "pca", group.by = "seurat_annotations", 
        label = TRUE) +
  NoLegend()
# Use the Annoy algorithm to find nearest neighbors between 3k and 10k datasets
n_neighbors <- 30  # Number of nearest neighbors to find k =30

# Create Annoy index for 10k PBMC dataset
annoy_index <- new(AnnoyAngular, ncol(cell_embeddings_10k)) ##use cosine distance Angular
for (i in 1:nrow(cell_embeddings_10k)) {
  annoy_index$addItem(i - 1, cell_embeddings_10k[i, ])
}
annoy_index$build(10)  # Build the index with 10 trees

# Find nearest neighbors for each cell in 3k dataset
nn_indices <- t(sapply(1:nrow(cell_embeddings_3k), function(i) {
  annoy_index$getNNsByVector(cell_embeddings_3k[i, ], n_neighbors)
}))

# nn_indices gives you the indices of nearest neighbors in the 10k PBMC dataset
# the rows are cells from 3k dataset, columns are the 30 nearest cells in the 10k dataset
dim(nn_indices)
head(nn_indices)
labels_10k<- as.character(pbmc10k$celltype)
# Transfer labels based on majority vote from nearest neighbors
transfer_labels <- apply(nn_indices, 1, function(neighbors) {
  # Get labels for the nearest neighbors
  neighbor_labels <- labels_10k[neighbors + 1]  # Add 1 for R's 1-based index
  
  # Return the most common label (majority vote)
  most_common_label <- names(sort(table(neighbor_labels), decreasing = TRUE))[1]
  return(most_common_label)
})

# Now, transfer_labels contains the predicted labels for the 3k PBMC dataset
head(transfer_labels)
pbmc3k$predicted<- transfer_labels

DimPlot(pbmc3k, reduction = "umap", group.by = "predicted", label = TRUE, repel=TRUE) +
  NoLegend()
# Step 1: Find transfer anchors
anchors <- FindTransferAnchors(
  reference = pbmc10k,     # The reference dataset
  query = pbmc3k,          # The query dataset
  dims = 1:100,            # The dimensions to use for anchor finding
  reduction = "pcaproject" # this is the default
)

# Step 2: Transfer labels
predictions <- TransferData(
  anchors = anchors,           # The anchors identified in the previous step
  refdata = pbmc10k$celltype, # Assuming 'label' is the metadata containing the true labels in seurat_10k
  dims = 1:30                  # Dimensions to use for transferring
)

# Step 3: Add predictions to the query dataset
pbmc3k <- AddMetaData(pbmc3k, metadata = predictions)

# predicted.id is from Seurat's wrapper function, predicted is from our naive implementation
table(pbmc3k$predicted, pbmc3k$predicted.id)
library(ComplexHeatmap)
table(pbmc3k$predicted, pbmc3k$predicted.id) %>%
  as.matrix() %>%
  scale() %>%
  Heatmap(cluster_rows = FALSE, cluster_columns= FALSE, name= "scaled\ncell number")
####################################multiple nearest neighbour#############################
library(RcppAnnoy)

# Number of nearest neighbors to find
n_neighbors <- 30

# Build an annoy index for the 10k dataset
#annoy_index_10k <- new(AnnoyEuclidean, ncol(cell_embeddings_10k))
annoy_index_10k <- new(AnnoyAngular, ncol(cell_embeddings_10k)) #use cosine distance instead

# Add each cell's PCA embeddings to the index
for (i in 1:nrow(cell_embeddings_10k)) {
  annoy_index_10k$addItem(i - 1, cell_embeddings_10k[i, ])  # 0-based index for Annoy
}

# Build the index for fast nearest neighbor search
annoy_index_10k$build(10)
nn_10k_for_3k <- t(sapply(1:nrow(cell_embeddings_3k), function(i) {
  annoy_index_10k$getNNsByVector(cell_embeddings_3k[i, ], n_neighbors)
}))

# Adjust for R's 1-based indexing
nn_10k_for_3k <- nn_10k_for_3k + 1  # convert to 1-based indexing for R

head(nn_10k_for_3k)

# annoy_index_3k <- new(AnnoyEuclidean, ncol(cell_embeddings_3k))
annoy_index_3k <- new(AnnoyAngular, ncol(cell_embeddings_3k)) 

for (i in 1:nrow(cell_embeddings_3k)) {
  annoy_index_3k$addItem(i - 1, cell_embeddings_3k[i, ])  # 0-based index for Annoy
}

annoy_index_3k$build(10)

nn_3k_for_10k <- t(sapply(1:nrow(cell_embeddings_10k), function(i) {
  annoy_index_3k$getNNsByVector(cell_embeddings_10k[i, ], n_neighbors)
}))

# Adjust for R's 1-based indexing
nn_3k_for_10k <- nn_3k_for_10k + 1  # convert to 1-based indexing for R

labels_10k <- as.character(labels_10k)

# Create empty vectors to store the scores and labels
pbmc3k_transferred_labels <- rep(NA, nrow(cell_embeddings_3k))
pbmc3k_transfer_scores <- rep(0, nrow(cell_embeddings_3k))

# Loop through each cell in the 3k dataset to find the mutual nearest neighbors
for (i in 1:nrow(cell_embeddings_3k)) {
  # Get nearest neighbors of the i-th 3k cell in 10k
  nn_in_10k <- nn_10k_for_3k[i, ]
  
  # Initialize count for mutual nearest neighbors
  mutual_count <- 0
  
  # Check mutual nearest neighbors
  for (nn in nn_in_10k) {
    # Check if i-th 3k cell is a nearest neighbor for the nn-th 10k cell
    if (i %in% nn_3k_for_10k[nn, ]) {  # Correct 1-based indexing
      mutual_count <- mutual_count + 1
      
      # Transfer the label from the 10k cell to the 3k cell
      pbmc3k_transferred_labels[i] <- labels_10k[nn]
    }
  }
  
  # Calculate the transfer score (mutual neighbor count / total neighbors)
  pbmc3k_transfer_scores[i] <- mutual_count / n_neighbors
}

# Fill in missing labels for cells without MNN based on nearest neighbor in 10k
for (i in 1:length(pbmc3k_transferred_labels)) {
  if (is.na(pbmc3k_transferred_labels[i])) {
    # Assign the label of the nearest 10k cell
    nearest_10k_cell <- nn_10k_for_3k[i, 1]  # First nearest neighbor
    pbmc3k_transferred_labels[i] <- labels_10k[nearest_10k_cell]
    
    # Assign a lower score for non-mutual neighbors
    pbmc3k_transfer_scores[i] <- 0.01  # assign a small score like 0.01 for non-mutual
  }
}

head(pbmc3k_transferred_labels)

head(pbmc3k_transfer_scores)

# Add predictions to the query dataset
pbmc3k$pbmc3k_transferred_labels<- pbmc3k_transferred_labels

# predicted.id is from Seurat's wrapper function, pbmc3k_transferred_labels is from our naive MNN implementation

table(pbmc3k$pbmc3k_transferred_labels, pbmc3k$predicted.id)
table(pbmc3k$pbmc3k_transferred_labels, pbmc3k$predicted.id) %>%
        as.matrix() %>%
        scale() %>%
        Heatmap(cluster_rows = FALSE, cluster_columns= FALSE, name= "scaled\ncell number")
