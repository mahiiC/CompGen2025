


########################################################################

# creating a new token on github for access to the required packages

if (!require("gitcreds", quietly = TRUE))
  install.packages("gitcreds")  

library(gitcreds)
gitcreds::gitcreds_set()

###########################################################################

# install necessary packages

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("BIMSBbioinfo/VoltRon")

if (!require("Biocmanager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")
BiocManager::install("RBioFormats")
install.packages("Seurat")
devtools::install_github("dmcable/spacexr")
BiocManager::install("ComplexHeatmap")
BiocManager::install("HDF5Array")
devtools::install_github("BIMSBbioinfo/ImageArray")
devtools::install_github("bnprks/BPCells/r@v0.3.0")
devtools::install_github("vitessce/vitessceR")
BiocManager::install("basilisk")
install.packages("ggnewscale")
devtools::install_github('immunogenomics/presto')

# Downloading data .zip file from terminal
#wget https://bimsbstatic.mdc-berlin.de/landthaler/VoltRon/workshop.zip
#unzip workshop.zip

# Calling necessary libraries
library(VoltRon)
library(rhdf5)

# import Breast cancer Visium data
Visium_breast_cancer <- importVisium("./workshop/data/BreastCancer/Visium/",
                                     sample_name = "Breast cancer")

# Omic profile clustering
# Processing and filtering

# Features
head(vrFeatures(Visium_breast_cancer))
length(vrFeatures(Visium_breast_cancer))

# Normalize and select features
Visium_breast_cancer <- normalizeData(Visium_breast_cancer)
Visium_breast_cancer <- getFeatures(Visium_breast_cancer, n = 2000)

# selected features
head(vrFeatureData(Visium_breast_cancer))
selected_features <- getVariableFeatures(Visium_breast_cancer)
head(selected_features, 20)

# Dimensional reduction
# Embedding
Visium_breast_cancer <- getPCA(Visium_breast_cancer, features = selected_features, dims = 30)
Visium_breast_cancer <- getUMAP(Visium_breast_cancer, dims = 1:30)
vrEmbeddingNames(Visium_breast_cancer)

# embedding visualization
vrEmbeddingPlot(Visium_breast_cancer, embedding = "umap")
vrSpatialPlot(Visium_breast_cancer)

## Clustering
# graph for neighbours
Visium_breast_cancer <- getProfileNeighbors(Visium_breast_cancer, dims = 1:30, k = 10, method = "SNN")
vrGraphNames(Visium_breast_cancer)

# clustering
Visium_breast_cancer <- getClusters(Visium_breast_cancer, resolution = 0.5, label = "Clusters", graph = "SNN")

# Visualization
# embedding
vrEmbeddingPlot(Visium_breast_cancer, embedding = "umap", group.by = "Clusters")
vrSpatialPlot(Visium_breast_cancer, group.by = "Clusters")
