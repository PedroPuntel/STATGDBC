# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: Pre-processamento dos conjuntos de dados
# Encoding: UTF-8

# Fontes
# . http://cs.uef.fi/sipu/datasets/
# . https://archive.ics.uci.edu/ml/datasets.php

# Imports
source("01_Scripts/02_MDSProjection.R")

# Pacotes
require("dplyr")
require("data.table")

#####################################
# Conjuntos de dados de Clusterização
#####################################

# Listagem dos arquivos de dados
all_cluster_files <- list.files("00_Data/Raw/Clustering/", full.names = T)

# Leitura dos arquivos de dados e conversão para matriz
all_cluster_data <- lapply(all_cluster_files, fread)
all_cluster_files <- lapply(all_cluster_data, as.matrix)

# Calcula o MDS
all_cluster_mds_projections <- lapply(all_cluster_data, function(i) MDSProjection(i))

# Nomes dos arquivos de saída
all_cluster_mds_destinations <- lapply(list.files("00_Data/Raw/Clustering/"), function(i) paste0("00_Data/Processed/Clustering/MDS_",i))

# Salva os resultados
lapply(1:length(all_cluster_files), function(i) fwrite(all_cluster_mds_projections[[i]], all_cluster_mds_destinations[[i]], sep=";"))

#####################################
# Conjuntos de dados de Classificação
#####################################

# Listagem dos arquivos de dados
all_classf_files <- list.files("00_Data/Raw/Classification/", full.names = T)

# Leitura dos arquivos de dados
all_classf_data <- lapply(all_classf_files, fread)

# Tratamento individual para cada base (separa o vetor de clusters do conjunto de dados original)
# . Manter a variável que identifica os clusters sempre como V3

all_classf_data[[3]]$V1 <- NULL
setnames(all_classf_data[[3]],c('V3','V9'),c('V9','V3'))
all_classf_data[[3]]

all_classf_data[[5]]$V1 <- NULL
all_classf_data[[5]]$V10 <- NULL
uniqueN(all_classf_data[[5]]$V11)
setnames(all_classf_data[[5]],c('V3','V11'),c('V11','V3'))
#all_classf_data[[5]]

setnames(all_classf_data[[6]],c('V4','V3'),c('V3','V4'))
#all_classf_data[[6]]

all_classf_data[[7]]$V2 <- NULL
all_classf_data[[7]]$V1 <- NULL
setnames(all_classf_data[[7]],c('V3','V35'),c('V35','V3'))
#all_classf_data[[7]]

setnames(all_classf_data[[8]],c('V3','V5'),c('V5','V3'))
#all_classf_data[[8]]

setnames(all_classf_data[[10]],c('V3','V1'),c('V1','V3'))
#all_classf_data[[10]]

setnames(all_classf_data[[14]],c('V8','V3'),c('V3','V8'))
#all_classf_data[[14]]

all_classf_data[[16]]$V1 <- NULL
setnames(all_classf_data[[16]],c('V2','V3'),c('V3','V2'))
#all_classf_data[[16]]

setnames(all_classf_data[[17]],c('V1','V3'),c('V3','V1'))
#all_classf_data[[17]]

all_classf_data[[18]]$V1 <- NULL
setnames(all_classf_data[[18]],c('V10','V3'),c('V3','V10'))
#all_classf_data[[18]]

# lapply(all_classf_data, function(i) uniqueN(i$V3)) # Inspeção do número de clusters

# Lista auxiliar contendo o vetor de clusters
all_classf_data_clusters <- lapply(all_classf_data, function(i) as.data.table(i$V3))

# Excluí o vetor de clusters da base original
all_classf_data <- lapply(all_classf_data, function(i) as.data.frame(i)) %>% lapply(function(i) i[,c(colnames(i) != 'V3')])

# Calcula o MDS
all_classf_mds_projections <- lapply(all_classf_data, function(i) MDSProjection(i))

# Nomes dos arquivos de saída
all_classf_mds_destinations <- lapply(list.files("00_Data/Raw/Classification/MDS/"), function(i) paste0("00_Data/Processed/Classification/MDS_",i))
all_classf_cvec_destinations <- lapply(list.files("00_Data/Raw/Classification/Clusters/"), function(i) paste0("00_Data/Processed/Classification/ClusterVec_",i))

# Salva os resultados
lapply(1:length(all_classf_files), function(i) fwrite(all_classf_mds_projections[[i]], all_classf_mds_destinations[[i]], sep=";"))
lapply(1:length(all_classf_files), function(i) fwrite(all_classf_data_clusters[[i]], all_classf_cvec_destinations[[i]], sep=";"))


