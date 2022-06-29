# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: Pre-processamento dos conjuntos de dados
# Encoding: UTF-8

# Fontes
# . http://cs.uef.fi/sipu/datasets/
# . https://archive.ics.uci.edu/ml/datasets.php

# Pacotes
require("dplyr")
require("data.table")

######################
# Subconjunto de Teste
######################

# Leitura dos conjuntos de dados 'as is'
all_test_subset_data <- list.files("00_Data/Raw/Subset/", full.names = T) %>% lapply(fread)

# Nomeação dos arquivos de saída
all_test_subset_data_destinations <- lapply(list.files("00_Data/Raw/Subset/"), function(i) paste0("00_Data/Processed/Subset/",i))

# Salva os resultados
lapply(1:length(all_test_subset_data), function(i) fwrite(all_test_subset_data[[i]], all_test_subset_data_destinations[[i]], sep=";"))

#####################################
# Conjuntos de dados de Clusterização
#####################################

# Leitura dos conjuntos de dados 'as is'
all_cluster_data <- list.files("00_Data/Raw/Clustering/", full.names = T) %>% lapply(fread)

# Nomeação dos arquivos de saída
all_cluster_data_destinations <- lapply(list.files("00_Data/Raw/Clustering/"), function(i) paste0("00_Data/Processed/Clustering/",i))

# Salva os resultados
lapply(1:length(all_cluster_data), function(i) fwrite(all_cluster_data[[i]], all_cluster_data_destinations[[i]], sep=";"))

#####################################
# Conjuntos de dados de Classificação
#####################################

# Leitura dos conjuntos de dados 'as is'
all_classf_data <- list.files("00_Data/Raw/Classification/", full.names = T) %>% lapply(fread)

# Tratamento individual para cada base (separação do vetor de clusters)
# --> Define o vetor de clusters sempre sob a variável V3 do conjunto de dados

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

# Verificação - Número de clusters
sapply(all_classf_data, function(i) uniqueN(i$V3))

# Lista auxiliar contendo o vetor de clusters
all_classf_data_clusters <- lapply(all_classf_data, function(i) as.data.table(i$V3))

# Excluí o vetor de clusters das bases originais
all_classf_data <- lapply(all_classf_data, function(i) as.data.frame(i)) %>% lapply(function(i) i[,c(colnames(i) != 'V3')])

# Nomeção dos arquivos de saída
all_classf_destinations <- lapply(list.files("00_Data/Raw/Classification/"), function(i) paste0("00_Data/Processed/Classification/",i))
all_classf_cvec_destinations <- lapply(list.files("00_Data/Raw/Classification/"), function(i) paste0("00_Data/Processed/Classification/Clusters/CVEC_",i))

# Salva os resultados
lapply(1:length(all_classf_data), function(i) fwrite(all_classf_data[[i]], all_classf_destinations[[i]], sep=";"))
lapply(1:length(all_classf_data), function(i) fwrite(all_classf_data_clusters[[i]], all_classf_cvec_destinations[[i]], sep=";"))


