# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: Script contendo as simulações/análises contempladas no capítulo 5 da monografia
# Encoding: UTF-8

# Notas de Desenvolvimento
# --> TODO: Investigar motivo do (ICS = NULL) para as instâncias: 200DATA, A1, BROKEN-RING, INDIAN, SONAR (ESG) | INDIAN (ASG)
# --> TODO: Investigar 'grid.method' = NULL

# Imports
source("01_Scripts/01_STATGDBC.R")
source("01_Scripts/02_PlotMDS.R")
source("01_Scripts/04_PlotGrid.R")
options(digits = 3)

# Pacotes
require("pdfCluster")

####################
# Rotinas auxiliares
####################

# TODO: Índice de Silhueta Médio
# TODO: Índice de Calinski-Harabasz

# Rotina auxiliar que calcula o Índice de Rand Ajustado
get_adj_rand_index <- function(statgdbc.obj, k.real) {
    return(pdfCluster::adj.rand.index(statgdbc.obj$cluster, k.real))
}

################################
# Leitura dos conjuntos de dados
################################

# Leitura das projeções EMD (dados de clusterização)
all_cluster_mds_data <- list.files("00_Data/Processed/Clustering/All/", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.matrix)

# Leitura das projeções EMD (dados de classificação)
all_classf_mds_data <- list.files("00_Data/Processed/Classification/MDS/", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.matrix)
all_classf_clusters <- list.files("00_Data/Processed/Classification/Clusters/", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.matrix)

# Subconjunto de bases de clusterização
subset_cluster_mds_data <- list.files("00_Data/Processed/Clustering/Subset/", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.matrix)
# lapply(subset_cluster_mds_data, function(i) PlotMDS(i))

##############################################
# 5.2.1 - Análise de Estabilidade do algoritmo
##############################################
# --> 10 replicações x 7 bases x 2 composições de grade (demais parâmetros default)
# --> Objetivo: Avaliar estabilidade em termos das diferentes composições de grade

# sim_ESG <- lapply(subset_cluster_mds_data, function(i) {
#     replicate(10, STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette"))
# }); saveRDS(sim_ESG, '00_Data/Results/Stability_Assessment/ESG/sim_ESG_v1.rds') # Approx 15 mints

# sim_ASG <- lapply(subset_cluster_mds_data, function(i) {
#     replicate(10, STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette"))
# }); saveRDS(sim_ASG, '00_Data/Results/Stability_Assessment/ASG/sim_ASG_v1.rds') # Approx 2hrs

# Resultados
sim_ESG <- readRDS('00_Data/Results/Stability_Assessment/ESG/sim_ESG_v1.rds')
sim_ASG <- readRDS('00_Data/Results/Stability_Assessment/ASG/sim_ASG_v1.rds')

# Inspeções iniciais
lapply(sim_ESG, function(i) i)
lapply(sim_ASG, function(i) i)

datasets <- c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2")

lapply(1:length(sim_ESG), function(i) {
    PlotMDS(cbind(sim_ESG[[i]][,1]$spatial.ppp.obj$x,sim_ESG[[i]][,1]$spatial.ppp.obj$y), sim_ESG[[i]][,1]$cluster, title=datasets[i])
})
lapply(1:length(sim_ASG), function(i) {
    PlotMDS(cbind(sim_ASG[[i]][,1]$spatial.ppp.obj$x,sim_ASG[[i]][,1]$spatial.ppp.obj$y), sim_ASG[[i]][,1]$cluster, title=datasets[i])
})

# Informações necessárias para construção da tabela
# datasets <- c('2-FACE','BROKEN-RING','BUPA','FORESTFIRES','GAUSS9','UNIFORM700','VOWEL2')

# ESG_min_k <- sapply(sim_ESG, function(i) min(apply(i, 2, function(j) j$k)) )
# ESG_max_k <- sapply(sim_ESG, function(i) max(apply(i, 2, function(j) j$k)) )
# ESG_med_k <- sapply(sim_ESG, function(i) mean(apply(i, 2, function(j) j$k)) ) %>% floor()
# ESG_cv_k <- sapply(sim_ESG, function(i) sd(apply(i, 2, function(j) j$k))/mean(apply(i, 2, function(j) j$k)) )
# ESG_min_is <-  sapply(sim_ESG, function(i) min(apply(i, 2, function(j) j$score)) )
# ESG_max_is <-  sapply(sim_ESG, function(i) max(apply(i, 2, function(j) j$score)) )
# ESG_med_is <-  sapply(sim_ESG, function(i) mean(apply(i, 2, function(j) j$score)) )
# ESG_cv_is <- sapply(sim_ESG, function(i) sd(apply(i, 2, function(j) j$score))/mean(apply(i, 2, function(j) j$score)) )
# ESG_min_exec <-  sapply(sim_ESG, function(i) min(apply(i, 2, function(j) j$exec)) )
# ESG_max_exec <-  sapply(sim_ESG, function(i) max(apply(i, 2, function(j) j$exec)) )
# ESG_med_exec <-  sapply(sim_ESG, function(i) mean(apply(i, 2, function(j) j$exec)) )
# ESG_cv_exec <- sapply(sim_ESG, function(i) sd(apply(i, 2, function(j) j$exec))/mean(apply(i, 2, function(j) j$exec)) )

# data.table(
#     "Base" = datasets,
#     "minK" = ESG_min_k,
#     "medK" = ESG_med_k,
#     "maxK" = ESG_max_k,
#     "CVk" = ESG_cv_k,
#     "minIS" = ESG_min_is,
#     "medIS" = ESG_med_is,
#     "maxIS" = ESG_max_is,
#     "CVIS" = ESG_cv_is,
#     "minEXEC" = ESG_min_exec,
#     "medEXEC" = ESG_med_exec,
#     "maxEXEC" = ESG_max_exec,
#     "CVEXEC" = ESG_cv_exec
# ) %>% fwrite('00_Data/Results/Stability_Assessment/ESG/ESG_Results_v1.csv', sep = ';')

# ASG_min_k <- sapply(sim_ASG, function(i) min(apply(i, 2, function(j) j$k)) )
# ASG_max_k <- sapply(sim_ASG, function(i) max(apply(i, 2, function(j) j$k)) )
# ASG_med_k <- sapply(sim_ASG, function(i) mean(apply(i, 2, function(j) j$k)) ) %>% floor()
# ASG_cv_k <- sapply(sim_ASG, function(i) sd(apply(i, 2, function(j) j$k))/mean(apply(i, 2, function(j) j$k)) )
# ASG_min_is <-  sapply(sim_ASG, function(i) min(apply(i, 2, function(j) j$score)) )
# ASG_max_is <-  sapply(sim_ASG, function(i) max(apply(i, 2, function(j) j$score)) )
# ASG_med_is <-  sapply(sim_ASG, function(i) mean(apply(i, 2, function(j) j$score)) )
# ASG_cv_is <- sapply(sim_ASG, function(i) sd(apply(i, 2, function(j) j$score))/mean(apply(i, 2, function(j) j$score)) )
# ASG_min_exec <-  sapply(sim_ASG, function(i) min(apply(i, 2, function(j) j$exec)) )
# ASG_max_exec <-  sapply(sim_ASG, function(i) max(apply(i, 2, function(j) j$exec)) )
# ASG_med_exec <-  sapply(sim_ASG, function(i) mean(apply(i, 2, function(j) j$exec)) )
# ASG_cv_exec <- sapply(sim_ASG, function(i) sd(apply(i, 2, function(j) j$exec))/mean(apply(i, 2, function(j) j$exec)) )

# data.table(
#     "Base" = datasets,
#     "minK" =ASG_min_k,
#     "medK" = ASG_med_k,
#     "maxK" = ASG_max_k,
#     "CVk" = ASG_cv_k,
#     "minIS" = ASG_min_is,
#     "medIS" = ASG_med_is,
#     "maxIS" = ASG_max_is,
#     "CVIS" = ASG_cv_is,
#     "minEXEC" = ASG_min_exec,
#     "medEXEC" = ASG_med_exec,
#     "maxEXEC" = ASG_max_exec,
#     "CVEXEC" = ASG_cv_exec
# ) %>% fwrite('00_Data/Results/Stability_Assessment/ASG/ASG_Results_v1.csv', sep = ';')

#################################################
# 5.2.2 - Estudo e Comparação entre os Parâmetros
#################################################
# TODO: Devido ao caráter estocástico, realizar mais de uma execução (melhor de 3 ou de 5, por exemplo)
# --> Execução do algoritmo para cada uma das bases de dados com parâmetros default
# --> Para as bases de classificação, comparar número de grupos obtidos + Índice Rand Ajustado

# datasets <- c(
#     "2-FACE"    
#    ,"400P3C"      
#    ,"200DATA"     
#    ,"A1"          
#    ,"BANKNOTE"    
#    ,"BREASTB"     
#    ,"BROKEN-RING"
#    ,"BUPA"       
#    ,"CHART"       
#    ,"CONCRETEDATA"
#    ,"DBLCA"       
#    ,"DBLCB"       
#    ,"DOWJONES"    
#    ,"FACE"        
#    ,"FORESTFIRES"
#    ,"GAMMA400"    
#    ,"GAUSS9"     
#    ,"HAYES-ROTH"  
#    ,"INDIAN"
#    ,"INDOCHINACOMBAT" 
#    ,"MARONNA"         
#    ,"MORESHAPES"      
#    ,"NORMAL300"       
#    ,"NUMBERS2"        
#    ,"OUTLIERS"        
#    ,"PARKINSONS"      
#    ,"PIB.MINAS"       
#    ,"PIB100"          
#    ,"RUSPINI"         
#    ,"SONAR"           
#    ,"SPHERICAL4D3C"   
#    ,"SPRDESP"         
#    ,"SYNTHETICCONTROL"
#    ,"TRIPADVISOR"     
#    ,"UNIFORM400"      
#    ,"UNIFORM700"     
#    ,"VOWEL2"         
#    ,"WAVEFORM21"              
# )

# clust_ESG_results <- lapply(c(all_cluster_mds_data, subset_cluster_mds_data), function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(clust_ESG_results, '00_Data/Results/Clustering/ESG/clust_ESG_v1.rds') # Approx. 17 mints

# clust_ESG_results <- readRDS('00_Data/Results/Clustering/ESG/clust_ESG_v1.rds')

# clust_esg_ics <- sapply(clust_ESG_results, function(i) i$grid.ics)
# clust_esg_ics[sapply(clust_esg_ics, is.null)] <- 0
# clust_esg_algo <- sapply(clust_ESG_results, function(i) i$grid.method)
# clust_esg_algo[sapply(clust_esg_algo, is.null)] <- 'ESGBRKGA'

# data.table(
#     "Base" = datasets,
#     "Algo" = unlist(clust_esg_algo),
#     "K" = sapply(clust_ESG_results, function(i) i$k),
#     "IS" = sapply(clust_ESG_results, function(i) i$score),
#     "ICS" = unlist(clust_esg_ics)
# ) %>% fwrite('00_Data/Results/Clustering/ESG/ESG_Results_v1.csv', sep = ';')

# clust_ASG_results <- lapply(c(all_cluster_mds_data, subset_cluster_mds_data), function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(clust_ASG_results, '00_Data/Results/Clustering/ASG/clust_ASG_v1.rds') # Approx. 3hrs

# clust_ASG_results <- readRDS('00_Data/Results/Clustering/ASG/clust_ASG_v1.rds')

# clust_asg_ics <- sapply(clust_ASG_results, function(i) i$best.grid.ics)
# clust_asg_ics[sapply(clust_asg_ics, is.null)] <- 0
# clust_asg_algo <- sapply(clust_ASG_results, function(i) i$grid.method)
# clust_asg_algo[sapply(clust_asg_algo, is.null)] <- 'ASGBRKGA'

# data.table(
#     "Base" = datasets,
#     "Algo" = unlist(clust_asg_algo),
#     "K" = sapply(clust_ASG_results, function(i) i$k),
#     "IS" = sapply(clust_ASG_results, function(i) i$score),
#     "ICS" = unlist(clust_asg_ics)
# ) %>% fwrite('00_Data/Results/Clustering/ASG/ASG_Results_v1.csv', sep = ';')

# classf_ESG_results <- lapply(all_classf_mds_data, function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(classf_ESG_results, '00_Data/Results/Classification/ESG/classf_ESG_v1.rds') # Approx. 10 mints

datasets <- c(
    "AGGREGATION"
    ,"COMPOUND"
    ,"GLASS"
    ,"ECOLI"
    ,"FLAME"
    ,"HABERMAN"
    ,"IRIS"
    ,"IONOSPHERE"
    ,"JAIN"
    ,"NEW-THYROID"
    ,"PATHBASED"
    ,"PIMA.INDIANS"
    ,"R15"
    ,"SEEDS"
    ,"SPIRAL"
    ,"WDBC"
    ,"WINE"
    ,"YEAST"     
)

classf_ESG_results <- readRDS('00_Data/Results/Classification/ESG/classf_ESG_v1.rds')
lapply(classf_ESG_results, function(i) {
    PlotMDS(cbind(i$spatial.ppp.obj$x,i$spatial.ppp.obj$y), i$cluster)
})

classf_esg_ics <- sapply(classf_ESG_results, function(i) i$grid.ics)
classf_esg_ics[sapply(classf_esg_ics, is.null)] <- 0
classf_esg_algo <- sapply(classf_ESG_results, function(i) i$grid.method)
classf_esg_algo[sapply(classf_esg_algo, is.null)] <- 'ESGBRKGA'

data.table(
    "Base" = datasets,
    "Algo" = unlist(classf_esg_algo),
    "K" = sapply(classf_ESG_results, function(i) i$k),
    "IS" = sapply(classf_ESG_results, function(i) i$score),
    "Rand" = sapply(1:length(classf_ESG_results), function(i) get_adj_rand_index(classf_ESG_results[[i]], all_classf_clusters[[i]])),
    "ICS" = unlist(classf_esg_ics)
) %>% fwrite('00_Data/Results/Classification/ESG/ESG_Results_v1.csv', sep = ';')

# classf_ASG_results <- lapply(all_classf_mds_data, function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(classf_ASG_results, '00_Data/Results/Classification/ASG/classf_ASG_v1.rds') # Approx. 40 mints

classf_ASG_results <- readRDS('00_Data/Results/Classification/ASG/classf_ASG_v1.rds')
lapply(classf_ASG_results, function(i) {
    PlotMDS(cbind(i$spatial.ppp.obj$x,i$spatial.ppp.obj$y), i$cluster)
})

classf_asg_ics <- sapply(classf_ASG_results, function(i) i$best.grid.ics)
classf_asg_ics[sapply(classf_asg_ics, is.null)] <- 0
classf_asg_algo <- sapply(classf_ASG_results, function(i) i$grid.method)
classf_asg_algo[sapply(classf_asg_algo, is.null)] <- 'ASGBRKGA'

data.table(
    "Base" = datasets,
    "Algo" = unlist(classf_asg_algo),
    "K" = sapply(classf_ASG_results, function(i) i$k),
    "IS" = sapply(classf_ASG_results, function(i) i$score),
    "Rand" = sapply(1:length(classf_ASG_results), function(i) get_adj_rand_index(classf_ASG_results[[i]], all_classf_clusters[[i]])),
    "ICS" = unlist(classf_asg_ics)
) %>% fwrite('00_Data/Results/Classification/ASG/ASG_Results_v1.csv', sep = ';')


