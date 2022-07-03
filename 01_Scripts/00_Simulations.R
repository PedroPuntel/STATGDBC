# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: Script contendo as simulações/análises contempladas no capítulo 5 da monografia
# Encoding: UTF-8

# Imports
source("01_Scripts/01_STATGDBC.R")
source("01_Scripts/02_PlotMDS.R")
source("01_Scripts/04_PlotGrid.R")
options(digits = 3)

# Pacotes
require("fpc")
require("cluster")
require("data.table")
require("pdfCluster")

####################
# Rotinas auxiliares
####################
# --> Rotinas externas que calculam os índices de qualidade para comparação dos resultados

# Rotina auxiliar que calcula o Índice de Silhueta Médio
get_ism <- function(dataset, statgdbc.obj) {
    
    cluster_vec <- statgdbc.obj$cluster
    
    if( !is.integer(which(cluster_vec == -1, arr.ind = T)) ) {
        to_remove <- which(cluster_vec == -1, arr.ind = T) 
        dataset <- dataset[-to_remove,]
        cluster_vec <- cluster_vec[cluster_vec != -1]
    }
    
    K <- uniqueN(cluster_vec)
    std_dist_mat <- apply(dataset, 2, function(j) {(j-mean(j))/sd(j)}) %>% dist() %>% as.matrix()
    
    if( K == 1)  {
        return(0)
    } else {
        score <- cluster::silhouette(cluster_vec, dmatrix=std_dist_mat) %>% as.data.frame()
        return(mean(score$sil_width))
    }
    
}

# Rotina auxiliar que calcula o Índice de Calinski-Harabasz
get_ich <- function(statgdbc.obj) {
    
    mds_coords <- cbind(statgdbc.obj$spatial.ppp.obj$x, statgdbc.obj$spatial.ppp.obj$y)
    cluster_vec <- statgdbc.obj$cluster
    
    if( !is.integer(which(cluster_vec == -1, arr.ind = T)) ) {
        to_remove <- which(cluster_vec == -1, arr.ind = T) 
        mds_coords <- mds_coords[-to_remove,]
        cluster_vec <- cluster_vec[cluster_vec != -1]
    } 
    
    K <- uniqueN(cluster_vec)
    
    if( K == 1)  {
        return(-Inf)
    } else {
        return(fpc::calinhara(x=mds_coords, clustering=cluster_vec, cn=K))
    }

}

# Rotina auxiliar que calcula o Índice de Rand Ajustado
get_adj_rand_index <- function(statgdbc.obj, k.real) {
    return(pdfCluster::adj.rand.index(statgdbc.obj$cluster, k.real))
}

################################
# Leitura dos conjuntos de dados
################################

# Leitura das projeções EMD (dados de clusterização)
all_cluster_data <- list.files("00_Data/Processed/Clustering/", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.data.frame)

# Leitura das projeções EMD (dados de classificação)
all_classf_data <- list.files("00_Data/Processed/Classification/", pattern="CLASSF*", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.data.frame)
all_classf_clusters <- list.files("00_Data/Processed/Classification/Clusters/", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.data.frame)

# Subconjunto de teste
test_subset_data <- list.files("00_Data/Processed/Subset/", full.names = T) %>% sort() %>% lapply(fread) %>% lapply(as.data.frame)

######################################
# 5.1 - Calibração de Parâmetros BRKGA
######################################
# --> 5 replicações x 7 bases de teste x 3 tamanhos de pop. x 3 tamanhos de ger. x 2 composições de grade: 630 execuções
# --> Objetivo: Por limitação de tempo, calibração dos parâmetros mais impactantes do BRKGA (demais seguem sugestão da literatura)

# Produto cartesiano dos parâmetros
# params <- expand.grid(
#     c(1:7),          # Conjuntos de dados
#     c(100,200,300),  # Tamanho de População
#     c(200,300,500),  # Número de gerações
#     c("esg","asg")   # Composição de grade
# ) %>% as.data.frame()
# colnames(params) <- c("datasets","p","iter","grid.type")

# Simulações (Approx 6hrs)
# cal_results <- lapply(1:nrow(params), function(i) {
#     replicate(5, STATGDBC(data = test_subset_data[params$datasets[i]], grid.type = params$grid.type[i], p = params$p[i], iter = params$iter[i]))
# }); saveRDS(cal_results, "00_Data/Results/Parameter_Calibration/Calibration_Results_v1.rds")

# Leitura dos resultados
# cal_results <- readRDS("00_Data/Results/Parameter_Calibration/Calibration_Results_v1.rds")

# Inclusão do nome dos datasets
# ds_names <- data.frame(datasets = c(1:7), names = c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2"))
# params <- params %>% left_join(ds_names, on='datasets') %>% select(c('names','p','iter','grid.type'))
# params$grid.type <- ifelse(params$grid.type == 'esg', 'ESGBRKGA', 'ASGBRKGA')
# params <- rename(params, grid.algo = grid.type)

# Resultados
# params$mean_ISM <- sapply(cal_results, function(i) mean(unlist(i['score',]))) %>% unlist()
# params$mean_K <- sapply(cal_results, function(i) mean(unlist(i['k',]))) %>% unlist()
# params$mean_T <- sapply(cal_results, function(i) mean(unlist(i['exec',]))) %>% unlist()

# Análise - Melhor combinação de parâmetros por algoritmo de grade
# cal_results_ESGBRKGA <- params %>%
#     group_by(grid.algo, p, iter) %>% 
#     summarise(Mean_ISM = mean(mean_ISM), Mean_T = mean(mean_T)) %>% 
#     arrange(desc(Mean_ISM)) %>% 
#     filter(grid.algo == 'ESGBRKGA')
# 
# cal_results_ASGBRKGA <- params %>%
#     group_by(grid.algo, p, iter) %>% 
#     summarise(Mean_ISM = mean(mean_ISM), Mean_T = mean(mean_T)) %>% 
#     arrange(desc(Mean_ISM)) %>% 
#     filter(grid.algo == 'ASGBRKGA')
# 
# cal_results_ESGBRKGA; cal_results_ASGBRKGA

# Exporta os resultados
# cal_results_ESGBRKGA %>% as.data.table() %>% fwrite("00_Data/Results/Parameter_Calibration/ESGBRKGA_Results_v1.csv", sep=';')
# cal_results_ASGBRKGA %>% as.data.table() %>% fwrite("00_Data/Results/Parameter_Calibration/ASGBRKGA_Results_v1.csv", sep=';')

# Comentários
# . Em geral, ainda persistem valores bem baixos de silhueta (independente da abordagem de grade)
# . Contudo, novamente valores superiores de silhueta para a abordagem assimétrica, porém não tão significantes levando-se em consideração o tempo de execução
# . ASGBRKGA mais sensível a escolha dos parâmetros 'iter' e 'p' (apesar de variabilidade pequena). Melhores valores: (p = 200, iter = 300)
# . ESGBRKGA praticamente invariante a escolha de parâmetros, inclusive em temros de tempo de execução. Assim, escolher os parâmetros mais 'custosos'
#   é justificável, simplesmente por uma questão de avaliar mais possibilidades (p = 500, iter = 300)

# Conclusão:
# . Parâmetros p/ ESGBRKGA: p = 500 | iter = 300
# . Parâmetros p/ ASGBRKGA: p = 200 | iter = 300

############################################
# 5.2 - Análise de Estabilidade do algoritmo
############################################
# --> 10 replicações x 7 bases x 2 composições de grade
# --> Objetivo: Avaliar estabilidade em termos das diferentes composições de grade

# Simulações
# sim_ESG <- lapply(test_subset_data, function(i) {
#     replicate(10, STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette"))
# }); saveRDS(sim_ESG, '00_Data/Results/Stability_Assessment/ESG/sim_ESG_v2.rds') # Approx 15 mints
# sim_ASG <- lapply(test_subset_data, function(i) {
#     replicate(10, STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette"))
# }); saveRDS(sim_ASG, '00_Data/Results/Stability_Assessment/ASG/sim_ASG_v2.rds') # Approx 1.5hrs

# Leitura dos resultados
# sim_ESG <- readRDS('00_Data/Results/Stability_Assessment/ESG/sim_ESG_v2.rds')
# sim_ASG <- readRDS('00_Data/Results/Stability_Assessment/ASG/sim_ASG_v2.rds')

# Inspeção visual dos agrupamentos formados
# test_datasets <- c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2")
# lapply(1:length(sim_ESG), function(i) {
#     PlotMDS(cbind(sim_ESG[[i]][,1]$spatial.ppp.obj$x,sim_ESG[[i]][,1]$spatial.ppp.obj$y), sim_ESG[[i]][,1]$cluster, title=test_datasets[i])
# })
# lapply(1:length(sim_ASG), function(i) {
#     PlotMDS(cbind(sim_ASG[[i]][,1]$spatial.ppp.obj$x,sim_ASG[[i]][,1]$spatial.ppp.obj$y), sim_ASG[[i]][,1]$cluster, title=test_datasets[i])
# })

# Tabela com resultados - Grade Reuglar (ESG)
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
#     "Base" = test_datasets,
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
# ) %>% fwrite('00_Data/Results/Stability_Assessment/ESG/ESG_Results_v2.csv', sep = ';')

# Tabela com resultados - Grade Irregular (ASG)
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
#     "Base" = test_datasets,
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
# ) %>% fwrite('00_Data/Results/Stability_Assessment/ASG/ASG_Results_v2.csv', sep = ';')

########################################
# 5.3.1 - Escalonamento Multidimensional
########################################
# --> Comparação visual das projeções EMD x coordenadas originais de algumas as instâncias bidimensionais



#################################################
# 5.3.2 - Estudo e Comparação entre os Parâmetros
#################################################
# --> Execução do algoritmo para cada uma das bases de dados com parâmetros default
# --> Para as bases de classificação, comparar número de grupos obtidos e Índice Rand Ajustado

# Nome das instâncias de teste
test_datasets <- c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2")

# Nomes das instâncias de clusterização 
clust_datasets <- c(
   "200DATA"
   ,"400P3C"  
   ,"A1"          
   ,"BANKNOTE"    
   ,"BREASTB"     
   ,"CHART"       
   ,"CONCRETEDATA"
   ,"DBLCA"       
   ,"DBLCB"       
   ,"DOWJONES"    
   ,"FACE"        
   ,"GAMMA400"    
   ,"HAYES-ROTH"  
   ,"INDIAN"
   ,"INDOCHINACOMBAT" 
   ,"MARONNA"         
   ,"NORMAL300"       
   ,"NUMBERS2"        
   ,"OUTLIERS"        
   ,"PARKINSONS"      
   ,"PIB.MINAS"       
   ,"PIB100"          
   ,"RUSPINI"         
   ,"SONAR"           
   ,"SPHERICAL4D3C"   
   ,"SPRDESP"         
   ,"SYNTHETICCONTROL"
   ,"TRIPADVISOR"     
   ,"UNIFORM400"      
   ,"WAVEFORM21"              
)

# Nomes das instâncias de classificação
classf_datasets <- c(
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

# Resultados - Bases de Clusterização + Grade Regular (ESG)
# clust_ESG_results <- lapply(c(all_cluster_data, test_subset_data), function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(clust_ESG_results, '00_Data/Results/Clustering/ESG/clust_ESG_v2.rds')
# 
# clust_ESG_results <- readRDS('00_Data/Results/Clustering/ESG/clust_ESG_v2.rds')
# 
# lapply(1:length(clust_ESG_results), function(i) {
#     PlotMDS(
#         cbind(clust_ESG_results[[i]]$spatial.ppp.obj$x, clust_ESG_results[[i]]$spatial.ppp.obj$y),
#         clust_ESG_results[[i]]$cluster,
#         title = paste0("ISM: ", clust_ESG_results[[i]]$score, "\n", "Dataset: ", c(clust_datasets, test_datasets)[i])
#     )
# })
# 
# clust_esg_ics <- unlist(sapply(clust_ESG_results, function(i) if(is.null(i$grid.ics)) return(i$best.grid.ics) else return(i$grid.ics)))
# clust_esg_algo <- unlist(sapply(clust_ESG_results, function(i) if(is.null(i$grid.method)) return(i$grids.method) else return(i$grid.method)))
# clust_esg_k <- unlist(sapply(clust_ESG_results, function(i) i$k))
# clust_esg_ism <- unlist(sapply(clust_ESG_results, function(i) i$score))
# clust_esg_ich <- unlist(sapply(clust_ESG_results, function(i) get_ich(i)))
# 
# data.table(
#     "Base" = c(clust_datasets, test_datasets),
#     "Algo" = clust_esg_algo,
#     "K" = clust_esg_k,
#     "ISM" = clust_esg_ism,
#     "ICH" = clust_esg_ich,
#     "ICS" = clust_esg_ics
# ) %>% fwrite('00_Data/Results/Clustering/ESG/ESG_Results_v2.csv', sep = ';')

# Resultados - Bases de Clusterização + Grade Irregular (ASG)
clust_ASG_results <- lapply(c(all_cluster_data, test_subset_data), function(i)
    STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette")
); saveRDS(clust_ASG_results, '00_Data/Results/Clustering/ASG/clust_ASG_v2.rds')

clust_ASG_results <- readRDS('00_Data/Results/Clustering/ASG/clust_ASG_v2.rds')

clust_asg_ics <- sapply(clust_ASG_results, function(i) i$best.grid.ics)
clust_asg_algo <- sapply(clust_ASG_results, function(i) i$grid.method)
clust_asg_algo[sapply(clust_asg_algo, is.null)] <- 'ASGBRKGA'

data.table(
    "Base" = datasets,
    "Algo" = unlist(clust_asg_algo),
    "K" = sapply(clust_ASG_results, function(i) i$k),
    "IS" = sapply(clust_ASG_results, function(i) i$score),
    "ICS" = unlist(clust_asg_ics)
) %>% fwrite('00_Data/Results/Clustering/ASG/ASG_Results_v2.csv', sep = ';')

# Resultados - Bases de Classificação + Grade Regular (ESG)
# classf_ESG_results <- lapply(all_classf_data, function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(classf_ESG_results, '00_Data/Results/Classification/ESG/classf_ESG_v2.rds')
# 
# classf_ESG_results <- readRDS('00_Data/Results/Classification/ESG/classf_ESG_v2.rds')
# 
# lapply(classf_ESG_results, function(i) {PlotMDS(cbind(i$spatial.ppp.obj$x,i$spatial.ppp.obj$y), i$cluster)})
# 
# classf_esg_ics <- unlist(sapply(classf_ESG_results, function(i) if(is.null(i$grid.ics)) return(i$best.grid.ics) else return(i$grid.ics)))
# classf_esg_algo <- unlist(sapply(classf_ESG_results, function(i) if(is.null(i$grid.method)) return(i$grids.method) else return(i$grid.method)))
# classf_esg_k <- unlist(sapply(classf_ESG_results, function(i) i$k))
# classf_esg_ism <- unlist(sapply(classf_ESG_results, function(i) i$score))
# classf_esg_ich <- unlist(sapply(classf_ESG_results, function(i) get_ich(i)))
# 
# data.table(
#     "Base" = classf_datasets,
#     "Algo" = classf_esg_algo,
#     "K" = classf_esg_k,
#     "ISM" = classf_esg_ism,
#     "ICH" = classf_esg_ich,
#     "Rand" = sapply(1:18, function(i) get_adj_rand_index(classf_ESG_results[[i]], all_classf_clusters[[i]]$V1)),
#     "ICS" = classf_esg_ics
# ) %>% fwrite('00_Data/Results/Classification/ESG/ESG_Results_v2.csv', sep = ';')

# Resultados - Bases de Classificação + Grade Irregular (ASG)
classf_ASG_results <- lapply(all_classf_data, function(i)
    STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette")
); saveRDS(classf_ASG_results, '00_Data/Results/Classification/ASG/classf_ASG_v2.rds')

classf_ASG_results <- readRDS('00_Data/Results/Classification/ASG/classf_ASG_v2.rds')

lapply(classf_ASG_results, function(i) {PlotMDS(cbind(i$spatial.ppp.obj$x,i$spatial.ppp.obj$y), i$cluster)})

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
) %>% fwrite('00_Data/Results/Classification/ASG/ASG_Results_v2.csv', sep = ';')

################################
# 5.3.3 - Critérios Estatísticos
################################

#####################################
# 5.4 - Comparações STATGDBC x DBSCAN
#####################################

# Pacote DBSCAN
require("dbscan")








