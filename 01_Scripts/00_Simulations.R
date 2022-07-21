# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: Script contendo as simulações/análises contempladas no capítulo 5 da monografia
# Encoding: UTF-8

# Imports
source("01_Scripts/01_STATGDBC.R")
source("01_Scripts/02_PlotMDS.R")
source("01_Scripts/02_MDSProjection.R")
source("01_Scripts/04_PlotGrid.R")
options(digits = 3)

# Pacotes
require("dbscan")
require('ggplot2')
require("cluster")
require("data.table")
require("pdfCluster")

####################
# Rotinas auxiliares
####################
# --> Rotinas externas que calculam os índices de qualidade para comparação dos resultados

# Rotina auxiliar que calcula a tendência de agrupamento de um conjunto de dados através da estatísticas de hopkins
hopkins <- function(D,m,a) {
    
    # D = Matriz com os Dados
    # m = Numero de Amostras
    # a = Tamanho das amostra em termos percentuais
    
    n <- dim(D)[1]
    na <- n*a
    obs <- t(replicate(m,sample(n,na)))
    d <- as.matrix(dist(D))
    vars_faixa <- t(apply(D,2,range))
    
    calcula_U <- function(a,D,vars_faixa) {
        uq <- apply(vars_faixa,1,function(x) runif(a,x[1],x[2]))
        uq <- sum(apply(uq,1,function(x) min(apply(t((t(D)-x)^2),1,function(z) sqrt(sum(z))))))
    }  
    
    calcula_W <- function(xi,d) {
        wq <- sum(apply(as.matrix(xi),1,function(xi) min(d[xi,-xi])))
        return(wq)
    }  
    
    U <- apply(as.matrix(rep(na,m)),1,function(ax) calcula_U(ax,D,vars_faixa))
    W <- apply(obs,1,function(xi) calcula_W(xi,d))
    H <- U/(U+W)
    
    return(mean(H))
    
}  

# Rotina auxiliar que calcula o Índice de Silhueta Médio
get_ism <- function(dataset, clust.obj, is.statgdbc = T) {
    
    # Debugging:
    # . dataset <- all_datasets[[1]]
    # . clust.obj <- dbscan_results[[1]]
    # . is.statgdbc = F
    
    cluster_vec <- clust.obj$cluster
    
    if( isTRUE(is.statgdbc) ) {
        
        if( any(cluster_vec == -1) ) {
            to_remove <- which(cluster_vec == -1, arr.ind = T) 
            dataset <- dataset[-to_remove,]
            cluster_vec <- cluster_vec[cluster_vec != -1]
        }

    } else {
        
        if( any(cluster_vec == 0) ) {
            to_remove <- which(cluster_vec == 0, arr.ind = T) 
            dataset <- dataset[-to_remove,]
            cluster_vec <- cluster_vec[cluster_vec != 0]
        }
    }
    
    K <- uniqueN(cluster_vec)
    std_dist_mat <- apply(dataset, 2, function(j) {(j-mean(j))/sd(j)}) %>% dist() %>% as.matrix()
    score <- cluster::silhouette(cluster_vec, dmatrix=std_dist_mat) %>% as.data.frame()
    return(mean(score$sil_width))
    
}

# Rotina auxiliar que calcula o Índice de Calinski-Harabasz
get_ich <- function(dataset, clust.obj, is.statgdbc = T) {
    
    # Debugging:
    # . dataset <- all_datasets[[1]]
    # . clust.obj <- dbscan_results[[1]]
    # . is.statgdbc = F
    
    calinhara <- function (x, clustering, cn = max(clustering))  {
        x <- as.matrix(x)
        p <- ncol(x)
        n <- nrow(x)
        cln <- rep(0, cn)
        W <- matrix(0, p, p)
        for (i in 1:cn) cln[i] <- sum(clustering == i)
        for (i in 1:cn) {
            clx <- x[clustering == i, ]
            cclx <- cov(as.matrix(clx))
            if (cln[i] < 2) 
                cclx <- 0
            W <- W + ((cln[i] - 1) * cclx)
        }
        S <- (n - 1) * cov(x)
        B <- S - W
        out <- (n - cn) * sum(diag(B), na.rm = T)/((cn - 1) * sum(diag(W), na.rm = T))
        out
    }
    
    cluster_vec <- clust.obj$cluster
    
    if( isTRUE(is.statgdbc) ) {
        
        if( any(cluster_vec == -1) ) {
            to_remove <- which(cluster_vec == -1, arr.ind = T) 
            dataset <- dataset[-to_remove,]
            cluster_vec <- cluster_vec[cluster_vec != -1]
        }
        
    } else {
        
        if( any(cluster_vec == 0) ) {
            to_remove <- which(cluster_vec == 0, arr.ind = T) 
            dataset <- dataset[-to_remove,]
            cluster_vec <- cluster_vec[cluster_vec != 0]
        }
    }
    
    K <- uniqueN(cluster_vec)
    std_dataset <- apply(dataset, 2, function(j) {(j-mean(j))/sd(j)})
    
    if( K == 1 )  {
        return(NA)
    } else {
        return(calinhara(x=std_dataset, clustering=cluster_vec, cn=K))
    }

}

# Rotina auxiliar que calcula o Índice de Rand Ajustado
get_adj_rand_index <- function(clust.obj, k.real) {
    
    # Debugging:
    # . clust.obj <- dbscan_results[26:43][[1]]
    # . k.real <- all_classf_clusters[[1]]
    
    return(pdfCluster::adj.rand.index(clust.obj$cluster, k.real))
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

# ---------------------------------
# Produto cartesiano dos parâmetros
# ---------------------------------
# params <- expand.grid(
#     c(1:7),          # Conjuntos de dados
#     c(100,200,300),  # Tamanho de População
#     c(200,300,500),  # Número de gerações
#     c("esg","asg")   # Composição de grade
# ) %>% as.data.frame()
# colnames(params) <- c("datasets","p","iter","grid.type")

# ------------------------
# Simulações (Approx 6hrs)
# ------------------------
# cal_results <- lapply(1:nrow(params), function(i) {
#     replicate(5, STATGDBC(data = test_subset_data[params$datasets[i]], grid.type = params$grid.type[i], p = params$p[i], iter = params$iter[i]))
# }); saveRDS(cal_results, "00_Data/Results/Parameter_Calibration/Calibration_Results_v1.rds")

# ----------------------
# Leitura dos resultados
# ----------------------
# cal_results <- readRDS("00_Data/Results/Parameter_Calibration/Calibration_Results_v1.rds")

# -----------------------------
# Inclusão do nome dos datasets
# -----------------------------
# ds_names <- data.frame(datasets = c(1:7), names = c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2"))
# params <- params %>% left_join(ds_names, on='datasets') %>% select(c('names','p','iter','grid.type'))
# params$grid.type <- ifelse(params$grid.type == 'esg', 'ESGBRKGA', 'ASGBRKGA')
# params <- rename(params, grid.algo = grid.type)

# ----------
# Resultados
# ----------
# params$mean_ISM <- sapply(cal_results, function(i) mean(unlist(i['score',]))) %>% unlist()
# params$mean_K <- sapply(cal_results, function(i) mean(unlist(i['k',]))) %>% unlist()
# params$mean_T <- sapply(cal_results, function(i) mean(unlist(i['exec',]))) %>% unlist()

# ----------------------------------------------------------------
# Análise - Melhor combinação de parâmetros por algoritmo de grade
# ----------------------------------------------------------------
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

# ---------------------
# Exporta os resultados
# ---------------------
# cal_results_ESGBRKGA %>% as.data.table() %>% fwrite("00_Data/Results/Parameter_Calibration/ESGBRKGA_Results_v1.csv", sep=';')
# cal_results_ASGBRKGA %>% as.data.table() %>% fwrite("00_Data/Results/Parameter_Calibration/ASGBRKGA_Results_v1.csv", sep=';')

# -----------
# Comentários
# -----------
# . Em geral, ainda persistem valores bem baixos de silhueta (independente da abordagem de grade)
# . Contudo, novamente valores superiores de silhueta para a abordagem assimétrica, porém não tão significantes levando-se em consideração o tempo de execução
# . ASGBRKGA mais sensível a escolha dos parâmetros 'iter' e 'p' (apesar de variabilidade pequena). Melhores valores: (p = 200, iter = 300)
# . ESGBRKGA praticamente invariante a escolha de parâmetros, inclusive em temros de tempo de execução. Assim, por simplicidade, serão escolhidos os mesmos
#   parâmetros que o ASGBRKGA (p = 200, iter = 300)

############################################
# 5.2 - Análise de Estabilidade do algoritmo    
############################################
# --> 10 replicações x 7 bases x 2 composições de grade
# --> Objetivo: Avaliar estabilidade em termos das diferentes composições de grade

# ----------
# Simulações
# ----------
# sim_ESG <- lapply(test_subset_data, function(i) {replicate(10, STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette"))})
# saveRDS(sim_ESG, '00_Data/Results/Stability_Assessment/ESG/sim_ESG_v3.rds')
# sim_ASG <- lapply(test_subset_data, function(i) {replicate(10, STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette"))})
# saveRDS(sim_ASG, '00_Data/Results/Stability_Assessment/ASG/sim_ASG_v3.rds')

# ----------------------
# Leitura dos resultados
# ----------------------
# sim_ESG <- readRDS('00_Data/Results/Stability_Assessment/ESG/sim_ESG_v3.rds')
# sim_ASG <- readRDS('00_Data/Results/Stability_Assessment/ASG/sim_ASG_v3.rds')

# -----------------------------
# Nomes das instâncias de teste
# -----------------------------
# test_datasets <- c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2")

# -----------------------------------------
# Inspeção visual dos agrupamentos formados
# -----------------------------------------
# lapply(1:length(sim_ESG), function(i) {
#     PlotMDS(cbind(sim_ESG[[i]][,1]$spatial.ppp.obj$x,sim_ESG[[i]][,1]$spatial.ppp.obj$y), sim_ESG[[i]][,1]$cluster, title=test_datasets[i])
# })
# lapply(1:length(sim_ASG), function(i) {
#     PlotMDS(cbind(sim_ASG[[i]][,1]$spatial.ppp.obj$x,sim_ASG[[i]][,1]$spatial.ppp.obj$y), sim_ASG[[i]][,1]$cluster, title=test_datasets[i])
# })

# -------------------------------------------
# Tabela com resultados - Grade Reuglar (ESG)
# -------------------------------------------
# ESG_min_k <- sapply(sim_ESG, function(i) min(unlist(i['k',])))
# ESG_max_k <- sapply(sim_ESG, function(i) max(unlist(i['k',])))
# ESG_med_k <- sapply(sim_ESG, function(i) mean(unlist(i['k',]))) %>% floor()
# ESG_cv_k <- sapply(sim_ESG, function(i) sd(unlist(i['k',]))/mean(unlist(i['k',])))
# ESG_min_is <-  sapply(sim_ESG, function(i) min(unlist(i['score',])))
# ESG_max_is <-  sapply(sim_ESG, function(i) max(unlist(i['score',])))
# ESG_med_is <-  sapply(sim_ESG, function(i) mean(unlist(i['score',])))
# ESG_cv_is <- sapply(sim_ESG, function(i) sd(unlist(i['score',]))/mean(unlist(i['score',])))
# ESG_min_exec <- sapply(sim_ESG, function(i) min(unlist(i['exec',])))
# ESG_max_exec <- sapply(sim_ESG, function(i) max(unlist(i['exec',])))
# ESG_med_exec <- sapply(sim_ESG, function(i) mean(unlist(i['exec',])))
# ESG_cv_exec <- sapply(sim_ESG, function(i) sd(unlist(i['exec',]))/mean(unlist(i['exec',])))
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
# ) %>% fwrite('00_Data/Results/Stability_Assessment/ESG/ESG_Results_v3.csv', sep = ';')

# ---------------------------------------------
# Tabela com resultados - Grade Irregular (ASG)
# ---------------------------------------------
# ASG_min_k <- sapply(sim_ASG, function(i) min(unlist(i['k',])))
# ASG_max_k <- sapply(sim_ASG, function(i) max(unlist(i['k',])))
# ASG_med_k <- sapply(sim_ASG, function(i) mean(unlist(i['k',]))) %>% floor()
# ASG_cv_k <- sapply(sim_ASG, function(i) sd(unlist(i['k',]))/mean(unlist(i['k',])))
# ASG_min_is <-  sapply(sim_ASG, function(i) min(unlist(i['score',]), na.rm = T))
# ASG_max_is <-  sapply(sim_ASG, function(i) max(unlist(i['score',]), na.rm = T))
# ASG_med_is <-  sapply(sim_ASG, function(i) mean(unlist(i['score',]), na.rm = T))
# ASG_cv_is <- sapply(sim_ASG, function(i) sd(unlist(i['score',]), na.rm = T)/mean(unlist(i['score',]), na.rm = T))
# ASG_min_exec <- sapply(sim_ASG, function(i) min(unlist(i['exec',])))
# ASG_max_exec <- sapply(sim_ASG, function(i) max(unlist(i['exec',])))
# ASG_med_exec <- sapply(sim_ASG, function(i) mean(unlist(i['exec',])))
# ASG_cv_exec <- sapply(sim_ASG, function(i) sd(unlist(i['exec',]))/mean(unlist(i['exec',])))
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
# ) %>% fwrite('00_Data/Results/Stability_Assessment/ASG/ASG_Results_v3.csv', sep = ';')

#######################################
# 5.3 - Estudo das Composições de Grade
#######################################
# --> Execução do algoritmo para cada uma das bases de dados com parâmetros default
# --> Para as bases de classificação, comparar número de grupos obtidos e Índice Rand Ajustado

# -------------------
# Nome das instâncias
# -------------------
# test_datasets <- c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2")
# 
# clust_datasets <- c(
#    "200DATA"
#    ,"400P3C"
#    ,"A1"
#    ,"BANKNOTE"
#    ,"BREASTB"
#    ,"CHART"
#    ,"CONCRETEDATA"
#    ,"DBLCA"
#    ,"DOWJONES"
#    ,"FACE"
#    ,"GAMMA400"
#    ,"HAYES-ROTH"
#    ,"INDIAN"
#    ,"MARONNA"
#    ,"NUMBERS2"
#    ,"OUTLIERS"
#    ,"PARKINSONS"
#    ,"RUSPINI"
#    ,"SONAR"
#    ,"SPHERICAL4D3C"
#    ,"SPRDESP"
#    ,"SYNTHETICCONTROL"
#    ,"TRIPADVISOR"
#    ,"UNIFORM400"
#    ,"WAVEFORM21"
# )
# 
# classf_datasets <- c(
#     "AGGREGATION"
#     ,"COMPOUND"
#     ,"GLASS"
#     ,"ECOLI"
#     ,"FLAME"
#     ,"HABERMAN"
#     ,"IRIS"
#     ,"IONOSPHERE"
#     ,"JAIN"
#     ,"NEW-THYROID"
#     ,"PATHBASED"
#     ,"PIMA.INDIANS"
#     ,"R15"
#     ,"SEEDS"
#     ,"SPIRAL"
#     ,"WDBC"
#     ,"WINE"
#     ,"YEAST"
# )
# 
# Instâncias de clusterização + subconjunto de teste
# clust_test_datasets <- c(all_cluster_data, test_subset_data)
 
# ---------------------------------------------------------
# Resultados - Bases de Clusterização + Grade Regular (ESG)
# ---------------------------------------------------------
# clust_ESG_results <- lapply(clust_test_datasets, function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(clust_ESG_results, '00_Data/Results/Clustering/ESG/clust_ESG_v3.rds')
#
# lapply(1:length(clust_ESG_results), function(i) {
#     PlotMDS(
#         cbind(clust_ESG_results[[i]]$spatial.ppp.obj$x, clust_ESG_results[[i]]$spatial.ppp.obj$y),
#         clust_ESG_results[[i]]$cluster,
#         title = paste0("ISM: ", clust_ESG_results[[i]]$score, "\n", "Dataset: ", c(clust_datasets, test_datasets)[i])
#     )
# })
#
# clust_esg_algo <- unlist(sapply(clust_ESG_results, function(i) if(is.null(i$grid.method)) return(i$grids.method) else return(i$grid.method)))
# clust_esg_k <- unlist(sapply(clust_ESG_results, function(i) i$k))
# clust_esg_ich <- unlist(sapply(1:length(clust_ESG_results), function(i) get_ich(clust_test_datasets[[i]], clust_ESG_results[[i]])))
# clust_esg_ism <- unlist(sapply(clust_ESG_results, function(i) i$score))
# clust_esg_ics <- unlist(sapply(clust_ESG_results, function(i) if(is.null(i$grid.ics)) return(i$best.grid.ics) else return(i$grid.ics)))
# 
# data.table(
#     "Base" = c(clust_datasets, test_datasets),
#     "Algo" = clust_esg_algo,
#     "K" = clust_esg_k,
#     "ISM" = clust_esg_ism,
#     "ICH" = clust_esg_ich,
#     "ICS" = clust_esg_ics
# ) %>% fwrite('00_Data/Results/Clustering/ESG/ESG_Results_v3.csv', sep = ';')

# -----------------------------------------------------------
# Resultados - Bases de Clusterização + Grade Irregular (ASG)
# -----------------------------------------------------------
# clust_ASG_results <- lapply(clust_test_datasets, function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(clust_ASG_results, '00_Data/Results/Clustering/ASG/clust_ASG_v3.rds')
# 
# lapply(1:length(clust_ASG_results), function(i) {
#     PlotMDS(
#         cbind(clust_ASG_results[[i]]$spatial.ppp.obj$x, clust_ASG_results[[i]]$spatial.ppp.obj$y),
#         clust_ASG_results[[i]]$cluster,
#         title = paste0("ISM: ", clust_ASG_results[[i]]$score, "\n", "Dataset: ", c(clust_datasets, test_datasets)[i])
#     )
# })
# 
# clust_asg_algo <- unlist(sapply(clust_ASG_results, function(i) i$grids.method))
# clust_asg_k <- unlist(sapply(clust_ASG_results, function(i) i$k))
# clust_asg_ism <- unlist(sapply(clust_ASG_results, function(i) i$score))
# clust_asg_ich <- unlist(sapply(1:35, function(i) get_ich(clust_test_datasets[[i]], clust_ASG_results[[i]])))
# clust_asg_ics <- unlist(sapply(clust_ASG_results, function(i) i$best.grid.ics))
# 
# data.table(
#     "Base" = c(clust_datasets, test_datasets),
#     "Algo" = clust_asg_algo,
#     "K" = clust_asg_k,
#     "ISM" = clust_asg_ism,
#     "ICH" = clust_asg_ich,
#     "ICS" = clust_asg_ics
# ) %>% fwrite('00_Data/Results/Clustering/ASG/ASG_Results_v3.csv', sep = ';')

# ---------------------------------------------------------
# Resultados - Bases de Classificação + Grade Regular (ESG)
# ---------------------------------------------------------
# classf_ESG_results <- lapply(all_classf_data, function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(classf_ESG_results, '00_Data/Results/Classification/ESG/classf_ESG_v3.rds')
# 
# classf_esg_algo <- unlist(sapply(classf_ESG_results, function(i) if(is.null(i$grid.method)) return(i$grids.method) else return(i$grid.method)))
# classf_esg_k <- unlist(sapply(classf_ESG_results, function(i) i$k))
# classf_esg_ism <- unlist(sapply(classf_ESG_results, function(i) i$score))
# classf_esg_ich <- unlist(sapply(1:length(classf_ESG_results), function(i) get_ich(all_classf_data[[i]], classf_ESG_results[[i]])))
# classf_esg_ira <- unlist(sapply(1:18, function(i) get_adj_rand_index(classf_ESG_results[[i]], all_classf_clusters[[i]]$V1)))
# classf_esg_ics <- unlist(sapply(classf_ESG_results, function(i) if(is.null(i$grid.ics)) return(i$best.grid.ics) else return(i$grid.ics)))
# 
# data.table(
#     "Base" = classf_datasets,
#     "Algo" = classf_esg_algo,
#     "K" = classf_esg_k,
#     "ISM" = classf_esg_ism,
#     "ICH" = classf_esg_ich,
#     "Rand" = classf_esg_ira,
#     "ICS" = classf_esg_ics
# ) %>% fwrite('00_Data/Results/Classification/ESG/ESG_Results_v3.csv', sep = ';')

# -----------------------------------------------------------
# Resultados - Bases de Classificação + Grade Irregular (ASG)
# -----------------------------------------------------------
# classf_ASG_results <- lapply(all_classf_data, function(i)
#     STATGDBC(i, alpha=.05, only.ics=0, grid.type="asg", density.test="clarkevans", clust.fobj="silhouette")
# ); saveRDS(classf_ASG_results, '00_Data/Results/Classification/ASG/classf_ASG_v3.rds')
# 
# classf_asg_algo <- unlist(sapply(classf_ASG_results, function(i) i$grids.method))
# classf_asg_k <- unlist(sapply(classf_ASG_results, function(i) i$k))
# classf_asg_ism <- unlist(sapply(classf_ASG_results, function(i) i$score))
# classf_asg_ich <- unlist(sapply(1:length(classf_ASG_results), function(i) get_ich(all_classf_data[[i]], classf_ASG_results[[i]])))
# classf_asg_ira <- unlist(sapply(1:18, function(i) get_adj_rand_index(classf_ASG_results[[i]], all_classf_clusters[[i]]$V1)))
# classf_asg_ics <- unlist(sapply(classf_ASG_results, function(i) i$best.grid.ics))
# 
# data.table(
#     "Base" = classf_datasets,
#     "Algo" = classf_asg_algo,
#     "K" = classf_asg_k,
#     "ISM" = classf_asg_ism,
#     "ICH" = classf_asg_ich,
#     "Rand" = classf_asg_ira,
#     "ICS" = classf_asg_ics
# ) %>% fwrite('00_Data/Results/Classification/ASG/ASG_Results_v3.csv', sep = ';')

# -----------------
# Pré-Processamento
# -----------------
# clust_ESG_results <- readRDS('00_Data/Results/Clustering/ESG/clust_ESG_v3.rds')
# clust_ASG_results <- readRDS('00_Data/Results/Clustering/ASG/clust_ASG_v3.rds')
# classf_ESG_results <- readRDS('00_Data/Results/Classification/ESG/classf_ESG_v3.rds')
# classf_ASG_results <- readRDS('00_Data/Results/Classification/ASG/classf_ASG_v3.rds')
# clust_esg_results <- fread('00_Data/Results/Clustering/ESG/ESG_Results_v3.csv') %>% as.data.frame()
# clust_asg_results <- fread('00_Data/Results/Clustering/ASG/ASG_Results_v3.csv') %>% as.data.frame()
# classf_esg_results <- fread('00_Data/Results/Classification/ESG/ESG_Results_v3.csv') %>% as.data.frame()
# classf_asg_results <- fread('00_Data/Results/Classification/ASG/ASG_Results_v3.csv')%>% as.data.frame()
# clust_esg_results <- clust_esg_results %>% filter(!Base %in% c('DBLCB','INDOCHINACOMBAT', 'NORMAL300')) # Instâncias descartadas p/ efeito de análise
# clust_asg_results <- clust_asg_results %>% filter(!Base %in% c('DBLCB','INDOCHINACOMBAT', 'NORMAL300')) # Instâncias descartadas p/ efeito de análise 
# clust_results <- rbind(clust_esg_results, clust_asg_results)
# classf_results <- rbind(classf_esg_results, classf_asg_results)
# clust_std_data <- lapply(1:35, function(i) apply(c(all_cluster_data, test_subset_data)[[i]], 2, function(j) {(j-mean(j))/sd(j)}))
# clust_tendency <- sapply(1:35, function(i) hopkins(clust_std_data[[i]], 10, 0.5)) # 10 amostras de tamanho 0.5 para avaliação da tendência de agrupamento
# clust_tendency[c(12,20,22)] <- NA # Descarte dos resultados associados as instâncias 'DBLCB','INDOCHINACOMBAT', 'NORMAL300'
# clust_tendency <- clust_tendency[!is.na(clust_tendency)]
# clust_results$clust_tendency <- rep(clust_tendency, 2)
# clust_results %>% as.data.table() %>% fwrite(file = '00_Data/Results/Parameter_Study_Comparisons/CLUST_Final.csv', sep=';')
# classf_std_data <- lapply(1:18, function(i) apply(all_classf_data[[i]], 2, function(j) {(j-mean(j))/sd(j)}))
# classf_tendency <- sapply(1:18, function(i) hopkins(classf_std_data[[i]], 10, 0.5)) # 10 amostras de tamanho 0.5 para avaliação da tendência de agrupamento
# classf_results$clust_tendency <- rep(classf_tendency, 2)
# classf_results %>% as.data.table() %>% fwrite(file = '00_Data/Results/Parameter_Study_Comparisons/CLASSF_Final.csv', sep=';')

# -------------------------------
# Conjuntos de dados para análise
# -------------------------------
clust_results <- fread('00_Data/Results/Grid_Assessment/CLUST_Final.csv') %>% as.data.frame()
classf_results <- fread('00_Data/Results/Grid_Assessment/CLASSF_Final.csv') %>% as.data.frame()

# -----------------------------------------------------------------------------------------
# Tabela auxiliar - Scores de tendência de agrupamento para as instâncias da base (sem EMD)
# -----------------------------------------------------------------------------------------
# rbind(clust_results[1:32,c('Base','clust_tendency')], classf_results[1:18,c('Base','clust_tendency')]) %>% 
#     as.data.table() %>% 
#     fwrite('00_Data/Results/Parameter_Study_Comparisons/HOPKINS_SCORES.csv', sep=';')

# -----------------------------------------------------------------------------------------
# Tabela auxiliar - Scores de tendência de agrupamento para as instâncias da base (com EMD)
# -----------------------------------------------------------------------------------------
# clust_mds_data <- lapply(1:35, function(i) MDSProjection(c(all_cluster_data, test_subset_data)[[i]]))
# clust_tendency <- sapply(1:35, function(i) hopkins(clust_mds_data[[i]], 10, 0.5)) # 10 amostras de tamanho 0.5 para avaliação da tendência de agrupamento
# clust_tendency[c(12,20,22)] <- NA # Descarte dos resultados associados as instâncias 'DBLCB','INDOCHINACOMBAT', 'NORMAL300'
# clust_tendency <- clust_tendency[!is.na(clust_tendency)]
# clust_hopkins_mds_ics <- clust_results %>% select(c('Base','Algo','ICS')) %>% cbind(clust_tendency = rep(clust_tendency,2))
# clust_hopkins_mds_ics$grid_type <- c(rep('ESG',32),rep('ASG',32))
# classf_mds_data <- lapply(1:18, function(i) MDSProjection(all_classf_data[[i]]))
# classf_tendency <- sapply(1:18, function(i) hopkins(classf_mds_data[[i]], 10, 0.5)) # 10 amostras de tamanho 0.5 para avaliação da tendência de agrupamento
# classf_hopkins_mds_ics <- classf_results %>% select(c('Base','Algo','ICS')) %>% cbind(clust_tendency = rep(classf_tendency,2))
# classf_hopkins_mds_ics$grid_type <- c(rep('ESG',18),rep('ASG',18))
# hopkins_mds_ics <- rbind(clust_hopkins_mds_ics, classf_hopkins_mds_ics)
# hopkins_mds_ics %>% as.data.table() %>% fwrite('00_Data/Results/Parameter_Study_Comparisons/HOPKINS_MDS_SCORES.csv', sep=';')

# ----------------------------------
# Análise - Performance geral do EMD
# ----------------------------------
# hopkins_raw <- fread('00_Data/Results/Grid_Assessment/HOPKINS_SCORES.csv') %>% as.data.frame() %>% select(Base,clust_tendency)
# hopkins_mds <- fread('00_Data/Results/Grid_Assessment/HOPKINS_MDS_SCORES.csv') %>% as.data.frame() %>% select(Base,clust_tendency)
# colnames(hopkins_mds) <- c('Base','clust_tendecy_emd')
# df <- left_join(hopkins_raw, hopkins_mds, on='Base') %>% unique() %>% arrange(Base)
# ggplot(data = df, aes(x=clust_tendency, y=clust_tendecy_emd)) +
#     geom_point(size=3, col='#f97306') +
#     geom_abline(slope=1, intercept=0, lwd=0.75, col='#e50000') +
#     labs(x="\nEstatística de Hopkins (Dados Originais)\n",  y="\nEstatística de Hopkins (Projeção EMD)\n") +
#     ggtitle(expression(atop('Estatística de Hopkins - Dados Originais x Projeção EMD (N = 50)'))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(0,1)) +
#     scale_x_continuous(limits=c(0,1)) +
#     theme(plot.title = element_text(hjust = 0.5))
# cor(df$clust_tendency, df$clust_tendecy_emd)

# -----------------------------------------------------
# Análise - Melhor composição de grade em termos de ISM
# -----------------------------------------------------
# ism_clust_analysis <- clust_results %>% select(c('Base','Algo','ISM','clust_tendency')); ism_clust_analysis$grid_type <- c(rep('ESG',32),rep('ASG',32))
# ism_clust_analysis <- ism_clust_analysis %>% filter(!Base %in% c('BUPA','OUTLIERS','BREASTB')) # Descarte das instâncias com (ISM = NA)
# ism_classf_analysis <- classf_results %>% select(c('Base','Algo','ISM','clust_tendency')); ism_classf_analysis$grid_type <- c(rep('ESG',18),rep('ASG',18))
# 
# boxplot(ism_clust_analysis$ISM ~ ism_clust_analysis$grid_type)$stats; boxplot(ism_clust_analysis$ISM ~ ism_clust_analysis$grid_type)$out
# boxplot(ism_classf_analysis$ISM ~ ism_classf_analysis$grid_type)$stats; boxplot(ism_classf_analysis$ISM ~ ism_classf_analysis$grid_type)$out
# mean(ism_clust_analysis$ISM); mean(ism_clust_analysis$clust_tendency)
# mean(ism_classf_analysis$ISM); mean(ism_classf_analysis$clust_tendency)
# 
# ggplot(data = ism_clust_analysis, aes(y=ISM, fill=factor(grid_type)))+
#     geom_boxplot() +
#     labs(x="\n\n",  y="\nÍndice de Silhueta Médio\n", fill='Grade') +
#     ggtitle(expression(atop('Distribuição do Índice de Silhueta Médio por Composição de Grade', atop('Instâncias clusterização (N = 29)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(-1,1)) +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# 
# ggplot(data = ism_classf_analysis, aes(y=ISM, fill=grid_type))+
#     geom_boxplot() +
#     labs(x="\n\n",  y="\nÍndice de Silhueta Médio\n", fill='Grade') +
#     ggtitle(expression(atop('Distribuição do Índice de Silhueta Médio por Composição de Grade', atop('Instâncias classificação (N = 18)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(-1,1)) +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# 
# ism_clust_analysis %>% as.data.table() %>% fwrite('00_Data/Results/Parameter_Study_Comparisons/ISM_CLSUT_Dataset.csv')
# ism_classf_analysis %>% as.data.table() %>% fwrite('00_Data/Results/Parameter_Study_Comparisons/ISM_CLASSF_Dataset.csv')

# -----------------------------------------------------
# Análise - Melhor composição de grade em termos de IRA
# -----------------------------------------------------
# ira_classf_analysis <- classf_results %>% select(c('Base','Algo','Rand','clust_tendency')); ira_classf_analysis$grid_type <- c(rep('ESG',18),rep('ASG',18))
# 
# boxplot(ira_classf_analysis$Rand ~ ira_classf_analysis$grid_type)$stats; boxplot(ira_classf_analysis$Rand ~ ira_classf_analysis$grid_type)$out
# mean(ira_classf_analysis$Rand); mean(ira_classf_analysis$clust_tendency)
# 
# ggplot(data = ira_classf_analysis, aes(y=Rand, fill=grid_type))+
#     geom_boxplot() +
#     labs(x="\n\n",  y="\nÍndice de Rand Ajustado\n", fill='Grade') +
#     ggtitle(expression(atop('Distribuição do Índice de Rand Ajustado por Composição de Grade', atop('Instâncias classificação (N = 18)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(-0.05,1)) +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# 
# ira_classf_analysis %>% as.data.table() %>% fwrite('00_Data/Results/Parameter_Study_Comparisons/IRA_CLASSF_Dataset.csv')

# ----------------------------------------------------------------------------
# Análise - Índices de Silhueta e Rand Ajustado x ICS, por composição de grade
# ----------------------------------------------------------------------------
# clust_ism_ics <- clust_results %>% select(c('Base','ISM','ICS','clust_tendency'))
# clust_ism_ics$grid_type <- c(rep('ESG',32),rep('ASG',32))
# clust_ism_ics <- clust_ism_ics %>% filter(!Base %in% c('BUPA','OUTLIERS','BREASTB'))
# classf_ism_ics <- classf_results %>% select(c('Base','Rand','ICS','clust_tendency'))
# classf_ism_ics$grid_type <- c(rep('ESG',18),rep('ASG',18))
# 
# ggplot(data = clust_ism_ics, aes(x=ICS, y=ISM, colour=grid_type)) +
#     geom_point(size=3) +
#     labs(x="\nÍndice de Tamanho de Cluster\n",  y="\nÍndice de Silhueta Médio\n", colour='Grade') +
#     ggtitle(expression(atop('Índice de Silhueta Médio x Índice de Tamanho de Cluster por Composição de Grade', atop('Instâncias clusterização (N = 29)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(-1,1)) +
#     scale_x_continuous(limits=c(0,max(clust_ism_ics$ICS))) +
#     theme(plot.title = element_text(hjust = 0.5))
# clust_ism_ics %>% group_by(grid_type) %>% summarise(Corr=cor(ISM,ICS))
# 
# ggplot(data = classf_ism_ics, aes(x=ICS, y=Rand, colour=grid_type)) +
#     geom_point(size=3) +
#     labs(x="\nÍndice de Tamanho de Cluster\n",  y="\nÍndice de Rand Ajustado\n", colour='Grade') +
#     ggtitle(expression(atop('Índice de Rand Ajustado x Índice de Tamanho de Cluster por Composição de Grade', atop('Instâncias classificação (N = 18)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(-0.1,1)) +
#     scale_x_continuous(limits=c(0,max(classf_ism_ics$ICS))) +
#     theme(plot.title = element_text(hjust = 0.5))
# classf_ism_ics %>% group_by(grid_type) %>% summarise(Corr=cor(Rand,ICS))

# ---------------------------------------------------------------------------
# Análise - Exemplos de composições de grade ESG/ASG (ilustrar o viés do ICS)
# ---------------------------------------------------------------------------
# esg_results <- readRDS('00_Data/Results/Clustering/ESG/clust_ESG_v3.rds')
# asg_results <- readRDS('00_Data/Results/Clustering/ASG/clust_ASG_v3.rds')
# 
# bases <- c("200DATA" ,"FACE" ,"MARONNA" ,"OUTLIERS")
# chosen <- c(1,11,17,19)
# 
# lapply(1:4, function(i) PlotGrid(esg_results[[chosen[i]]]$spatial.ppp.obj, esg_results[[chosen[i]]]$grid, F, bases[i]))
# lapply(1:4, function(i) PlotGrid(asg_results[[chosen[i]]]$spatial.ppp.obj, asg_results[[chosen[i]]]$best.grid, T, bases[i]))

# -----------------------------------------------------------------------
# (DEP) Análise - Tendência de Agrupamento x ICS, por composição de grade
# -----------------------------------------------------------------------
# hopkins_mds_ics <- fread('00_Data/Results/Grid_Assessment/HOPKINS_MDS_SCORES.csv') %>% as.data.frame()
# 
# ggplot(data = hopkins_mds_ics, aes(x=ICS, y=clust_tendency, colour=grid_type)) +
#     geom_point(size=3) +
#     labs(x="\nÍndice de Tamanho de Cluster\n",  y="\nEstatística de Hopkins\n", colour='Grade') +
#     ggtitle(expression(atop('Estatística de Hopkins x Índice de Tamanho de Cluster por Composição de Grade',
#                             atop('Considera projeções obtidas via EMD (N = 50)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(0,1)) +
#     scale_x_continuous(limits=c(0,max(hopkins_mds_ics$ICS))) +
#     theme(plot.title = element_text(hjust = 0.5))
# hopkins_mds_ics %>% group_by(grid_type) %>% summarise(Corr=cor(clust_tendency,ICS))

##############################
# 5.4 - Critérios Estatísticos
##############################
# --> Objetivo: Definir o melhor conjunto de parâmetros associados
# aos critérios estatísticos do STATGDBC, condicional à composição
# de grade ESG.

# --------------------------------------------------------------------------------------
# Análise - (K, Hopkins, ISM, ICH, ICS) x only.ics p/ instâncias do subconjunto de teste
# --------------------------------------------------------------------------------------
# sim_only_ics_yes <- lapply(test_subset_data, function(i)  STATGDBC(i, only.ics = 1))
# saveRDS(sim_only_ics_yes, '00_Data/Results/Statistical_Parameters_Analysis/SIM_ONLYICS_1.rds')
# sim_only_ics_no <- lapply(test_subset_data, function(i)  STATGDBC(i, only.ics = 0))
# saveRDS(sim_only_ics_no, '00_Data/Results/Statistical_Parameters_Analysis/SIM_ONLYICS_0.rds')
#
# sim_only_ics_yes <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/SIM_ONLYICS_1.rds')
# sim_only_ics_no <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/SIM_ONLYICS_0.rds')
# 
# test_datasets <- c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2")
# 
# sim_only_ics_yes <- data.frame(
#     Base = test_datasets,
#     K = unlist(sapply(1:7, function(i) sim_only_ics_yes[[i]]$k )),
#     Hopkins =  fread('00_Data/Results/Grid_Assessment/HOPKINS_SCORES.csv') %>%
#         as.data.frame() %>%
#         select(Base,clust_tendency) %>%
#         filter(Base %in% test_datasets) %>%
#         select(clust_tendency) %>%
#         unlist() %>%
#         as.numeric(),
#     ISM = unlist(sapply(1:7, function(i) sim_only_ics_yes[[i]]$score )),
#     ICH = unlist(sapply(1:7, function(i) get_ich(test_subset_data[[i]], sim_only_ics_yes[[i]]))),
#     ICS = unlist(sapply(1:7, function(i) sim_only_ics_yes[[i]]$grid.ics)),
#     only.ics = rep(1, 7)
# )
# 
# sim_only_ics_no <- data.frame(
#     Base = test_datasets,
#     K = unlist(sapply(1:7, function(i) sim_only_ics_no[[i]]$k )),
#     Hopkins =  fread('00_Data/Results/Grid_Assessment/HOPKINS_SCORES.csv') %>%
#         as.data.frame() %>%
#         select(Base,clust_tendency) %>%
#         filter(Base %in% test_datasets) %>%
#         select(clust_tendency) %>%
#         unlist() %>%
#         as.numeric(),
#     ISM = unlist(sapply(1:7, function(i) sim_only_ics_no[[i]]$score )),
#     ICH = unlist(sapply(1:7, function(i) get_ich(test_subset_data[[i]], sim_only_ics_no[[i]]))),
#     ICS = unlist(sapply(1:7, function(i) sim_only_ics_no[[i]]$grid.ics)),
#     only.ics = rep(0, 7)
# )
# 
# sim_only_ics <- rbind(sim_only_ics_yes, sim_only_ics_no)
# sim_only_ics %>% as.data.table() %>% fwrite('00_Data/Results/Statistical_Parameters_Analysis/SIM_ONLYICS_FINAL.csv', sep=';')

# --------------------------------------------------------------------------------
# Análise - (K,ISM, ICH, ICS) x density.test p/ instâncias do subconjunto de teste
# --------------------------------------------------------------------------------
# sim_density_test_hopkins <- lapply(test_subset_data, function(i)  STATGDBC(i, only.ics = 1, density.test = 'hopkins'))
# saveRDS(sim_density_test_hopkins, '00_Data/Results/Statistical_Parameters_Analysis/SIM_DENSITYTEST_HOPKINS.rds')
# sim_density_test_clarkevans <- lapply(test_subset_data, function(i)  STATGDBC(i, only.ics = 1, density.test = 'clarkevans'))
# saveRDS(sim_density_test_clarkevans, '00_Data/Results/Statistical_Parameters_Analysis/SIM_DENSITYTEST_CLARKEVANS.rds')
# 
# sim_density_test_hopkins <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/SIM_DENSITYTEST_HOPKINS.rds')
# sim_density_test_clarkevans <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/SIM_DENSITYTEST_CLARKEVANS.rds')
# 
# test_datasets <- c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2")
# 
# sim_density_test_hopkins <- data.frame(
#     Base = test_datasets,
#     K = unlist(sapply(1:7, function(i) sim_density_test_hopkins[[i]]$k )),
#     ISM = unlist(sapply(1:7, function(i) sim_density_test_hopkins[[i]]$score )),
#     ICH = unlist(sapply(1:7, function(i) get_ich(test_subset_data[[i]], sim_density_test_hopkins[[i]]))),
#     ICS = unlist(sapply(1:7, function(i) sim_density_test_hopkins[[i]]$grid.ics)),
#     density.test = rep('Hopkins', 7)
# )
# 
# sim_density_test_clarkevans <- data.frame(
#     Base = test_datasets,
#     K = unlist(sapply(1:7, function(i) sim_density_test_clarkevans[[i]]$k )),
#     ISM = unlist(sapply(1:7, function(i) sim_density_test_clarkevans[[i]]$score )),
#     ICH = unlist(sapply(1:7, function(i) get_ich(test_subset_data[[i]], sim_density_test_clarkevans[[i]]))),
#     ICS = unlist(sapply(1:7, function(i) sim_density_test_clarkevans[[i]]$grid.ics)),
#     density.test = rep('Clark-Evans', 7)
# )
# 
# sim_density_test <- rbind(sim_density_test_hopkins, sim_density_test_clarkevans)
# sim_density_test %>% as.data.table() %>% fwrite('00_Data/Results/Statistical_Parameters_Analysis/SIM_DENSITYTEST_FINAL.csv', sep=';')

# ------------------------------------------------------------------
# Análise - (K', K, IRA) x clust.fobj p/ instâncias de classificação
# ------------------------------------------------------------------
# sim_clust_fobj_ich <- lapply(all_classf_data, function(i)  STATGDBC(i, only.ics = 1, density.test = 'hopkins', clust.fobj = 'calinski'))
# saveRDS(sim_clust_fobj_ich, '00_Data/Results/Statistical_Parameters_Analysis/SIM_CLUSTFOBJ_ICH.rds')
# sim_clust_fobj_ism <- lapply(all_classf_data, function(i)  STATGDBC(i, only.ics = 1, density.test = 'hopkins', clust.fobj = 'silhouette'))
# saveRDS(sim_clust_fobj_ism, '00_Data/Results/Statistical_Parameters_Analysis/SIM_CLUSTFOBJ_ISM.rds')
# 
# sim_clust_fobj_ich <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/SIM_CLUSTFOBJ_ICH.rds')
# sim_clust_fobj_ism <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/SIM_CLUSTFOBJ_ISM.rds')
# 
# classf_datasets <- c(
#     "AGGREGATION"
#     ,"COMPOUND"
#     ,"GLASS"
#     ,"ECOLI"
#     ,"FLAME"
#     ,"HABERMAN"
#     ,"IRIS"
#     ,"IONOSPHERE"
#     ,"JAIN"
#     ,"NEW-THYROID"
#     ,"PATHBASED"
#     ,"PIMA.INDIANS"
#     ,"R15"
#     ,"SEEDS"
#     ,"SPIRAL"
#     ,"WDBC"
#     ,"WINE"
#     ,"YEAST"
# )
# 
# sim_clust_fobj_ich <- data.frame(
#     Base = classf_datasets,
#     K = unlist(sapply(1:18, function(i) sim_clust_fobj_ich[[i]]$k )),
#     IRA = unlist(sapply(1:18, function(i) get_adj_rand_index(sim_clust_fobj_ich[[i]], all_classf_clusters[[i]]$V1))),
#     clust.fobj = rep('Calinski-Harabasz', 18)
# )
# 
# sim_clust_fobj_ism <- data.frame(
#     Base = classf_datasets,
#     K = unlist(sapply(1:18, function(i) sim_clust_fobj_ism[[i]]$k )),
#     IRA = unlist(sapply(1:18, function(i) get_adj_rand_index(sim_clust_fobj_ism[[i]], all_classf_clusters[[i]]$V1))),
#     clust.fobj = rep('Silhueta Média', 18)
# )
# 
# sim_clust_fobj <- rbind(sim_clust_fobj_ich, sim_clust_fobj_ism)
# sim_clust_fobj %>% as.data.table() %>% fwrite('00_Data/Results/Statistical_Parameters_Analysis/SIM_CLUSTFOBJ_FINAL.csv', sep=';')

# -----------------------------------------------------------------------
# Análise - Tabela final STATGDBC (parâmetros ajustados) para as 50 bases
# -----------------------------------------------------------------------
# all_datasets <- c(all_cluster_data,all_classf_data,test_subset_data)
# statgdbc_final <- lapply(all_datasets, function(i) STATGDBC(i, only.ics = 1, grid.type = 'esg', density.test = 'hopkins', clust.fobj = 'silhouette'))
# statgdbc_final %>% saveRDS('00_Data/Results/Statistical_Parameters_Analysis/STATGDBC_FINAL.rds')
#
# statgdbc_final <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/STATGDBC_FINAL.rds')
# 
# datasets_names <- c(
#     "200DATA"
#     ,"400P3C"
#     ,"A1"
#     ,"BANKNOTE"
#     ,"BREASTB"
#     ,"CHART"
#     ,"CONCRETEDATA"
#     ,"DBLCA"
#     ,"DOWJONES"
#     ,"FACE"
#     ,"GAMMA400"
#     ,"HAYES-ROTH"
#     ,"INDIAN"
#     ,"MARONNA"
#     ,"NUMBERS2"
#     ,"OUTLIERS"
#     ,"PARKINSONS"
#     ,"RUSPINI"
#     ,"SONAR"
#     ,"SPHERICAL4D3C"
#     ,"SPRDESP"
#     ,"SYNTHETICCONTROL"
#     ,"TRIPADVISOR"
#     ,"UNIFORM400"
#     ,"WAVEFORM21"
#     ,"AGGREGATION"
#     ,"COMPOUND"
#     ,"GLASS"
#     ,"ECOLI"
#     ,"FLAME"
#     ,"HABERMAN"
#     ,"IONOSPHERE"
#     ,"IRIS"
#     ,"JAIN"
#     ,"NEW-THYROID"
#     ,"PATHBASED"
#     ,"PIMA.INDIANS"
#     ,"R15"
#     ,"SEEDS"
#     ,"SPIRAL"
#     ,"WDBC"
#     ,"WINE"
#     ,"YEAST"
#     ,"2-FACE"
#     ,"BROKEN-RING"
#     ,"BUPA"
#     ,"FORESTFIRES"
#     ,"GAUSS9"
#     ,"UNIFORM700"
#     ,"VOWEL2"
# )
# 
# statgdbc_final <- data.frame(
#     Base = datasets_names,
#     K = unlist(sapply(statgdbc_final, function(i) i$k)),
#     ISM = unlist(sapply(statgdbc_final, function(i) i$score))
# )
# 
# statgdbc_final %>% as.data.table() %>% fwrite('00_Data/Results/STATGDBC_DBSCAN_Comparison/STATGDBC_FINAL.csv', sep=';')

# ---------
# Conclusão
# ---------
# --> only.ics = 1
# --> density.test = Hopkins
# --> clust.fobj = Silhueta (diferença pouco significante, varia por instância, porém tem interpretação mais simples)


#####################################
# 5.5 - Comparações STATGDBC x DBSCAN
#####################################
# --> Objetivo: Comparar as soluções produzidas pelos algoritmos em termos
# do Índice de Rand Ajustado, Índice de Silhueta e Número de grupos, além
# da inspeção visual de algumas das soluções de agrupamento produzidas.

# ---------------------------------------------------------------
# Estimação do parâmetro 'epsilon' (raio de vizinhança) do DBSCAN
# ---------------------------------------------------------------
# all_datasets_std <- c(all_cluster_data,all_classf_data,test_subset_data) %>% lapply(function(i) apply(i, 2, function(j) {(j-mean(j))/sd(j)}))
# all_datasets_emd <- c(all_cluster_data,all_classf_data,test_subset_data) %>% lapply(function(i) MDSProjection(i))
# eps_hat <- lapply(all_datasets_emd, function(i) mean(dbscan::kNNdist(i, k=3))) %>% unlist()
# eps_hat %>% saveRDS('00_Data/Results/STATGDBC_DBSCAN_Comparison/Epsilon_Estimates_V3.rds')
# eps_hat <- readRDS('00_Data/Results/STATGDBC_DBSCAN_Comparison/Epsilon_Estimates_V3.rds')

# -------------------
# Resultados - DBSCAN
# -------------------
# dbscan_results <- lapply(1:50, function(i) dbscan::dbscan(MDSProjection(all_datasets_emd[[i]]), eps_hat[i], minPts=5))
# dbscan_results %>% saveRDS('00_Data/Results/STATGDBC_DBSCAN_Comparison/DBSCAN_Results_V3.rds')
# dbscan_results <- readRDS('00_Data/Results/STATGDBC_DBSCAN_Comparison/DBSCAN_Results_V3.rds')

# -----------------------------------------------
# Análise - Tabela Geral com resultados do DBSCAN
# -----------------------------------------------
# datasets_names <- c(
#     "200DATA"
#     ,"400P3C"
#     ,"A1"
#     ,"BANKNOTE"
#     ,"BREASTB"
#     ,"CHART"
#     ,"CONCRETEDATA"
#     ,"DBLCA"
#     ,"DOWJONES"
#     ,"FACE"
#     ,"GAMMA400"
#     ,"HAYES-ROTH"
#     ,"INDIAN"
#     ,"MARONNA"
#     ,"NUMBERS2"
#     ,"OUTLIERS"
#     ,"PARKINSONS"
#     ,"RUSPINI"
#     ,"SONAR"
#     ,"SPHERICAL4D3C"
#     ,"SPRDESP"
#     ,"SYNTHETICCONTROL"
#     ,"TRIPADVISOR"
#     ,"UNIFORM400"
#     ,"WAVEFORM21"
#     ,"AGGREGATION"
#     ,"COMPOUND"
#     ,"GLASS"
#     ,"ECOLI"
#     ,"FLAME"
#     ,"HABERMAN"
#     ,"IONOSPHERE"
#     ,"IRIS"
#     ,"JAIN"
#     ,"NEW-THYROID"
#     ,"PATHBASED"
#     ,"PIMA.INDIANS"
#     ,"R15"
#     ,"SEEDS"
#     ,"SPIRAL"
#     ,"WDBC"
#     ,"WINE"
#     ,"YEAST"
#     ,"2-FACE"
#     ,"BROKEN-RING"
#     ,"BUPA"
#     ,"FORESTFIRES"
#     ,"GAUSS9"
#     ,"UNIFORM700"
#     ,"VOWEL2"
# )
# 
# dbscan_final <- data.frame(
#     Base = datasets_names,
#     K = unlist(sapply(dbscan_results, function(i) uniqueN(i$cluster[i$cluster != 0]) )),
#     ISM = unlist(sapply(1:50, function(i) get_ism(all_datasets_emd[[i]], dbscan_results[[i]], is.statgdbc = F)))
# )
# 
# dbscan_final %>% as.data.table() %>% fwrite('00_Data/Results/STATGDBC_DBSCAN_Comparison/DBSCAN_Results.csv', sep=';')

# -------------------------------------------------------------------------------
# Análise - Boxplot comparativo STATGDBC x DBSCAN em termos de ISM de forma geral
# -------------------------------------------------------------------------------
# dbscan_final <- fread('00_Data/Results/STATGDBC_DBSCAN_Comparison/DBSCAN_Results.csv') %>% as.data.frame()
# statgdbc_final <- fread('00_Data/Results/STATGDBC_DBSCAN_Comparison/STATGDBC_FINAL.csv') %>% as.data.frame()
# dbscan_final$algo <- rep('DBSCAN',50); statgdbc_final$algo <- rep('STATGDBC',50)
# dbscan_clean <- dbscan_final %>% filter(!Base %in% c("BREASTB","CHART","DBLCA","IONOSPHERE","NEW-THYROID","PARKINSONS","SONAR"))
# statgdbc_clean <- statgdbc_final %>% filter(!Base %in% c('BREASTB')) # Calibração dos parâmetros estatísticos ajudou
# df <- rbind(dbscan_clean, statgdbc_clean)
# 
# boxplot(df$ISM ~ df$algo)$stats; boxplot(df$ISM ~ df$algo)$out
# df %>% group_by(algo) %>% summarise(MeanISM = mean(ISM, na.rm=T))
# 
# ggplot(data = df, aes(y=ISM, fill=factor(algo)))+
#     geom_boxplot() +
#     labs(x="\nAlgoritmo\n",  y="\nÍndice de Silhueta Médio\n", fill='Algoritmo') +
#     ggtitle(expression(atop('Comparação STATGBDC x DBSCAN (Índice de Silhueta Médio)', atop('DBSCAN (N = 43) | STATGDBC (N = 49)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(-1,1)) +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# ------------------------------------------------------------------------------------------
# Análise - Boxplot comparativo STATGDBC x DBSCAN em termos de IRA p/ bases de classificação
# ------------------------------------------------------------------------------------------
# dbscan_results <- readRDS('00_Data/Results/STATGDBC_DBSCAN_Comparison/DBSCAN_Results_V3.rds')
# statgdbc_results <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/STATGDBC_FINAL.rds')
# 
# dbscan_final <- fread('00_Data/Results/STATGDBC_DBSCAN_Comparison/DBSCAN_Results.csv') %>% as.data.frame()
# statgdbc_final <- fread('00_Data/Results/STATGDBC_DBSCAN_Comparison/STATGDBC_FINAL.csv') %>% as.data.frame()
# 
# dbscan_final$algo <- rep('DBSCAN',50); statgdbc_final$algo <- rep('STATGDBC',50)
# dbscan_final$data_type <- c(rep('CLUST',25),rep('CLASSF',18),rep('CLUST',7))
# statgdbc_final$data_type <- c(rep('CLUST',25),rep('CLASSF',18),rep('CLUST',7))
# df_classf <- rbind(dbscan_final, statgdbc_final) %>% filter(data_type == 'CLASSF') %>% arrange(Base)
# dbscan_rand_results <- sapply(1:18, function(i) get_adj_rand_index(dbscan_results[26:43][[i]], all_classf_clusters[[i]]$V1))
# statgdbc_rand_results <- sapply(1:18, function(i) get_adj_rand_index(statgdbc_results[26:43][[i]], all_classf_clusters[[i]]$V1))
# df_classf$IRA <- c(dbscan_rand_results,statgdbc_rand_results)
# 
# boxplot(df_classf$IRA ~ df_classf$algo)$stats; boxplot(df_classf$IRA ~ df_classf$algo)$out
# df_classf %>% group_by(algo) %>% summarise(MeanIRA = mean(IRA))
# 
# ggplot(data = df_classf, aes(y=IRA, fill=factor(algo)))+
#     geom_boxplot() +
#     labs(x="\nAlgoritmo\n",  y="\nÍndice de Rand Ajustado\n", fill='Algoritmo') +
#     ggtitle(expression(atop('Comparação STATGBDC x DBSCAN (Índice de Rand Ajustado)', atop('(N = 18)')))) +
#     theme_classic() +
#     scale_y_continuous(limits=c(-0.15,1)) +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# --------------------------------------------------------------------------------------------
# Análise - Comparativo STATGDBC x DBSCAN em termos do número de grupos p/ bases classificação
# --------------------------------------------------------------------------------------------
# df_classf$K_prime <- rep(sapply(all_classf_clusters, uniqueN), each=2)
# df_classf %>% group_by(algo) %>% summarise(MeanK = mean(K), MAE_K = mean(abs(K-K_prime)), SD_K = sd(K))
# df_classf %>% as.data.table() %>% fwrite('00_Data/Results/STATGDBC_DBSCAN_Comparison/CLASSF_Results.csv')

# ------------------------------------------------------------------------
# Análise - Inspeção visual das soluções de agrupamento do DBSCAN/STATGDBC
# ------------------------------------------------------------------------
dbscan_results <- readRDS('00_Data/Results/STATGDBC_DBSCAN_Comparison/DBSCAN_Results_V3.rds')
statgdbc_results <- readRDS('00_Data/Results/Statistical_Parameters_Analysis/STATGDBC_FINAL.rds')

all_datasets <- c(all_cluster_data, all_classf_data, test_subset_data)
all_names <- c(
    # Instâncias de Clusterização (25)
    "200DATA"
    ,"400P3C"
    ,"A1"
    ,"BANKNOTE"
    ,"BREASTB"
    ,"CHART"
    ,"CONCRETEDATA"
    ,"DBLCA"
    ,"DOWJONES"
    ,"FACE"
    ,"GAMMA400"
    ,"HAYES-ROTH"
    ,"INDIAN"
    ,"MARONNA"
    ,"NUMBERS2"
    ,"OUTLIERS"
    ,"PARKINSONS"
    ,"RUSPINI"
    ,"SONAR"
    ,"SPHERICAL4D3C"
    ,"SPRDESP"
    ,"SYNTHETICCONTROL"
    ,"TRIPADVISOR"
    ,"UNIFORM400"
    ,"WAVEFORM21"
    # Instâncias de Classificação (18)
    ,"AGGREGATION"
    ,"COMPOUND"
    ,"GLASS"
    ,"ECOLI"
    ,"FLAME"
    ,"HABERMAN"
    ,"IONOSPHERE"
    ,"IRIS"
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
    # Instâncias de Teste (7)
    ,"2-FACE"
    ,"BROKEN-RING"
    ,"BUPA"
    ,"FORESTFIRES"
    ,"GAUSS9"
    ,"UNIFORM700"
    ,"VOWEL2"
)

# -------------------------------
# STATGDBC - Todos os resultados:
# -------------------------------
lapply(1:50, function(i) {
    PlotMDS(
        cbind(statgdbc_results[[i]]$spatial.ppp.obj$x, statgdbc_results[[i]]$spatial.ppp.obj$y),
        statgdbc_results[[i]]$cluster,
        title = paste0("STATGDBC\n ISM: ", statgdbc_results[[i]]$score, "\nInstância: ", all_names[i])
    )
})

# -----------------------------
# DBSCAN - Todos os resultados:
# ----------------------------
lapply(1:50, function(i) {
    PlotMDS(
        cbind(statgdbc_results[[i]]$spatial.ppp.obj$x, statgdbc_results[[i]]$spatial.ppp.obj$y),
        dbscan_results[[i]]$cluster,
        title = paste0("DBSCAN\n ISM: ", dbscan_results[[i]]$score, "\nInstância: ", all_names[i])
    )
})

# -------------------------------------------
# Detalhamento - Instâncias de classificação:
# -------------------------------------------

# AGGREGATION
PlotMDS(
    cbind(statgdbc_results[[26]]$spatial.ppp.obj$x, statgdbc_results[[26]]$spatial.ppp.obj$y),
    statgdbc_results[[26]]$cluster,
    title = paste0("STATGDBC\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[26])
)
PlotMDS(
    cbind(statgdbc_results[[26]]$spatial.ppp.obj$x, statgdbc_results[[26]]$spatial.ppp.obj$y),
    dbscan_results[[26]]$cluster,
    title = paste0("DBSCAN\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[26])
)

# FLAME
PlotMDS(
    cbind(statgdbc_results[[30]]$spatial.ppp.obj$x, statgdbc_results[[30]]$spatial.ppp.obj$y),
    statgdbc_results[[30]]$cluster,
    title = paste0("STATGDBC\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[30])
)
PlotMDS(
    cbind(statgdbc_results[[30]]$spatial.ppp.obj$x, statgdbc_results[[30]]$spatial.ppp.obj$y),
    dbscan_results[[30]]$cluster,
    title = paste0("DBSCAN\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[30])
)

# HABERMAN
PlotMDS(
    cbind(statgdbc_results[[31]]$spatial.ppp.obj$x, statgdbc_results[[31]]$spatial.ppp.obj$y),
    statgdbc_results[[31]]$cluster,
    title = paste0("STATGDBC\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[31])
)
PlotMDS(
    cbind(statgdbc_results[[31]]$spatial.ppp.obj$x, statgdbc_results[[31]]$spatial.ppp.obj$y),
    dbscan_results[[31]]$cluster,
    title = paste0("DBSCAN\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[31])
)

# YEAST
PlotMDS(
    cbind(statgdbc_results[[43]]$spatial.ppp.obj$x, statgdbc_results[[43]]$spatial.ppp.obj$y),
    statgdbc_results[[43]]$cluster,
    title = paste0("STATGDBC\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[43])
)
PlotMDS(
    cbind(statgdbc_results[[43]]$spatial.ppp.obj$x, statgdbc_results[[43]]$spatial.ppp.obj$y),
    dbscan_results[[43]]$cluster,
    title = paste0("DBSCAN\n Índice de Rand Ajustado: ", "\nInstância: ", all_names[43])
)

