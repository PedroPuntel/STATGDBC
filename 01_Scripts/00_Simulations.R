# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: Script contendo as simulações/análises contempladas no capítulo 5 da monografia
# Encoding: UTF-8

# Imports
source("01_Scripts/01_STATGDBC.R")
source("01_Scripts/02_PlotMDS.R")
source("01_Scripts/04_PlotGrid.R")
options(digits = 3)

# Pacotes
require("data.table")
require("pdfCluster")

####################
# Rotinas auxiliares
####################

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
# --> 3 replicações x 7 bases de teste x 3 tamanhos de pop. x 3 tamanhos de ger. x 2 composições de grade: 378 execuções
# --> Objetivo: Por limitação de tempo, calibração dos parâmetros mais impactantes do BRKGA (demais seguem sugestão da literatura)

# Produto cartesiano dos parâmetros
params <- expand.grid(
    c(1:7),          # Conjuntos de dados
    c(100,200,300),  # Tamanho de População
    c(200,300,500),  # Número de gerações
    c("esg","asg")   # Composição de grade
) %>% as.data.frame()
colnames(params) <- c("datasets","p","iter","grid.type")

# Simulações (Approx 6hrs)
# cal_results <- lapply(1:nrow(params), function(i) {
#     replicate(3, STATGDBC(data = test_subset_data[params$datasets[i]], grid.type = params$grid.type[i], p = params$p[i], iter = params$iter[i]))
# }); saveRDS(cal_results, "00_Data/Results/Parameter_Calibration/Calibration_Results_v1.rds")
cal_results <- readRDS("00_Data/Results/Parameter_Calibration/Calibration_Results_v1.rds")

# Inclusão do nome dos datasets
ds_names <- data.frame(datasets = c(1:7), names = c("2-FACE","BROKEN-RING","BUPA","FORESTFIRES","GAUSS9","UNIFORM700","VOWEL2"))
params <- params %>% left_join(ds_names, on='datasets') %>% select(c('names','p','iter','grid.type'))
params$grid.type <- ifelse(params$grid.type == 'esg', 'ESGBRKGA', 'ASGBRKGA')
params <- rename(params, grid.algo = grid.type)

# Resultados
params$mean_ISM <- sapply(cal_results, function(i) mean(unlist(i['score',]))) %>% unlist()
params$mean_K <- sapply(cal_results, function(i) mean(unlist(i['k',]))) %>% unlist()
params$mean_T <- sapply(cal_results, function(i) mean(unlist(i['exec',]))) %>% unlist()

# Análise - Melhor combinação de parâmetros por algoritmo de grade
cal_results_ESGBRKGA <- params %>%
    group_by(grid.algo, p, iter) %>% 
    summarise(Mean_ISM = mean(mean_ISM), Mean_T = mean(mean_T)) %>% 
    arrange(desc(Mean_ISM)) %>% 
    filter(grid.algo == 'ESGBRKGA')

cal_results_ASGBRKGA <- params %>%
    group_by(grid.algo, p, iter) %>% 
    summarise(Mean_ISM = mean(mean_ISM), Mean_T = mean(mean_T)) %>% 
    arrange(desc(Mean_ISM)) %>% 
    filter(grid.algo == 'ASGBRKGA')

cal_results_ESGBRKGA; cal_results_ASGBRKGA

# Comentários
# . Em geral, ainda persistem valores bem baixos de silhueta (independente da abordagem de grade)
# . Contudo, novamente valores superiores de silhueta para a abordagem assimétrica, porém não tão significantes levando-se em consideração o tempo de execução
# . ASGBRKGA mais sensível a escolha dos parâmetros 'iter' e 'p' (apesar de variabilidade pequena). Melhores valores: (p = 200, iter = 300)
# . ESGBRKGA praticamente invariante a escolha de parâmetros, de forma que a escolha dos parâmetros 'menos custosos' é justificável (p = 100, iter = 200)

# Conclusão:
# . Parâmetros p/ ESGBRKGA: p = 100 | iter = 200
# . Parâmetros p/ ASGBRKGA: p = 200 | iter = 300

# TODO - Gráficos comparativos

##############################################
# 5.2.1 - Análise de Estabilidade do algoritmo
##############################################
# --> 10 replicações x 7 bases x 2 composições de grade)
# --> Objetivo: Avaliar estabilidade em termos das diferentes composições de grade

# sim_ESG <- lapply(test_subset_data, function(i) {
#     replicate(10, STATGDBC(i, alpha=.05, only.ics=0, grid.type="esg", density.test="clarkevans", clust.fobj="silhouette"))
# }); saveRDS(sim_ESG, '00_Data/Results/Stability_Assessment/ESG/sim_ESG_v1.rds') # Approx 15 mints

# sim_ASG <- lapply(test_subset_data, function(i) {
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

# clust_ESG_results <- lapply(c(all_cluster_data, test_subset_data), function(i)
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

# clust_ASG_results <- lapply(c(all_cluster_data, test_subset_data), function(i)
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

# classf_ESG_results <- lapply(all_classf_data, function(i)
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

# classf_ASG_results <- lapply(all_classf_data, function(i)
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


