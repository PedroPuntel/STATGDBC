# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: DGBClust's main script
# Encoding: UTF-8

# Imports
source("01_Scripts/03_AxisInfo.R")

# Pacotes & Configurações
require("fpc")
require("dplyr")    
require("hopkins")
require("cluster")
require("data.table")
require("spatstat.core")
require("spatstat.geom")
options(scipen = 999)

# Rotina principal que implementa a abordagem proposta do DGBClust como um todo
DGBClust_main <- function(data, ppp.obj, elite.grids, elite.grids.scores, which.method, alpha=.05, density.test='clarkevans', clust.fobj="silhouette") {
    
    # --> Entradas
    # . data <data.table | data.frame>: Conjunto de dados
    # . ppp.obj <ppp>: Objeto ppp (retornado pelos algoritmos na etapa de grade) associado a instância avaliada
    # . elite.grids <list>: Conjunto elite de indivíduos propagados da etapa de grade
    # . elite.grids.scores <list>: Scores associados ao indivíduos do conjunto elite
    # . which.method <bool>: Identifica de qual algoritmo 'BFESGA', 'ESGBRKGA' ou 'ASGBRKGA' a solução pertence.
    # . alpha <float>: Nível de significância considerado no testes estatísticos. Default é 5%.
    # . density.test <str>: Teste estatístico a ser utilizado na etapa de densidade. Um dentre 'hopkins' e 'clarkevans'. Default é 'clarkevans'.
    # . clust.fobj <str>: Função objetivo para avaliação dos clusters na etapa de densidade. Um dentre 'silhouette' e 'calinski'. Default é 'silhouette'.
    
    ####################
    # Rotinas Auxiliares
    ####################
    
    # Rotina auxiliar que identifica a qual célula da grade cada ponto de dados pertence
    get_cell_mapping <- function(ppp.obj, which.method, this.grid) {
        
        # this.grid <- elite.grids[[1]]
        
        # Obtêm as coordenadas resultantes do EMD
        mds_coords <- as.matrix(cbind(ppp.obj$x, ppp.obj$y))
        colnames(mds_coords) <- c("Xcoord","Ycoord")
        
        # Caso o algoritmo de grade aplicado seja aquele da grade assimétrica
        if(which.method == 'ASGBRKGA') {
            
            # Extraí diretamente os intervalos onde estão definidas as células
            xbreaks <- this.grid[1] %>% unlist() %>% as.numeric()
            ybreaks <- this.grid[2] %>% unlist() %>% as.numeric()
            
        } else {
            
            # Do contrário, identifica os intervalos onde estão definidas as células
            nx <- this.grid[1] %>% as.numeric()
            ny <- this.grid[2] %>% as.numeric()
            max_x_coord <- max(mds_coords[,1])
            min_x_coord <- min(mds_coords[,1])
            max_y_coord <- max(mds_coords[,2])
            min_y_coord <- min(mds_coords[,2])
            delta_x <- (max_x_coord - min_x_coord)/nx
            delta_y <- (max_y_coord - min_y_coord)/ny
            xbreaks <- seq(min_x_coord, max_x_coord, by=delta_x)
            ybreaks <- seq(min_y_coord, max_y_coord, by=delta_y)
            
        }
        
        # Adiciona uma pequena porção de ruído aleatório as 'pontas' dos grids para garantir
        # que todos os prontos serão compreendidos no intervalo [Xmin,Xmax] e [Ymin,Ymax]. É
        # preciso testar se as pontas são positivas ou negativas pois, a depender do caso,
        # estaríamos 'encurtando' a grade ao invés de 'expandi-la', como queremos.
        xbreaks[1] <- ifelse(xbreaks[1] >= 0, xbreaks[1] + runif(1,0,1e-3), xbreaks[1] - runif(1,0,1e-3))
        xbreaks[length(xbreaks)] <- ifelse(xbreaks[length(xbreaks)] >= 0, xbreaks[length(xbreaks)] + runif(1,0,1e-3), xbreaks[length(xbreaks)] - runif(1,0,1e-3))
        ybreaks[1] <- ifelse(ybreaks[1] >= 0, ybreaks[1] + runif(1,0,1e-3), ybreaks[1] - runif(1,0,1e-3))
        ybreaks[length(ybreaks)] <- ifelse(ybreaks[length(ybreaks)] >= 0, ybreaks[length(ybreaks)] + runif(1,0,1e-3), ybreaks[length(ybreaks)] - runif(1,0,1e-3))
        
        # Mapeia os pontos as células (excluindo células vazias)
        grid_mapped <- apply(mds_coords, 1, function(i) {
            c(findInterval(i[1], xbreaks, rightmost.closed=F, left.open=F),findInterval(i[2], ybreaks, rightmost.closed=F, left.open=F))
        }) %>% t() %>% as.data.frame() 
        # A depender do algoritmo de grade aplicado, obtêm as contagens de pontos por celula manipulando o objeto .ppp associado
        if(which.method == 'ASGBRKGA') {
            grid_qcount <- quadratcount(ppp.obj, xbreaks=this.grid[[1]], ybreaks=this.grid[[2]])
        } else {
            grid_qcount <- quadratcount(ppp.obj, nx=nx, ny=ny)                
        }
        
        # Adiciona os índices as células
        grid_labels <- cbind(as.vector(apply(grid_qcount,2,rev)), 1:(nrow(grid_qcount)*ncol(grid_qcount))) %>% as.data.frame()
        colnames(grid_labels) <- c('Count','Index')
        grid_mapped_grouped <- grid_mapped %>% group_by(V1,V2) %>% summarise(Count=n())
        grid_mapped_grouped$Index <- 1:nrow(grid_mapped_grouped)
        grid_mapped_grouped <- grid_mapped_grouped %>% left_join(grid_labels, on=c('Count','Index'))
        grid_mapped <- left_join(grid_mapped, grid_mapped_grouped, on=c('V1','V2'))
        colnames(grid_mapped) <- c("Xcell","Ycell","Count","Index")
        
        # Verifica se para esta composição de grade, existem células vazias
        if( 0 %in% grid_qcount ) {
            
            # Instancia 'grid_mapped' com zeros inicialmente
            grid_index_mat <- matrix(data=0, nrow=nrow(grid_qcount), ncol=ncol(grid_qcount))
            
            # Variável auxiliar
            counter <- 1
            
            # Indexa as células não-nulas
            for(col in 1:ncol(grid_qcount)) {
                for(row in nrow(grid_qcount):1) {
                    if(grid_qcount[row,col] != 0) {
                        grid_index_mat[row,col] <- counter
                        counter <- counter + 1
                    } 
                }
            }
            
            # Indexa as células nulas
            grid_index_mat[grid_index_mat == 0] <- setdiff(grid_labels$Index, c(1:(counter-1)))
            empty_cells <- setdiff(grid_labels$Index, c(1:(counter-1)))
            
        } else {
            
            # OBS: Numeração crescente de baixo para cima, da esquerda para a direita. Para obter
            # tal numeração, precisamos inverter a ordem, por coluna, dos índices
            grid_index_mat <- matrix(1:(nrow(grid_qcount)*ncol(grid_qcount)), nrow=length(ybreaks)-1, ncol=length(xbreaks)-1) %>% apply(2, rev)
            empty_cells <- NULL
        }
        
        # Retorna as informações
        return(list(
            "points_cell_mapping" = as.data.frame(cbind(mds_coords,grid_mapped[,c("Xcell","Ycell","Index")])),
            "cells_index" = grid_index_mat,
            "which_empty" = empty_cells
        ))    
        
    }
    
    # Rotina auxiliar que retorna as células vizinhas a uma célula especifica
    get_neighbor_cells <- function(grid_index_mat, empty_cells) {
        
        n <- nrow(grid_index_mat)
        k <- ncol(grid_index_mat)
        cell_positions <- cbind(expand.grid(1:n,1:k), as.vector(apply(grid_index_mat, 2, function(i) i))) %>% as.matrix()
        
        get_adjacent_cells <- function(grid_index_mat, n, k, empty_cells, this.cell) {
            
            if( grid_index_mat[this.cell[1],this.cell[2]] %in% empty_cells ) {
                
                return(NULL)
                
            } else {
                
                neighbors <- NULL
                
                for(row in c(-1,0,1)) {
                    for(col in c(-1,0,1)) {
                        row_check <- ((this.cell[1]+row) >= 1) & ((this.cell[1]+row) <= n)
                        col_check <- ((this.cell[2]+col) >= 1) & ((this.cell[2]+col) <= k)
                        if( isTRUE(row_check) && isTRUE(col_check) ) {
                            neighbors <- c(neighbors, grid_index_mat[this.cell[1]+row, this.cell[2]+col])
                        }
                    }
                }
                
                return(sort(neighbors, decreasing=F))                
            }
        }
        
        out <- apply(cell_positions, 1, function(i) get_adjacent_cells(grid_index_mat, n, k, empty_cells, i))
        names(out) <- cell_positions[,3]
        
        return(out)
        
    }    
    
    # Rotina auxiliar que implementa os testes estatísticos de Hopkins/Clark-Evans
    assess_cluster_tendency <- function(mapped_matrix, cells_to_merge, density.test, alpha) {
        
        # Filtra somente os pontos compreendidos pelas células em 'cells_to_merge'
        aux <- mapped_matrix[mapped_matrix$Index %in% cells_to_merge, c('Xcoord','Ycoord')]
        
        # Aplica o teste estatístico conforme escolha do usuário
        if(density.test == 'hopkins') {
            
            # Caso a célula possua menos de 5 pontos, não avalia esta célula (não dá para ter uma ideia da variabilidade)
            if(nrow(aux) <= 5) {
                
                flag_merge <- -1
                
            } else {
                
                # Do contrário, aplica o teste de Hopkins | ?hopkins::hopkins
                # H0: PP ~ CSR
                # H1: PP != CSR
                m <- nrow(aux)-1 # Utiliza os |M|-1 pontos contidos na célula para compor a amostra da estatística de teste
                hopkins_stat <- hopkins(X=aux, m=m)
                flag_merge <- ifelse(hopkins.pval(hopkins_stat, n=m) > alpha, TRUE, FALSE)
                
                # Comentário:
                # Para avaliarmos as junções entre as células, queremos juntá-las somente caso tenhamos
                # evidência estatística suficiente de que a distribuição dos objetos compreendidos pelas
                # células avaliadas seja homogênea. Alternativamente, queremos que os objetos pertencentes
                # a um mesmo cluster, sejam regularmente distribuídos.
                
            }
            
        } else {
            
            # Caso a célula possua menos de 5 pontos, não avalia esta célula (não dá para ter uma ideia da variabilidade)
            if(nrow(aux) <= 5) {
                
                flag_merge <- -1
                
            } else {
            
                # Define o objeto ppp referente às células mescladas
                x.dim.info <- AxisInfo(aux, ax=1)
                y.dim.info <- AxisInfo(aux, ax=2)
                grid.dims <- list("xrange" = c(x.dim.info[1],x.dim.info[2]), "yrange" = c(y.dim.info[1],y.dim.info[2]))
                aux_ppp_obj <- as.ppp(aux, owin(xrange=grid.dims$xrange, yrange=grid.dims$yrange))
                
                # Aplica o teste Clark-Evans | ?spatstat.core::clarkevans
                # OBS: Não utiliza correção dos efeitos de borda.
                regular_pp_test <- ifelse(clarkevans.test(aux_ppp_obj, correction="none", alternative="regular", nsim=1000)$p.val <= alpha, TRUE, FALSE)
                random_pp_test <- ifelse(clarkevans.test(aux_ppp_obj, correction="none", alternative="two.sided", nsim=1000)$p.val > alpha, TRUE, FALSE)
                
                if( isTRUE(regular_pp_test) || isTRUE(random_pp_test) ) {
                    flag_merge <- TRUE
                } else {
                    flag_merge <- FALSE
                }
                
                # Comentário:
                # Para avaliarmos as junções entre as células, queremos juntá-las somente caso tenhamos
                # evidência estatística suficiente de que a distribuição dos objetos compreendidos pelas
                # células avaliadas seja homogênea/aleatória. Alternativamente, queremos que os objetos pertencentes
                # a um mesmo cluster, sejam regularmente/aleatoriamente distribuídos.
            
            }
            
        }
        
        # Retorna a decisão
        return(flag_merge)
        
    }
    
    # Rotina auxiliar que implementa a lógica de mesclagem das células
    assess_grid_density <- function(this.grid, cell.neighborhood, density.test, alpha) {
        
        # BFESGA
        # --> this.grid <- grid_cell_mapping; cell.neighborhood <- grid_neighbor_cells
        
        # ASGBRKGA/ESGBRKGA
        # --> this.grid <- grids_cell_mapping[[3]]; cell.neighborhood <- grids_neighbor_cells[[3]]
        
        # Vetor auxiliar contendo as células não visitadas
        to_visit <- as.vector(this.grid$cells_index)
        
        # Descarta as células vazias
        to_visit <- sort(to_visit[!to_visit %in% this.grid$which_empty])
        
        # Lista auxiliar contendo as mesclagens realizadas
        all_merges <- vector(mode="list", length=length(to_visit))
        
        # Para a i-ésima célula da grade
        for(i in 1:length(to_visit)) {
            
            # Lista as suas células vizinhas (incluindo ela mesma)
            neighbors <- as.vector(unlist(cell.neighborhood[toString(to_visit[i])]))
            
            # Descarta as células vizinhas nulas (caso existam)
            neighbors <- neighbors[!neighbors %in% this.grid$which_empty]
                
            # Caso a lista de vizinhos fique vazia, ou seja, caso não existam junções a serem avaliadas
            if(is.empty(neighbors)) {
                
                # Avalia se a célula em questão define, por si só, um cluster
                # OBS: Note que neste caso, queremos testar se existem evidências de agrupamento na célula em questão
                # o que seria o contrário (negação) do resultado retornado pelos Testes Estatísticos ('p-valor é baixo?')
                eval <- suppressWarnings(assess_cluster_tendency(this.grid$points_cell_mapping, to_visit[i], density.test, alpha))
                flag_merge <- ifelse(eval == -1, -1, !eval)
                
                # Caso defina
                if(flag_merge != -1 && isTRUE(flag_merge)) {
                    
                    # Trata esta célula individualmente como um cluster
                    all_merges[[toString(i)]] <- to_visit[i]
                    
                } else {
                    
                    # Do contrário, trata esta como ruído
                    all_merges[["-1"]] <- to_visit[i]
                    
                }
                
            } else {
                    
                # Do contrário, define uma variável auxiliar para guardar a junções feitas (caso existam)
                merges <- NULL
                
                # Para a j-ésima célula vizinha
                for(j in 1:length(neighbors)) {
                    
                    # Avalia a mesclagem
                    flag_merge <- suppressWarnings(assess_cluster_tendency(this.grid$points_cell_mapping, c(to_visit[i],neighbors[j]), density.test, alpha))
                    
                    # Se for estatisticamente relevante
                    if(flag_merge) {
                        
                        # Realiza-a. Do contrário, não precisa fazer nada
                        merges <- unique(c(merges,neighbors[j]))
                        
                    }
                    
                }
                
                # Atualiza a lista com todas as mesclagens
                # OBS: Caso não tenha realizado nenhuma mesclagem, retorna o próprio índice da célula
                all_merges[[toString(to_visit[i])]] <- ifelse(is.null(merges),to_visit[i],merges)
                    
            }
                
        }
            
        # Excluí as entradas nulas da lista de mesclagens
        all_merges[sapply(all_merges, is.null)] <- NULL
        
        # Vetor de clusters (referente às grades)
        cluster_vec <- sort(unlist(all_merges))
        cluster_vec <- factor(cluster_vec, levels=unique(cluster_vec), labels=1:length(unique(cluster_vec)))
        cluster_vec <- cbind(as.numeric(names(cluster_vec)),as.numeric(cluster_vec)) %>% as.data.frame()
        colnames(cluster_vec) <- c("Index","Cluster")
        
        # Lógica de reprocessamento
        # --> Avalia novamente cada um dos clusters formados (células mesclados e/ou células individuais)
        # buscando por evidências de agrupamento (hipótese contrária aquela dos testes estatísticos). Caso
        # os clusters formados não tenham evidências de agrupamento, estes são marcados como ruído.
        cluster_vec$Clustered <- suppressWarnings(
                sapply(cluster_vec$Index, function(i) assess_cluster_tendency(this.grid$points_cell_mapping, i, density.test, alpha))
        )
        
        # OBS: Note que neste caso, queremos testar se existem evidências de agrupamento na célula em questão
        # o que seria o contrário (negação) do resultado retornado pelos Testes Estatísticos ('p-valor é baixo?')
        cluster_vec$Clustered <- ifelse(cluster_vec$Clustered == -1, -1, !cluster_vec$Clustered)
        
        # Atualiza o vetor de clusters com o resultado do reprocessamento
        cluster_vec$Cluster <- ifelse(cluster_vec$Clustered == -1, -1, cluster_vec$Cluster)
        cluster_vec$Clustered <- NULL
        
        # Retorna as informações
        return(cluster_vec)

    }
    
    # Rotina auxiliar que dado um vetor de clusters, calcula os Índices de Silhueta/Calinksi-Harabasz
    assess_cluster_quality <- function(data, this.cluster, clust.fobj) {
        
        # this.cluster <- points_cluster_mapping[[1]]

        # Remove os objetos do tipo ruído
        if(!is.integer(which(this.cluster$Cluster == -1, arr.ind = T))) {
            to_remove <- which(this.cluster$Cluster == -1, arr.ind = T) 
            data <- data[-to_remove,] # Também remove do conjunto de dados original (sem projeção EMD)
            this.cluster <- this.cluster[-to_rmeove, c("Xcoord","Ycoord","Cluster")]
            
        } else {
            this.cluster <- this.cluster[, c("Xcoord","Ycoord","Cluster")]            
        }
    
        # Padronização das variáveis do conjunto de dados
        std_dist_mat <- apply(data, 2, function(j) {(j-mean(j))/sd(j)}) %>% dist() %>% as.matrix()
        
        if(clust.fobj == 'silhouette') {
            
            # Calcula o Índice de Silhoueta Médio
            if(uniqueN(this.cluster$Cluster) == 1) {
                score <- 0
            } else {
                score <- cluster::silhouette(this.cluster$Cluster, dmatrix=std_dist_mat) %>% as.data.frame()
                score <- mean(score$sil_width)
            }
        
        } else {
            
            # Calcula o Índice de Calinksi-Harabasz
            if(uniqueN(this.cluster$Cluster) == 1) {
                score <- -Inf
            } else {
                score <- fpc::calinhara(x=this.cluster, clustering=this.cluster$Cluster, cn=uniqueN(this.cluster$Cluster))
            }
        
        }
        
        # Retorna os resultados
        return(score)
        
    }

    ##############
    # Main Program
    ##############
    
    # Mapeia os pontos as células da grade e estrutura de vizinhança, diferenciando a abordagem a depender do algoritmo de grade aplicado
    if(which.method == "BFESGA") {
        
        # Seleciona o melhor indivíduo obtido via força bruta
        elite.grids <- elite.grids[1,]
        
        # Mapeia os objetos as células da grade, indexando as mesmas de forma sistemática
        grid_cell_mapping <- suppressMessages(get_cell_mapping(ppp.obj, which.method, elite.grids))
        
        # Enumeração da estrutura de vizinhança entre as células
        grid_neighbor_cells <- suppressMessages(get_neighbor_cells(grid_cell_mapping$cells_index, grid_cell_mapping$which_empty))
        
        # Análise de densidade: junção das células da grade
        grid_cluster_mapping <- suppressMessages(assess_grid_density(grid_cell_mapping, grid_neighbor_cells, density.test, alpha))
        
        # Mapeamento dos objetos aos clusters formados
        points_cluster_mapping <- suppressMessages(left_join(grid_cell_mapping$points_cell_mapping, grid_cluster_mapping, on="Index"))
        
        # Avaliação dos clusters com base no índice de qualidade escolhido
        cluster_score <- suppressMessages(assess_cluster_quality(data, points_cluster_mapping, clust.fobj))
        
        message('..... Finalizado DGBClust')
        
        # Retorna a melhor solução de agrupamento e suas respectivas informações
        return(list(
            "cluster" = points_cluster_mapping$Cluster,
            "score" = cluster_score,
            "k" = ifelse(-1 %in% unique(points_cluster_mapping$Cluster), uniqueN(points_cluster_mapping$Cluster)-1, uniqueN(points_cluster_mapping$Cluster)),
            "dist" = table(points_cluster_mapping$Cluster),
            "grid" = elite.grids,
            "grid.ics" = elite.grids.scores[1],
            "grid.method" = which.method,
            "spatial.ppp.obj" = ppp.obj
        ))
        
    } else {
        
        # Mapeia os objetos as células da grade, indexando as mesmas de forma sistemática
        grids_cell_mapping <- suppressMessages(lapply(elite.grids, function(i) get_cell_mapping(ppp.obj, which.method, i))) 
        
        # Enumeração da estrutura de vizinhança entre as células
        grids_neighbor_cells <- suppressMessages(lapply(grids_cell_mapping, function(i) get_neighbor_cells(i$cells_index, i$which_empty)))
        
        # Análise de densidade: junção das células da grade
        grids_cluster_mapping <- suppressMessages(lapply(1:length(grids_cell_mapping), function(i) {
            assess_grid_density(grids_cell_mapping[[i]], grids_neighbor_cells[[i]], density.test, alpha)
        }))
        
        # Mapeamento dos objetos as clusters formados
        points_cluster_mapping <- suppressMessages(lapply(1:length(grids_cell_mapping), function(i) {
            left_join(grids_cell_mapping[[i]]$points_cell_mapping, grids_cluster_mapping[[i]], on="Index") 
        }))
        
        # Avaliação dos clusters com base no índice de qualidade escolhido
        cluster_scores <- sapply(points_cluster_mapping, function(i) assess_cluster_quality(data, i, clust.fobj))
        best_cluster <- which.max(cluster_scores)[1] # Independente do índice escolhido, queremos maximizar ambos. Em caso de empate, seleciona o primeiro.
        
        message('..... Finalizado DGBClust')
        
        # Retorna a melhor solução de agrupamento e suas respectivas informações
        return(list(
            "cluster" = points_cluster_mapping[[best_cluster]]$Cluster,
            "score" = cluster_scores[best_cluster],
            "k" = ifelse(-1 %in% unique(points_cluster_mapping[[best_cluster]]$Cluster),
                         uniqueN(points_cluster_mapping[[best_cluster]]$Cluster)-1,
                         uniqueN(points_cluster_mapping[[best_cluster]]$Cluster)),
            "dist" = table(points_cluster_mapping[[best_cluster]]$Cluster),
            "best.grid" = elite.grids[[best_cluster]],
            "best.grid.ics" = elite.grids.scores[best_cluster],
            "grids" = elite.grids,
            "grids.method" = which.method,
            "spatial.ppp.obj" = ppp.obj
        ))
        
    }
    
}

############
# Validações
############
# --> Debugging
# data <- datasets[[2]]; ppp.obj <- esg_outputs[[2]]$ppp.obj; elite.grids <- esg_outputs[[2]]$all.indv
# elite.grids.scores <- esg_outputs[[2]]$fit.scores; which.method <- esg_outputs[[2]]$grid.method
# alpha=.05; density.test='hopkins'; clust.fobj="silhouette"

# --> Validação
# datasets <- all_cluster_data[c(1,3,14,24)]
# 
# esg_clust <- lapply(1:length(datasets), function(i) {
#     DGBClust_main(datasets[[i]], esg_outputs[[i]]$ppp.obj, esg_outputs[[i]]$all.indv, esg_outputs[[i]]$fit.scores, esg_outputs[[i]]$grid.method)
# })
# 
# lapply(1:length(esg_clust), function(i) {
#     PlotMDS(cbind(esg_clust[[i]]$spatial.ppp.obj$x, esg_clust[[i]]$spatial.ppp.obj$y), esg_clust[[i]]$cluster,
#             title = paste0("ISM: ", esg_clust[[i]]$score, "\n", "ICS: ", esg_clust[[i]]$grid.ics) )
# })
# 
# asg_clust <- lapply(1:length(datasets), function(i) {
#     DGBClust_main(datasets[[i]], asg_outputs[[i]]$ppp.obj, asg_outputs[[i]]$all.indv, asg_outputs[[i]]$fit.scores, asg_outputs[[i]]$grid.method)
# })
# 
# lapply(1:length(asg_clust), function(i) {
#     PlotMDS(cbind(asg_clust[[i]]$spatial.ppp.obj$x, asg_clust[[i]]$spatial.ppp.obj$y), asg_clust[[i]]$cluster,
#             title = paste0("ISM: ", asg_clust[[i]]$score, "\n", "ICS: ", asg_clust[[i]]$best.grid.ics))
# })

# --> Investigação pontual (instância OUTLIERS)
# dataset <- all_cluster_data[[19]]; ppp.obj <- esg_output$ppp.obj; elite.grids <- esg_output$all.indv
# elite.grids.scores <- esg_output$fit.scores; which.method <- esg_output$grid.method; alpha=.05; density.test='clarkevans'; clust.fobj="silhouette"
#
# outliers_esg <- DGBClust_main(dataset, esg_output$ppp.obj, esg_output$all.indv, esg_output$fit.scores, esg_output$grid.method)
# PlotMDS(cbind(outliers_esg$spatial.ppp.obj$x, outliers_esg$spatial.ppp.obj$y), outliers_esg$cluster,
#     title = paste0("ISM: ", outliers_esg$score, "\n", "ICS: ", outliers_esg$grid.ics)) 
#
# outliers_asg <- DGBClust_main(dataset, asg_output$ppp.obj, asg_output$all.indv, asg_output$fit.scores, asg_output$grid.method)
# PlotMDS(cbind(outliers_asg$spatial.ppp.obj$x, outliers_asg$spatial.ppp.obj$y), outliers_asg$cluster,
#      title = paste0("ISM: ", outliers_asg$score, "\n", "ICS: ", outliers_asg$best.grid.ics)) 
