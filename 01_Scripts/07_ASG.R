# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: ASG's main script
# Encoding: UTF-8

# Imports
source("01_Scripts/03_AxisInfo.R")
source("01_Scripts/04_CheckChiSqr.R")
source("01_Scripts/05_ComputeICS.R")

# Pacotes & Configuraçãoes
require("dplyr")
require("data.table")
require("spatstat.core")
require("spatstat.geom")
options(scipen = 999, digits = 10)

# Rotina principal que implementa a abordagem proposta de ASG como um todo
ASG_main <- function(mds.proj, alpha=.05, only.ics=0, iter=300, tol=0.2, p=200, pe=0.25, pm=0.3, rho=0.8, verbose=F) {
    
    # --> Entradas
    # . mds.proj <matrix>: Projeção resultante da rotina MDSProjection()
    # . alpha <float>: Nível de significância para o Teste dos Quadrats. Default é 5%.
    # . only.ics <int>: Se '1' (TRUE), avalia os grids somente com base no ICS. Se '0' (FALSE), com base no ICS e Teste dos Quadrats.
    # . iter <int>: Número de iterações do algoritmo. Default é 300 (definido com base em estudo de calibração de parâmetros)
    # . tol <int>: Proporção de gerações sem melhoria na função objetivo permitido. Default 20%.
    # . p <int>: Tamanho da poplução do algoritmo. Default é 200 (definido com base em estudo de calibração de parâmetros)
    # . pe <int>: Proporção de indivíduos que irão compor o conjunto elite. Default é 0.25.
    # . pm <int>: Proporção de indivíduos mutantes. Default é 0.3.
    # . rho <float>: Probabilidade de cruzamento dos indivíduos. Defaul é 0.8.
    # . verbose <bool>: Se a rotina deve imprimir os resultados da função objetivo ao longo das iteracoes. Default é FALSE.
    
    # --> Saídas:
    # Objeto tipo lista, contendo:
    # . fit.best: Melhor valor encontrado da função objetivo em cada iteração
    # . best.indv: Melhores indivíduos de cada iteração
    # . fit.scores: Scores dos indivíduos (distintos) do conjunto elite
    # . all.indv: Indivíduos (distintos) do conjunto elite
    # . ppp.obj: Objeto ppp associado ao grid
    # . iter: Número de iterações necessárias
    
    ####################
    # Rotinas Auxiliares
    ####################
    
    # Define a criação dos indivíduos da população
    population <- function(mds.proj, x.dim.info, y.dim.info, p) {
        
        # Informações de interesse sobre cada dimensão
        # . Menor coordenada, Maior coordenada e Tamanho mínimo dos intervalos
        min_x <- x.dim.info[1]; max_x <- x.dim.info[2]; min_delta_x <- x.dim.info[4]
        min_y <- y.dim.info[1]; max_y <- y.dim.info[2]; min_delta_y <- y.dim.info[4]
        
        # Enumera todas as escolhas de intervalos possíveis para cada eixo (considerando tamanho mínimo)
        all_intervals_in_x <- seq(min_x + min_delta_x, max_x - min_delta_x, min_delta_x)
        all_intervals_in_y <- seq(min_y + min_delta_y, max_y - min_delta_y, min_delta_y)
        
        # Determina o total de intervalos possíveis (considerando a restrição [2;N_X],[2;N_Y])
        n <- nrow(mds.proj); size_x <- length(all_intervals_in_x); size_y <- length(all_intervals_in_y)
        N_X <- ceiling(min(sqrt(n),size_x)); N_Y <- ceiling(min(sqrt(n),size_y))
        
        # Sorteio dos chaves aleatórias e composição da população
        # . Amostra sem repetição os intervalos em cada eixo
        pop <- lapply(1:p, function(i) list(
            c(min_x, sort( sample(all_intervals_in_x, sample(1:N_X,1), replace=F) ), max_x),
            c(min_y, sort( sample(all_intervals_in_y, sample(1:N_Y,1), replace=F) ), max_y)
        ))
        
        # Retorna a população
        return(pop)
        
    }
    
    # Implementa o decodificador
    decoder <- function(individual, grid.ppp, alpha, only.ics) {
        
        # Caso o tipo de avaliação desejada seja somente através do ICS
        if( only.ics == 1 ) { return(ComputeICS(grid.ppp, individual, T)) }
        
        # Caso contrário, considera tanto o Teste dos Quadrats e ICS como forma de avaliação
        else {
            
            # Verifica a conformidade em relação a restrição da estatística de teste (Chi-Quadrado) utilizada no Teste Dos Quadrats
            chisqr_restriction <- CheckChiSqr(grid.ppp, individual, T)

            if(is.na(chisqr_restriction)) {

                # Caso o grid não atenda a restrição chi-quadrado, retorna NA
                return(NA)
                
            } else {
                
                # Porém, caso tenha atendido, aplica o teste
                tess.obj <- as.tess(quadrats(grid.ppp, xbreaks=individual[[1]], ybreaks=individual[[2]]))
                quadrat_test <- suppressWarnings(quadrat.test(X=grid.ppp, tess = tess.obj))
                
                # Mas caso não tenhamos evidências estatísticas significantes de agrupamento, retorna seu ICS como zero
                if ( quadrat_test$p.value > alpha ) {
                    return(0)
                }
                
                # Do contrário, calculamos o seu ICS
                else { 
                    return(ComputeICS(grid.ppp, individual, T)) 
                    
                }
            }
            
        }
        
    }
    
    # Define o processo de cruzamento
    crossover <- function(pop.elite, pop.mutant, elite.size, mutant.size, crossover.size, rho) {
        
        # Pré-aloca a lista que guardará os indivíduos resultantes do cruzamento
        pop.crossover <- vector(mode = "list", length = crossover.size)
        
        for(i in 1:crossover.size) {
            
            # Sorteia o pai da população elite
            parent.elite <- pop.elite[[sample(1:elite.size,1)]]
            
            # Sorteia o pai da população mutante
            parent.mutant <- pop.mutant[[sample(1:mutant.size,1)]]
            
            # Probabilidade de herança de cada alelo (isto é, dos intervalos ao longo de cada eixo)
            alele_probs <- c( runif(1), runif(1) )
            
            # Constrói o filho
            pop.crossover[[i]] <- c(
                ifelse(alele_probs[1] <= rho, parent.elite[1], parent.mutant[1]),
                ifelse(alele_probs[2] <= rho, parent.elite[2], parent.mutant[2]) 
            )
        }
        
        # Retorna os resultados
        return(pop.crossover)
        
    }
    
    # ASGBRKGA = Asymetric Spaced Grid Biased Random Key Genetic Algorithm
    ASGBRKGA <- function(mds.proj, x.dim.info, y.dim.info, grid.ppp, alpha, only.ics, iter, tol, p, pe, pm, rho, verbose) {
        
        # Calcula o tamanho das populações elite, mutante e crossover
        elite.size <- ceiling(pe*p)
        mutant.size <- ceiling(pm*p)
        crossover.size <- (p - elite.size - mutant.size)
        
        # Pré-aloca os objetos que guardarão as informações das iterações
        elite <- vector(mode = "list", length = iter)
        fit.best <- vector(length = iter)
        best.indv <- vector(mode = "list", length = iter)
        
        # População inicial
        # . Verificacao: lapply(pop, function(i) c(length(i[[1]]), length(i[[2]])))
        pop <- population(mds.proj, x.dim.info, y.dim.info, p)
        
        # Contadores
        current.tol <- 0
        gen <- 1
        flag_forced_ics_execution <- FALSE
        
        # Enquanto o critério de parada não for atingido
        while(gen <= iter && current.tol <= (tol*iter)) {
            
            # Aplica o decodificador
            fit.scores <- lapply(pop, function(i) decoder(i, grid.ppp, alpha, only.ics)) %>% unlist()
            
            # Safeguard
            # . Apesar de ser raro, observou-se através da instância OUTLIERS, a possibilidade de
            # nenhuma composição de grade atenda as restrições do teste chi-quadrado. Assim, por
            # questões de compatibilidade com os demais procedimentos do algoritmo, nesta situação
            # específica, será emitido um warning ao usuário e o ICS associado a grade será calculado
            if( all(is.na(fit.scores)) ) {
                flag_forced_ics_execution = TRUE
                fit.scores <- lapply(pop, function(i) ComputeICS(grid.ppp, i, asymetrical = T)) %>% unlist() # Força o cálculo do ICS
            }            
            
            # Ordena os indivíduos com base no seu score de ICS
            pop <- pop[order(fit.scores, decreasing = T)]
            
            # Reordena os scores de ICS
            fit.scores <- fit.scores[order(fit.scores, decreasing = T)]
            
            # Atualiza e imprime o melhor valor da função objetivo até o momento
            if(gen == 1) {
                
                # Exclusivamente para a primeira iteração, irá imprimir fit.best e 
                # popular as listas de melhores indivíduos e conjunto elite
                
                fit.best[gen] <- fit.scores[1]
                best.indv[[gen]] <- pop[[1]]
                elite[[gen]] <- pop[1:elite.size]
                if (verbose == T) message(paste0("Fbest: ",fit.scores[1]," | Iteracao: ",gen,"\n"))
                
            } else {
                
                # Caso o valor de fit.best da iteração anterior seja superior
                # ao da iteração atual, não atualiza fit.best e nem o propaga
                # o mesmo conjunto elite da iteração passada para a próxima
                # iteração.
                
                fit.best[gen] <- fit.scores[1] # funciona pois 'fit.scores' é ordenado de forma decrescente
                best.indv[[gen]] <- pop[[1]]
                test <- fit.best[gen-1] >= fit.best[gen] # Problema de maximização
                
                if(isTRUE(test) || is.na(test)) { # is.na(test) existe pois é possível que o ICS seja NA (raro, porém possível)
                    
                    elite[[gen]] <- elite[[gen-1]]
                    if (verbose == T) message(paste0("Sem melhoria nesta geracao (tol = ",current.tol,")","\n"))
                    current.tol <- current.tol + 1 # ou seja, não houve melhoria nesta geração
                    
                } else {
                    
                    elite[[gen]] <- pop[1:elite.size]
                    if (verbose == T) message(paste0("Fbest: ",fit.scores[1]," | Ieracao: ",gen,"\n"))
                    current.tol <- 0 # ou seja, houve melhoria nesta geração. Logo, reseta o contador
                
                }
            }
            
            # Conjuntos elite, mutante e crossover de indivíduos
            indv.mutants <- population(mds.proj, x.dim.info, y.dim.info, mutant.size)
            indv.crossover <- crossover(elite[[gen]], indv.mutants, elite.size, mutant.size, crossover.size, rho)
            
            # Define a população da próxima iteração
            pop <- c(elite[[gen]], indv.crossover, indv.mutants)
            
            # Incrementa os contadores
            gen <- gen + 1
            
        }
        
        # Remove entradas nulas da lista de indivíduos (caso existam)
        best.indv[sapply(best.indv, is.null)] <- NULL
        
        # Safeguard - Warning
        if(flag_forced_ics_execution == T)
            message("WARNING: Nenhuma composicao de grade avaliada atendeu a restricao Chi-Quadrado. Executado considerando only.ics = 1...")
        
        # Retorna os resultados
        return(list(
            "fit.best" = max(unique(fit.best[fit.best != 0])),
            "best.indv" = unique(best.indv[fit.best != 0])[which.max(unique(fit.best[fit.best != 0]))],
            "fit.scores" = unique(fit.best[fit.best != 0]),
            "all.indv" = unique(best.indv[fit.best != 0]),
            "ppp.obj" = grid.ppp,
            "grid.method" = "ASGBRKGA",
            "iter" = gen
        ))
        
    }
    
    ##############
    # Main Program
    ##############
    
    # Calcula apriori as informações sobre as dimensões do grid
    x.dim.info <- AxisInfo(mds.proj, ax=1)
    y.dim.info <- AxisInfo(mds.proj, ax=2)
    
    # Objetos auxiliares para a rotina spatstat.core::quadrat.test()
    grid.dims <- list("xrange" = c(x.dim.info[1],x.dim.info[2]), "yrange" = c(y.dim.info[1],y.dim.info[2]))
    grid.ppp <- as.ppp(mds.proj, owin(xrange=grid.dims$xrange, yrange=grid.dims$yrange))
    
    # Aplica o algoritmo
    # if(verbose == T) message('Abordagem de Grade: ASGBRKGA')
    message('Abordagem de Grade: ASGBRKGA')
    out <- ASGBRKGA(mds.proj, x.dim.info, y.dim.info, grid.ppp, alpha, only.ics, iter, tol, p, pe, pm, rho, verbose)
    
    # Retorna as informações
    return(out)
    
}

############
# Validações
############
# asg_outputs <- lapply(mds_projections, function(i) ASG_main(i, verbose=T)) # Approx. 2min
# lapply(asg_outputs, function(i) PlotGrid(i$ppp.obj, i$best.indv[[1]], T))

# --> Investigação pontual (instância OUTLIERS)
# mds.proj <- mds_projections; alpha=.05; only.ics=0; iter=300; tol=0.2; p=200; pe=0.25; pm=0.3; rho=0.8; verbose=T
# asg_output <- ASG_main(mds_projections, verbose=T)
# PlotGrid(asg_output$ppp.obj, asg_output$best.indv[[1]], T)
