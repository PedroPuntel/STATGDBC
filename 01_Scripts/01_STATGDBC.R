# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: STATGDBC main script
# Encoding: UTF-8

# Notas de Desenvolvimento
# --> TODO: Implementar lógica que permita ao usuário fornecer valores aos parâmetros do BRKGA

# Imports
source("01_Scripts/06_ESG.R")
source("01_Scripts/07_ASG.R")
source("01_Scripts/08_DGBClust.R")

# Rotina que implementa o algoritmo STATGDBC por completo
STATGDBC <- function(mds.proj, alpha=.05, only.ics=0, grid.type="esg", density.test="hopkins", clust.fobj="silhouette") {
    
    # --> Entradas:
    # . mds.proj <matrix>: Projeção bidimensional N x 2 resultante do Escalonamento Multidimensional
    # . alpha <float>: Nível de significância considerado no testes estatísticos. Default é 5%.
    # . only.ics <int>: Se '1' (TRUE), avalia os grids somente com base no ICS. Se '0' (FALSE), com utiliza além do ICS o Teste dos Quadrats.
    # . grid.type <str>: Tipo de composição de grade a ser considerada, sendo um dentre 'esg' ou 'asg'. Default é 'asg'.
    # . density.test <str>: Teste estatístico a ser utilizado na etapa de densidade. Um dentre 'hopkins' e 'clarkevans'. Default é 'clarkevans'.
    # . clust.fobj <str>: Função objetivo para avaliação dos clusters na etapa de densidade. Um dentre 'silhouette' e 'calinski'. Default é 'clarkevans'.

    # --> Saídas (BFESGA)
    # . cluster <vector> : 
    # . score <float> : 
    # . k <int> : 
    # . dist <table> :
    # . grid <vector> :
    # . grid.ics <float> :
    # . spatial.ppp.obj <spatial.ppp> : 
    # . exec <float> : 
    
    # --> Saída (ESGBRKGA/ASGBRKGA)
    # . cluster <vector> : 
    # . score <float> : 
    # . k <int> : 
    # . dist <table> :
    # . best.grid <list>: 
    # . best.grid.ics <float> :
    # . grids <list> :
    # . spatial.ppp.obj <spatial.ppp> : 
    # . exec <float> : 
    
    # Inicia o cronômetro
    tic <- Sys.time()
    
        #########################
        # Fase 1 - Etapa de Grade
        #########################
        
        if( grid.type ==  "esg" ) {
            
            # Abordagem de grid simétrico (regular)
            grid.results <- ESG_main(
                mds.proj = mds.proj,
                alpha = alpha,
                only.ics = only.ics
            )
            
        } else {
            
            # Abordagem de grid assimétrico (irregular)
            grid.results <- ASG_main(
                mds.proj = mds.proj,
                alpha = alpha,
                only.ics = only.ics
            )
            
        }
    
        #############################
        # Fase 2 - Etapa de Densidade
        #############################
        
        density.results <- DGBClust_main(
            ppp.obj = grid.results$ppp.obj,
            elite.grids = grid.results$best.indv,
            elite.grids.scores = grid.results$fit.best,
            which.method = grid.results$grid.method,
            alpha = alpha,
            density.test = density.test,
            clust.fobj = clust.fobj
        )

    ########################
    # Retorna as informações
    ########################
    
    # Encerra o cronômetro e calcula o tempo decorrido
    toc <- Sys.time()
    elapsed <- round(as.numeric(difftime(time1 = toc, time2 = tic, units = "secs")),2)
    density.results$exec = elapsed
    
    # Retorna os objetos de saída
    return(density.results)
    
}


