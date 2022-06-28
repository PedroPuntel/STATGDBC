# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Description: Script responsável pelos procedimentos de Pré-Processamento dos conjuntos de dados
# Encoding: UTF-8

# Rotina auxiliar que implementa os procedimentos de Padronização (Z-Score) e Escalonamento Multidimensional
MDSProjection <- function(data) {
    
    # --> Entradas:
    # . data <data.table, data.frame, tibble ...>: Conjunto de dados convertível para matrix.
    
    # --> Saídas:
    # . data <matrix>: Projeção de dimensção N x 2 resultante do Escalonamento Multidimensional

    # Padronização Z-score (se desejado)
    data <- apply(X=data, MARGIN=2, function(j) { (j-mean(j))/sd(j) })
    
    # Matriz de distâncias
    dist_mat <- as.matrix(dist(data, method = "euclidean"))
    
    # Escalonamento Multidimensional
    data <- cmdscale(d=dist_mat, k=2)
    
    # Retorna as resulados
    return(data)

}

