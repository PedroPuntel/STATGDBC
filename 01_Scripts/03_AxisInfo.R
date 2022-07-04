# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Encoding: UTF-8

# Rotina auxiliar calcula as estatísticas necessárias sobre as projeções obtidas vias EMD
AxisInfo <- function(mds.proj, ax, tol=0.1) {
    
    # Menor coordenada do eixo
    min_ax_coord <- min(mds.proj[,ax])
    
    # Índice da menor coordenada do eixo
    index_min_ax_coord <- which.min(mds.proj[,ax])
    
    # Maior coordenada do eixo
    max_ax_coord <- max(mds.proj[,ax])
    
    # Maior variação ao longo do eixo
    max_ax_delta <- max_ax_coord - min_ax_coord
    
    # Menor variação ao longo do eixo
    # dist_to_min_ax_coord <- as.matrix(dist(rbind(mds.proj, mds.proj[index_min_ax_coord,ax])))[nrow(mds.proj) + 1, 1:nrow(mds.proj)]
    # nearest_point_from_min_ax_coord <- dist_to_min_ax_coord[which.min(dist_to_min_ax_coord)]
    # index_of_nearest_point <- names(nearest_point_from_min_ax_coord) %>% as.numeric()
    # min_ax_delta <- unique(sort(mds.proj[,ax], decreasing = F))[2] - min_ax_coord
    min_ax_delta <- min(abs(diff(mds.proj[,ax])))
    
    # Safeguard
    # . Observou-se que para conjuntos de daods onde existem muitos pontos "coexistentes", a menor variação
    # ao longo do eixo pode ser extremamente pequena (e.g 1e-19), o que impede e enumerção de todos os pontos
    # de corte para consrução os individuos no ASGBRKGA, por exemplo. Assim, como solução paliativa, caso o
    # valor de min_ax_delta seja inferior ao parâmetro 'tol', considera-se a menor variação ao longo do eixo
    # como o 'tol'.
    if (min_ax_delta < tol) {
        min_ax_delta <- tol
    }
    
    # Retorna as informações
    return(c(min_ax_coord, max_ax_coord, max_ax_delta, min_ax_delta))
    
}