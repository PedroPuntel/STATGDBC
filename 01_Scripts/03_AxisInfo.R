# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Encoding: UTF-8

# Rotina auxiliar calcula as estatísticas necessárias sobre as projeções obtidas vias EMD
AxisInfo <- function(mds.proj, ax) {
    
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
    min_ax_delta <- unique(sort(mds.proj[,ax], decreasing = F))[2] - min_ax_coord
    
    # Retorna as informações
    return(c(min_ax_coord, max_ax_coord, max_ax_delta, min_ax_delta))
    
}