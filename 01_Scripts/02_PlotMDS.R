# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Encoding: UTF-8

# Pacotes
require('data.table')
require("ggplot2")

# Rotina auxiliar que plota umna projeção obtida via EMD
PlotMDS <- function(mds.proj, cluster.vec=NULL, title="", col="#7e1e9c") {
    
    # --> Entradas:
    # . mds.proj <matrix> : Conjunto de dados resultante da projeção obtida via EMD
    # . cluster.vec <matrix>  : Vetor de clusters
    
    # --> Saídas:
    # . Gráfico resultante
    
    if(is.null(cluster.vec)) {
        
        # Gráfico
        mds.plot <- ggplot(as.data.frame(mds.proj), mapping=aes(x=V1, y=V2)) +
            geom_point(color=col, size=3) +
            labs(x="\nX\n",  y="\nY\n", title=title) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5))
        
    } else {
        
        df <- as.data.frame(cbind(mds.proj,cluster.vec))
        colnames(df) <- c('X','Y','Cluster')
        
        if( -1 %in% unique(df$Cluster) ) {
            K <- uniqueN(df$Cluster)-1
            clust_labels <- data.frame(Cluster = sort(unique(df$Cluster), F), Cluster_Fixed_Labels = c(-1,1:K))
            df <- df %>% left_join(clust_labels, on = 'Cluster') %>% select(c('X','Y','Cluster_Fixed_Labels'))
            colnames(df) <- c('X','Y','Cluster')
        } else {
            K <- uniqueN(df$Cluster)
            clust_labels <- data.frame(Cluster = sort(unique(df$Cluster), F), Cluster_Fixed_Labels = 1:K)
            df <- df %>% left_join(clust_labels, on = 'Cluster') %>% select(c('X','Y','Cluster_Fixed_Labels'))
            colnames(df) <- c('X','Y','Cluster')            
        }
        
        df$Cluster <- factor(df$Cluster, levels = sort(unique(df$Cluster)))
        #f <- df %>% arrange(Cluster)
        #df$Cluster <- as.character(df$Cluster)
        #legend_order <- sort(unique((df$Cluster)))
        
        # Gráfico
        mds.plot <- ggplot(df, mapping=aes(x=X, y=Y, colour=Cluster)) +
            geom_point(size=3) +
            labs(x="\nX\n",  y="\nY\n", title=title) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5))  
        
    }
    
    # Retorna
    return(mds.plot)
}
