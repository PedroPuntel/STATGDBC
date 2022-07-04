# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Encoding: UTF-8

# Pacotes
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
        df$Cluster <- as.character(df$Cluster)
        
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
