# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Encoding: UTF-8

# Pacotes
require("spatstat.geom")

# Rotina auxiliar que verifica a adequação de um grid às restriçõs do Teste Chi-Quadrado
CheckChiSqr <- function(ppp.object, this.grid, asymetrical = F) {
    
    if (asymetrical == F) {
        
        this.grid <- unlist(this.grid)
        nx <- as.numeric(this.grid[1])
        ny <- as.numeric(this.grid[2])
        
        # Identifica o total de células da grade
        m <- prod(unlist(this.grid))
        
        # Pontos por célula
        pp_grid <- as.data.frame(quadratcount(ppp.object, nx=nx, ny=ny))$Freq
        
        # Média de pontos por célula
        Xbarra <- sum(pp_grid)/m
        
        # Retorna o resultado
        if(m > 6 && Xbarra > 1) {
            return(0)
        } else {
            return(NA)
        }        
        
    } else {
        
        xbreaks <- as.numeric(this.grid[[1]])
        ybreaks <- as.numeric(this.grid[[2]])
        
        # Identifica o total de células da grade
        m <- (length(xbreaks)-1) * (length(ybreaks)-1)
        
        # Pontos por célula
        pp_grid <- as.data.frame(quadratcount(ppp.object, xbreaks=xbreaks, ybreaks=ybreaks))$Freq
        
        # Média de pontos por célula
        Xbarra <- sum(pp_grid)/m
        
        # Retorna o resultado
        if(m > 6 && Xbarra > 1) {
            return(0)
        } else {
            return(NA)
        }  
        
    }

}