# Author: Pedro Puntel (pedro.puntel@gmail.com)
# Encoding: UTF-8

# Imports
source("01_Scripts/04_CheckChiSqr.R")
source("01_Scripts/05_ComputeICS.R")

# Pacotes
require("spatstat.core")
require("spatstat.geom")

# Rotina auxiliar para plotagem dos grids
PlotGrid <- function(ppp.object, this.grid, asymetrical = F, data_name="") {
    
    if( asymetrical == F ) {
        
        this.grid <- unlist(this.grid)
        nx <- as.numeric(this.grid[1])
        ny <- as.numeric(this.grid[2])

        # grid.title <- paste0(
        #     "\n Grid (", nx," x ", ny ,")",
        #     "\n Tipo: Regular",
        #     "\n ICS: ", round(ComputeICS(ppp.object, this.grid, F),2),
        #     "\n Restricao Chi-Sqr: ", ifelse(CheckChiSqr(ppp.object, this.grid, F) == 0, 'Sim', 'Nao'),
        #     "\n Teste dos Quadrats (p < alpha): ", ifelse(suppressWarnings(quadrat.test(X=ppp.object, nx=nx, ny=ny)$p.value) <= 0.05, "Sim", "Nao")
        # )
        
        grid.title <- paste0(
            "\n Grid (", nx," x ", ny ,")",
            "\n Tipo de Grade: Regular (ESG)",
            "\n ICS: ", round(ComputeICS(ppp.object, this.grid, F),2),
            "\n Instância: ", data_name
        )
        
        qcount <- quadratcount(ppp.object, nx=nx, ny=ny)
        plot(qcount, main=grid.title, col="red", lwd=2)
        plot(ppp.object, add=TRUE, pch=19, col="#005DB3")        
        
    } else {
        
        xbreaks <- as.numeric(this.grid[[1]])
        ybreaks <- as.numeric(this.grid[[2]])
        
        # grid.title <- paste0(
        #     "\n Grid (", length(xbreaks)-1," x ", length(ybreaks)-1,")",
        #     "\n Tipo: Irregular",
        #     "\n ICS: ", round(ComputeICS(ppp.object, this.grid, T),2),
        #     "\n Restricao Chi-Sqr: ", ifelse(CheckChiSqr(ppp.object, this.grid, T) == 0, 'Sim', 'Nao'),
        #     "\n Teste dos Quadrats (p < alpha): ", ifelse(suppressWarnings(quadrat.test(X=ppp.object, xbreaks=xbreaks, ybreaks=ybreaks)$p.value) <= 0.05, "Sim", "Nao")
        # )
        
        grid.title <- paste0(
            "\n Grid (", length(xbreaks)-1," x ", length(ybreaks)-1,")",
            "\n Tipo de Grade: Irregular (ASG)",
            "\n ICS: ", round(ComputeICS(ppp.object, this.grid, T),2),
            "\n Instância: ", data_name
        )
        
        qcount <- quadratcount(ppp.object, xbreaks=xbreaks, ybreaks=ybreaks)
        plot(qcount, main=grid.title, col="red", lwd=2)
        plot(ppp.object, add=TRUE, pch=19, col="#F8766D")          

    }
    
}