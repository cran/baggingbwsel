".sm.Options" <-
    list(hmult = 1, h.weights = NA, period = NA,
         add = FALSE, band = NA, props = c(75, 50, 25), nbins = NA,
         positive = FALSE, delta = NA, display = NA, 
         hscale = 1, vscale = 1,
         xlab = NA, ylab = NA, zlab = NA, 
         xlim = NA, ylim = NA, zlim = NA, yht = NA,
         model = "none", reference = "none",
         panel = FALSE, panel.plot = NA,
         ngrid = NA, eval.points = NA, rugplot = TRUE, 
         col = NA, col.band = "cyan", col.mesh = "black", 
         col.points = "black",
         col.palette = topo.colors(12), col.palette.fn = topo.colors, 
         superimpose = NA,
         se = NA, se.breaks = c(-3, -2, 2, 3), lty = 1, lwd = 1, pch = 1, cex = NA,
         theta = -30, phi = 40, size = 2, scaling = NULL, 
         alpha = 0.7, alpha.mesh = 1, lit = FALSE,
         poly.index = 1, diff.ord = 2, test = NA, hull = TRUE, verbose = 1, 
         df = NA, method = NA, structure.2d = "scaled", nboot = 100, 
         describe = TRUE, show.script = TRUE, eval.grid = TRUE,
         mask.method = "hull", partial.residuals = TRUE, nlevels = 20)

.onAttach <- function(library, pkg)
{
    pkg.info <- drop( read.dcf(file = system.file("DESCRIPTION", package = "baggingbwsel"),
                               fields = c("Title", "Version", "Date") ))
    packageStartupMessage( 
      # paste0(" Package baggingbwsel: ", pkg.info["Title"], ",\n"),
      "\n Package baggingbwsel: Bagging Bandwidth Selection\n",
      " in Kernel Density and Regression Estimation,\n",
      paste0(" version ", pkg.info["Version"], " (built on ", pkg.info["Date"], ").\n"))
}


isMatrix <- function(x) length(dim(x)) == 2

isInteger <- function(x) all(x == round(x))


