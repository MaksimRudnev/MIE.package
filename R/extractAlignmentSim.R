#' Extracts summaries of alignment simulations from Mplus output file
#' 
#' @param sim.outputs a character vector of file names containing simulation results of alignment ran in Mplus. 
#' @param silent Logical. Used for debugging.
#' @details Best used as part of \code{\link[MIE]{runAlignment}}
#' 
#' @examples 
#' \dontrun{  
#'      align.sim.summ <- extractAlignmentSim (c("sim500.out", "sim100.out", "sim1000.out")) 
#'      }
#' 
#' @return Invisibly returns a summary table.
#' @seealso \code{\link[MIE]{runAlignment}}  and \code{\link[MIE]{extractAlignment}} 
#' @export
extractAlignmentSim <- function(sim.outputs = c("sim500.out", "sim100.out", "sim1000.out"), silent = FALSE) {
  
  
  otp<- lapply(sim.outputs, function(x) {
    f <- paste(readLines(x), collapse="\n")
    cor.tabs <- sub(".*CORRELATIONS AND MEAN SQUARE ERROR OF POPULATION AND ESTIMATE VALUES *(.*?) *QUALITY OF NUMERICAL RESULTS.*", "\\1", f)
    cor.tabs <- strsplit(cor.tabs, "CORRELATION AND MEAN SQUARE ERROR OF THE AVERAGE ESTIMATES")[[1]]
    cor.tabs1 <- strsplit(cor.tabs[1], "\n")[[1]][-c(1:4)]
    cor.tabs1 <- cor.tabs1[!cor.tabs1==""]
    
    dt1 <-sapply(seq(1, length(cor.tabs1), by=3), function(f.id) unlist(read.table(text=cor.tabs1[f.id+1:2])[-1]))
    colnames(dt1) = gsub(" ", "", cor.tabs1[seq(1, length(cor.tabs1), by=3)])
    row.names(dt1)<-c("Correlations Mean Average", "Correlations Variance Average", "Correlations Mean SD", "Correlations Variance SD",
                      "MSE Mean Average", "MSE Variance Average", "MSE Mean SD", "MSE Variance SD")
    
    
    
    
    
    cor.tabs2 <- strsplit(cor.tabs[2], "\n")[[1]]
    cor.tabs2 <- cor.tabs2[!cor.tabs2==""]
    dt2 <- read.table(text=gsub(".*(Mean|Variance) *", "", cor.tabs2))
    dt2 <- sapply(seq(1, nrow(dt2), by=2), function(y) unlist(dt2[y:(y+1),]) )
    colnames(dt2) <- gsub(" *", "", gsub("(Mean|Variance).* *", "",   cor.tabs2)) [seq(1, length(cor.tabs2), by=2)]
    rownames(dt2) <- c("Correlation of average means with true", "Correlation of average variances with true",
                       "MSE of average means with true","MSE of average variances with true")
    
    list('Correlations and mean square error of population and estimate values' = dt1,
         'Correlation and mean square error of the average estimates' = dt2)
    
    
  })
  
  names(otp)<- sim.outputs
  
  if(!silent) { 
    
    for(i in colnames(otp[[1]][[1]])  ) {
      cat("\n", "⎯⎯⎯⎯⎯⎯⎯⎯⎯ ", i," ", rep("⎯", getOption("width", 80)-nchar(i)-2),  "\n", sep="") 
      print(sapply(otp, function(x) x[[1]][,i] ), digits = 2)
      print(sapply(otp, function(x) x[[2]][,i] ), digits = 2)
    }
  }
  
  invisible(otp)
  
}