#' Run the same cfa within every group and extract fit indices
#'
#' @param model Model in lavaan syntax
#' @param data The data
#' @param group Character. Grouping variable.
#' @param out Character. "fit" returns some fit indices for every group; "models" list of entire fittted models.
#' @param ... Other arguments passed to lavaan::cfa
#'
#' @export


groupwiseMI <- function(model,  data, group, ..., out = c("fit", "models")) {
  fit.list <- lapply(unique(data[, group]), function(gr) {
    print(gr)
    lavaan::cfa(model, data = data[data[, group]==gr, ], ...)
  })
  
  names(fit.list)<- unique(data[, group])
  
  tb.countrywise <- lapply(fit.list, function(x) data.frame(
    converged = x@optim$converged, 
    CFI=ifelse (x@optim$converged, fitMeasures(x)[c("cfi")],  NA),
    RMSEA=ifelse (x@optim$converged, fitMeasures(x)[c("rmsea")],NA),
    mod.ind=ifelse (x@optim$converged, 
                    paste(modindices(x, sort = T)[1,1:3], collapse = ""), ""),
    mod.ind.v=ifelse (x@optim$converged,  round(modindices(x, sort = T)[1,4]), ""),
    stringsAsFactors = F))
  tb.countrywise1 <- Reduce("rbind", tb.countrywise)
  rownames(tb.countrywise1) <- names(tb.countrywise)
  tb.countrywise1
  tb.countrywise1$mod.ind <-  gsub("~~",  " W ",  tb.countrywise1$mod.ind )
  
  
  if("models" %in% out) return(fit.list)
  if("fit" %in% out) { 
    df_to_viewer(tb.countrywise1[order(tb.countrywise1$CFI, decreasing=T), c("CFI", "RMSEA")])
    invisible(tb.countrywise1)
  }
}