#' Run the same cfa within every group and extract fit indices
#'
#' @param model Model in lavaan syntax
#' @param data The data
#' @param group Character. Grouping variable.
#' @param out Character. "fit" returns some fit indices for every group; "models" list of entire fittted models.
#' @param what Fit indices to print in case out="fit". Possible values: "cfi", "rmsea", "chisq", "mod" (parameter with the largest modification index and its value).
#' @param ... Other arguments passed to lavaan::cfa
#'
#' @export


groupwiseCFA <- function(model,  data, group, ..., out = c("fit", "models"),
                         what  = c("cfi", "rmsea", "chisq", "mod")) {
  
  if(length(data[, group])==1) {
    group.names = na.omit(unique(data[, group][[1]])) 
  } else {
  group.names = na.omit(unique(data[, group]))
  }
  
  fit.list <- lapply(group.names, function(gr) {
    print(gr)
    lavaan::cfa(model, data = data[data[, group]==gr, ], ...)
  })
  
  names(fit.list)<- group.names
  
  if("fit" %in% out) { 
  tb.countrywise <- lapply(fit.list, function(x) {
    if(x@optim$converged) fm = fitMeasures(x)
    if(x@optim$converged & "mod" %in% what) mi = modindices(x, sort = T)
    
    data.frame(
      converged = x@optim$converged, 
      CFI=      ifelse (x@optim$converged & "cfi" %in% what, fm["cfi"],  NA),
      RMSEA=    ifelse (x@optim$converged & "rmsea" %in% what, fm["rmsea"],NA),
      CHI.sq =  ifelse (x@optim$converged & "chisq" %in% what, fm["chisq"],NA),
      Pvalue =  ifelse (x@optim$converged & "chisq" %in% what, fm["pvalue"],NA),
      mod.ind=  ifelse (x@optim$converged & "mod" %in% what,  paste(mi[1,1:3], collapse = ""), ""),
      mod.ind.v=ifelse (x@optim$converged & "mod" %in% what,  round(mi[1,4]), ""),
      N = nobs(x),
      stringsAsFactors = F)
    
    })
  tb.countrywise1 <- Reduce("rbind", tb.countrywise)
  rownames(tb.countrywise1) <- names(tb.countrywise)
  tb.countrywise1
  tb.countrywise1$mod.ind <-  gsub("~~",  " W ",  tb.countrywise1$mod.ind )
  }
  
  if("models" %in% out) return(fit.list)
  if("fit" %in% out) { 
    b=tb.countrywise1[order(tb.countrywise1$CFI, decreasing=T), 
                      c("CFI", "RMSEA", "CHI.sq", "Pvalue")]
    b$CFI <- round(b$CFI, 2)
    b$RMSEA <- round(b$RMSEA, 2)
    b$CHI.sq <- round(b$CHI.sq, 2)
    b$Pvalue <- round(b$Pvalue, 3)
    print(b)
    invisible(tb.countrywise1)
  }
}

#' Get more comprehensible output from lavTestScore
#'
#' Simply applies \code{\link{lavTestScore}} and attaches group labels and parameter names to the output.
#'
#' @param lavaan.fit Model fitted with lavaan
#' @param ... Arguments passed to lavTestScore
#'
#' @export

lavTestScore_clean <- function(lavaan.fit,  ...) {
  require("lavaan")
  lvts <- lavTestScore(lavaan.fit, ...)
  
  for(lvts.part in names(lvts)[names(lvts) %in% c("uni", "cumulative")]) {
    
    partab.a<- partable(lavaan.fit)[,c(c("lhs", "op",  "rhs", "group",  "plabel"))]  %>%
      dplyr::filter(., plabel!="")
    
    names(partab.a)[1:3] <- c("one", "two", "three")
    
    out<- merge(as.data.frame(lvts[[lvts.part]]),
                partab.a,
                by.x=c("lhs"), by.y=c("plabel"),
                all.x=T)
    out2<- merge(out,
                 partab.a,
                 by.x=c("rhs"), by.y=c("plabel"),
                 all.x=T, suffixes = c(".lhs", ".rhs"))
    
    out2$group.lhs <- factor(out2$group.lhs, levels=1:length(lavInspect(lavaan.fit, "group.label")), labels=lavInspect(lavaan.fit, "group.label"))
    out2$group.rhs <- factor(out2$group.rhs, levels=1:length(lavInspect(lavaan.fit, "group.label")), labels=lavInspect(lavaan.fit, "group.label"))
    
    out3 <- data.frame(Term = paste(out2$one.lhs, out2$two.lhs, out2$three.lhs, sep=""),
                       Group1 = out2$group.lhs,
                       Group2 = out2$group.rhs,
                       Chi.square=round(out2$X2, 3), df=out2$df, p.value=round(out2$p.value,3),
                       "."=format(as.character(sapply(out2$p.value, function(x) ifelse(x>0.05, "", ifelse(x>0.01, "*", ifelse(x>0.001, "**", "***"))))), justify = "left")
    )
    
    lvts[[lvts.part]]<-out3
    if(lvts.part=="uni") attr(lvts[[lvts.part]], "header") <- "Chi-square test of releasing single constraints, equivalent to modification indices"
    if(lvts.part=="cumulative") attr(lvts[[lvts.part]], "header") <- "Chi-square test of releasing multiple constraints at the same time"
    class(lvts[[lvts.part]]) <- c("lavaan.data.frame","data.frame")
    
  }
  
  
  if(any(names(lvts) == c("epc"))) {
    lvts[["epc"]]$group <- factor(lvts[["epc"]]$group,
                                  levels=1:length(lavInspect(lavaan.fit, "group.label")),
                                  labels=lavInspect(lavaan.fit, "group.label"))
  }
  
  
  return(lvts)
  
  
}

#' Quick diagnostics for multiple-group lavaan models
#' @param lavaan.model Fitted lavaan multiple group model.
#' @param output A character list of what diagnostics should be computed. 
#' #' \describe{
#'   \item{_overall_}{Prints chi-sq, CFI, RMSEA, and SRMR }
#'   \item{_neg.var_}{Checks if there are negative variances of latent variables, if yes, identifies the groups and prints them.}
#'   \item{_mi_}{Aggregates modification indices across groups, finding the most impactful. Also prints the top part of sorted table of modification indices.}
#'   \item{_constraints_}{ Optional. Applies and prints \code{\link{lavTestScore_clean}}.}
#' }
#'
#' @details Doesn't return anything.
#'
#'
#' @md
#'
#' @export

mgcfa_diagnose <- function(lavaan.model, output=c("overall", "neg.var", "mi")) {
  require(magrittr)
  if("overall" %in% output) {
    # Print general info
    cat(paste("\n\n  Model fitted to ", lavInspect(lavaan.model, "ngroups"), "groups, and ",
              sum(lavInspect(lavaan.model, "nobs")), "observations.\n"))
    
    # Print overall fit measures
    cat("\n\n  Overall fit measures: \n")
    a<-fitmeasures(lavaan.model)
    fit.mes.index <- c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr")
    if( all( c("cfi.scaled", "rmsea.scaled")  %in% names(a))) 
      fit.mes.index <- append(fit.mes.index, c("cfi.scaled", "rmsea.scaled"))
    
        
    
    print(data.frame(Index=round(a[fit.mes.index],3),
                     df=c(round(a["df"]),  rep("", length(fit.mes.index)-1)),
                     p.value=c(round(a["pvalue"],4),  rep("", length(fit.mes.index)-1))
    ), na.print=NULL, quote=F, print.gap=4)
  }
  # Print negative LV variances if there are any
  if("neg.var" %in% output) {
    
    a<-lavInspect(lavaan.model,"cov.lv") %>% sapply(diag) %>% t %>% round(5)
    if(any(a<0)) {
      cat("\n\n  Negative LV variances found: \n")
      if(sum(rowSums(a<0))==1) {
        data.frame(a[rowSums(a<0)>0,]) %>% t %>%
          `row.names<-`(row.names(a)[rowSums(a<0)>0]) %>%
          format(digits=4) %>% print(quote=F, print.gap = 4)
      } else {
        a[rowSums(a<0)>0,] %>%
          format(digits=4) %>% print(quote=F, print.gap = 4, row.names=T)
      }
    } else {
      cat("\n\n  All latent variables' variances are positive. \n")
    }
    rm(a)
    
    b<-lavInspect(lavaan.model,"cov.ov") %>% sapply(diag) %>% t %>% round(5)
    if(any(b<0)) {
      cat("\n\n  Negative LV variances found: \n")
      if(sum(rowSums(b<0))==1) {
        data.frame(b[rowSums(b<0)>0,]) %>% t %>%
          `row.names<-`(row.names(b)[rowSums(b<0)>0]) %>%
          format(digits=4) %>% print(quote=F, print.gap = 4)
      } else {
        b[rowSums(b<0)>0,] %>%
          format(digits=4) %>% print(quote=F, print.gap = 4, row.names=T)
      }
    } else {
      cat("\n\n  All residuals' variances are positive. \n")
    }
    rm(b)
    
  }
  # Find most common misspecisfications
  if("mi" %in% output) {
    a<-try(modindices(lavaan.model, sort.=T))
    
    if(any(class(a)=="try-error")) {
      cat("\n\n  *modification indices cannot be computed.")
    } else {
      cat("\n\n  Aggregate modification indices \n")
      aggregate(a$mi, list(paste(a$lhs,a$op, a$rhs)), sum) %>%
        `names<-`(c("Term", "Sum of modif.indices across groups")) %>%
        `[`(., order(`[`(., ,2), decreasing = T),) %>%
        print(right = F, row.names=F)
      # Head of modification indices list
      cat("Head of modification indices list \n")
      print(as.data.frame(a)[1:10,], row.names=F)
      cat("* Use 'modindices()' to see the full list.")
    }
    rm(a)
  }
  # Find problems with constraints
  if("constraints" %in% output) {
    cat("\n\n What if some constraints are relaxed?  \n")
    try(lavTestScore_clean(lavaan.model))
  }
}