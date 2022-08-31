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

#' Get the number of groups which show invariance with each of the groups
#' @param fit Output of \code{\link{incrementalFit}}.
#' @param sort  Logical. Should the results be sorted?
#' @param fit.index Character. What measure should be used?
#' @param cutoff The cutoff of the difference in `fit.index` to decide whether invariance holds.
#' @param drop Character list of the groups that should be excluded from calculations.
#'
#' @returns  Returns a data frame with two columns: group name and the number of groups with invariant parameters.
#'
#' @export

n_invariant <- function(fit, sort=F, fit.index="cfi", cutoff=.01, drop=NA) {
  
  a=cbind(attr(fit$bunch,"pairs.of.groups"), fit$detailed[[fit.index]])
  if(!is.na(drop)) a = a[!( a$V1 %in% drop |  a$V2 %in% drop),]
  a$invariant <- a$fit.decrease<cutoff
  a1=reshape2::melt(a[, c("V1", "V2", "invariant")], id.vars = "invariant")
  a2 = aggregate(a1$invariant, list(a1$value), sum )
  
  if(sort) 
    return(a2[order(a2$x),])
  else 
    return(a2)
  
}


#' Get the tree of 
#' @param ... See \code{\link{n_invariant}} except `sort` and `drop` arguments.
#' @returns  Returns a data frame with N-1 columns, where N is number of groups.
#' @details Runs \code{\link{n_invariant}} iteratively, by dropping the least invariant group at each iteration until there are 2 groups left.
#' @export
#' 
n_invariant_matrix <- function(...) {
  
  n.inv <- MIE:::n_invariant(..., sort=T, drop=NA)
  dropped <-NA
  n.inv.list <- n.inv
  
  for(i in 1:(length(n.inv[,1])-2) ) {
    
    dropped <- na.omit(c(dropped, n.inv[1,1]))
    
    n.inv <- MIE:::n_invariant(..., sort=T, drop=dropped)
    
    n.inv.list <- merge(n.inv.list, n.inv, by = "Group.1", all.x=T)
    
  }
  
  
  sorted.n.inv.list <- n.inv.list[order(rowSums(is.na(n.inv.list[,-1])), decreasing = T),]
  
  print( as.matrix(sorted.n.inv.list) , na.print = "" , quote = FALSE )
  invisible(n.inv.list)
}

#' Append lavaan syntax with group-specific covariances
#'
#' @param model character, lavaan syntax model
#' @param group character, grouping variable
#' @param data data frame
#' @param cov character, covariance to add, e.g. "variable1 ~~ variable2"
#' @param focal.groups Character vector for the groups to add the cov.
#' @examples cov.model <-  "F =~ v1 + v2 v3 + v4 + v5"
#' cov.model.custom.covs <-
#'    lav.mod %>%
#'      add_custom_covs("country", dat1,
#'                      "v1 ~~ v3", c("China", "Indonesia")) %>%
#'
#'      add_custom_covs("country.f", dat1,
#'                      "v2 ~~ v3", c("Israel"))
#'
#' @export
add_custom_covs <- function(model, group, data, cov, focal.groups) {
  
  vector.zeros.nas <- paste(collapse="",
                            capture.output(
                              dput(
                                sapply(unique(data[,group]),function(x) ifelse(x %in% focal.groups , NA, 0)))
                            ))
  new.synt <- paste("\n  ",
                    strsplit(cov, "~~")[[1]][1],
                    " ~~ ",
                    vector.zeros.nas,
                    "*",
                    strsplit(cov, "~~")[[1]][2])
  
  paste(model,  new.synt, collapse = ";\n  ")
  
}


#' Run clustered measurement invariance tests
#'
#' @param model character, lavaan syntax model
#' @param group character, grouping variable
#' @param data data frame
#' @param strata A list of character vectors of the group names to create strata.
#' @param parameters character vector, "all", "loadings", "thresholds", or "intercepts". strata are applied to this subset of parameters.
#' @examples clusteredMI("F =~ v1 + v2 + v3 + v4", 
#'          group = "country", 
#'          data = Dat1, 
#'          strata = list(North = c("Norway", "Denmark", "Finland"), 
#'                          South = c("Spain", "Portugal", "Italy")
#'                          )
#'           )
#'
#' @export

stratifiedMI <- function(model, group, data, strata, parameters = c("loadings", "intercepts"), ref = "configural", additional = c("scalar"), ...) {
  
  config = cfa(model=model, data=data,  group=group,  do.fit = F, ...)
  model2 <- lavaanify(model, ngroups = lavInspect(config, "ngroups"), meanstructure = T, auto=T)
  op <- c("loadings" = "=~","intercepts" = "~1", "thresholds"= "|")[parameters]
  labs <- lavInspect(config, "group.label")
  
  for(m in names(strata) ) {
    if(length(strata[[m]])>1) {
      
      groupings = combn(match(strata[[m]], labs), 2)
      for(i in 1:ncol(groupings)) {
        labs1 = model2[model2$group==groupings[1,i] & model2$op %in% op & model2$free!=0, "plabel" ]
        labs2 = model2[model2$group==groupings[2,i] & model2$op %in% op & model2$free!=0, "plabel" ]
        nd<- data.frame(id = (max(model2$id)+1):(max(model2$id)+length(labs1)),
                        lhs = labs1, 
                        op = rep("==",length(labs1)),
                        rhs = labs2, 
                        user = rep(2,length(labs1)),
                        block=rep(0,length(labs1)), 
                        group=rep(0,length(labs1)),
                        free=rep(0,length(labs1)), 
                        ustart=rep(NA,length(labs1)),exo=rep(0,length(labs1)),
                        label=rep("",length(labs1)), plabel=rep("",length(labs1)))
        model2 <- rbind(model2, nd)
      }}}
  
  if(ref == "metric") {
    metrc.model <- lavaanify(model, ngroups = lavInspect(config, "ngroups"), 
                             meanstructure = T, auto=T, group.equal = "loadings")
    model2 <- rbind(model2, metrc.model[metrc.model$op=="==",])
    
    # # ! Fix means in one group within each cluster
    # # remove constraints on means
    # means.except.1.group <- 
    #   model2$lhs %in% names(lavInspect(config, "mean.lv")[[1]]) &
    #   model2$op == "~1" & 
    #   model2$group != 1
    
    # remove constraints on means
    first.groups.in.strata <- sapply(strata, function(x) x[[1]])
    first.group.not.in.cluster <- labs[ !labs %in% as.vector(unlist(strata))]
    if(length(first.group.not.in.cluster)>0) {
      first.groups <- c(first.groups.in.strata, first.group.not.in.cluster)
    } else {
      first.groups <- first.groups.in.strata
    }
    
    means.except.1.group <- 
        model2$lhs %in% names(lavInspect(config, "mean.lv")[[1]]) &
        model2$op == "~1" &
        !model2$group %in% match(labs, first.groups)
  

    model2[means.except.1.group,"free"] <- (max(model2$free)+1):
                                            max(model2$free+sum(means.except.1.group))
    
  }
  
  fits <- list(
    
    clustered  = cfa(model=model2, data=data,  group=group, ...),
    configural = cfa(model=model, data=data,  group=group, ...),
    metric     = cfa(model=model, data=data,  group=group, group.equal = "loadings", ...),
    scalar     = cfa(model=model, data=data,  group=group, group.equal = c("loadings", "intercepts", "thresholds"), ...)
    )
  
  if(length(additional)==0) {
    
      LittleHelpers::lav_compare(fits[[ref]], fits[["clustered"]])
  
    } else if(length(additional)==1) {
    
      LittleHelpers::lav_compare(fits[[ref]], fits[["clustered"]], fits[[additional]])
    } else if(length(additional)==2) {
      LittleHelpers::lav_compare(fits[[ref]], fits[["clustered"]], 
                                 fits[[additional[1]]], fits[[additional[2]]])
    }
  
  invisible(fits)
  
}

