#' Run the same cfa within every group and extract fit indices
#'
#' @param model Model in lavaan syntax
#' @param data The data
#' @param group Character. Grouping variable.
#' @param out Legacy argument, keeps only `fit` or `models` parts of the model.
#' @param what Fit indices to print and return. Possible values are any fit measures returned by `fitmeasures()`, a special value `mod` adds a parameter with the largest modification index and its MI value.
#' @param ... Other arguments passed to lavaan::cfa
#'
#' @export


groupwiseCFA <- function(model,  data, group, ..., out = NULL,
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
  
  tb.countrywise <- lapply(fit.list, function(x) {
    conv = lavInspect(x, "converged")
    
    if(conv) fm.all = fitMeasures(x)
    
    fm = lapply(setNames(nm=what[what !="mod"]), 
                function(w) {
                  if(conv)
                    ifelse(w %in% names(fm.all), fm.all[[w]],  NA) 
                  else
                    return(NA)
                  })
    
    if("mod" %in% what) { 
        if(conv) mi = modindices(x, sort = T)
        fm$mod.ind =  ifelse(conv, gsub("~~",  " W ",  paste(mi[1,1:3], collapse = "")), NA)
        fm$mod.ind.v = ifelse(conv, round(mi[1,4]), NA)
      }

    fm = append(fm, setNames(conv, nm = "converged"), after = 0)
    fm$N = nobs(x)
    
    return(as.data.frame(fm))
  
    # data.frame(
    #     converged = conv, #x@optim$converged, 
    #     CFI=      ifelse (conv & "cfi" %in% what, fm["cfi"],  NA),
    #     RMSEA=    ifelse (conv & "rmsea" %in% what, fm["rmsea"],NA),
    #     CHI.sq =  ifelse (conv & "chisq" %in% what, fm["chisq"],NA),
    #     Pvalue =  ifelse (conv & "chisq" %in% what, fm["pvalue"],NA),
    #     mod.ind=  ifelse (conv & "mod" %in% what,  paste(mi[1,1:3], collapse = ""), ""),
    #     mod.ind.v=ifelse (conv & "mod" %in% what,  round(mi[1,4]), ""),
    #     N = nobs(x)
    #     )
    
    })
  #print(tb.countrywise)
  tb.countrywise1 <- Reduce("rbind", tb.countrywise)
  rownames(tb.countrywise1) <- names(tb.countrywise)

  # print some fits
  b=tb.countrywise1[order(tb.countrywise1[,grepl("cfi", colnames(tb.countrywise1))], 
                          decreasing=T), 
                    grepl("cfi|rmsea|pvalue|srmr|mod|N", colnames(tb.countrywise1))]
  # b$CFI    <- round(b$CFI, 3)
  # b$RMSEA  <- round(b$RMSEA, 3)
  # b$CHI.sq <- round(b$CHI.sq, 2)
  # b$Pvalue <- round(b$Pvalue, 3)
  print(b, digits=2, print.gap = 3)
  
  class(fit.list) <- append(class(fit.list) , "groupwiseCFA")
  
  if(!is.null(out)) {
    if(out == "fit") 
      invisible(tb.countrywise1)
      else if(out == "models")
        invisible(fit.list)
    
  } else  {
    invisible(list(fit = tb.countrywise1, 
                   models = fit.list))
  }
        
}

#' Summarize a list of CFAs produced by lavaan
#'
#' @param model.list A list of models in lavaan syntax. Can be an object of class groupwise.CFA
#'
#' @export
summary.groupwiseCFA <- function(model.list, what = c("cfi", "rmsea", "srmr"), diagnose = T, corrs = F, mod = F, n = F) { # list of fitted CFAs
 
  
  
  # fix duplicated names
  
  if(any(duplicated(names(model.list)))) {
    
    old.names = names(model.list)
    new.nmz = old.names
    
    for(x in old.names) 
      new.nmz[old.names== x] = paste0(old.names[old.names == x], ".",
                                      1:sum(old.names == x))
    
    names(model.list) <- new.nmz
  }
  
  
  
  fit.idx = sapply(model.list, function(x) {
    if(lavInspect(x, "converged"))
      fitMeasures(x)[what]
    else 
      sapply(what, function(yy) NA)
  })
  
  
  if(diagnose) {
    diagnostics = data.frame(
      converged=ifelse(sapply(model.list, lavInspect, "converged"), "", "NO"),
      corr.greater.1 =
        
        # sapply(model.list, function(x) 
        # ifelse(any(abs(unlist(lavInspect(x,"cor.lv")))>1), "YES", "")),
      
      ifelse(any(
        sapply(model.list, 
               function(x) 
                 any(abs(unlist(lavInspect(x, "cor.lv", drop.list.single.group = F))) > 1))
               
        ), "YES", ""),
      
      negat.lv.variance = 
        # sapply(model.list, function(x)
        # ifelse(any(diag(lavInspect(x, "cov.lv"))<0), "YES", "")),
      
      ifelse(any(
        sapply(model.list, 
               function(x) 
                 any(unlist(sapply(lavInspect(x, "cov.lv", drop.list.single.group = F), diag)) < 0)
               
        )), "YES", ""),
      
      
      negat.ov.variance = 
        # sapply(model.list, function(x)
        # ifelse(any(diag(lavInspect(x, "cov.ov"))<0), "YES", ""))
      
        ifelse(any(
          sapply(model.list, 
                 function(x) 
                   any(sapply(lavInspect(x, "cov.ov", drop.list.single.group = F), diag) < 0)
                 
        )), "YES", "")
    ) 
    
    out  = cbind(diagnostics, t(fit.idx))
    
    
    
  } else {
      out = t(fit.idx)
  }
  
  if(corrs) {
    
    corrs <- lapply(model.list, function(x) {
      corrs.mx = lavInspect(x, "cor.lv", drop.list.single.group = F) 
      corrs.mx.mlt = reshape2::melt(lapply(corrs.mx, function(x) { x[upper.tri(x, diag = T)] <- NA; x}))
      
      if(all(is.na(corrs.mx.mlt$value))) {
        data.frame(factors = as.character(NA),
                   cor.factors = NA)
      } else {
        corrs.mx = filter(corrs.mx.mlt, !is.na(value)) %>%
          filter(value == max(value))
        data.frame(factors = paste(ifelse(is.na(corrs.mx[[4]]), "", paste0(corrs.mx[[4]], ":")), 
                                   corrs.mx[[1]], corrs.mx[[2]]),
                   cor.factors = corrs.mx[[3]])
      }
    })
    
    max.corrs = melt(corrs, id.vars = c("factors")) %>% select(1,3) %>% set_colnames(c("Max_Cor_Factors", "Max_Corr"))
    

    
    out = cbind(out, max.corrs)
    
  }
  
  if(mod) {
    
    max.mi <- sapply(model.list, function(x) {
      if(lavInspect(x, "converged")) { 
        mi = modindices(x, sort = T)
        max.mi = mi[1,1:4]
        max.mi.name =  gsub("~~",  " W ",  paste(max.mi[1:3], collapse = ""))
        max.mi.value = max.mi[[4]]
        return(data.frame(max.mi.name = max.mi.name, max.mi.value = max.mi.value))
      } else {
        return(data.frame(max.mi.name = NA, max.mi.value = NA))
      }
      
    })
    
    out = cbind(out, t(max.mi))
    
    }
    
  if(n) {
    
    n <- sapply(model.list, function(x) {
      if(lavInspect(x, "converged")) { 
        sum(lavInspect(x, "nobs"))
      } else {
        NA
      }
      
    })
    
    out = cbind(out, n = n)
    
  }

  
  return(out)
     
}


#' @rdname summary.groupwiseCFA
#' @export
gsumCFA <- summary.groupwiseCFA



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

mgcfaDiagnose <- function(lavaan.model, output=c("overall", "neg.var", "mi")) {
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
#' @noRd
nInvariant_ <- function(fit, sort=F, fit.index="cfi", cutoff=.01, drop=NA) {
  
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


#' Get the number of groups which show invariance with each of the groups
#' @param fit Output of \code{\link{incrementalFit}}.
#' @param sort  Logical. Should the results be sorted?
#' @param fit.index Character. What measure should be used?
#' @param cutoff The cutoff of the difference in `fit.index` to decide whether invariance holds.
#' @param drop Character list of the groups that should be excluded from calculations.
#' @param df Logical, if the data.frame of the comparison tree should be returned. If TRUE, `sort` and `drop` are ignored.
#'
#' @returns  If df = FALSE, returns a data frame with two columns: group name and the number of groups with invariant parameters. IF df = TRUE, calculates the number of gorups iteratively, by dropping the least invariant group at each iteration until there are 2 groups left. Returns a data frame with N-1 columns, where N is number of groups.

#' @export
nInvariant <- function(fit, sort=F, fit.index="cfi", cutoff=.01, drop=NA, df=FALSE) {
  
  if(!df) {
    MIE:::nInvariant_(fit = fit, fit.index = fit.index, cutoff = cutoff, sort=sort, drop=drop)
  } else {
    
  n.inv <- MIE:::nInvariant_(fit = fit, fit.index = fit.index, cutoff = cutoff, sort=T, drop=NA)
  dropped <-NA
  n.inv.list <- n.inv
  
  for(i in 1:(length(n.inv[,1])-2) ) {
    
    dropped <- na.omit(c(dropped, n.inv[1,1]))
    
    n.inv <- MIE:::nInvariant_(fit = fit, fit.index = fit.index, cutoff = cutoff, sort=T, drop=dropped)
    
    n.inv.list <- merge(n.inv.list, n.inv, by = "Group.1", all.x=T)
    
  }
  
  sorted.n.inv.list <- n.inv.list[order(rowSums(is.na(n.inv.list[,-1])), decreasing = T),]
  print( as.matrix(sorted.n.inv.list) , na.print = "" , quote = FALSE )
  invisible(n.inv.list)
}}

#' Append lavaan syntax with group-specific covariances
#'
#' @param model character, lavaan syntax model
#' @param group character, grouping variable
#' @param data data frame
#' @param cov character, covariance to add, e.g. "variable1 ~~ variable2"
#' @param focal.groups Character vector for the groups to add the cov.
#' @examples
#' \dontrun{
#' cov.model <-  "F =~ v1 + v2 v3 + v4 + v5"
#' cov.model.custom.covs <-
#'    lav.mod %>%
#'      add_custom_covs("country", dat1,
#'                      "v1 ~~ v3", c("China", "Indonesia")) %>%
#'
#'      add_custom_covs("country.f", dat1,
#'                      "v2 ~~ v3", c("Israel"))
#' }
#'
#' @noRd
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
#' @details
#' This function builds a single model with constraints applied to subsets of groups, and compares it to the reference model (less constrained) as well as to global invariance tests.
#' 
#' @examples  
#' \dontrun{
#' stratifiedMI("F =~ v1 + v2 + v3 + v4", 
#'          group = "country", 
#'          data = Dat1, 
#'          strata = list(North = c("Norway", "Denmark", "Finland"), 
#'                          South = c("Spain", "Portugal", "Italy")
#'                          )
#'           )
#' }
#' @export

stratifiedMI <- function(model, group, data, strata, parameters = c("loadings", "intercepts"), ref = "configural",  ...) {
  
  #additional = c("scalar"),
  
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
  
  # if(length(additional)==0) {
  #   
       MIE:::lav_compare(fits[[ref]], fits[["clustered"]])
  # 
  #   } else if(length(additional)==1) {
  #   
  #     LittleHelpers::lav_compare(fits[[ref]], fits[["clustered"]], fits[[additional]])
  #   } else if(length(additional)==2) {
  #     LittleHelpers::lav_compare(fits[[ref]], fits[["clustered"]], 
  #                                fits[[additional[1]]], fits[[additional[2]]])
  #   }
  
  invisible(fits)
  
}


lav_compare <- function(..., what = c("cfi", "tli", "rmsea", "srmr", "bic", "df"), LRT = F, print  = T) {
  
  
  if(length(what)<1) warning("Please choose at least one fit statistic to report.")
  if(length(list(...))==1 & inherits(list(...)[[1]], "list") ) {
    modellist = list(...)[[1]]
    modelnames <- names(modellist)
  } else {
    modellist = list(...)
    if(is.null(names(modellist)))
      modelnames <- as.character(substitute(...()) )
    else
      modelnames <- names(modellist)
  }
  
  
  out2<- t(sapply(modellist,  function(x) {
    all.fit= fitmeasures(x)
    sapply(what, function(y) ifelse(any(names(all.fit) == y), all.fit[[y]], NA))
  }))
  diffs <- apply(out2, 2, function(v) v - c(NA, v[-length(v)]))
  out2 <- cbind(out2, diffs)
  
  
  
  rownames(out2)<-modelnames
  
  out = out2
  
  if(LRT) {
    # LRT <- do.call(lavaan::lavTestLRT, append(modellist, list(model.names = modelnames)))
    # above does not work properly
    args.for.lav.LRT = append(modellist, list(modelnames))
    names(args.for.lav.LRT) <- c("object", rep("", length(modellist)-1), "model.names")
    LRT <- do.call(lavaan::lavTestLRT, args.for.lav.LRT)
    
    if(print) print(LRT)
    out = list(fit=out2, LRT = LRT)
  }
  
  if(print) {
    cat("\n")
    print(round(out2, 3), digits=3, row.names = TRUE, na.print = "" )
  }
  invisible(out)
}

