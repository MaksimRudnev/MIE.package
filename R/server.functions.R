## Global functions #####

#' Test for measureent invariance in three steps
#' @description Computes configural, metric, and scala invariance models and compares them.
#' 
#' @param ... Formula, group, and all the other eligible arguments of \code{\link[lavaan]{cfa}} function.
#'
#' @export
globalMI <- function(..., chi.sq=FALSE, omit = "") {
  
  if(!"configural" %in% omit) {
    print("Running configural model")
    configural<-lavaan::cfa(...)
    
    # Refit with different method if models didn't converge with defaults
    if(!configural@optim$converged) {
      print("Refitting configural")
      configural  <-lavaan::cfa(..., optim.method = "L-BFGS-B")
    }
    
  }
  
  thresholds.in.the.model <- any(lavaan::parameterTable(configural)$op =="|")
  
  if(thresholds.in.the.model & "thresholds" %in% omit) {
    print("Running thresholds model")
    thresholds <- lavaan::cfa(..., group.equal = c("thresholds"))
    
    if(!thresholds@optim$converged) {
      print("Refitting thresholds") 
      thresholds <- lavaan::cfa(..., group.equal = c("thresholds"), optim.method = "L-BFGS-B")
    }
    
  } else if(!"metric" %in% omit){
    
    print("Running metric model")
    metric<-lavaan::cfa(..., group.equal = "loadings")
    
    if(!metric@optim$converged) {
      print("Refitting metric")
      metric<-lavaan::cfa(..., optim.method = "L-BFGS-B")
    }
  }
  if(!"scalar" %in% omit) {
    print("Running scalar model")
    scalar<-lavaan::cfa(..., group.equal = c("loadings", "intercepts", "thresholds"))
    
    if(!scalar@optim$converged) {
      print("Refitting scalar")
      scalar<-lavaan::cfa(..., optim.method = "L-BFGS-B")
    }
  }

  # Summarize the results

 # if(thresholds.in.the.model)  {
    
  mdls =  c("configural", "metric", "scalar")[(!c("configural", "metric", "scalar") %in% omit)]
  if(thresholds.in.the.model) mdls = ifelse("configural" %in% mdls, 
                                            append(mdls, "thresholds", after = 1), 
                                            append(mdls, "thresholds"))
  
  out <- try(do.call(lavaan::lavTestLRT, lapply(mdls, as.name)))
  
     
    
  fit.mes.index = c("cfi", "tli", "rmsea", "srmr", "chisq", "df")
  if( all( c("cfi.scaled", "tli.scaled", "rmsea.scaled")  %in% names( do.call(lavaan::fitMeasures, list(as.name(mdls[1])))))) 
    fit.mes.index <- append(fit.mes.index, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "chisq.scaled"))
  

  
  
  model.list <- do.call(list, lapply(mdls, as.name))
  
  out2<-t(sapply(model.list, function(x) 
    if(x@optim$converged) 
      lavaan::fitMeasures(x)[fit.mes.index] 
    else 
      rep(NA, length(fit.mes.index))
    ))
  
  
  # out2.2 <- apply(out2, 2, function(x) (c(NA, x[2]-x[1], x[3]-x[2]  )))
  # colnames(out2.2)<- paste("diff.", toupper(colnames(out2.2)), sep="")
  # out2.3 <- t(Reduce("rbind", lapply(1:ncol(out2), function(x) rbind(out2[,x], out2.2[,x]))))
  # colnames(out2.3)<-toupper(as.vector(sapply(1:length(colnames(out2)), function(i) c(colnames(out2)[i], colnames(out2.2)[i]) )))

  out2.3 = Reduce("cbind", lapply(1:ncol(out2), function(column) {
    a= data.frame(
      out2[,column], 
      diff=sapply(1:length(out2[,column]), 
                  function(x) if(x==1) NA else out2[x,column] - out2[x-1,column])
    )
    
    names(a) <- c(colnames(out2)[column], paste0("diff.", colnames(out2)[column] ))
    a
  }))
  
  #rownames(out2.3)<-c("Configural", "Metric", "Scalar")
  #if(any(lavaan::parameterTable(r.conf)$op =="|")) rownames(out2.3)<-c("Configural", "Thresholds", "Scalar")
  rownames(out2.3) <- mdls
  
  if(chi.sq) {
    print(fit.mes.index)
    chi=ifelse(any(fit.mes.index=="chisq.scaled"), "chisq.scaled", "chisq")
    out2.3$chisq.df <- paste0(round(out2.3[,chi],1), "(", out2.3[,"df"], ")")
    out2.3[,chi] <- NULL        
    out2.3[,"df"]  <- NULL
    out2.3
  } else {
    
    out2.3[,"chisq.scaled"]<- NULL
    out2.3[,"chisq"]       <- NULL
    out2.3[,"df"]          <- NULL
    
  }
    
  print(out)
  cat("\n")
  print(#round(out2.3, 3), 
        out2.3,
        digits=3, row.names = TRUE, na.print = "" )
  
  invisible(out2.3)
}


#' Run Shiny MIE application
#' @param model Character, lavaan formula. Special value "demo" to use simulated dataset and model.
#' @param data Data.frame containing only the variables to be used in calculations. May contain only numeric variables. If `group` argument is not specified, the group variable should be the first in the data.
#' @param group Character, name of the variable in data.
#' 
#' @examples if (interactive()) { 
#' runMIE("f =~ a1 + a2 + a3", mydr, "municipality") 
#' })
#' 
#' @export
runMIE <- function(model = NULL, data = NULL,  group = NULL, verbose = FALSE) {
  
  
  #if(!is.null(data) & is.null(group)) stop("The grouping variable should be specified.")
  
  .GlobalEnv$.data <- data 
  .GlobalEnv$.model <- model
  .GlobalEnv$.group <- group
  .GlobalEnv$.verbose <- verbose
  on.exit(rm(list=c(".data", ".model", ".group"), envir=.GlobalEnv))                       
  shiny::runApp(appDir = system.file("application", package = "MIE"))
  #shiny::runApp(appDir = "inst/application")
}

#' Returns all possible pairs of countries without duplicates based on variable 
#' 
#' @param variable Grouping variable to find the set of all possible unique pairs of values.
#' 
pairs_of_groups <- function(variable) {
  as.data.frame(t(utils::combn(unique(as.character(variable)), 2)), stringsAsFactors = F)
}

#.. Compute the MGCFA model for a given pair of groups #####
#' Computes MGCFA model for a given pair of groups
#' @param model Model in lavaan syntax
#' @param data Data object
#' @param group Character. Group variable name.
#' @param constraints Can be "", "loadings", or c("loadings", "intercepts").
#' @param pairs.of.groups Full list of pairs of groups, used in shiny only.
#' @param message Notification in shiny.
#' @param shiny If it is executed in shiny environment.
#' @param signed Only for DMACS, if the signed version should be used. See Nye et al 2019
#' @param ... Arguments passed to lavaan 'cfa' function.
# @example # pairwiseFit(model="F =~ impfree + iphlppl + ipsuces + ipstrgv", data = ess6, group = "cntry", "loadings")
#' 
#' @return  The function returns matrix of fit indices for multiple group CFA models fitted to each possible pair of groups.
#' @details Mostly for internal use  within \code{\link[MIE]{incrementalFit}}
#' 
#' @export
pairwiseFit <- function(model,
                        data, 
                        group = "cntry",
                        constraints=c(""), 
                        pairs.of.groups=NULL, 
                        message="Fitting pairwise lavaan models",
                        shiny=FALSE,
                        signed = F
                        , ...
                        
) {
  
  if(is.null(pairs.of.groups)) 
    pairs.of.groups <- MIE:::pairs_of_groups(data[[group]])
  
  # rearrange columns so that the group variable goes first
  data=data[,c(group, colnames(data)[colnames(data)!=group])]
  
  colnames(data)[1]<-"cntry"
  
  runPairwiseModels <- function(...) {
    
    model.lav<- lavaan::cfa(model, 
                            data=data[data$cntry==pairs.of.groups[1,1] | 
                                       data$cntry==pairs.of.groups[2,2],],
                    group="cntry",
                    #group=group, 
                    group.equal=constraints
                    , ...
    )
    
    print("LAV CALL")
    print(model.lav@call)
    
    # FN is a number of fit indices currently provided by lavaan  (currently 41) + 1 for dmacs
    FN = length(lavaan::fitmeasures(model.lav)) + 1
    
    if(lavaan::lavInspect(model.lav, "converged")) {
      mod<- lavaan::fitmeasures(model.lav) 
      # add dmacs
      mod<- c(mod, average.dmacs =  mean(MIE:::dmacs_lavaan(model.lav, signed = F), na.rm=T))
      
      dmacs.signed.list<-list()
      if(signed) dmacs.signed.list[[1]] <- MIE:::dmacs_lavaan(model.lav, signed = T)
      
      
      
    } else { 
      mod <- rep(999, FN)
      mod <- c(mod, average.dmacs = 999)
    }
    
    
    mod<-matrix( c(mod, rep(0,FN*(nrow(pairs.of.groups)-1))), nrow=FN, dimnames=list(names(mod), NULL))
    
    #Non-positive definite?
    non.positive <- FALSE
    
    # Uses 'for' instead of 'sapply' in order to show the progress bar
    for(x in 2:nrow(pairs.of.groups)) {
      model.lav<- lavaan::cfa(model, 
                              data=data[data$cntry==pairs.of.groups[x,1] | 
                                         data$cntry==pairs.of.groups[x,2],],
                      group="cntry", 
                      #group = group,
                      group.equal=constraints
                      , ... # extra.options
      )
      
      #If converged, record fitmeasure; if not converged add a missing sign 999.
      if(lavaan::lavInspect(model.lav, "converged")) {
        mod[,x]<- c(lavaan::fitmeasures(model.lav),
                    average.dmacs = 
                      mean(MIE:::dmacs_lavaan(model.lav, signed = F),na.rm=T))
        
        
        if(signed) dmacs.signed.list[[x]] <- MIE:::dmacs_lavaan(model.lav, signed = T)
        
        
        
      } else {
        mod[,x]<- rep(999, FN)
      }
        
      
      #Save non-positive definite status
      
      non.positive[[x]]<-lavaan::lavInspect(model.lav, "post.check")==FALSE
      
      # # add dmacs
      # if(!non.positive[[x]]) {
      #   mod[,x]<- cbind(mod[,x], 
      #                   dmacs = 
      #                     mean(
      #                       dmacs::lavaan_dmacs(model.lav, 
      #                         RefGroup=lavaan::lavInspect(model.lav, "group.label")[1],  
      #                         "pooled" )$DMACS,
      #                       na.rm=T)
      #                   )
      # } else {  
      #   mod[,x]<- cbind(mod[,x], dmacs = NA)
      # }
      
      
      
     if(shiny) {
       incProgress(1/nrow(pairs.of.groups), 
                  detail = paste("Compute for pair of", pairs.of.groups[x,1], "and", pairs.of.groups[x,2]))
     } else {
       utils::setTxtProgressBar(pb, x/nrow(pairs.of.groups))
     }
    }
    attr(mod, "pairs.of.groups")<-pairs.of.groups
    #attr(mod, "model.formula")<-model
    
    # Show the status of non.positive
    if(sum(non.positive)>0 & shiny) { 
      showNotification(
      "Non-positive definite matrix for ",
      action=a(href = paste("javascript:alert('",
                            paste(apply(pairs.of.groups[non.positive,], 1, paste, collapse=" & "), 
                                  collapse=",\n"),
                            "');"), paste( sum(non.positive),  "models (paired groups subsamples).") ), type="warning", duration=NULL  )
   
       } else if(sum(non.positive)>0 & !shiny) {
    
      warning("Non-positive definite matrix for ", 
              paste(apply(pairs.of.groups[non.positive,], 1, paste, collapse=" & "), collapse=",\n"), paste(sum(non.positive),  "models (paired groups subsamples)."))
      
       }
    
    colnames(mod)<-apply(pairs.of.groups, 1, paste, collapse="_")
    attr(mod, "model.formula") <- model
    if(signed) attr(mod, "dmacs.signed.list")<-  dmacs.signed.list
     
    return(mod)
    
  }
  
    if(shiny) {
  withProgress(message = message, value = 0, {   #Create progress bar
    mod <- runPairwiseModels(...)
    mod
  }) #close progress bar 
    } else {
      pb <- utils::txtProgressBar(title="Fitting pairwise models", style=3)
      mod <- runPairwiseModels(...)
      mod
    }
}

# .. Extracts attribute ------

#' Extract attribute 'pairs.of.groups'
#' 
#' @param df Any object with an attribute "pairs.of.groups"
#' 
#' @export
get_pairs <- function(df)  attr(df, "pairs.of.groups")





# Computes MGCFA models for all available groups in the data (or their subset) -----
#' Computes MGCFA models for all available groups in the data
#' 
#' @param model Model in lavaan syntax
#' @param data The data
#' @param group Character. Grouping variable.
#' @param parameters Character. If "loadings" then configural model is fitted and loadings are returned. If "intercepts" then metric invariance model is fitted and intercepts are returned.
#' @param extra.options Currently not used
#' @param shiny Logical. If it is evaluated in a shiny context. Default is TRUE.
#' 
#'  
#' @export 
MGCFAparameters <- function(model=NULL,
                            data, 
                            group = "cntry", 
                            parameters = "loadings",  
                            extra.options=NULL, 
                            shiny = FALSE) {
  data=data[,c(group, colnames(data)[colnames(data)!=group])]
  if(is.null(model)) {
    if(shiny) {
    showNotification("The model is not specified", type="error")
    } else {
      warning("Please specify a model")
    }
    # Errors$nomodel = TRUE
    # c("error")
  } else {
    
    constraint <- ifelse(parameters == "intercepts", c("loadings"), "")
    operator   <- ifelse(parameters == "intercepts", "~1", "=~")
    #left.or.right <- ifelse(configural.or.metric=="metric", "lhs", "rhs")
    
    extractPars <- function() {
      #mod<-cfa(model, data, group="cntry", group.equal=constraint)
      old_option <- getOption("show.error.messages")
      options(show.error.messages = FALSE)
      
      cfa.argument.list <- c(extra.options, list(model=model, data=data, group=group, group.equal=constraint))
      
      mod<-try(do.call("cfa",  cfa.argument.list, quote = FALSE), silent=TRUE)
      options(show.error.messages = TRUE) 
      # Show error message or extract the parameters
      if(class(mod)=="try-error") {
        
        er1 <- paste(strsplit(mod, "ERROR: ")[[1]][2])
        if(length(er1>299)) er1 <- paste(strtrim( er1, 300), "...", sep="")
        
      if(shiny)  showNotification("Error in lavaan:", action=  er1, type="error", duration = NULL)
        stop("Stopped due to error in lavaan")
        rm(er1)
        
      } else {
        
        #print("lavCALL"); print(mod@call[-3])
        
        #If non-positive definite, show notification
        if(lavaan::lavInspect(mod, "post.check")==FALSE & shiny) {
          showNotification("The model produced non-positive definite matrix.",
                           action = a(href = "javascript:location.reload();", "Reload page"))
        } else if(lavInspect(mod, "post.check")==FALSE & !shiny) {
          warning("The model produced non-positive definite matrix.")
        }
        
        if(shiny) incProgress(1/2)
        # parameters<-cbind(mod@ParTable[["group"]][mod@ParTable$op==operator & mod@ParTable$free!=0],
        #                   mod@ParTable[[left.or.right]][mod@ParTable$op==operator & mod@ParTable$free!=0],
        #                   mod@ParTable[["est"]][mod@ParTable$op==operator & mod@ParTable$free!=0]
        # )
        # parameters.t<- parameters %>%  as.data.frame(.) %>% melt(., c("V1", "V2")) %>% dcast(., V1 ~ V2) 
        requireNamespace("magrittr")
        parameters.t <- lavaan::parTable(mod) %>%
          base::subset(., select=c("group", "lhs", "op", "rhs", "est"), subset = op==operator & free!=0 ) %>%
          dplyr::mutate(par=paste(lhs, ifelse(operator=="=~", "_by_", ""), rhs, sep="")) %>%
          reshape2::melt(., c("group", "par"), measure.vars="est" ) %>% reshape2::dcast(., group ~ par) 
        
        
        se.t <- lavaan::parTable(mod) %>%
          base::subset(., select=c("group", "lhs", "op", "rhs", "se"), subset = op==operator & free!=0 ) %>%
          dplyr::mutate(par=paste(lhs, ifelse(operator=="=~", "_by_", ""), rhs, sep="")) %>%
          reshape2::melt(., c("group", "par"), measure.vars="se" ) %>% reshape2::dcast(., group ~ par) 
        
        #print(paste("Fit model with", paste(mod@Data@group.label, collapse=",")))
        #print("parameters.t"); print(parameters.t)
        
        parameters<-matrix(as.numeric(unlist(parameters.t[,-1])), 
                           ncol=ncol(parameters.t)-1,
                           nrow=nrow(parameters.t),
                           dimnames=list(lavInspect(mod, "group.label"), #mod@Data@group.label, #parameters.t[,1], 
                                         names(parameters.t[,-1]) ))
        se<-matrix(as.numeric(unlist(se.t[,-1])), 
                           ncol=ncol(se.t)-1,
                           nrow=nrow(se.t),
                           dimnames=list(lavInspect(mod, "group.label"), #mod@Data@group.label, #se.t[,1], 
                                         names(se.t[,-1]) ))
        attr(parameters, "fit")<-fitmeasures(mod)
        attr(parameters, "se")<-fitmeasures(mod)
        class(parameters)<- c("MGCFAparameters")
        return(parameters)
      }
    }
    
    
    if(shiny) {
      withProgress(message = paste("Computing", parameters, "MGCFA model for ALL the selected groups."), value = 0, { 
      incProgress(1/2)
      
      mod <- extractPars()
      return(mod)
      
    })
    } else {
      mod<-extractPars()
      return(mod)
    }
    
  }
}

#.. Computes and formats covariance matrices ----
#' Computes and formats covariance matrices
#' @param data Data containing group variable
#' @param group Character, group variable
#' 
#' @export
computeCovariance <- function(data, group) {
  
  #    if(input$weights=="noweight") {
  # message("Computing covariance matrix")
  
  #Split dataset and compute variance-covariance for each group separately
  if(is.null(group)) {
    dat.split<-split(data[,2:ncol(data)], data[,1], drop = T)
  } else {
    dat.split<-split(data[,colnames(data)!=group], data[,group], drop = T)
  }
  
  
  #                                  #covariance matrix        #without duplications
  tab<-sapply(dat.split, function(x) cov(x, use="complete.obs")[lower.tri(var(x, use="complete.obs"), diag = F)])
  
  if(is.null(group)) {

    rownames(tab)<- apply(combn(names(data[,2:ncol(data)]), 2), 2, paste, collapse="_")
    
  } else {

    rownames(tab)<- apply(combn(names(data[,colnames(data)!=group]), 2), 2, paste, collapse="_")
    
  }
  
  tab <- t(tab)
  class(tab) <- "covariances"
  tab
  #New version
  #sapply(unique(dt$dat$cntry), function(x) 
  #  var(dt$dat[dt$dat$cntry==x,-1], use="complete.obs")[lower.tri(var(dt$dat[dt$dat$cntry==x,-1], use="complete.obs"),diag = F)]
  #)
}

#' Computes and formats Fisher's transformed correlation matrices
#' 
#' @param data Individual data
#' @param group Character. Grouping variable
#' 
#' @export
computeCorrelation <- function(data, group) {
  
  #    if(input$weights=="noweight") {
  # message("Computing correlation matrix")
  
  #Split dataset and compute variance-covariance for each group separately
  #
  if(is.null(group)) {
    dat.split<-split(data[,2:ncol(data)], data[,1], drop = T)
  } else {
  dat.split<-split(data[,colnames(data)!=group], data[,group], drop = T)
  }
 
  tab<-sapply(dat.split, function(x) {
    rho = cor(x, use="complete.obs")[lower.tri(var(x, use="complete.obs"), diag = F)]
    0.5 * log((1 + rho)/(1 - rho)) #Fisher z transformation
    })
  if(is.null(group)) {
    
    rownames(tab)<- apply(combn(names(data[,2:ncol(data)]), 2), 2, paste, collapse="_")
    
  } else {
    
    rownames(tab)<- apply(combn(names(data[,colnames(data)!=group]), 2), 2, paste, collapse="_")
    
  }
  tab <- t(tab)
  class(tab) <- "correlations"
  tab
  
  #New version
  #sapply(unique(dt$dat$cntry), function(x) 
  #  var(dt$dat[dt$dat$cntry==x,-1], use="complete.obs")[lower.tri(var(dt$dat[dt$dat$cntry==x,-1], use="complete.obs"),diag = F)]
  #)
}



#' Run pairwise models and compute decrease in fit
#'
#' @description Runs two MGCFA models for each possible pair of groups and computes change in fit.
#' @param ... The arguments passed to \code{\link[MIE]{pairwiseFit}}. Required arguments are \code{'model'}, \code{'data'}, and \code{'group'}.
#' @param level Character. A model set to be computed. The function will compute a decrease in fit between two models:
#' \describe{
#' \item{metric}{(default) between configural and metric models.}
#' \item{scalar}{between metric and scalar models.}
#' \item{at.once}{between configural and scalar models.}
#' \item{intercepts.first}{between configural and model with equal intercepts/thresholds, but free loadings (Wu & Estabrook, 2017's step one). Appropriate when all the indicators are ordlinal.}
#' \item{intercepts.scalar}{between a model with equal intercepts/thresholds, but free loadings and a full scalar model (Wu & Estabrook, 2017's step two). Appropriate when all the indicators are ordlinal.}
#' 
#' }
#' @return  Returns a list with detailed output on every available fit index, and a large matrix used for plotting with \code{\link[MIE]{plotDistances}}
#'
#'@export
incrementalFit <- function(..., level="metric") {
  
  
  if(level=="metric") {
  #fit.decrease <- abs(pairwiseFit(..., constraints = c("")) - pairwiseFit(..., constraints = c("loadings")))
    cat("Fitting configural models\n")
  configural = pairwiseFit(..., constraints = c(""))
    cat("\nFitting metric models\n")
  metric = pairwiseFit(..., constraints = c("loadings"))
  fit.decrease <- abs(configural - metric)
  detailed<-lapply(rownames(fit.decrease), 
                   function(f) 
                     cbind(configural=configural[f,], 
                           metric=metric[f,], 
                           fit.decrease=fit.decrease[f, ])  
                   )
  names(detailed)<-rownames(fit.decrease)
  
  } else  if(level=="scalar") {
      cat("Fitting metric models\n")
    metric = pairwiseFit(..., constraints = c("loadings"))
      cat("\nFitting scalar models\n")
    scalar = pairwiseFit(..., constraints = c("loadings", "intercepts", "thresholds"))
    fit.decrease <- abs(metric - scalar)
    detailed<-lapply(rownames(fit.decrease), function(f) cbind(metric=metric[f,], scalar=scalar[f,], fit.decrease=fit.decrease[f, ])  )
    names(detailed)<-rownames(fit.decrease)
    
  } else  if(level=="at.once") {
    cat("Fitting configural models\n")
    configural = pairwiseFit(..., constraints = c(""))
    cat("\nFitting scalar models\n")
    scalar = pairwiseFit(..., constraints = c("loadings", "intercepts", "thresholds"))
    fit.decrease <- abs(configural - scalar)
    detailed<-lapply(rownames(fit.decrease), function(f) cbind(configural=configural[f,], scalar=scalar[f,], fit.decrease=fit.decrease[f, ])  )
    names(detailed)<-rownames(fit.decrease)
  
    }   else  if(level=="intercepts.first") {
    cat("Fitting configural models\n")
      configural = pairwiseFit(..., constraints = c(""))
    cat("\nFitting equal-intercepts models\n")
    int.only = pairwiseFit(..., constraints = c("intercepts", "thresholds"))
    fit.decrease <- abs(configural - int.only)
    detailed<-lapply(rownames(fit.decrease), function(f) cbind(configural=configural[f,], int.only=int.only[f,], fit.decrease=fit.decrease[f, ])  )
    names(detailed)<-rownames(fit.decrease)
    
    }    else  if(level=="intercepts.scalar") {
      
      cat("Fitting equal-intercepts models\n")
      int.only = pairwiseFit(..., constraints = c("intercepts", "thresholds"))
      cat("\nFitting scalar models\n")
      scalar = pairwiseFit(..., constraints = c("loadings", "intercepts", "thresholds"))
      fit.decrease <- abs(int.only - scalar)
      detailed<-lapply(rownames(fit.decrease), function(f) cbind(int.only=int.only[f,], scalar=scalar[f,], fit.decrease=fit.decrease[f, ])  )
      names(detailed)<-rownames(fit.decrease)
      
    } 
  
  
 out <- list(detailed=detailed, bunch=fit.decrease)

 class(out)<-c("incrementalFit")
  
return(out)
 
}

#' Plot distances using computed measures of distance
#' 
#' @param measures Can be result of `MGCFAparameters`, `computeCovariance`, `computeCorrelation`, `incrementalFit`.
#' @param n.clusters Number of clusters.
#' @param fit.index Index to be used to in represneting measurement invariance distances. Only if the `measures` argument is an output of `incrementalFit`. Can be anything returned by `lavaan::fitMeasures`
#' @param drop Vector of group names to be dropped from the plot.
#' 
#' @return Computes distances and performs multidimensional scaling (two-dimensional projection). Returns ggplot-based plot. 
#' @seealso \code{\link{plotCutoff}}
#' @export
plotDistances <- function(measures, n.clusters = "auto", fit.index="cfi", drop = NULL, dist.method = NULL, shiny = FALSE) {
  
  pam1 = function(x, k){list(cluster = cluster::pam(x,k, diss = T, cluster.only=TRUE))}
  
  

    if("incrementalFit" %in% class(measures) | "average.dmacs" %in% fit.index) {
      
      if(fit.index=="average.dmacs") measures <- list(bunch=measures)

      dist1 <- reshape2::dcast(rbind( cbind(get_pairs(measures$bunch), measures$bunch[fit.index,]),
                            cbind(`colnames<-`(get_pairs(measures$bunch)[,2:1], c("V1", "V2")), measures$bunch[fit.index,])
      ), V2 ~ V1, value.var = "measures$bunch[fit.index, ]")
      row.names(dist1)<- dist1$V2
      dist1$V2 <- NULL
      dist1<-as.matrix(dist1)
      diag(dist1)<-rep(0, nrow(dist1))
      
      # remove dropped groups
      if(!is.null(drop)) dist1 <- dist1[!rownames(dist1) %in% drop, !(colnames(dist1) %in% drop) ]
      
      # new line below converting raw fit measures to distances with varying method
      dist2 <- stats::dist(dist1, method = ifelse(is.null(dist.method), "maximum", dist.method   ))
      
      if(n.clusters == "auto") {

      set.seed(1234)
      gskmn = cluster::clusGap(as.matrix(dist2), FUN=pam1, K.max = attr(dist2, "Size")-1, B = 50, verbose = F)
      n.clusters <- cluster::maxSE(f = gskmn$Tab[, "gap"], SE.f = gskmn$Tab[, "SE.sim"], method = "Tibs2001SEmax", SE.factor = 1) 
        if(shiny) {
          #updateSliderInput(session, "nclusters", max = groups - 1, value = n.clusters)
        } else {
          cat("\nOptimal number of clusters is ", n.clusters)
        }
      }
      
        
        #clusters <- stats::kmeans(dist1, n.clusters)$cluster
      clusters <- cluster::pam(dist2, diss =T, k = n.clusters)$clustering
      
      mds1 <- stats::cmdscale(dist2 * 1, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
      
      
    } else if( any(c("covariances", "correlations", "MGCFAparameters")  %in% class(measures)) ) {
       
      # remove dropped groups
      if(!is.null(drop)) measures <- measures[!rownames(measures) %in% drop, ]
      
      dist1 <- stats::dist(measures, 
                           method = ifelse(is.null(dist.method), "euclidean", dist.method ))
      mds1 <- stats::cmdscale(dist1 * 1, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
      
      if(n.clusters == "auto") {
        
        set.seed(1234)
        gskmn = cluster::clusGap(as.matrix(dist1), FUN=pam1, K.max = attr(dist1, "Size")-1, B = 50, verbose = F)
        n.clusters <- cluster::maxSE(f = gskmn$Tab[, "gap"], SE.f = gskmn$Tab[, "SE.sim"], method = "Tibs2001SEmax", SE.factor = 1) 
        
        if(shiny) {
          #shiny::updateSliderInput(session, "nclusters", value = n.clusters)
        } else {
          cat("\nOptimal number of clusters is ", n.clusters)
        }
      }
      
      #clusters <- stats::kmeans(measures, n.clusters)$cluster
      clusters <- cluster::pam(dist1, diss =T, k = n.clusters)$clustering
      
    } else if (fit.index == "signed.dmacs"){
      
      dmacs.signed.list <- attr(measures, "dmacs.signed.list")
      
      dmacs.per.parameter <- lapply(1:nrow(dmacs.signed.list[[1]]), function(x) {
        par.across.groups <- sapply(dmacs.signed.list, function(y) na.omit(y[x,]) )
      })
      
            
   distList <- lapply(1:length(dmacs.per.parameter), function(i) {
     
         dist1 <- reshape2::dcast(rbind( 
                                   cbind(get_pairs(measures), 
                                          x=dmacs.per.parameter[[i]]),
                                    cbind(`colnames<-`(get_pairs(measures)[,2:1], c("V1", "V2")), 
                                          x=dmacs.per.parameter[[i]])
                                   ), 
                                 V2 ~ V1, value.var = "x")
          row.names(dist1)<- dist1$V2
          dist1$V2 <- NULL
          dist1<-as.matrix(dist1)
          diag(dist1)<-rep(0, nrow(dist1))
  
        # remove dropped groups
         if(!is.null(drop)) dist1 <- dist1[!rownames(dist1) %in% drop, !(colnames(dist1) %in% drop) ]
        
        # new line below converting raw fit measures to distances with varying method
         dist2 <- stats::dist(dist1, method = ifelse(is.null(dist.method), "maximum", dist.method   ))
   })
         
   dist2 = Reduce("+", distList)/length(distList)
         
          
         if(n.clusters == "auto") {
           
           set.seed(1234)
           gskmn = cluster::clusGap(as.matrix(dist2), FUN=pam1, K.max = attr(dist2, "Size")-1, B = 50, verbose = F)
          n.clusters <- cluster::maxSE(f = gskmn$Tab[, "gap"], SE.f = gskmn$Tab[, "SE.sim"], method = "Tibs2001SEmax", SE.factor = 1)
          if(shiny) {
            #updateSliderInput(session, "nclusters", max = groups - 1, value = n.clusters)
          } else {
            cat("\nOptimal number of clusters is ", n.clusters)
          }
        }
      #   
      #   
      #   #clusters <- stats::kmeans(dist1, n.clusters)$cluster
        clusters <- cluster::pam(dist2, diss =T, k = n.clusters)$clustering
      #   
         mds1 <- stats::cmdscale(dist2 * 1, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
      #   
      #   
      #   
      #   
      # })
      
      #Reduce("+", dmacs.signed.list) / length(dmacs.signed.list)
      
    } else {
      
      stop("This function can accept only objects of classes 'incrementalFit', 'MGCFAparameters', 'covariances', and 'correlations', or fit.index='signed.dmacs' ")
      
    }
    
    colnames(mds1) <- c("dim1", "dim2")
    d <- as.data.frame(mds1)
    d$group <- rownames(d)
    d$cluster <- clusters
   
  
   # library(ggplot2); library(ggrepel)
    requireNamespace("ggplot2")
    requireNamespace("ggforce")
    requireNamespace("ggrepel")
    find_hull <- function(df) df[chull(df$dim1, df$dim2), ]
    hulls <- plyr::ddply(d, "cluster", find_hull)
    #rad.and.exp = 10*(max(unlist(hulls[, 1:2]))-min(unlist(hulls[, 1:2])))/nrow(d)
    
 g<-  ggplot(d, aes(dim1, dim2,  col=as.factor(cluster)))+
    geom_point( size=5, show.legend = F)+labs(x="", y="", col="")+
    geom_text_repel(aes(label=group), point.padding = unit(.3, "lines"), show.legend=F)+
    theme_minimal()+
    coord_fixed()+
    scale_colour_hue(l = 50, c = 120)+
    theme(panel.grid = element_blank(), axis.line=element_line(size=.5),axis.ticks=element_line(size=.5), plot.title=element_text(face="bold", size=18))+
   geom_shape(data = hulls, alpha = .4, linetype="blank", aes(fill = as.factor(cluster)), 
               expand = unit(10, "points"), radius = unit(10, "points"),
               show.legend = F)+
    labs(caption=paste("Distances represent measurement non-invariance"))
 

 
 
 
 # Add the circle
 if(any(class(measures)=="incrementalFit") & (fit.index == "cfi" | fit.index == "rmsea")) {
   
  
   r = ifelse(fit.index == "cfi", .005, .0075)
   xc = min(d$dim1) +  sqrt(r^2/2)
   yc = min(d$dim2) +   sqrt(r^2/2)
   
   circle.data <- data.frame(x = xc+r*cos(seq(0,2*pi,length.out=100)), 
                             y=  yc+r*sin(seq(0,2*pi,length.out=100)))
   
   g<-g+geom_path(data=circle.data, aes(x,y#, label=NULL
   ), col="lightgray", linetype="dashed")
   g<-g+labs(caption="The circle has diameter .01, meaning the increment \nin the fit index is within interval recommended by Chen")
   
 } else if ( "MGCFAparameters" %in% class(measures)) {
   
   fits<-attr(measures, "fit")[c("cfi", "rmsea", "srmr")]
   
   g<-g+labs(caption=paste(paste(c("CFI=", "RMSEA=", "SRMR="), sep=""),
                  paste(round(fits,3), sep=""), collapse=", "))
     }
   
   

g
 
 #invisible( d[,3:4])
 
}



#' Plots network based on pairwise fit indices reduced to Chen's cutoffs
#'  
#' @param measures The result of \code{\link[MIE]{incrementalFit}}
#' @param fit.index Index to be used to in representing measurement invariance distances. Only if the `measures` argument is an output of \code{\link[MIE]{incrementalFit}}. Should be "cfi", "rmsea", or "srmr", because only for these indices the cutoffs were suggested by Chen (2007).
#' @param drop Vector of group names to be dropped from the plot.
#' @param weighted Logical. If weighted graph should be created. See \code{\link[igraph]{graph_from_adjacency_matrix}} for details.
#' 
#' @details The function extracts a given fit indeces from pairwise fitted MGCFAs, and uses cutoff of .01 to identify edges between groups (nodes), so that the groups for whom  invariance is supported, are connected on the plot. The results are plotted using \code{\link[igraph]{cluster_label_prop}}.
#' @seealso \code{\link[MIE]{plotDistances}}
#' @export         
plotCutoff <- function(measures, fit.index = "cfi", cutoff = NULL, weighted = TRUE, drop = NULL, shiny = F) {
  
  
  if(any(!fit.index %in% c("cfi", "rmsea", "srmr"))) {
    if (shiny == T) {
      showNotification("Cutoffs for this fit measure are not available. Using convenient .01 (unrealiable!).\n Consider switching off 'Use cutoffs' option.", type = "warning", duration = NULL, id = "nocutoffs")  
    } else {
      warning("Cutoffs for this fit measure are not available. Using convenient .01 (unrealiable!).\n Interpret cautiously.")
     }
    
    } 
  
  # remove dropped groups
  #if(!is.null(drop)) measures <- measures[!rownames(measures) %in% drop, ]
  
  #abs.thrshld <- switch(fit.index, cfi=function(x) `>`(x, .90), rmsea = function(x) `<`(x, .05))
  if(is.null(cutoff)) {
    #Chen's
   # thrshld <- switch(fit.index, 
   #                       cfi   = function(x) `<`(x, .01),
   #                       rmsea = function(x) `<`(x, .01),
   #                       srmr  = function(x) `<`(x, .01),
   #                       nnfi  = function(x) `<`(x, .01) 
   # 
   #                   )
    
    thrshld <- function(x) `<`(x, .01)
  
  } else {
    thrshld <- function(x) `<`(x, cutoff)
  }
  
  # mtrx <- measures$detailed[[fit.index]]
  # mtrx.df <- data.frame(
  #   i = gsub("^.*_", "",  rownames(mtrx)), 
  #   j = gsub("_.*$", "",  rownames(mtrx)),
  #   # absolute
  #   #tie = abs.thrshld(unname(mtrx[,2]))
  #   
  #   # incremental
  #   tie = thrshld(unname(mtrx[,"fit.decrease"]))
  # )
  
  dist1 <- reshape2::dcast(rbind( cbind(get_pairs(measures$bunch), measures$bunch[fit.index,]),
                                  cbind(`colnames<-`(get_pairs(measures$bunch)[,2:1], c("V1", "V2")), measures$bunch[fit.index,])
  ), V2 ~ V1, value.var = "measures$bunch[fit.index, ]")
  row.names(dist1)<- dist1$V2
  dist1$V2 <- NULL
  dist1<-as.matrix(dist1)
  # remove dropped groups
  if(!is.null(drop)) dist1 <- dist1[!rownames(dist1) %in% drop, !(colnames(dist1) %in% drop) ]
  
  thrshld.dist1 = 
    
  dist2<-thrshld(dist1)*1 
  dist1[!thrshld(dist1)] <- 0 
  #diag(dist1)<-rep(0, nrow(dist1))
 
  if(weighted) {
  
    dist3 <- dist1
    diag(dist1)<-0
    dist1[!thrshld(dist1)] <-0
    dist1[ thrshld(dist1)] <- 1/(dist1[thrshld(dist1)]+1)

    diag(dist2)<-0
    dist1[dist2==1]<-1/(dist1[dist2==1]*100+1)
    dist1[dist2==0]<-0
    # 
    
    net <- igraph::graph_from_adjacency_matrix(dist1, diag = F, mode = "lower", weighted = TRUE)
  } else {
    
    dist1 <- dist2
    net <- igraph::graph_from_adjacency_matrix(dist2, diag = F, mode = "lower", weighted = NULL)
  }
  

  
  
  
  clp <- igraph::cluster_label_prop(net)
  # set.seed(123)
  # igraph:::plot.communities(clp, net, edge.color = "darkgrey", layout = layout_with_fr)
  set.seed(123)
   coords <- layout_with_fr(net, dim = 2, niter = 500)
   colnames(coords)<- c("dim1", "dim2")
   coords<- as.data.frame.matrix(coords)
   coords$group <- rownames(dist1)
   clusters <- Reduce("rbind",  lapply(1:length(clp), function(x) data.frame( 
     group = clp[[x]], 
     cluster = rep(x, length(clp[[x]])), stringsAsFactors = F)))
    coords <- merge(coords, clusters, by = "group")
    find_hull <- function(df) df[chull(df$dim1, df$dim2), ]
    hulls <- plyr::ddply(coords, "cluster", find_hull)
    
    d1 <- reshape2::melt(dist1)
    d1 <- d1[d1$value!=0,]
    d2 <- merge(d1, coords, by.x = "Var1", by.y = "group", all.x = T)
    d3 <- merge(d2, coords, by.x = "Var2", by.y = "group", all.x = T)
    
    requireNamespace("ggplot2")
    requireNamespace("ggforce")

    g<-  ggplot(coords, aes(dim1, dim2,  col=as.factor(cluster)))+
      geom_segment(data = d3, aes(x = dim1.x, xend = dim1.y, y = dim2.x, yend = dim2.y), col = "black", alpha = .5)+
      #geom_text_repel(aes(label=group), point.padding = unit(.3, "lines"), show.legend=F)+
      
      geom_shape(data = hulls, aes(fill = as.factor(cluster) ),
                 alpha = ifelse(length(unique(coords$cluster))>1, .4, 0), linetype="blank",
                 expand = unit(10, "points"), radius = unit(10, "points"),
                 show.legend = F)+
      geom_point( size=5, show.legend = F)+labs(x="", y="", col="")+
      geom_label(aes(label=group), show.legend=F, alpha = 1)+
      theme_void()+
      coord_fixed()+
      scale_colour_hue(l = 50, c = 120)+
      theme(panel.grid = element_blank(), axis.line=element_line(size=.5),axis.ticks=element_line(size=.5), plot.title=element_text(face="bold", size=18))+
      labs(caption=paste("Lines represent measurement invariance"))
    
g
  # clp <- igraph::cluster_label_prop(igraph::graph_from_edgelist(as.matrix(mtrx.df[mtrx.df$tie,-3]), directed = F))
  # net <- igraph::graph_from_edgelist(as.matrix(mtrx.df[mtrx.df$tie,-3]), directed = F)
  # igraph:::plot.communities(clp, net)
  
}


# igraph:::plot.communities(clp, net, edge.color = "darkgrey", layout = layout_with_fr)
# p <- recordPlot()
# plot.new() ## clean up device
# p


# Effects and DMACs #####


#'Dmacs 
#' @description  Simply a wrapper arounnd dmacs::lavaan_dmacs() function
#' @param lav.model fitted lavaan model
dmacs <- function(lav.model) {
  as.matrix(
  dmacs::lavaan_dmacs(configural, 
                      RefGroup=lavaan::lavInspect(configural, "group.label")[1], 
                      "pooled" )$DMACS
  )
}


dmacs.unsigned <- function(LambdaR, ThreshR, LambdaF, ThreshF, MeanF, VarF, SD) {
  
  expected_value <- function(Lambda, Thresh, Theta) { Thresh + Lambda * Theta }
  z <- seq(-5, 5, .001)
  
  Y_j1 = expected_value(LambdaF, ThreshF, MeanF + z * sqrt(VarF))
  Y_j2 = expected_value(LambdaR, ThreshR, MeanF + z * sqrt(VarF))
  
  sqrt(sum((Y_j1- Y_j2)^2 * dnorm(z) * .001 * sqrt(VarF)))/SD
}

dmacs.signed <- function(LambdaR, ThreshR, LambdaF, ThreshF, MeanF, VarF, SD) {
  
  expected_value <- function(Lambda, Thresh, Theta) { Thresh + Lambda * Theta }
  z <- seq(-5, 5, .001)
  
  Y_j1 = expected_value(LambdaF, ThreshF, MeanF + z * sqrt(VarF))
  Y_j2 = expected_value(LambdaR, ThreshR, MeanF + z * sqrt(VarF))
  
  sum((Y_j1- Y_j2) * dnorm(z) * .001 * sqrt(VarF))/SD
}





dmacs_lavaan <- function(fit, signed = F, pooled.sd = T, RefGroup = 1) {
  
  if(class(fit)!='lavaan') warning("'fit' argument should be lavaan fitted object!")
  est <- lapply(lavaan::lavInspect(fit, "est"), function(x) x[ c("lambda","alpha", "psi", "nu" )])

  Groups <- names(lavaan::lavInspect(fit, "est"))
   #which(Groups == RefGroup)
  
  if(pooled.sd) {
    
    data.ref.group = lavaan::lavInspect(fit, "data")[[RefGroup]]
    
    refsd <- apply(data.ref.group, 2, sd, na.rm = TRUE)
    refn <- colSums(!is.na(data.ref.group))
    focsd <- apply(data.ref.group, 2, sd, na.rm = TRUE)
    focn <- colSums(!is.na(data.ref.group))
    
    SDs <-  ((focn - 1) * focsd + (refn - 1) * refsd) / ((focn -  1) + (refn - 1))
    
  } else {
    
    SDs <- apply(lavaan::lavInspect(fit, "data")[[-RefGroup]], 2, sd,  na.rm = TRUE)
    
  }
  
  DMACS <- mapply(function(...) if (signed) dmacs.signed(...) else dmacs.unsigned(...),
                  est[[RefGroup]]$lambda,
                  est[[RefGroup]]$nu,
                  est[[-RefGroup]]$lambda,
                  est[[-RefGroup]]$nu,
                  est[[-RefGroup]]$alpha,
                  diag(est[[-RefGroup]]$psi),
                  SDs)
  
  DMACS= matrix(DMACS, nrow= dim(est[[RefGroup]]$lambda)[1], ncol =  dim(est[[RefGroup]]$lambda)[2])
  rownames(DMACS) <- rownames(est[[RefGroup]]$lambda)
  colnames(DMACS) <- colnames(est[[RefGroup]]$lambda)
  DMACS[est[[RefGroup]]$lambda==0] <- NA
  DMACS
  
}




#' Pairwise DMACS
#' @param ... arguments passed to pairwiseFit()
#'
#' @export
pairwiseDmacs <- function(..., signed = T) {
  
  
  cat("Fitting configural models\n")
  configural = pairwiseFit(..., constraints = c(""))
  fit.decrease <- configural
  detailed<-configural
  
}

