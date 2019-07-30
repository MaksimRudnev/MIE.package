## Global functions #####

## ..return all possible pairs of countries without duplicates based on variable #####
#' Produce all possible pairs of groups
#' 
#' @param variable Grouping variable to find the set of all possible unique pairs of values.
#' 
#' @export
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
# @example # pairwiseFit(model="F =~ impfree + iphlppl + ipsuces + ipstrgv", data = ess6, group = "cntry", "loadings")
#' 
#' @return  The function returns matrix of fit indices for multiple group CFA models fitted to each possible pair of groups.
#' 
#' @export
pairwiseFit <- function(model,
                        data, 
                        group = "cntry",
                        constraints=c(""), 
                        pairs.of.groups=NULL, 
                        message="Fitting pairwise lavaan models",
                        shiny=FALSE
                        #, extra.options
                        
) {
  
  if(is.null(pairs.of.groups)) pairs.of.groups <- pairs_of_groups(data[[group]])
  data=data[,c(group, colnames(data)[colnames(data)!=group])]
  
  colnames(data)[1]<-"cntry"
  
  runPairwiseModels <- function() {
    model.lav<- lavaan::cfa(model, data=data[data$cntry==pairs.of.groups[1,1] | 
                                       data$cntry==pairs.of.groups[2,2],],
                    group="cntry", group.equal=constraints#, extra.options
    )
    
    if(lavaan::lavInspect(model.lav, "converged")) mod<- lavaan::fitmeasures(model.lav) else mod <- rep(999, 41)
    
    # 41 is a number of fit indices currently provided by lavaan  
    mod<-matrix( c(mod, rep(0,41*(nrow(pairs.of.groups)-1))), nrow=41, dimnames=list(names(mod), NULL))
    
    #Non-positive definite?
    non.positive <- FALSE
    
    # Uses 'for' instead of 'sapply' in order to show the progress bar
    for(x in 2:nrow(pairs.of.groups)) {
      model.lav<- lavaan::cfa(model, data=data[data$cntry==pairs.of.groups[x,1] | 
                                         data$cntry==pairs.of.groups[x,2],],
                      group="cntry", group.equal=constraints#, extra.options
      )
      
      #If converged, record fitmeasure; if not converged add a missing sign 999.
      if(lavaan::lavInspect(model.lav, "converged")) mod[,x]<- lavaan::fitmeasures(model.lav) else mod[,x]<- rep(999, 41)
      
      #Save non-positive definite status
      
      non.positive[[x]]<-lavaan::lavInspect(model.lav, "post.check")==FALSE
      
      
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
    return(mod)
    
  }
  
    if(shiny) {
  withProgress(message = message, value = 0, {   #Create progress bar
    mod <- runPairwiseModels()
    mod
    print(mod)
  }) #close progress bar 
    } else {
      pb <- utils::txtProgressBar(title="Fitting pairwise models", style=3)
      mod <- runPairwiseModels()
      colnames(mod)<-apply(pairs.of.groups, 1, paste, collapse="_")
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
        
        print("lavCALL"); print(mod@call[-3])
        
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
        
        print(paste("Fit model with", paste(mod@Data@group.label, collapse=",")))
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
        class(mod)<- "MGCFAparameters"
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
  message("Computing covariance matrix")
  
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
  message("Computing correlation matrix")
  
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
  class(tab) <- "correlations"
  tab
  
  #New version
  #sapply(unique(dt$dat$cntry), function(x) 
  #  var(dt$dat[dt$dat$cntry==x,-1], use="complete.obs")[lower.tri(var(dt$dat[dt$dat$cntry==x,-1], use="complete.obs"),diag = F)]
  #)
}

#' Test for measureent invariance in three steps
#' @description Computes configural, metric, and scala invariance models and compares them.
#' 
#' @param ... Formula, group, and all the other eligible arguments of `lavaan::cfa` function.
#'
#' @export
globalMI <- function(...) {
  
  r.conf<-lavaan::cfa(...)
  r.metric<-lavaan::cfa(..., group.equal = "loadings")
  r.scalar<-lavaan::cfa(..., group.equal = c("loadings", "intercepts"))
  
  #print(r.conf@call[-3])
  
  out <- lavaan::lavTestLRT(r.conf, r.metric, r.scalar, model.names = c("Configural", "Metric", "Scalar"))
  #out <- out[,names(out)[c(4, 1, 5, 6)]]
  out2<-t(sapply(list(r.conf, r.metric, r.scalar),   fitmeasures)[c("cfi", "tli", "rmsea", "srmr"),])
  out2.2 <- apply(out2, 2, function(x) (c(NA, x[2]-x[1], x[3]-x[2]  )))
  colnames(out2.2)<- paste("diff.", toupper(colnames(out2.2)), sep="")
  out2.3 <- t(Reduce("rbind", lapply(1:ncol(out2), function(x) rbind(out2[,x], out2.2[,x]))))
  colnames(out2.3)<-toupper(as.vector(sapply(1:length(colnames(out2)), function(i) c(colnames(out2)[i], colnames(out2.2)[i]) )))
  rownames(out2.3)<-c("Configural", "Metric", "Scalar")
  print(out)
  cat("\n")
  print(round(out2.3, 3), digits=3, row.names = TRUE, na.print = "" )
}


#' Run Shiny MIE application
#' 
#' @param data Data.frame containing only the variables to be used in calculations. May contain only numeric variables. If `group` argument is not specified, the group variable should be the first in the data.
#' @param model Character, lavaan formula.
#' @param group Character, name of the variable in data.
#' 
#' @export
runMIE <- function(model = NULL, data = NULL,  group = NULL) {

  .GlobalEnv$.data <- data
  .GlobalEnv$.model <- model
  .GlobalEnv$.group <- group
  on.exit(rm(list=c(".data", ".model", ".group"), envir=.GlobalEnv))                       
  shiny::runApp(appDir = system.file("application", package = "MIE"))
  #shiny::runApp(appDir = "inst/application")
}

# foo <- function(x) {   
#   e <- environment() # current environment
#   p <- parent.env(e)
#   cat(paste(ls(envir=parent.frame(n=1)), collapse="\n"))
#   
#   }
# 
# runMIE(.data = ess_trimmed)

# foo <- function(.data) { obj.name = deparse(substitute(.data));  get(obj.name, envir = .GlobalEnv)  }
# foo(.data = ess_trimmed)

#' Run pairwise models and compute decrease in fit
#'
#' @param ... The arguments passed to `pairwiseFit``
#' @param level Character. If "metric" (default) the decrease of fit between configural and metric model is computed, if "scalar", metric is compared to scalar model fit.
#' @return  Returns a list with detailed output on every available fit index, and a large matrix used for plotting with `plotDistances`
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
  otpt<-lapply(rownames(fit.decrease), function(f) cbind(configural=configural[f,], metric=metric[f,], fit.decrease=fit.decrease[f, ])  )
  names(otpt)<-rownames(fit.decrease)
  
  } else  if(level=="scalar") {
      cat("Fitting metric models\n")
    metric = pairwiseFit(..., constraints = c("loadings"))
      cat("\nFitting scalar models\n")
    scalar = pairwiseFit(..., constraints = c("loadings", "intercepts"))
    fit.decrease <- abs(metric - scalar)
    otpt<-lapply(rownames(fit.decrease), function(f) cbind(metric=metric[f,], scalar=scalar[f,], fit.decrease=fit.decrease[f, ])  )
    names(otpt)<-rownames(fit.decrease)
  }
  
  
 out <- list(detailed=otpt, bunch=fit.decrease)
 
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
#' @export
plotDistances <- function(measures, n.clusters = 4, fit.index="cfi", drop = NULL) {
  if(class(measures) == "covariances" | class(measures) == "correlations") measures <- t(measures)
    
    if(class(measures)=="incrementalFit") {

      dist1 <- reshape2::dcast(rbind( cbind(get_pairs(measures$bunch), measures$bunch[fit.index,]),
                            cbind(`colnames<-`(get_pairs(measures$bunch)[,2:1], c("V1", "V2")), measures$bunch[fit.index,])
      ), V2 ~ V1, value.var = "measures$bunch[fit.index, ]")
      row.names(dist1)<- dist1$V2
      dist1$V2 <- NULL
      dist1<-as.matrix(dist1)
      diag(dist1)<-rep(0, nrow(dist1)) 
      mds1 <- stats::cmdscale(dist1 * 1, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
      clusters <- stats::kmeans(dist1, n.clusters)$cluster
      
    } else {
      
      mds1 <- stats::cmdscale(dist(measures) * 1, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
      clusters <- stats::kmeans(measures, n.clusters)$cluster
      
    }
    
    colnames(mds1) <- c("dim1", "dim2")
    d <- as.data.frame(mds1)
    d$group <- rownames(d)
    d$cluster <- clusters
    if(!is.null(drop)) d <- d[!(d$group %in% drop), !(colnames(d) %in% drop) ]
  
   # library(ggplot2); library(ggrepel)
    requireNamespace(ggplot2)
 g<-  ggplot(d, aes(dim1, dim2,  col=as.factor(cluster)))+
    geom_point( size=5, show.legend = F)+labs(x="", y="", col="")+
    geom_text_repel(aes(label=group), point.padding = unit(.3, "lines"), show.legend=F)+
    theme_minimal()+
    coord_fixed()+
    scale_colour_hue(l = 50, c = 120)+
    theme(panel.grid = element_blank(), axis.line=element_line(size=.5),axis.ticks=element_line(size=.5), plot.title=element_text(face="bold", size=18))+
    labs(caption=paste("Distances represent measurement non-invariance"))
 
 
 # Add the circle
 if(class(measures)=="incrementalFit" & (fit.index == "cfi" | fit.index == "rmsea")) {
   
  
   r = ifelse(fit.index == "cfi", .005, .0075)
   xc = min(d$dim1) +  sqrt(r^2/2)
   yc = min(d$dim2) +   sqrt(r^2/2)
   
   circle.data <- data.frame(x = xc+r*cos(seq(0,2*pi,length.out=100)), 
                             y=  yc+r*sin(seq(0,2*pi,length.out=100)))
   
   g<-g+geom_path(data=circle.data, aes(x,y#, label=NULL
   ), col="lightgray", linetype="dashed")
   g<-g+labs(caption="The circle has diameter .01, meaning the increment \nin the fit index is within interval recommended by Chen")
   
 } else if (class(measures)=="MGCFAparameters") {
   
   fits<-attr(measures, "fit")[c("cfi", "rmsea", "srmr")]
   
   g<-g+labs(caption=paste(paste(c("CFI=", "RMSEA=", "SRMR="), sep=""),
                  paste(round(fits,3), sep=""), collapse=", "))
     }
   
   

 g
 
 
 
}
  