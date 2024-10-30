
# ADD IF conditions to pairwise, and add extra.options to pairwise.fit function


#setwd("~/Dropbox/STAT/European Social Survey/ShinyValues/MeasurementInvarianceExplorer")


# TO develop: errors handlers; 
# Speed of pairwise (think of what can be done with non-pairwise)

#dat<-simMGCFA

require("shiny", quietly = T)
require("lavaan", quietly = T)
require("magrittr", quietly = T)
require("reshape2", quietly = T)
require("markdown", quietly = T)
require("DT", quietly = T)
require("ggplot2", quietly = T)
require("ggrepel", quietly = T)
require("dplyr", quietly = T)
options(shiny.maxRequestSize=100*1024^2) 


#source("server.functions.R")


# Begin the SHINY code #####
shinyServer(function(session, input, output) {
 

### Values ####
  dt<-reactiveValues ( 
    dat = NULL,
    model = NULL,
    extra.options = NULL
  )  
  
  temp <- reactiveValues (
    old.model.configural.MGCFA = NULL,
    old.extra.options.configural.MGCFA = NULL,
    old.model.metric.MGCFA = NULL,
    old.extra.options.metric.MGCFA = NULL
  )
  

  vals <- reactiveValues( 
    
    keeprows = NULL,
    excluded = NULL
  )
  
  modelStorage <- reactiveValues(
    
    covariance=NULL, #lower triangle of covariance matrix of all but first columns in the datafile
    correlation=NULL, # same
    loadings=NULL, #loadings from MGCFA configural model
    intercepts=NULL, #intercepts from MGCFA metric model
    
    conf=NULL,   #
    metric=NULL,
    scalar=NULL
    
  )
  
  Errors <- reactiveValues(
    nomodel = FALSE,
    nodata = FALSE,
    current = NULL
  )
  
  
  # Reading in the data from globalEnv ####
# reactive({ 
  if(!is.null(.data) & !is.null(.group) ) { 
    local.data <- .data[, c(colnames(.data)[colnames(.data) == .group],
                       colnames(.data)[colnames(.data) != .group])]
    print(colnames(local.data))
    colnames(local.data)[1] <- "cntry"
    print("changed the data order")
  } else if (!is.null(.data) & is.null(.group) ) {
    local.data <- .data
    
  }
  
  
  if(!is.null(.model)) { 
    dt$model = .model
    updateCheckboxInput(session, "use.formula", value=TRUE)
  }
  
  
  if(!is.null(.data)) {
     dt$dat = local.data
  #print("unique(dt$dat$cntry)"); print(unique(dt$dat$cntry))
  vals$keeprows = unique(isolate(dt$dat$cntry))
  vals$excluded <- NULL
  modelStorage$covariance <- t(computeCovariance(isolate(dt$dat), group = "cntry"))
  
  showNotification("Using data from the R object.", type="message", duration=10)
  }
# }) 
  
  
  
  # .. Button for using simulated data #####
  observeEvent(input$useSimulated, {
    dt$dat = read.csv("simulated2.csv")
    dt$model = "#By default a simulated data - model mimics Schwartz values ESS scale
  person=~ ipcrtiv +impfree +impfun +ipgdtim +impdiff +ipadvnt+ imprich +iprspot +ipshabt +ipsuces;
    social=~ impenv +ipeqopt +ipudrst +iplylfr +iphlppl +impsafe +ipstrgv +ipfrule +ipbhprp +ipmodst +imptrad;"
    updateCheckboxInput(session, "use.formula", value=TRUE)
    #print("unique(dt$dat$cntry)"); print(unique(dt$dat$cntry))
    vals$keeprows = unique(dt$dat$cntry)
    vals$excluded <- NULL
    print(paste("Button 'play with fake data' has been used."))
    
    modelStorage$covariance <- computeCovariance(dt$dat, group = "cntry")
    
    showNotification("Using fake data for testing the tool.", type="warning", duration=10)
  })
  
##Event input new data file #### 
  #It resets the settings and computes covariance
  
  observeEvent(input$file1, {
 
    #Force the measure selection to covariance 
    #input$measure<-"covariance"
    
    message("New file was uploaded.")
    
    updateRadioButtons(session, "measure", selected="covariance")
    print("input$measure"); print(input$measure)
    #Set all model results (including covariance) to NULL
    for(x in names(modelStorage)) modelStorage[[x]]<-NULL
    
    #Read the data
    showNotification("Reading data...")
    d<-read.csv(input$file1$datapath, header = T)
    names(d)[1]<-c("cntry")
    d$cntry<-as.factor(d$cntry) # Hmm...
    dt$dat <- d
    print(paste(colnames(d), collapse="; "))
    rm(d)
    #Set subsets to NULL
        vals$keeprows <- unique(dt$dat$cntry)
        vals$excluded = NULL
        
    #Set all model results (beside covariance) to NULL
        #for(x in names(modelStorage)[-1]) modelStorage[[x]]<-NULL

        ##modelStorage <- lapply(names(modelStorage), function(x) modelStorage[[x]]<-NULL)
    
        

   #     showNotification("Computing covarances for a new data...", type="warning", duration=3)    
    
    #Compute covariance matrix   
    #Split dataset and compute variance-covariance for each group separately
        modelStorage$covariance<-computeCovariance(isolate(dt$dat), group = "cntry")
        
        #print( str(modelStorage$covariance))
       # print(head(dt$dat))
  })

 
  
  
#Subset the data #####  
  selectedData <- reactive({
    
    
    print(paste(" selectedData() subset the raw data for", paste(vals$keeprows, collapse=",")))
      dat<-dt$dat

      dat<-dat[dat$cntry %in% vals$keeprows, ] 
      if(is.factor(dat$cntry)) dat$cntry <- droplevels(dat$cntry) 
   
    dat
    
    })


  
  

  
#Subset measures that saved to modelStorage. Reactive subsettngMatrices() ####
subsettingMatrices <- reactive ({
  
  ###### Covariance ----------------------------------------------------------
  
  if(input$measure =="covariance") {   # This reacts to changes in modelStorage$covariance AND vals$keeprows and drops the excluded columns
    #message("Subsetting covariance matrix")
    #message(paste("Subsetting for vals$keeprows = ", paste(vals$keeprows, collapse=", "), ".", sep=""))
    #print(head(modelStorage$covariance))
    tab<-modelStorage$covariance[, vals$keeprows]
    
    # Compute distances between groups based on a subset of covariances
    dist<-dist(t(tab))
    
    # Formats the table to show
    additional<-data.frame(group=rownames(t(tab)), round(t(tab), 3))
    
    # Export
    list(dist=dist, additional=additional)
    
   ###### Correlation ----------------------------------------------------------
    
  } else if ( input$measure =="correlation") {
    
    
    modelStorage$correlation<-t(computeCorrelation(isolate(dt$dat), group="cntry"))[, vals$keeprows]
    
   #print("modelStorage$correlation"); print(class(modelStorage$correlation))
    
    tab<-modelStorage$correlation
    
    #print("tab"); print(tab)
    dist<-dist(t(tab))
    additional<-data.frame(group=rownames(t(tab)), round(t(tab), 3))
    list(dist=dist, additional=additional)
    
    
    
    ##### Fitting configural vs metric pairwise ---------------------------------
    

  } else if (input$measure=="fitincrement.metric") { 
    
    # This takes saved model results from modelStorage 
    #  and subsets it to fit the list of included countries and chosen fit measure
    # or, if there are extra pairs not computed yet, computes it
    
    # Fitting configural pairwise models  
    if(is.null(modelStorage$conf) |
       sum(vals$keeprows %in% unique(unlist(attr(modelStorage$conf, "pairs.of.groups")))==FALSE) >0 |
       ifelse( !is.null(dt$model),
               ifelse(!is.null(attr(modelStorage$conf, "model.formula")),
                      dt$model!=attr(modelStorage$conf, "model.formula"), TRUE), FALSE))     {
      #Logging
      print("There are some uncomputed pairs of configural models, or model has changed, or compute for the first time.")
      print("Recorded model is:"); print(attr(modelStorage$conf, "model.formula"));
      print("User updated model is:"); print(dt$model)
      
      if(!is.null(modelStorage$conf) & !ifelse( !is.null(dt$model),
                                                ifelse(!is.null(attr(modelStorage$conf, "model.formula")),
                                                       dt$model!=attr(modelStorage$conf, "model.formula"), TRUE), FALSE)) {
        extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$conf, "pairs.of.groups")))]
        pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$conf, "pairs.of.groups"))), stringsAsFactors = F)
        names(pairs.c)<-names(attr(modelStorage$conf, "pairs.of.groups"))
        print("Configural models were already computed, but there are some new pairs to compute.")
        
        
      } else {
        print("No configural model in storage, or new formula, so compute the whole thing for all selected groups")
        pairs.c <- MIE:::pairs_of_groups(as.character(vals$keeprows))
      }
      
      print("pairs.c"); print(pairs.c)
      # Compute lacking pairs of conf models
      conf.pairwise<- pairwiseFit(model = dt$model,
                                  data  = dt$dat, 
                                  group = "cntry",
                                  constraints = c(""),
                                  pairs.of.groups = pairs.c, 
                                  message = 'Fitting pairwise configural models by lavaan',
                                  shiny = TRUE#,
                                  #extra.options = dt$extra.options
                                  )
      
      
      # Merge with previous fits (if any)
      temp<- cbind(modelStorage$conf, conf.pairwise)
      attr(temp, "pairs.of.groups")<- rbind(attr(modelStorage$conf, "pairs.of.groups"),
                                               attr(conf.pairwise, "pairs.of.groups"))
      attr(temp, "model.formula") <- dt$model
      
      modelStorage$conf <- temp
      
    }
    
    # Fitting metric pairwise models
    
    
    if(is.null(modelStorage$metric) |
       sum(vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))==FALSE) >0|
       ifelse( !is.null(dt$model),
               ifelse(!is.null(attr(modelStorage$metric, "model.formula")),
                      dt$model!=attr(modelStorage$metric, "model.formula"), TRUE), FALSE) )  {
      
      #Logging
      print("There are some uncomputed pairs of metric models, or model has changed, or compute for the first time.")
      print("Recorded model is:"); print(attr(modelStorage$metric, "model.formula"));
      print("User updated model is:"); print(dt$model)
      
      if(!is.null(modelStorage$metric) & !ifelse( !is.null(dt$model),
                                                  ifelse(!is.null(attr(modelStorage$metric, "model.formula")),
                                                         dt$model!=attr(modelStorage$metric, "model.formula"), TRUE), FALSE)) {
        print("Making a subset of existing metric models")
        
        extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))]
        pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$metric, "pairs.of.groups"))), stringsAsFactors = F)
        names(pairs.c)<-names(attr(modelStorage$metric, "pairs.of.groups"))
        
      } else {
        
        print("None of metric models were computed")
        pairs.c<-MIE:::pairs_of_groups(as.character(vals$keeprows))
      }
      
      ##Compute lacking pairs of metric models
      metric.additional<- pairwiseFit(model = dt$model, 
                                      data = dt$dat, 
                                      group="cntry",
                                      constraints = c("loadings"),
                                      pairs.of.groups = pairs.c, 
                                      message='Fitting pairwise metric models by lavaan',
                                      shiny=TRUE#,
                                      #extra.options = dt$extra.options
                                      )
      
      
      # Export 
      temp<- cbind(modelStorage$metric, metric.additional)
      attr(temp, "pairs.of.groups")<- rbind(attr(modelStorage$metric, "pairs.of.groups"),
                                               attr(metric.additional, "pairs.of.groups"))
      attr(temp, "model.formula") <- dt$model
      modelStorage$metric <- temp
      
    }
  

    
# Formatting fit indices for export
   
   list.included.pairs.conf   <- modelStorage$conf   %>% get_pairs %>% apply(2, is_in, vals$keeprows) %>% rowSums(.)==2
   list.included.pairs.metric <- modelStorage$metric %>% get_pairs %>% apply(2, is_in, vals$keeprows) %>% rowSums(.)==2
   
   # if(sum(list.included.pairs.conf %in% list.included.pairs.metric==FALSE)>0) { 
   #   showNotification("Something went wrong... Switching to covariance view", type="error")
   #   updateRadioButtons(session, "measure", selected="covariance")
   #   }
   
   # Subset tables of fit indices
   conf_subset<-   modelStorage$conf   %>%  extract(,list.included.pairs.conf) %>% extract(rownames(.)  == input$fitincrement.chosen,)
   metric_subset<- modelStorage$metric %>%  extract(,list.included.pairs.conf) %>% extract(rownames(.)  == input$fitincrement.chosen,)
   fit.decrease <- abs(conf_subset - metric_subset)
   
   # Get a subset of group pairs names
   pair.names.subset <- get_pairs(modelStorage$conf)[list.included.pairs.conf,]
   

   # Create a distance matrix for MDS   
   dist<- sapply(as.character(vals$keeprows), function(colname) sapply(as.character(vals$keeprows), function(rowname) {
     
     fit.decrease[pair.names.subset[,1]==colname & pair.names.subset[,2]==rowname |
                   pair.names.subset[,1]==rowname & pair.names.subset[,2]==colname]
     
   } )) %>% inset(., sapply(., length)==0, 0)
   
   dist <- matrix(unlist(dist), nrow=dim(dist)[1], ncol=dim(dist)[2], dimnames = dimnames(dist))
   
  
   # Export 
   additional<-data.frame( "Group 1"=pair.names.subset[,1],
                            "Group 2"=pair.names.subset[,2],
                            configural=round(conf_subset,3) ,
                            metric=round(metric_subset,3),
                            difference=round(fit.decrease, 3))
    
    list(dist=dist, additional=additional)
    
 
 ######### Fitting scalar vs metric pairwise -----------------------------------------------------    
    
    
  } else if(input$measure=="fitincrement.scalar") {
  
    message("Fitting increment scalar/metric...")
    
    # Fitting metric models
    if(is.null(modelStorage$metric) | sum(vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))==FALSE) >0 ) 
    {
      
      if(!is.null(modelStorage$metric)) {
        extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))]
        
        pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$metric, "pairs.of.groups"))), stringsAsFactors = F)
        names(pairs.c)<-names(attr(modelStorage$metric, "pairs.of.groups"))
      } else {
        pairs.c <- MIE:::pairs_of_groups(as.character(vals$keeprows))
      }
      
      print("pairs.c"); print(str(pairs.c))
      ##Compute lacking pairs of metric models
      metric.pairwise<- pairwiseFit(model = dt$model, 
                                    data = dt$dat, 
                                    group="cntry",
                                    constraints = c("loadings"), 
                                    pairs.of.groups = pairs.c, 
                                    message='Fitting extra pairwise scalar models by lavaan',
                                    shiny=TRUE#,
                                    #extra.options = dt$extra.options
                                    )
      
      temp<- cbind(modelStorage$metric, metric.pairwise)
      attr(temp, "pairs.of.groups")<- rbind(attr(modelStorage$metric, "pairs.of.groups"),
                                               attr(metric.pairwise,      "pairs.of.groups"))
      modelStorage$metric <- temp
      
    }
    
    
    # Fitting scalar models
    if(is.null(modelStorage$scalar) | sum(vals$keeprows %in% unique(unlist(attr(modelStorage$scalar, "pairs.of.groups")))==FALSE) >0 ) 
    {
      if(!is.null(modelStorage$scalar)) {
      extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$scalar, "pairs.of.groups")))]
      
      pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$scalar, "pairs.of.groups"))), stringsAsFactors = F)
      names(pairs.c)<-names(attr(modelStorage$scalar, "pairs.of.groups"))
      
      } else {
        pairs.c <- MIE:::pairs_of_groups(as.character(vals$keeprows))
      }
      
      
      ##Compute lacking pairs of scalar models
      scalar.pairwise<- pairwiseFit(model = dt$model, 
                                    data = dt$dat, 
                                    group="cntry",
                                    constraints = c("loadings", "intercepts"), 
                                    pairs.of.groups = pairs.c,
                                    message='Fitting extra pairwise scalar models by lavaan',
                                    shiny=TRUE#,
                                    #extra.options = dt$extra.options
                                    )
      
      temp<- cbind(modelStorage$scalar, scalar.pairwise)
      
      attr(temp, "pairs.of.groups") <- rbind(attr(modelStorage$scalar, "pairs.of.groups"),
                                                              attr(scalar.pairwise, "pairs.of.groups"))
      modelStorage$scalar <- temp
      
    }


  # Formatting fit indices for export
    
    list.included.pairs.scalar <- modelStorage$scalar %>% get_pairs %>% apply(2, is_in, vals$keeprows) %>% rowSums(.)==2
    list.included.pairs.metric <- modelStorage$metric %>% get_pairs %>% apply(2, is_in, vals$keeprows) %>% rowSums(.)==2
    
    # if(sum(list.included.pairs.scalar %in% list.included.pairs.metric==FALSE)>0) { 
    #   showNotification("Something went wrong... Switching to covariance view", type="error")
    #   updateRadioButtons(session, "measure", selected="covariance")
    #   }
    
    
    # Subset tables of fit indices
    scalar_subset<- modelStorage$scalar %>%  extract(,list.included.pairs.scalar) %>% extract(rownames(.)  == input$fitincrement.chosen,)
    metric_subset<- modelStorage$metric %>%  extract(,list.included.pairs.metric) %>% extract(rownames(.)  == input$fitincrement.chosen,)
    fit.decrease <- abs(scalar_subset - metric_subset)
    
    # Get a subset of group pairs names
    pair.names.subset <- get_pairs(modelStorage$scalar)[list.included.pairs.scalar,]
    
    
    # Create a distance matrix for MDS   
    dist<- sapply(as.character(vals$keeprows), function(colname) sapply(as.character(vals$keeprows), function(rowname) {
      
      fit.decrease[pair.names.subset[,1]==colname & pair.names.subset[,2]==rowname |
                     pair.names.subset[,1]==rowname & pair.names.subset[,2]==colname]
      
    } )) %>% inset(., sapply(., length)==0, 0)
    
    dist <- matrix(unlist(dist), nrow=dim(dist)[1], ncol=dim(dist)[2], dimnames = dimnames(dist))
    
    
    # Export 
    additional<-data.frame( "Group 1"=pair.names.subset[,1],
                            "Group 2"=pair.names.subset[,2],
                            metric=round(metric_subset,3)  ,
                            scalar=round(scalar_subset, 3),
                            difference=round(fit.decrease,3) )
    
    list(dist=dist, additional=additional)
    
    
    
    
######### Compute configural MGCFA ------------------------------------------------------------    
    
  } else  if(isolate(input$measure == "parameters.loadings")) {
      
    isolate(isolated.modelStorage.loadings <- modelStorage$loadings)

    
    
    print("Begin computing loadings...")
    
    
    print("CHECK:"); print(is.null(temp$old.model.configural.MGCFA))#; print(old.model.configural.MGCFA)
    
    if(is.null(isolated.modelStorage.loadings) |
       sum((vals$keeprows %in% dimnames(isolated.modelStorage.loadings)[[1]])==FALSE) >0 |
        ifelse(is.null(temp$old.model.configural.MGCFA), TRUE, temp$old.model.configural.MGCFA != dt$model) |
       
        ifelse(is.null(temp$old.extra.options.configural.MGCFA), TRUE,
            paste(deparse(temp$old.extra.options.configural.MGCFA, 
                          control=c("quoteExpressions")), collapse="") !=
            paste(deparse(dt$extra.options,
                          control=c("quoteExpressions")), collapse="")
            )

       ) { 
    
          
    print("Subset is longer than computed model -> I am going to recalculate the whole model")
       #a<-Sys.time()
       
       subset.loadings <- MGCFAparameters(model = dt$model,
                                          data = selectedData(),
                                          parameters = "loadings", 
                                          extra.options=dt$extra.options, 
                                          shiny=TRUE)
       subset.loadings <- unclass(subset.loadings)
    
      #print(paste("Computed in", round(Sys.time()-a), "seconds."))
      
      modelStorage$loadings <- subset.loadings
      
      temp$old.model.configural.MGCFA <- dt$model
      temp$old.extra.options.configural.MGCFA <- dt$extra.options 
      
    } else {
    
    print("Subset existing configural model")
      subset.loadings <- isolated.modelStorage.loadings[vals$keeprows,]
    }
      

      dist<-stats::dist(subset.loadings)
      additional<-data.frame(group=rownames(subset.loadings), round(subset.loadings,3))
      list(dist=dist, additional=additional)
      
      
     
######### Compute metric MGCFA ------------------------------------------------------------
      
  } else  if(input$measure == "parameters.intercepts") {
    
    print("Getting parameter intercepts...")
    
    isolated.modelStorage.intercepts <-isolate(modelStorage$intercepts)
    
    # Old condition
    # if(is.null(isolated.modelStorage.intercepts) |
    #    sum((vals$keeprows %in% dimnames(isolated.modelStorage.intercepts)[[1]])==FALSE) >0 ) {
    
    # New condition
    if(is.null(isolated.modelStorage.intercepts) |
       sum((vals$keeprows %in% dimnames(isolated.modelStorage.intercepts)[[1]])==FALSE) >0 |
       ifelse(is.null(temp$old.model.metric.MGCFA), TRUE, temp$old.model.metric.MGCFA != dt$model) |
       
       ifelse(is.null(temp$old.extra.options.metric.MGCFA), TRUE,
              paste(deparse(temp$old.extra.options.metric.MGCFA, control=c("quoteExpressions")), collapse="") !=
              paste(deparse(dt$extra.options,                   control=c("quoteExpressions")), collapse="")
       )
       
    ) {
    
    print("Subset is longer than computed model -> I am going to recalculate the whole model")
    
      subset.intercepts<- MGCFAparameters(model = dt$model,
                                          data = selectedData(), 
                                          parameters = "intercepts", 
                                          extra.options=dt$extra.options, 
                                          shiny=TRUE)
      subset.intercepts <- unclass(subset.intercepts)
      modelStorage$intercepts <- subset.intercepts
      
      temp$old.model.metric.MGCFA <- dt$model
      temp$old.extra.options.metric.MGCFA <- dt$extra.options 
      
    } else {
      
      print("Subset existing configural model")
      subset.intercepts<- isolated.modelStorage.intercepts[vals$keeprows,]
      
    }
        #validate(need(!is.null(subset.intercepts), "Something went wrong when extracting intrcepts"))
        
        dist<-dist(subset.intercepts)
        additional<-data.frame(group=rownames(subset.intercepts), round(subset.intercepts,3))
        list(dist=dist, additional=additional)        
        
        
        
  } else {  
        warning("Problems with input$measure")
    }
      
   
   
})  


  # Run MDS using output of subsettingMatrices() ####
  mds1 <- reactive({ 
    
    print("Trying to compute MDS.")

   # If there is some data loaded and number of countries is not more than 2
      if( !is.null(dt$dat) &  length(vals$keeprows) < 3  ) { 
        
        showNotification("The number of groups should be more than 2. Resetting to the initial number of groups.",
                         type="warning", duration=5)
        
        vals$keeprows <- unique(dt$dat[,1])
        vals$excluded<-NULL
        
        } else if (is.null(dt$dat)) {
          print("Didn't compute MDS, because dt$dat is null.")
        } else {
    
        print(paste("Computimg MDS for vals$keeprows=",paste(vals$keeprows, collapse=",")))
          
        mds.res<- cmdscale(subsettingMatrices()$dist, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
        #print(mds.res)
        
        if(ncol(mds.res)==1) {
          
          print("Problem in computing MDS, trying to fix it by multuplication by 10.")
 #         showNotification("There was a problem in computing MDS. Trying to fix by multiplying the distance matrix by 10.", type="warning", duration=5)
          
          mds.res<-cmdscale(subsettingMatrices()$dist * 1, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
#          if(ncol(mds.res)==1)  showNotification("Problem wasn't solved", type="error")

        }
        
        mds.res
  
        } 
        
    
    })
  

  # k-means cluster analysis using output of subsettingMatrices() ####
  clusters <- reactive({ 
   # groups<-ifelse(class(covMatrices()$dist)=="dist", attr(covMatrices()$dist, "Size"), nrow(covMatrices()$dist))
    print("Clustering...")
      groups <- length(unique(vals$keeprows))
      updateSliderInput(session, "nclusters", max = groups - 1)
  
      kmeans(isolate(subsettingMatrices()$dist), input$nclusters)$cluster
      
    })
  
  
  

  # CREATE OUTPUTS #####
  
  #..Event formula area and controls ####
  output$formulaArea <- renderUI({
    textAreaInput("model", "Your model lavaan syntax", rows=7, cols=25, placeholder="Paste you model code here"        #,
                  # value = 
                  #   "#By default a simulated data - model mimics Schwartz values ESS scale
                  # person=~ ipcrtiv +impfree +impfun +ipgdtim +impdiff +ipadvnt+ imprich +iprspot +ipshabt +ipsuces;
                  # social=~ impenv +ipeqopt +ipudrst +iplylfr +iphlppl +impsafe +ipstrgv +ipfrule +ipbhprp +ipmodst +imptrad;")
  )})
  
  observeEvent(input$use.formula, {
    
    if(input$use.formula==TRUE) {
      #print("input$model:"); print(input$model)
      print("Updating formula area.")
      if(input$model=="")  {
        showNotification("The model formula should be specified!", type="error")
        updateCheckboxInput(session, "use.formula", value=FALSE)
        #dt$model <- ""
      } else {
      
      dt$model <- input$model
      output$formulaArea <- renderUI( pre(dt$model) )
      print(dt$model)
      #current.status<-input$measure
      #updateRadioButtons(session, "measure", selected="covariance")
      #updateRadioButtons(session, "measure", selected=current.status)
      }
    } else {
      
      output$formulaArea <- renderUI({
        textAreaInput("model", "Your model lavaan syntax",
                      rows=7, cols=25, placeholder="Paste you model code here",
                                                     value = dt$model)})
      
    }
    
  })
  
 
  
  # ..Options button and dialog #####
  
  observeEvent(input$lavaan.options, {
    
    extra.options.string <- paste(deparse(dt$extra.options, control=c("quoteExpressions")), collapse="")
  
    if(extra.options.string=="NULL") {
      extra.options.string<-NULL 
      } else {
        extra.options.string<-gsub("^list\\(|\\)$", "", extra.options.string)
        }
    
    showModal(modalDialog(
      title = "lavaan options",
      "Override defaults of `lavaan::cfa()` function. Type here any argument beside formula, data, and group. This is an EXPERIMENTAL option, beware!",
      textAreaInput("model.options", "",
                    rows=7, cols=25, 
                    placeholder=
"ordered = c('impfree', 'impfun'),
orthogonal = TRUE,
group.partial = c('person =~ impfree') ",
                    value = extra.options.string ),
      actionButton("save.model.options", "Save", class="buttonHighlighted"),
      modalButton("Close", icon=icon("times")),
      footer = NULL,
      size = "s", easyClose=T
    ))
  })
  
  observeEvent(input$save.model.options, {
    l <- gsub("\"|“|”", "'", input$model.options)
    l <- paste( "list(", l, ")" )
    
    dt$extra.options <- eval(parse(text=l))
    
    # group.partial = c("person =~ impfree", "person=~ impfun")
    removeModal()
  } )
  
  
  
  #..verbatim text ( mostly for globalMI) ####

  observeEvent(input$semTools, {
 if(input$semTools==TRUE) {
   
   output$verbatimText <- renderUI({
     tagList(
       verbatimTextOutput("verb"))})
   
  output$verb <- renderText ({
    
    withProgress(message = 'Computing global MI test', value = 0, {
      incProgress(1/2, detail = "")
      #library("semTools")
      d=selectedData()
      cfa.argument.list <- c(dt$extra.options, list(model=dt$model, data=d, group="cntry"))
      r<-capture.output(do.call("globalMI",  cfa.argument.list, quote = FALSE))
      
      # r<-capture.output(globalMI(dt$model, data=d, group="cntry", dt$extra.options))

      paste("Global MI output:","\n",
            paste(r, collapse="\n"))
    
  })
  }) } else {
     output$verbatimText <- renderUI({
       tagList(
         #helpText("To run three omnibus tests click the checkbox at the left.")
         )
       
       })
  }
  })  

  
  #..Header for  table of computed measures ####
  output$table_header <- renderUI({
  
    validate(need(!is.null(dt$dat), message="No data"),
             need(input$measure=="covariance"|input$measure=="correlation" |  !is.null(dt$model), 
                  message="No model"))
    
    a<-data.frame(a=c("covariance", "correlation", "parameters.loadings", "parameters.intercepts", "fitincrement.metric", "fitincrement.scalar"),
                  b=c("Covariances between all the variables in the dataset",
                      "Correlations between all the variables in the dataset",
        "Loadings from configural MGCFA model",
        "Intercepts from metric MGCFA model",
        paste(toupper(input$fitincrement.chosen),"difference between configural and metric models"),
        paste(toupper(input$fitincrement.chosen),"difference between metric and scalar models"))) 
    one<-

    h5("Computed measures:", a[a$a==input$measure,"b"])
    
  
  })
  
  #..Table of computed measures ####
  
  output$tabl.DT <- DT::renderDataTable({
    
    validate(need(!is.null(dt$dat), message="Data need to be loaded"),
             need(input$measure=="covariance" | input$measure=="correlation"| !is.null(dt$model), 
                  message="Model needs to be specified"))
      
    DT::datatable(
      subsettingMatrices()$additional, rownames=F, options = list(paging = FALSE, searching = FALSE, info=FALSE)
      )
    })
  
  
  #..Panel of excluded groups ####
  output$excluded <- renderUI({
    if(!is.null(vals$excluded)) {
      
      if(!length(vals$excluded)==0) {
        fluidRow(id="excludePanel",
                 
                 checkboxGroupInput("incl", "Excluded groups (click to include): ", vals$excluded, inline=T),
                 div(actionLink("resetExcluded", "", icon =icon("remove-circle", lib="glyphicon")), align="right")
        )
      } else {
        helpText("Click points on the plot to exclude from analysis")
      }
      
    } else {  
        helpText("Click points on the plot to exclude from analysis") 
      }
  })
  
  #....click observer for panel of excluded groups ####
  observeEvent(input$incl, {
    showNotification(paste("Group", input$incl, "was returned to the plot."),
                     type="message", duration=3)
    
    vals$keeprows<-c(vals$keeprows, input$incl)   
    vals$excluded<-vals$excluded[! vals$excluded %in% input$incl]
  })
  
  observeEvent(input$resetExcluded, {
    vals$keeprows<-unique(dt$dat$cntry)   
    vals$excluded<-NULL
  })
  

  
 #..The Plot ####
  output$distPlot <- renderPlot({
    print("Attempting to plot")
    #print("dt$model"); print(dt$model)
    
    validate(need(!is.null(dt$dat), message="Data need to be loaded"),
             need(input$measure=="covariance"|input$measure=="correlation" |  !is.null(dt$model), 
                  message="Model needs to be specified"),
             need(  ifelse(is.matrix(mds1()), ncol(mds1())==2, TRUE),
             message="Can't create two-dimensional representation, because got negative eigenvalue. \nTry to include/exclude groups or use another measure. It's also possible that you have already found a set of invariant groups.")
             )
    
    
    
    d<-mds1() %>% set_colnames(c("dim1", "dim2")) %>% as.data.frame %>%
      mutate(group=rownames(.), cluster=clusters())
    
    g<-ggplot(d, aes(dim1, dim2,  col=as.factor(cluster)))+
      geom_point( size=5, show.legend = F)+labs(x="", y="", col="")+
      geom_text_repel(aes(label=group),point.padding = unit(.3, "lines"), show.legend=F)+
      
      #lims(x=c(min(d$dim1), ifelse(max(d$dim1)-min(d$dim1) < 0.01, min(d$dim1)+0.015, max(d$dim1))),
      #     y=c(min(d$dim2), ifelse(max(d$dim2)-min(d$dim2) < 0.01, min(d$dim2)+0.015, max(d$dim2))))+
      theme_minimal()+
      #coord_fixed()+
      scale_colour_hue(l = 50, c = 120)+
      theme(panel.grid = element_blank(), axis.line=element_line(size=.5),axis.ticks=element_line(size=.5), plot.title=element_text(face="bold", size=18))+
      ggtitle(paste("Clustering based on", isolate(input$measure) ))
    # d<-cmdscale(dist(cov(data.Benelux[, 10:15], use="complete.obs")), 2)%>% set_colnames(c("dim1", "dim2")) %>% as.data.frame %>% dplyr::mutate(group=rownames(.), cluster=clustr)
    
    
    if(isolate(input$measure=="fitincrement.scalar"|input$measure=="fitincrement.metric")) {
           
      
   r=   switch(isolate(input$fitincrement.chosen), 
             cfi=.005,
             rmsea=.0075,
             srmr=0,
             nnfi=0,
             gfi=0,
             rmsea.ci.upper = 0
             )
   
      print("this is r")
      print(isolate(input$measure))
      
      #r=.005
      xc = min(d$dim1) +  sqrt(r^2/2)
      yc = min(d$dim2) +   sqrt(r^2/2)
      
      circle.data <- data.frame(x = xc+r*cos(seq(0,2*pi,length.out=100)), 
                                y=  yc+r*sin(seq(0,2*pi,length.out=100)))
      
      g<-g+geom_path(data=circle.data, aes(x,y#, label=NULL
                                           ), col="lightgray", linetype="dashed")
      if(r!=0) {
      g<-g+labs(caption=paste("The circle has diameter", round(r*2, 3), "meaning the increment \nin the fit index is within interval recommended by Chen"))
      }
      
    } else if (isolate(input$measure == "parameters.loadings")) {
      
        fits<-attr(modelStorage$loadings, "fit")[c("cfi", "rmsea", "srmr")]
    
        caption<-paste(paste(c("CFI=", "RMSEA=", "SRMR="), sep=""),
                       paste(round(fits,3), sep=""), collapse=", ")
    
    
    
        if(  nrow(modelStorage$loadings)>length(vals$keeprows)) {
          caption <- paste("Fit for ", nrow(modelStorage$loadings), "groups:", caption, ". \nThe graph is based on a subset of parameters from larger model.")
    
          #output$forceFitLink <- renderUI({ actionLink("forceFitting", "Refit for current groups") })
    
        }
        
        g<-g+labs(caption=caption)
      
    } else if (isolate(input$measure == "parameters.intercepts")) {
      
      fits<-attr(modelStorage$intercepts, "fit")[c("cfi", "rmsea", "srmr")]
      
      caption<-paste(paste(c("CFI=", "RMSEA=", "SRMR="), sep=""),
                     paste(round(fits,3), sep=""), collapse=", ")
      
      
      
      if(  nrow(modelStorage$intercepts)>length(vals$keeprows)) {
        caption <- paste("Fit for ", nrow(modelStorage$intercepts), "groups:", caption, ". \nThe graph is based on a subset of parameters from larger model.")
        
        #output$forceFitLink <- renderUI({ actionLink("forceFitting", "Refit for current groups") })
        
        }
      
      g<-g+labs(caption=caption)
      
      
    } else {}
    
    g
  })
  
 
  
  #...click observer for the plot #####
  observeEvent(input$plot1_click, {
    #covMatrices()
    
    d<-as.data.frame(mds1(), row.names=rownames(mds1()))
    res <- nearPoints(d, xvar=names(d)[1], yvar=names(d)[2], coordinfo=input$plot1_click, maxpoint=1, allRows=F, threshold=10)
    vals$keeprows<-rownames(mds1())[! rownames(mds1()) %in% rownames(res) ]
    vals$excluded<-unique(dt$dat[,1])[!unique(dt$dat[,1]) %in% vals$keeprows]
    #output$vals$excluded<- unique(dt$dat[,1])[!unique(dt$dat[,1]) %in% vals$keeprows]
    #groups<-vals$keeprows
    
  })
  
  plot.size<- reactiveVal(value=550)
  
  observeEvent(input$zoomIn, {
    new <- plot.size() + 50
    plot.size(new)
    #print(plot.size())
  })
  observeEvent(input$zoomOut, {
    new <- plot.size() - 50
    plot.size(new)
    #print(plot.size())
  })

  output$plot<- renderUI({  plotOutput("distPlot", height=paste(plot.size(), "px", sep=""),
             click = "plot1_click")
             #,brush = brushOpts(id = "plot1_brush")
  })
  

  
})


