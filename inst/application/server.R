
# ADD IF conditions to pairwise, and add extra.options to pairwise.fit function


#setwd("~/Dropbox/STAT/European Social Survey/ShinyValues/MeasurementInvarianceExplorer")


# TO develop: errors handlers; 
# Speed of pairwise (think of what can be done with non-pairwise)



requireNamespace("shiny", quietly = T)
requireNamespace("lavaan", quietly = T)
requireNamespace("magrittr", quietly = T)
requireNamespace("reshape2", quietly = T)
requireNamespace("markdown", quietly = T)
requireNamespace("DT", quietly = T)
requireNamespace("ggplot2", quietly = T)
requireNamespace("ggrepel", quietly = T)
requireNamespace("dplyr", quietly = T)
requireNamespace("shinyjs", quietly = T)
requireNamespace("shinyWidgets", quietly = T)
options(shiny.maxRequestSize=100*1024^2) 


#source("server.functions.R")


# Begin the SHINY code #####
shinyServer(function(session, input, output) {

verb <- function(...) { if(.verbose) cat("\n", ...)  }

if(.verbose) cli::cat_rule("Begin log")

### Values ####
  dt<-reactiveValues ( 
    dat = NULL,
    model = NULL,
    extra.options = NULL,
    plot = NULL
     
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
  if( ifelse(length(.model == "demo") == 0, FALSE, .model == "demo") ) { 
    shinyjs::click("useSimulated")
    
    } else {
      
      
  if(!is.null(.data) & !is.null(.group) ) { 
    if(any(class(.data)=="tbl")) .data <- as.data.frame.list(lapply(unclass(.data), unclass))
    local.data <- .data[, c(colnames(.data)[colnames(.data) == .group],
                       colnames(.data)[colnames(.data) != .group])]
    verb(colnames(local.data))
    colnames(local.data)[1] <- "grp"
    verb("changed the data order")
    
    if(any(!sapply(local.data[-1], is.numeric))) 
      stopApp(stop("All the variables in the dataset (beside group) must be numeric"))
    
    
  } else if (!is.null(.data) & is.null(.group) ) {
    #local.data <- .data
    stop("The grouping variable must be specified.")
  }
  
  
  if(!is.null(.model)) { 
    dt$model = .model
    updateCheckboxInput(session, "use.formula", value=TRUE)
  }
  
  
  if(!is.null(.data)) {
     dt$dat = local.data
  #print("unique(dt$dat$grp)"); print(unique(dt$dat$grp))
  vals$keeprows = unique(isolate(dt$dat$grp))
  vals$excluded <- NULL
  modelStorage$covariance <- t(computeCovariance(isolate(dt$dat), group = "grp"))
  
  showNotification("Using data from the R object.", type="message", duration=10)
  }
}
# }) 
  
  
  
  # .. Button for using simulated data #####
  observeEvent(input$useSimulated, {
    dt$dat = read.csv("simulated2.csv")
    dt$model = "#By default a simulated data - model mimics Schwartz values ESS scale
  person=~ ipcrtiv +impfree +impfun +ipgdtim +impdiff +ipadvnt+ imprich +iprspot +ipshabt +ipsuces;
    social=~ impenv +ipeqopt +ipudrst +iplylfr +iphlppl +impsafe +ipstrgv +ipfrule +ipbhprp +ipmodst +imptrad;"
    updateCheckboxInput(session, "use.formula", value=TRUE)
    #print("unique(dt$dat$grp)"); print(unique(dt$dat$grp))
    vals$keeprows = unique(dt$dat$grp)
    vals$excluded <- NULL
    print(paste("Button 'play with fake data' has been used."))
    
    modelStorage$covariance <- computeCovariance(dt$dat, group = "grp")
    
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
    names(d)[1]<-c("grp")
    d$grp<-as.factor(d$grp) # Hmm...
    dt$dat <- d
    print(paste(colnames(d), collapse="; "))
    rm(d)
    #Set subsets to NULL
        vals$keeprows <- unique(dt$dat$grp)
        vals$excluded = NULL
        
    #Set all model results (beside covariance) to NULL
        #for(x in names(modelStorage)[-1]) modelStorage[[x]]<-NULL

        ##modelStorage <- lapply(names(modelStorage), function(x) modelStorage[[x]]<-NULL)
    
        

   #     showNotification("Computing covarances for a new data...", type="warning", duration=3)    
    
    #Compute covariance matrix   
    #Split dataset and compute variance-covariance for each group separately
        modelStorage$covariance<-computeCovariance(isolate(dt$dat), group = "grp")
        
        #print( str(modelStorage$covariance))
       # print(head(dt$dat))
  })

 
  
  
#Subset the data #####  
  selectedData <- reactive({
    
    
    print(paste(" selectedData() subset the raw data for", paste(vals$keeprows, collapse=",")))
      dat<-dt$dat

      dat<-dat[dat$grp %in% vals$keeprows, ] 
      if(is.factor(dat$grp)) dat$grp <- droplevels(dat$grp) 
   
    dat
    
    })

 
#Subset measures that are saved to modelStorage. Reactive subsettngMatrices() ####
subsettingMatrices <- reactive ({
  
  if(input$measure=="covariance") {
    
  computeCovariance(data  = dt$dat, group = "grp")
  
    } else if (input$measure=="correlation"){
    
  computeCorrelation(data  = dt$dat, group = "grp")

      ######### Fitting conf vs metric pairwise -----------------------------------------------------    
    } else if (input$measure=="fitincrement.metric") { 
    
    # This takes saved model results from modelStorage 
    #  and subsets it to fit the list of included countries and chosen fit measure
    # or, if there are extra pairs not computed yet, computes it
    
    # Fitting configural pairwise models  
      
      it.is.required.to.refit.model.to.new.formula <- 
        ifelse( !is.null(dt$model),  # if there is model in dt storage
                ifelse(!is.null(attr(modelStorage$conf, "model.formula")), # and there is model in configural fit
                       dt$model!=attr(modelStorage$conf, "model.formula"), # but they differ
                       TRUE), # or there is no model in configural fit but there is some in dt storage
                FALSE) # there is no model in dt storage
      
  if(is.null(modelStorage$conf) | #no configural fit
     any(! vals$keeprows %in% unique(unlist(attr(modelStorage$conf, "pairs.of.groups")))) | # new groups added
     it.is.required.to.refit.model.to.new.formula
     )     {
     
      
      
      
      if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$conf) ) {
        
        verb("Recorded model is:", attr(modelStorage$conf, "model.formula"))
        verb("Updated model is:", paste(dt$model))
        pairs.c <- MIE:::pairs_of_groups(as.character(vals$keeprows))
        
      } else {
        
        verb("There are some uncomputed pairs of configural models.")
        extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$conf, "pairs.of.groups")))]
        pairs.c <- MIE:::pairs_of_groups(as.character(extra.countries))
        
        #pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$conf, "pairs.of.groups"))), stringsAsFactors = F)
        #names(pairs.c)<-names(attr(modelStorage$conf, "pairs.of.groups"))
        
      }
      
      
      verb("Computing missing pairs of conf models");

        
      # conf.pairwise<- pairwiseFit(model = dt$model,
      #                             data  = dt$dat, 
      #                             group = "grp",
      #                             constraints = c(""),
      #                             pairs.of.groups = pairs.c, 
      #                             message = 'Fitting pairwise configural models by lavaan',
      #                             shiny = TRUE#,
      #                             #extra.options = dt$extra.options
      #                             )
                                  
      # attempt to include extra options
      conf.pairwise<-    do.call("pairwiseFit", append(list(
        model = dt$model,
        data  = dt$dat, 
        group = "grp",
        constraints = c(""),
        pairs.of.groups = pairs.c, 
        message = 'Fitting pairwise configural models by lavaan',
        shiny = TRUE
      ), dt$extra.options))
      
      
      
      if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$conf))  { #rewrite
        
        attr(conf.pairwise, "model.formula") <- dt$model
        modelStorage$conf <- conf.pairwise
        rm(conf.pairwise)
        
      } else { # Merge with previous fits 
        
        temp<- cbind(modelStorage$conf, conf.pairwise)
        attr(temp, "pairs.of.groups")<- rbind(attr(modelStorage$conf, "pairs.of.groups"),
                                              attr(conf.pairwise, "pairs.of.groups"))
        attr(temp, "model.formula") <- dt$model
        modelStorage$conf <- temp
        rm(temp)
      }

   
    }
    
    # Fitting metric pairwise models
      
      it.is.required.to.refit.model.to.new.formula <- 
        ifelse( !is.null(dt$model),  # if there is model in dt storage
                ifelse(!is.null(attr(modelStorage$metric, "model.formula")), # and there is model in metric fit
                       dt$model!=attr(modelStorage$metric, "model.formula"), # but they differ
                       TRUE), # or there is no model in metric fit but there is some in dt storage
                FALSE) # there is no model in dt storage
      
      if(is.null(modelStorage$metric) | #no metric fit
         any(! vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))) | # new groups added
         it.is.required.to.refit.model.to.new.formula
      )     {
        
        
        
        
        if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$metric) ) {
          
          verb("Recorded model is:", attr(modelStorage$metric, "model.formula"))
          verb("Updated model is:", paste(dt$model))
          pairs.c <- MIE:::pairs_of_groups(as.character(vals$keeprows))
          
        } else {
          
          verb("There are some uncomputed pairs of metric models.")
          extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))]
          pairs.c <- MIE:::pairs_of_groups(as.character(extra.countries))
          
          #pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$metric, "pairs.of.groups"))), stringsAsFactors = F)
          #names(pairs.c)<-names(attr(modelStorage$metric, "pairs.of.groups"))
          
        }
        
        
        verb("Computing metric models for "); print(pairs.c);
        
        # metric.pairwise<- pairwiseFit(model = dt$model,
        #                             data  = dt$dat, 
        #                             group = "grp",
        #                             constraints = c("loadings"),
        #                             pairs.of.groups = pairs.c, 
        #                             message = 'Fitting pairwise metric models by lavaan',
        #                             shiny = TRUE#,
        #                             #extra.options = dt$extra.options
        # )
        
        # attempt to include extra options
        metric.pairwise<-   do.call("pairwiseFit", append(list(
          model = dt$model,
          data  = dt$dat, 
          group = "grp",
          constraints = c("loadings"),
          pairs.of.groups = pairs.c, 
          message = 'Fitting pairwise metric models by lavaan',
          shiny = TRUE
        ), dt$extra.options))
        
        
        
        
        if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$metric))  { #rewrite
          
          modelStorage$metric <- metric.pairwise
          
        } else { # Merge with previous fits 
          
          temp<- cbind(modelStorage$metric, metric.pairwise)
          attr(temp, "pairs.of.groups")<- rbind(attr(modelStorage$metric, "pairs.of.groups"),
                                                attr(metric.pairwise, "pairs.of.groups"))
          attr(temp, "model.formula") <- dt$model
          modelStorage$metric <- temp
          rm(temp)
        }
        
        rm(metric.pairwise)
        
      }
  

    
# Formatting fit indices for export
   verb("Formatting fit indices for export")
   conf <- isolate(modelStorage$conf )
   metric <- isolate(modelStorage$metric)
   fit.decrease <- abs(conf - metric)
    
    detailed<-lapply(rownames(fit.decrease), function(f) {
        cbind(configural= conf[f,], 
             metric= metric[f,], 
             fit.decrease=fit.decrease[f, ])
      })
  
    
  names(detailed)<-rownames(fit.decrease)
  out <- list(detailed=detailed, bunch=fit.decrease)
  class(out)<-c("incrementalFit")
  out
  
 
  ###### Fitting scalar vs metric pairwise -----------
  } else if(input$measure=="fitincrement.scalar") {
  
    
    # This takes saved model results from modelStorage 
    #  and subsets it to fit the list of included countries and chosen fit measure
    # or, if there are extra pairs not computed yet, computes it
    
    
    # Fitting metric pairwise models
    
    it.is.required.to.refit.model.to.new.formula <- 
      ifelse( !is.null(dt$model),  # if there is model in dt storage
              ifelse(!is.null(attr(modelStorage$metric, "model.formula")), # and there is model in metric fit
                     dt$model!=attr(modelStorage$metric, "model.formula"), # but they differ
                     TRUE), # or there is no model in metric fit but there is some in dt storage
              FALSE) # there is no model in dt storage
    
    if(is.null(modelStorage$metric) | #no metric fit
       any(! vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))) | # new groups added
       it.is.required.to.refit.model.to.new.formula
    )     {
      
      
      
      
      if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$metric) ) {
        
        verb("Recorded model is:", attr(modelStorage$metric, "model.formula"))
        verb("Updated model is:", paste(dt$model))
        pairs.c <- MIE:::pairs_of_groups(as.character(vals$keeprows))
        
      } else {
        
        verb("There are some uncomputed pairs of metric models.")
        extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$metric, "pairs.of.groups")))]
        pairs.c <- MIE:::pairs_of_groups(as.character(extra.countries))
        
        #pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$metric, "pairs.of.groups"))), stringsAsFactors = F)
        #names(pairs.c)<-names(attr(modelStorage$metric, "pairs.of.groups"))
        
      }
      
      
      verb("Computing metric models for "); print(pairs.c);
      
      # metric.pairwise<- pairwiseFit(model = dt$model,
      #                               data  = dt$dat, 
      #                               group = "grp",
      #                               constraints = c("loadings"),
      #                               pairs.of.groups = pairs.c, 
      #                               message = 'Fitting pairwise metric models by lavaan',
      #                               shiny = TRUE#,
      #                               #extra.options = dt$extra.options
      # )
      
      # attempt to include extra options
      metric.pairwise<-   do.call("pairwiseFit", append(list(
        model = dt$model,
        data  = dt$dat, 
        group = "grp",
        constraints = c("loadings"),
        pairs.of.groups = pairs.c, 
        message = 'Fitting pairwise metric models by lavaan',
        shiny = TRUE
      ), dt$extra.options))
      
      
      if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$metric))  { #rewrite
        
        modelStorage$metric <- metric.pairwise
        
      } else { # Merge with previous fits 
        
        temp<- cbind(modelStorage$metric, metric.pairwise)
        attr(temp, "pairs.of.groups")<- rbind(attr(modelStorage$metric, "pairs.of.groups"),
                                              attr(metric.pairwise, "pairs.of.groups"))
        attr(temp, "model.formula") <- dt$model
        modelStorage$metric <- temp
        rm(temp)
      }
      
      rm(metric.pairwise)
      
    }
    
    # Fitting scalar pairwise models  
    
    it.is.required.to.refit.model.to.new.formula <- 
      ifelse( !is.null(dt$model),  # if there is model in dt storage
              ifelse(!is.null(attr(modelStorage$scalar, "model.formula")), # and there is model in scalar fit
                     dt$model!=attr(modelStorage$scalar, "model.formula"), # but they differ
                     TRUE), # or there is no model in scalar fit but there is some in dt storage
              FALSE) # there is no model in dt storage
    
    if(is.null(modelStorage$scalar) | #no scalar fit
       any(! vals$keeprows %in% unique(unlist(attr(modelStorage$scalar, "pairs.of.groups")))) | # new groups added
       it.is.required.to.refit.model.to.new.formula
    )     {
      
      
      
      
      if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$scalar) ) {
        
        verb("Recorded model is:", attr(modelStorage$scalar, "model.formula"))
        verb("Updated model is:", paste(dt$model))
        pairs.c <- MIE:::pairs_of_groups(as.character(vals$keeprows))
        
      } else {
        
        verb("There are some uncomputed pairs of scalar models.")
        extra.countries <- vals$keeprows[!vals$keeprows %in% unique(unlist(attr(modelStorage$scalar, "pairs.of.groups")))]
        pairs.c <- MIE:::pairs_of_groups(as.character(extra.countries))
        
        #pairs.c <- expand.grid(extra.countries, unique(unlist(attr(modelStorage$scalar, "pairs.of.groups"))), stringsAsFactors = F)
        #names(pairs.c)<-names(attr(modelStorage$scalar, "pairs.of.groups"))
        
      }
      
      
      verb("Computing missing pairs of scalar models");
      # scalar.pairwise<- pairwiseFit(model = dt$model,
      #                             data  = dt$dat, 
      #                             group = "grp",
      #                             constraints = c("loadings", "intercepts"),
      #                             pairs.of.groups = pairs.c, 
      #                             message = 'Fitting pairwise scalar models by lavaan',
      #                             shiny = TRUE#,
      #                             #extra.options = dt$extra.options
      # )
      
      # attempt to include extra options
      scalar.pairwise<-   do.call("pairwiseFit", append(list(
        model = dt$model,
        data  = dt$dat, 
        group = "grp",
        constraints = c("loadings", "intercepts"),
        pairs.of.groups = pairs.c, 
        message = 'Fitting pairwise scalar models by lavaan',
        shiny = TRUE
      ), dt$extra.options))
      
      if(it.is.required.to.refit.model.to.new.formula | is.null(modelStorage$scalar))  { #rewrite
        
        attr(scalar.pairwise, "model.formula") <- dt$model
        modelStorage$scalar <- scalar.pairwise
        rm(scalar.pairwise)
        
      } else { # Merge with previous fits 
        
        temp<- cbind(modelStorage$scalar, scalar.pairwise)
        attr(temp, "pairs.of.groups")<- rbind(attr(modelStorage$scalar, "pairs.of.groups"),
                                              attr(scalar.pairwise, "pairs.of.groups"))
        attr(temp, "model.formula") <- dt$model
        modelStorage$scalar <- temp
        rm(temp)
      }
      
      
    }
    
    
    
    # Formatting fit indices for export
    verb("Formatting fit indices for export")
    metric <- isolate(modelStorage$metric)
    scalar <- isolate(modelStorage$scalar )
    fit.decrease <- abs(scalar - metric)
    
    detailed<-lapply(rownames(fit.decrease), function(f) {
      cbind(metric= metric[f,], 
            scalar= scalar[f,], 
            fit.decrease=fit.decrease[f, ])
    })
    
    
    names(detailed)<-rownames(fit.decrease)
    out <- list(detailed=detailed, bunch=fit.decrease)
    class(out)<-c("incrementalFit")
    out
    
    
    
    
    
    
    
######### Compute configural MGCFA ------------------------------------------------------------    
    
  } else  if(isolate(input$measure == "parameters.loadings")) {
      
    isolate(isolated.modelStorage.loadings <- modelStorage$loadings)

    
    
    verb("Begin computing loadings...")

    
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
    
          
    verb("Subset is longer than computed model -> I am going to recalculate the whole model")
       #a<-Sys.time()
       
       subset.loadings <- MGCFAparameters(model = dt$model,
                                          data = selectedData(),
                                          parameters = "loadings", 
                                          group = "grp",
                                          extra.options=dt$extra.options, 
                                          shiny=TRUE)

      
      modelStorage$loadings <- subset.loadings
      
      
    } 
      
    isolate(modelStorage$loadings)
     
######### Compute metric MGCFA ------------------------------------------------------------
      
  } else  if(input$measure == "parameters.intercepts") {
    
    verb("Getting parameter intercepts...")
    
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
    
    verb("Subset is longer than computed model -> I am going to recalculate the whole model")
    
      subset.intercepts<- MGCFAparameters(model = dt$model,
                                          data = selectedData(), 
                                          parameters = "intercepts", 
                                          group = "grp",
                                          extra.options=dt$extra.options, 
                                          shiny=TRUE)
      
      modelStorage$intercepts <- subset.intercepts
      
     
      
    } 
    
    
        isolate(modelStorage$intercepts)
        
        
  } else {  
        warning("Problems with input$measure")
    }
      
   
   
})  


  
  
  

  # CREATE OUTPUTS #####
  
  #..Event 'formula area' and controls ####
  output$formulaArea <- renderUI({
    textAreaInput("model", "Your model lavaan syntax", rows=7, cols=25, placeholder="Paste you model code here" )})
  
  observeEvent(input$use.formula, {
    
    if(input$use.formula==TRUE) {

      verb("Updating formula area.")
      if(input$model=="")  {
        showNotification("The model formula should be specified!", type="error")
        updateCheckboxInput(session, "use.formula", value=FALSE)
      } else {
      
      dt$model <- input$model
      output$formulaArea <- renderUI( pre(dt$model) )
      verb(dt$model)

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
    
    #extra.options.string <- paste(deparse(dt$extra.options, control=c("quoteExpressions")), collapse="")
    # if(extra.options.string=="NULL") {
    #   extra.options.string<-NULL 
    # } else {
    #   extra.options.string<-gsub("^list\\(|\\)$", "", extra.options.string)
    # }
    
    if( length(dt$extra.options)!=0) {
    extra.options.string <- paste(sapply(1:length(dt$extra.options), 
                     function(x) paste(names(dt$extra.options)[x], "=", 
                                       ifelse(is.character(dt$extra.options[[x]]),
                                          paste0("'", dt$extra.options[[x]], "'"),
                                          dt$extra.options[[x]])
                                       
                                       )), collapse = ",\n")
    
    } else {
      extra.options.string = dt$extra.options
    }
  
  
    
    showModal(modalDialog(
      title = "lavaan options",
      "Override defaults. Type here any argument available in lavaan's `cfa` and `lavOptions`, except formula, data, and group. This is an EXPERIMENTAL option, beware!",
      textAreaInput("model.options", "",
                    rows=7, cols=25, 
                    placeholder=
"group.partial = c('person =~ impfree'),
 orthogonal = TRUE",
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

  observeEvent(input$globalMI, {
 if(input$globalMI==TRUE) {
   
   output$verbatimText <- renderUI({
     tagList(
       verbatimTextOutput("verb"))})
   
  output$verb <- renderText ({
    
    withProgress(message = 'Computing global MI test', value = 0, {
      incProgress(1/2, detail = "")

      d=selectedData()
      cfa.argument.list <- c(dt$extra.options, list(model=dt$model, data=d, group="grp"))
      r<-capture.output(do.call("globalMI",  cfa.argument.list, quote = FALSE))
      
      # r<-capture.output(globalMI(dt$model, data=d, group="grp", dt$extra.options))

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
  table.header <- reactive({
  
    validate(need(!is.null(dt$dat), message="No data"),
             need(input$measure=="covariance"|input$measure=="correlation" |  !is.null(dt$model), 
                  message="No model"))
    
    a<-data.frame(a=c("covariance",
                      "correlation",
                      "parameters.loadings",
                      "parameters.intercepts",
                      "fitincrement.metric",
                      "fitincrement.scalar"),
                  b=c("Covariances between all the variables in the dataset",
                      "Correlations between all the variables in the dataset",
                      "Loadings from configural MGCFA model",
                      "Intercepts from metric MGCFA model",
        paste(toupper(input$fitincrement.chosen),"difference between configural and metric models"),
        paste(toupper(input$fitincrement.chosen),"difference between metric and scalar models"))) 
    one<-

    paste("Computed measures:", a[a$a==input$measure,"b"])
    
  
  })
  
  #..Table of computed measures ####
  
  output$tabl.DT <- DT::renderDataTable({
    
    validate(need(!is.null(dt$dat), message="Data need to be loaded"),
             need(input$measure=="covariance" | input$measure=="correlation"| !is.null(dt$model), 
                  message="Model needs to be specified"))
    
      dd <- isolate(subsettingMatrices())
      dd1 <- switch(class(dd),
                   correlations =    format(round(unclass(dd[vals$keeprows,]), 2), 
                                            digits = 2),
                   covariances =     format(round(unclass(dd[vals$keeprows,]), 2), 
                                            digits = 3, nsmall = 1, scientific = F),
                   MGCFAparameters = format(round(unclass(dd[vals$keeprows,]), 2), 
                                            digits = 3, nsmall = 1, scientific = F),
                   incrementalFit =  {
                     ind.row <- rowSums(apply(attr(dd$bunch,"pairs.of.groups"), 2,  
                                              `%in%`, vals$keeprows))==2
                     format(round(dd$detailed[[input$fitincrement.chosen]][ind.row,], 2), 
                            digits = 3, nsmall = 1, scientific = F)
                     
                 })

    DT::datatable(dd1, rownames=T, options = list(paging = FALSE, searching = FALSE, info=FALSE),
                  caption=table.header() )
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
    vals$keeprows<-unique(dt$dat$grp)   
    vals$excluded<-NULL
  })
  
  
#... The switch for plot #####
  observeEvent(input$netSwitch, {
    
    if( !any(input$measure %in% c("fitincrement.metric", "fitincrement.scalar"))) {
      
      updateMaterialSwitch(session, "netSwitch", value = FALSE)
    }
    if(!input$netSwitch) removeNotification("nocutoffs")
    
  })
  
  
 #..The Plot ####
  output$distPlot <- renderPlot({
    verb("Attempting to plot")
    
    
    # validate(need(!is.null(dt$dat), message="Data need to be loaded"),
    #          need(input$measure=="covariance"|input$measure=="correlation" |  !is.null(dt$model),
    #               message="Model needs to be specified"),
    #          need(  ifelse(is.matrix(mds1()), ncol(mds1())==2, TRUE),
    #          message="Can't create two-dimensional representation, because of negative eigenvalue. \nTry to include/exclude groups or use another measure. It's also possible that you have already found a set of invariant groups.")
    #          )

      # If there is some data loaded and number of countries is not more than 2
      if( !is.null(dt$dat) &  length(vals$keeprows) < 3  ) {

        showNotification("The number of groups should be more than 2. Resetting to the initial number of groups.", type="warning", duration=5)

        vals$keeprows <- unique(dt$dat[,1])
        vals$excluded<-NULL

      } else if (is.null(dt$dat)) {
        verb("Didn't compute MDS, because dt$dat is null.")

      } else {

        subsettingMatrices()

         }
 

    
d <- isolate(subsettingMatrices())
verb("The data type is ", class(d), "\n")

if(!input$netSwitch) {
  
  g<-  plotDistances(d,
                  n.clusters = "auto",
                  fit.index = input$fitincrement.chosen,
                  drop = vals$excluded,
                  shiny = TRUE)

  dt$plot <- g$data

  g+ggtitle(table.header())
  
} else {
  
  g<-  plotCutoff(d, 
             fit.index = input$fitincrement.chosen, 
             drop = vals$excluded,
             weighted = T,
             shiny = T)
  dt$plot <- g$data
  g+ggtitle(table.header())
  
}
      })
  
 
  
  #...click observer for the plot #####
  observeEvent(input$plot1_click, {
    d<- isolate(dt$plot)
    #d<-as.data.frame(mds1(), row.names=rownames(mds1()))
    res <- nearPoints(d, coordinfo=input$plot1_click, 
                      xvar="dim1", yvar="dim2",  
                      #x="x", y="y",  
                      maxpoints=1, allRows=F, threshold=5)
    verb("res is "); print(res)
    vals$keeprows<-d$group[! d$group %in% res$group ]
    vals$excluded<-unique(dt$dat[,1])[!unique(dt$dat[,1]) %in% vals$keeprows]

  })
  


  output$plot<- renderUI({  plotOutput("distPlot", height=paste(plot.size(), "px", sep=""),
             #click = if(input$netSwitch) NULL else "plot1_click"
             click = "plot1_click"
             )#, brush = brushOpts(id = "plot1_brush")
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
  
})


