#' Run alignment procedure using Mplus
#' 
#' Facilitates running frequentist alignmnet procedure in Mplus. It creates Mplus input and data from R data.frame, runs the input using Mplus, and provides its summaries back in R.
#' 
#' @param model Character. Formula in Mplus format, e.g. "Factor1 BY item1 item2 item3 item4; item4 WITH item3;", see example.
#' @param group Character, name of the grouping variable.
#' @param dat Data frame containing data. 
#' @param categorical Character vector of variable names. Indicators that are binary or ordinal. Default is NULL.
#' @param sim.samples Vector of integers. Group sample sizes for simulation,  the length of this vector also determines  a number of simulation studies. Default is `c(100, 500, 1000)`. May take a substantial amount of time. Use NULL to avoid running simulations.
#' @param sim.reps A number of simulated datasets in each simulation. Default is 500.
#' @param Mplus_com  Sometimes you don't have a direct access to Mplus, so this argument specifies what to send to a system command line. Default value is "mplus".
#' @param path Where all the .inp, .out, and .dat files should be stored?
#' @param summaries If the \code{\link[MIE]{extractAlignment}}     and \code{\link[MIE]{extractAlignmentSim}}    should be run after all the Mplus work is done. Default is FALSE.
#' @details The function runs in four steps to facilitate setting up and running alignment:
#' \enumerate{
#'     \item Converts data.drame into a dat file compatible with Mplus. Saves it on the disk.
#'     \item Creates Mplus input file for free alignemnent. Saves it on the disk.
#'     \item Runs this input using the saved data, which produces an .out file.
#'     \item Reads in the output, selects the appropriate reference group, creates a new input file for fixed alignment, and runs it.
#'     \item (Optional) Runs simulations to check the reliability of the fixed alignment results.
#'     \item (Optional) Using \code{\link[MIE]{extractAlignment}} and/or \code{\link[MIE]{extractAlignmentSim}} reads in the output file of the fixed alignment and summarizies them in a convenient, publishable tables.
#' }
#' The sequence of free, fixed alignment and simulations follows recommendations of Muthen & Asparouhov (2014).
#' 
#' All the files created during these steps stay on the disk and can be used independently for reporting. For a detailed tutorial on alignment in general see \url{https://maksimrudnev.com/2019/05/01/alignment-tutorial/}
#' 
#' 
#' @examples 
#' \dontrun{ aling.res = runAlignment(model = "Moral1 BY prostit homosex abortion divorce;
#' Moral2 BY benefits taxes bribes;", 
#' group = "country",
#'  dat = wvs.s1,
#'  Mplus_com = "mplus",
#'  summaries = T
#'  )
#'  }
#'  @seealso \code{\link[MIE]{extractAlignment}}, \code{\link[MIE]{extractAlignmentSim}}, \code{\link[MIE]{extractAlignmentSim}}, \code{\link[MIE]{measurementInvarianceMplus}}
#'  
#' @export
runAlignment <- function(
  model, 
  group,
  dat, 
  categorical=NULL,
  sim.samples = c(100, 500, 1000), 
  sim.reps = 500,
  Mplus_com = "mplus",
  path = getwd(),
  summaries = FALSE
) {
  
  oldwd <- getwd()
  setwd(path)
  
  message("Creating input for free alignment.\n")
  
  var.list <- strsplit(model, ";|\n") [[1]]
  var.list <- gsub("\\s+", " ",  var.list)
  var.list <-   var.list[!var.list==""]
  var.list <-   var.list[!var.list==" "]
  var.list <- gsub("^.*((by)|(BY))", "", var.list)
  #var.list <-   unlist(strsplit(var.list, "(?i)(by)", perl=TRUE))
  #var.list <-   unlist(strsplit(var.list[seq(2, length(var.list), by=2)], " "))
  var.list <-   unlist(strsplit(var.list, " "))
  var.list <-   var.list[!var.list=="WITH"]
  var.list <-   var.list[!var.list==""]
  var.list <- unique(var.list)
  #var.list <- paste(unique(unlist(var.list)), collapse=" ")
  #var.list <- strsplit(var.list, " ")[[1]]
  #var.list <-   var.list[!var.list==""]
  
  

  
  # var.list <- paste0("; ", model, " ;")
  # var.list<- gsub("\n", ";", var.list)
  # var.list <- paste(sapply(var.list, function(i) sub(".*BY *(.*?) *;.*", "\\1", i)), collapse=" ")
  # var.list <- strsplit(gsub(" BY | +|;", " ", var.list), " ")[[1]]
  # var.list <- var.list[!var.list ==""]
  
  d <- dat[c(group, var.list)]
  for(i in colnames(d)) d[,i] <- unclass(d[,i])
  rm(i)
  
  if(!is.numeric(d[,group])) {
    #d[,group] <- gsub(" ", "_", as.character( d[,group] )  )
    message("The group variable must be numeric!")
    
  }
  
  #require(MplusAutomation)
  #inp <- capture.output(prepareMplusData(d,  "mplus_temp.tab"))
  
  write.table(d, "mplus_temp.tab", quote=F, sep="\t", row.names=F, col.names=F, na=".")
  
  #var.list <- gsub("\\.", "_", var.list)
  
  list.of.groups = unique(as.matrix(d[,1]))
  ngroups = length(list.of.groups)
  
  inp <- c("DATA:","\n",
           "   file = 'mplus_temp.tab';", "\n",
           " VARIABLE:", "\n",
           "   names =", gsub("\\.", "_", group), " ", paste(gsub("\\.", "_", var.list), collapse="\n"), ";\n",
           "   missing = .;", "\n",
           ifelse(any(is.null(categorical)),
                  "\n",
                  paste("   categorical = ", paste(categorical, collapse = " "), ";\n")
           ),
           
           "   classes = c(", ngroups, ");\n",
           "   knownclass = c(", paste0(gsub("\\.", "_", group), " = ", list.of.groups, " \n    ", collapse=""),
           
           ");\n\n",
           
           "ANALYSIS:\n",
           "  type = mixture;\n",
           "  estimator = ml;\n",
           "  alignment =", kind = "", ";\n", 
           ifelse(any(is.null(categorical)),
                  "\n",  
                  "  algorithm = integration;\n\n"),
           
           "MODEL:\n",
           "  %OVERALL%\n",
           model, 
           "\n\n",
           
           "OUTPUT: align tech8 SVALUES;", 
           "\n\n",
           
           "SAVEDATA: ", "\n",
           "  RANKING = ranking.dat; "
           
  )
  
  
  
  inp["kind"]<-"FREE"
  cat(inp, file = "free.inp", sep="")
  message("Run free in Mplus.")
  trash <- system(paste(Mplus_com, "free.inp"))
  
  
  outFree <- paste(readLines("free.out"), collapse = "\n") 
  if(grepl("TO AVOID MISSPECIFICATION USE THE GROUP WITH VALUE", outFree)) {
    refGroup <- sub(".*TO AVOID MISSPECIFICATION USE THE GROUP WITH VALUE *(.*?) *AS THE BASELINE GROUP.*", "\\1", outFree)
  } else {
    
    free.tab.means <- sub(".*FACTOR MEAN COMPARISON AT THE 5% SIGNIFICANCE LEVEL IN DESCENDING ORDER *(.*?) *QUALITY OF NUMERICAL RESULTS.*", "\\1", outFree)
    refGroup <- as.character(read.table(text=sub(".*\n *(.*?) *\n\n\n\n\n.*", "\\1", free.tab.means))[3])
    
    
  }
  
  inp["kind"]<-paste0("FIXED(", refGroup, ")")
  cat(inp, file = "fixed.inp", sep="")
  message("Run fixed in Mplus.")
  trash <- system(paste(Mplus_com, "fixed.inp"))
  
  # Creating simulations
  if(!is.null(sim.samples)) {
    
    outFixed <- paste(readLines("fixed.out"), collapse = "\n") 
    
    
    stValues <- sub(".*MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES *(.*?) *\n\n\n\n.*", "\\1", outFixed)
    stValues <- gsub("%C#", "%g#", stValues)
    stValues <- gsub("c#", "g#", stValues)
    
    corrupt.code <- sub(".*%OVERALL% *(.*?) *%g#1%.*", "\\1", stValues)
    correction <-strsplit(corrupt.code, "\n")[[1]]
    correction <- correction[grep(" BY ",  correction)]
    correction <- gsub(";", "*1;", correction)
    
    stValues <- paste(paste(correction, collapse="\n"), "\n", substr(stValues, regexpr("%g#1%", stValues), nchar(stValues)))
    
    if(!any(is.null(categorical))) {
      g1 <- sub(".*%g#1% *(.*?) *%g#2%.*", "\\1", stValues)
      g1 <- strsplit(g1, "\n")[[1]]
      g1 <- g1[grep("\\[", g1)]
      g1 <- g1[grep("\\$", g1)]
      g1 <- sapply(g1 , function(x)   sub(" *\\[ *(.*?) *\\$.*", "\\1", x))
      gen.cat <- paste0(names(table(g1)), " (", table(g1), ")")
    }
    
    
    
    
    for(x in sim.samples) { 
      code <- c("MONTECARLO:",
                " NAMES = ", paste(gsub("\\.", "_", var.list), collapse = " "), ";\n",
                " ngroups = ", ngroups, ";\n", 
                " NOBSERVATIONS =", ngroups, "(", x, ");\n", 
                " NREPS =", sim.reps, ";\n\n",
                ifelse(any(is.null(categorical)),
                       "\n",  
                       paste(
                         " CATEGORICAL =", paste(categorical, collapse = " "), ";\n", 
                         " GENERATE = ", paste(gen.cat, collapse = " "),
                         ";\n\n"  )),
                
                
                
                "ANALYSIS:",
                " TYPE = MIXTURE;",
                " ESTIMATOR = ml;",
                " alignment = fixed;\n",
                ifelse(any(is.null(categorical)),
                       "\n",  
                       " algorithm = integration;\n\n"),
                
                "MODEL POPULATION:",
                " %OVERALL%\n",
                paste(stValues, collapse="\n"),
                "\nMODEL:",
                " %OVERALL%\n",
                paste(stValues, collapse="\n")
      )
      cat(code, sep="", file = paste0("sim", x , ".inp"))
    }
    
    for (x in sim.samples) {
      message("Run simulation", x, "in Mplus.\n")
      trash <- system(paste(Mplus_com, paste0("sim", x, ".inp")))
      
    }
  }
  
  # Return summaries
  
  if(summaries) {
    
    if(!is.null(sim.samples)) {
      
      otpt <- list(fixed= extractAlignment("fixed.out", silent = TRUE),
                   free = extractAlignment("free.out", silent = TRUE),
                   simulations = extractAlignmentSim(sapply(sim.samples, function(x) paste0("sim", x, ".out")), silent = TRUE)
      )
      cat("\n", "⎯⎯⎯⎯⎯⎯⎯⎯⎯ ", "Results of Free alignemnt", rep("⎯", getOption("width", 80)-20),  "\n", sep="")
      print(otpt$free$summary)
      
      cat("\n", "⎯⎯⎯⎯⎯⎯⎯⎯⎯ ", "Results of Fixed alignemnt", rep("⎯", getOption("width", 80)-20),  "\n", sep="") 
      print(otpt$fixed$summary)
      
      cat("\n", "⎯⎯⎯⎯⎯⎯⎯⎯⎯ ", "Results of simulations", rep("⎯", getOption("width", 80)-20),  "\n", sep="") 
      print(otpt$simulations)
      
      
      
    } else {
      otpt <- list(fixed = extractAlignment("fixed.out", silent = TRUE),
                   free =  extractAlignment("free.out", silent = TRUE))
      
      
      cat("\n", "⎯⎯⎯⎯⎯⎯⎯⎯⎯ ", "Results of Free alignemnt", rep("⎯", getOption("width", 80)-20),  "\n", sep="")
      print(otpt$free$summary)
      
      cat("\n", "⎯⎯⎯⎯⎯⎯⎯⎯⎯ ", "Results of Fixed alignemnt", rep("⎯", getOption("width", 80)-20),  "\n", sep="") 
      print(otpt$fixed$summary)
      
      
    }
    
  } else {
    message("Done running models. Refer to the free.out, fixed.out, ranking.dat and some sim###.out files.\nConsider  using `extractAlignment()` and `extractAlignmentSim()` to extract important parts.")
  }
  
  setwd(oldwd)
  if(summaries) invisible(otpt)
}
