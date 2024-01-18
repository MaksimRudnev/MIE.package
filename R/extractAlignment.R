#' Extracts summaries from Mplus alignment output file
#' 
#' @param file filename of the Mplus alignment output file
#' @param nice.tables Logical. If tables should be send to RStudio Viewer pane using `kable` and `kableExtra` packages.
#' @param silent Logical. Used for debugging.
#' @param what Character vector. What exactly to extract. Possible values: "summary", "ranking", "comparisons", "contributions", "savedata".
#' @examples
#' \dontrun{  
#'    align.summ <- extractAlignment("fixed.out") 
#'    }
#' @return A list of summary tables.
#' @seealso \code{\link[MIE]{runAlignment}}  and \code{\link[MIE]{extractAlignmentSim}} 
#' @export
extractAlignment <-  function(file = "fixed.out", 
                              nice.tables = FALSE, 
                              silent = FALSE, 
                              what = c("summary", "ranking", "comparisons", "contributions", "savedata") ) {
  
  if(!file.exists(file)) warning("File not found")
  
  # Basic extraction function  
  extractBetween <- function(begin, end, string) {
    mapply(function(a, b) substr(string, a, b),
           gregexpr(begin, string)[[1]]+nchar(begin),
           gregexpr(end, string)[[1]]-1
    )  
  }
  
  
 

  # Read file
  b.vector = readLines(file)
  b.vector <- gsub("!.*", "", b.vector)
  b.string <-  paste(b.vector, collapse="\n")
 
  
  # Extract estimator
  sum.of.analysis.s <-  extractBetween( "SUMMARY OF ANALYSIS",  "Input data format  ", b.string)
  sum.of.analysis.s <-  strsplit(sum.of.analysis.s,"\n")[[1]]
  estimator <- sum.of.analysis.s[grep("Estimator", sum.of.analysis.s)]
  estimator <- sub("^(Estimator)\\s*", "", estimator)
  
  # TYPE of analysis
   type.analysis <- sub(".*ANALYSIS\\s*(.*?) *MODEL.*", "\\1", b.string, ignore.case = T)
   
   type.analysis.s <-  strsplit(type.analysis,"\n")[[1]]
   type.analysis <- type.analysis.s[grep("type", type.analysis.s, ignore.case = T)]
   type.analysis = ifelse(length(type.analysis)>0, 
                          sub("^(type)\\s*", "", type.analysis),
                          "MG")
  
  # Version of Mplus
 mplus.version <-  MIE:::MplusVersion(out.string = b.string) 

  
  # Var list
 names.begin.index = gregexpr("names", b.string, ignore.case = T)[[1]]
 names.end.index = gregexpr(";", b.string, ignore.case = T)
 names.end.index <- names.end.index[[1]][names.end.index[[1]] > names.begin.index][[1]]
 string.varnames = substr(b.string, names.begin.index + attr(names.begin.index, "match.length"), names.end.index-1)
 
  

  expand.varlist <- function(string.varnames) {
    
      var.list = scan(text=string.varnames, 
                      what="character", quiet = T)[-1]
      var.list <- var.list[!var.list %in% c("ARE", "=")]
      if(any(var.list == "ALL")) {
        
        var.list <- "ALL"
        
      } else {
      var.list <- unname(unlist(
                  sapply(var.list, 
                         function(x) 
                           if(grepl("-", x)) { 
        list.ends = strsplit(x, "-")[[1]]
        list.margins = as.numeric(gsub("\\D", "", list.ends))
        vars.at.the.end = gsub("\\d", "",list.ends)
        if(vars.at.the.end[[1]]==vars.at.the.end[[2]]) {
               paste0(vars.at.the.end[[1]], list.margins[[1]]:list.margins[[2]])
        } else {
              # message("Error expanding list of variables")
          return(list.ends)
        }
        
        } else {
          x
        })))}
      return(var.list)
      }
  
  var.list = expand.varlist(string.varnames)
  
  
  # usevariables
  names.begin.index = gregexpr("USEVARIABLES", b.string, ignore.case = T)[[1]]
  names.end.index = gregexpr(";", b.string, ignore.case = T)
  names.end.index <- names.end.index[[1]][names.end.index[[1]] > names.begin.index][[1]]

  string.used.vars = substr(b.string, names.begin.index + attr(names.begin.index, "match.length"), names.end.index-1)
  
  used.vars.list = expand.varlist(string.used.vars)
  
if(is.matrix(used.vars.list)) {
  used.vars.list = var.list[which(used.vars.list[1,1] == var.list):which(used.vars.list[2,1] == var.list)]
} else if(used.vars.list[[1]] == "ALL") {
  used.vars.list = var.list
} 
  
  
  if(grepl("categorical", b.string)) {
    
    categorical.var.string = sub(".*categorical\\s*(.*?) *;.*", "\\1", b.string, ignore.case = T)
    
    categorical.var.list = expand.varlist(categorical.var.string)
    
    if(is.matrix(categorical.var.list)) {
      categorical.var.list = var.list[which(categorical.var.list[1,1] == var.list):which(categorical.var.list[2,1] == var.list)]
    } else if(categorical.var.list[[1]] == "ALL") {
      categorical.var.list = used.vars.list
    } 
    
  } else {
    categorical.var.list = NULL
  }
  
  
 
  
  # Parameterization
  parameterization = ifelse(any(grepl("Parameterization", sum.of.analysis.s)),
                            trimws(gsub("Parameterization" ,"",sum.of.analysis.s[grepl("Parameterization", sum.of.analysis.s)])),
                            NA)

  if(!mplus.version %in% c("8.8", "8.9", "8.10"))
    warning("Mplus versions 8.7 and earlier are not officially supported,
    but trying to extract summary anyways.")
  
  # output list to be filled
  output <- list()
  
  # Extract mean comparison ######
  if(!mplus.version %in% c("8.9", "8.10") & 
     estimator!="BAYES") {
    
    mean.comparison <- extractBetween(
      "FACTOR MEAN/INTERCEPT COMPARISON AT THE 5% SIGNIFICANCE LEVEL IN DESCENDING ORDER", 
      "\n\n\n\n\n", b.string)
    
  } else if(mplus.version %in% c("8.9", "8.10") & estimator!="BAYES" ) {
    
    mean.comparison <- extractBetween(
      "FACTOR INTERCEPT COMPARISON AT THE 5% SIGNIFICANCE LEVEL IN DESCENDING ORDER", 
      "\n\n\n\n\n", b.string)
    
   } else {
    
    mean.comparison <- extractBetween(
      "FACTOR MEAN COMPARISON AT THE 5% SIGNIFICANCE LEVEL IN DESCENDING ORDER", 
      
      "\n\n\n\n\n", b.string)
  }
  
  
  mean.comparison<-mean.comparison[!mean.comparison==""]
  mean.comparison<- strsplit(mean.comparison,"(Results for Factor)")[[1]][-1]
  
  mean.comparison<- gsub("\n\n$", "", mean.comparison)
  
  names(mean.comparison) <- sapply(mean.comparison, function(x)  substr(x, 2, regexpr("\n", x)-1))
  
  mean.comp <- lapply(mean.comparison, function(x) {
    read.fwf(file=textConnection(x), skip=4, widths = c(7, 10, 10, 12, 1000),
             col.names = c("Ranking", "Latent class", "Group value", "Factor mean", "Groups With Significantly Smaller Factor Mean")) 
    
  })
  
mean.comp <- lapply(mean.comp, function(x) {
  x.new <- x[which(!is.na(x$Ranking)),]
  x.new$Groups.With.Significantly.Smaller.Factor.Mean <- 
      c(sapply(1:(nrow(x.new)-1), function(y) {
       index.of.groups <- which(!is.na(x$Ranking))[[y]]:(which(!is.na(x$Ranking))[[y+1]]-1)
       paste(x$Groups.With.Significantly.Smaller.Factor.Mean[index.of.groups], collapse="")
      }), "NA")
 

  Groups.With.Significantly.Smaller.Factor.Mean <- 
    strsplit(x.new$Groups.With.Significantly.Smaller.Factor.Mean, " ")
  Groups.With.Significantly.Smaller.Factor.Mean <- 
    lapply(Groups.With.Significantly.Smaller.Factor.Mean, function(b) b[b!="" & b!="NA"])
  x.new$NGroups.With.Smaller.Factor.Mean <- 
    sapply(Groups.With.Significantly.Smaller.Factor.Mean, length)
  x.new$Groups.With.Significantly.Smaller.Factor.Mean <- sapply(Groups.With.Significantly.Smaller.Factor.Mean, paste, collapse=", ")
  x.new
  })
  
  output[["mean.comparison"]] <- mean.comp
  
  
  
  if("ranking" %in% what) {
  # Extract ranking table #####
  if(grepl("Factor Mean Ranking Tables", b.string)) {
    
    rankingFile <- sub(".*Factor Mean Ranking Tables *(.*?) *Save format.*", "\\1", b.string)
    rankingFile <- gsub("Save file|\n| ", "", rankingFile)
    rankingFileName <- paste0(substr(file, 1, regexpr("/.*$", file)), rankingFile)
    if(file.exists(rankingFileName)) {
     ranking.tabs <- paste(readLines(rankingFileName), collapse="\n")
     ranking.tabs <- strsplit(ranking.tabs, "Ranking table for ")[[1]][-1]
     names(ranking.tabs)<- sapply(ranking.tabs, function(x) substr(x, 1, regexpr("\n", x)-1))
     ranking.tab <- lapply(ranking.tabs, function(x) {
       rt <- read.csv(text=x, skip=2, row.names = 1)
       colnames(rt)<- rownames(rt)
       rt[,-ncol(rt)]
     })
    
    
    output[["ranking.table"]] <-ranking.tab
    } else {
      warning("Could not find ranking file.")
    }
  }
  }
  
  if(any(c("summary", "comparisons") %in% what)) {
  # Extract pairwise comparisons #####
  
  # extract alignment part
  align.outp <- extractBetween("ALIGNMENT OUTPUT", "Average Invariance index", b.string)
  #separate intercepts/thresholds and loadings
  
  # correction for version 8.10
  align.outp <- gsub("Approximate Invariance Was Not Found For This Parameter.", 
       "Approximate Measurement Invariance Holds For Groups:", align.outp)

  if(mplus.version %in% c("8.8", "8.9", "8.10") & estimator!="BAYES") { 
    
      align.outp1 <- strsplit(align.outp, "\n\n\n Loadings for ", fixed=T)[[1]]
     
      
      # intercepts.thresholds.comparison
      al.i.th <- align.outp1[[1]]
      al.i.th <- sub("\n\nINVARIANCE ANALYSIS\n\n Intercepts/Thresholds\n ", " ",  al.i.th)
      al.pw.i <- strsplit(al.i.th, "(Intercept for )|(Threshold )")[[1]]
    
      al.loads <-align.outp1[-1]
      names(al.loads)<-names(mean.comp)
      al.loads <-lapply(seq_along(al.loads), 
                       function(x) 
                        gsub("Loadings for ",
                             paste("Loadings for ", names(al.loads)[[x]], "by "),
                             al.loads[[x]][[1]])
                       )
      
      al.loads <-lapply(al.loads,   function(x) sub("^\\S*", "", x) )
      
      # strsplit(al.loads[[2]], "\n")[[1]][1:15]
      # strsplit(sub("^\\S*", "", x), "\n")[[1]][1:15]
      
      al.pw.l <- unlist(strsplit(paste(al.loads, collapse="\n"), "Loadings for "))
      #str(al.pw.l)
      
      #names(al.pw.l)
      # strsplit(al.pw.l[[9]], "\n")[[1]][1:15]
      
  } else {
    
  align.outp <- strsplit(align.outp, "Loadings\n")[[1]]
  align.outp <- sub("\n\nINVARIANCE ANALYSIS\n\n Intercepts/Thresholds\n ", " ", align.outp)
  align.outp <- sub("\n\nINVARIANCE ANALYSIS\n\n Intercepts\n ", " ", align.outp)
  al.pw.i <- strsplit(align.outp[1], "Intercept for |Threshold ")[[1]]
  al.pw.l <- unlist(strsplit(align.outp[2:length(align.outp)], "Loadings for "))
  }
  
  # drop first header
if(estimator=="MLR") {
  # this is required to name the fit contribution
  loading.names.by.factor <- {
    a <- strsplit(align.outp[2:length(align.outp)], "Loadings for ")
    lapply(a, function(y) {
      m <- sapply(y, function(b) substr( b, 1, regexpr("\n", b)-1))
      m <- m[!m==""]
      unname(m)
    })
  }}
  
  
  al.pw.i <-al.pw.i[!al.pw.i %in% c(" ",NA,"")]
  al.pw.l <-al.pw.l[!al.pw.l %in% c(" ",NA,"", "\n ")]
  
  al.pw.i.names1 <- sapply(1:length(al.pw.i), function(b) substr( al.pw.i[b], 1, regexpr("\n", al.pw.i[b])-1))
  al.pw.i.names2 <- sapply(al.pw.i.names1, function(nmz) ifelse( grepl("\\$", nmz), "Threshold", "Intercept"))
  
  al.pw <- c(paste(al.pw.i.names2, al.pw.i), paste("Loadings", al.pw.l))
  al.pw.names <- sapply(1:length(al.pw), function(b) substr( al.pw[b], 1, regexpr("\n", al.pw[b])-1))
  var.list.model <- gsub("\\$.*$", "", gsub("^.* ", "", al.pw.names))
  var.list.model <- var.list.model[!duplicated(var.list.model)]
  
  if(estimator == "BAYES") {
    al.pw <- gsub("Approximate Invariance \\(Noninvariance\\) Holds For Groups:", 
                  "Approximate Measurement Invariance Holds For Groups:", 
                  al.pw)
  }
  al.pw.list <-strsplit(al.pw, " Approximate Measurement Invariance Holds For Groups:")
  
  names(al.pw.list)<- al.pw.names
  
  align.outp <- lapply(setNames(nm=al.pw.names), function(x) { 
    #print(x)
    #x=al.pw.names[[5]]
    x <- al.pw.list[[x]]
    
    
    pairwise.tab <- read.table(text=x[1], stringsAsFactors = FALSE, skip=2)
    if(estimator=="BAYES") {
      colnames(pairwise.tab) <- c("Group1", "Group2", "Est_in_G1", "Est_in_G2", "Difference", "SE", "P_value",
                                  "lowerCI", "upperCI")
      
    } else {
      colnames(pairwise.tab) <- c("Group1", "Group2", "Est_in_G1", "Est_in_G2", "Difference", "SE", "P_value" )
    }
    
    invariant.groups <- substr(x[2], 2, regexpr("Weighted Average Value Across Invariant Groups:",x[2])-1)
    invariant.groups <- unlist(strsplit(readLines(textConnection(invariant.groups)), " "))
    invariant.groups <- invariant.groups[!invariant.groups==""]
    
    all.groups <- unique(unlist(pairwise.tab[,1:2]))
    non.invariant.groups <- all.groups[!all.groups %in% invariant.groups]
    
    
    AlignedPar = as.numeric(sub(".*Weighted Average Value Across Invariant Groups: *(.*?) *\n.*", "\\1", x[2]))
    R2 = as.numeric(sub(".*R-square/Explained variance/Invariance index: *(.*?) *\n.*", "\\1", x[2]))
    
    if(grepl("Invariant Group Values", x[2])) {
      
    
    inv.comparison <- substr(x[2], regexpr("Invariant Group Values, Difference to Average and Significance", x[2]), nchar(x[2]))
    inv.comparison <- read.table(text = inv.comparison, skip=2)
    } else { 
      inv.comparison <- t(rep(NA,5))
      }
    
    if(estimator=="BAYES") {
      colnames(inv.comparison) <- c("Group","Value", "Difference","SE","P.value", "lowerCI", "upperCI")
    } else {
      colnames(inv.comparison) <- c("Group","Value", "Difference","SE","P.value")
    }
    
    list(
      "Pairwise comparison" = pairwise.tab,
      "Aligned parameter" = AlignedPar,
      "R2" = R2,
      "Comparison of aligned pars to average" = inv.comparison,
      "N_noninvariant" = length(non.invariant.groups),
      "N_invariant" = length(invariant.groups),
      "Invariant groups" = invariant.groups,
      "Non-invariant groups" = non.invariant.groups
    )
    
  })

  output[["alignment.output"]] <- align.outp
  }
  
  # Summaries
  
  
  # invariant.groups.pars.tab
  
  all.groups <- c(align.outp[[1]]$`Invariant groups`, align.outp[[1]]$`Non-invariant groups`)
  non.invariant.groups.pars.tab <- sapply(align.outp, function(x) { 
    
    temp<-rep("", length(all.groups))
    temp[!all.groups %in% x[["Invariant groups"]]]<-"X"
    temp
  })
  rownames(non.invariant.groups.pars.tab)<-all.groups
  
  output[["non.invariant.pars"]] <- t(non.invariant.groups.pars.tab)
  
  
  
  # Extract summaries: R2, aligned parameters, list of invariant and non-invariant groups #####
  
  summ <- 
    lapply(align.outp, function(x)  data.frame(AlignedParameter = x[["Aligned parameter"]],
                                        R2 = x[["R2"]], 
                                        N_invariant = x[["N_invariant"]],
                                        N_noninvariant = x[["N_noninvariant"]],
                                        invariant.gr = paste(x[["Invariant groups"]], collapse=" "), 
                                        non.invar.gr = paste(x[["Non-invariant groups"]], collapse = " ")))
  
  summ1 <- Reduce("rbind", summ)
  rownames(summ1) <- names(summ)
  
  output[["summary"]] <- summ1
  
  
  # Fit contribution ######
  if(estimator %in% c("MLR", "ML") & "contributions" %in% what) {
    
    if(grepl("TECHNICAL 8 OUTPUT", b.string))  {
      
      tech8 <-  substr(b.string, regexpr("TECHNICAL 8 OUTPUT", b.string), nchar(b.string))
      tech8 <-  sub("TECHNICAL 8 OUTPUT\n\n\n", "", tech8)
      tech8 <-  substr(tech8, 1, regexpr("(\n\n\n)", tech8))
      tech8.align <- strsplit(tech8, "ALIGNMENT RESULTS FOR ")[[1]][-1]
      
      # turn string output to a table
      f.names <- sapply(tech8.align, function(x) substr(x, 1, regexpr("\n", x)-1))
      #
      #fit.contrib <- lapply(1:length(tech8.align), function(x) {
      
      f.contrib.l <-  sub(".* Fit Function Loadings Contribution By Variable *(.*?) *Fit Function Loadings Contribution By Group.*", "\\1", tech8.align)
      f.contrib.i <-  sub(".* Fit Function Intercepts Contribution By Variable *(.*?) *Fit Function Intercepts Contribution By Group.*", "\\1", tech8.align)
      
    if(length(f.contrib.i)>0) {
      
      contrib.i.tab <- read.table(text=f.contrib.i)
      contrib.l.tab <- read.table(text=f.contrib.l)
      contrib <- unname(c(unlist(contrib.i.tab), unlist(contrib.l.tab)))
    
      # names(contrib) <- c( 
      #   paste("Intercept", unlist(loading.names.by.factor[[x]])),
      #   paste("Loadings", loading.names.by.factor[[x]])
      # )
      
      nmz.th.int <- paste(al.pw.i.names2, al.pw.i.names1)
      nmz.loadings <- al.pw.names[!al.pw.names %in% nmz.th.int]
      
      
      # because Mplus prints duplicates for factor loadings when thresholds are present, we need to drop duplicates
      if(nrow(contrib.l.tab)>length(nmz.loadings)) {
        # delete duplicated contributions of factor loadings
        contrib.l.tab <- contrib.l.tab[!duplicated(gsub("\\$.*", "", nmz.th.int)),]
        contrib <- unname(c(unlist(contrib.i.tab), unlist(contrib.l.tab)))
        
      }
      
      names(contrib) <- c(nmz.th.int, nmz.loadings)
    
      fit.contrib = data.frame(Fit.contribution = contrib, 
                               #             Factor = rep(f.names[1], length(contrib)), 
                               #             row.names = names(contrib),
                               stringsAsFactors = F)
      
      # })
      
      # f.names <- sapply(tech8.align, function(x) substr(x, 1, regexpr("\n", x)-1))
      # 
      # fit <- data.frame(Fit.contribution = unlist(fit.contrib),
      #                   Factor = rep(unname(f.names), each=length(f.contrib.l)),
      #                   stringsAsFactors = FALSE)
      # 
      # fit <- aggregate(fit,  list(names(unlist(fit.contrib))), function(x) x[1])
      # rownames(fit) <- fit$Group.1
      
      #fit.contrib <- Reduce("rbind", fit.contrib)
      output$summary <- merge(output[["summary"]], 
                              fit.contrib, 
                              by = "row.names")
      rownames(output$summary)<-output$summary$Row.names
      output$summary <- output$summary[,-1]
    }
    } else {
      warning("Extracting fit contributions failed. Please check if Tech8 is included in the input file.")
    } 
  }
  
  # Savedata ####
  if( "savedata" %in% what) {
    
    if(!grepl("SAVEDATA INFORMATION", b.string)) {
      
      warning("Couldn't find the Savedata information.")
      
    } else {
  str.out <- strsplit(b.string, "\n")[[1]]
  save.dat.str <- str.out[ 
    grep("SAVEDATA INFORMATION", str.out):
      grep("Save missing symbol", str.out)]
  
  filename <- trimws(save.dat.str[grep("Save file$", save.dat.str)+1])
  fileformat <- save.dat.str[ 
    grep("Order and format of variables", save.dat.str):
      grep("Save file format", save.dat.str)]
  
  fileformat <- fileformat[3:(length(fileformat)-2)]
  fileformat <- read.table(text = paste(fileformat, sep = "\n"), col.names = c("varname", "type"))
  saveddata <- read.fwf(filename, 
                        widths = as.numeric(gsub("F|I|\\.\\d", "", fileformat$type)), 
                        col.names = fileformat$varname, 
                        strip.white = T, 
                        na.strings = trimws(gsub("Save missing symbol", "", 
                                                 save.dat.str[grepl("Save missing symbol", 
                                                                    save.dat.str)])),
                        colClasses = "numeric")
  
  
  output$savedata <- saveddata
  }}
  
  # Final message ####
  sum.noninv = sum(output$summary[,"N_noninvariant"])
  n.params.total = sum(output$summary[,c("N_invariant", "N_noninvariant")])
  
  final.message = paste("There are", sum.noninv, 
        "non-invariant parameters out of",
        n.params.total, 
        "which is", 
        scales::percent(sum.noninv/n.params.total),
        ".\n This is", ifelse(sum.noninv/n.params.total < .25,
                              "smaller than the recommended cutoff of 25%,\n thereby the results support approximate invariance. \nRunning simuations is still recommended.\n",
                              "greater than the recommended cutoff of 25%,\nthereby the results do NOT support approximate invariance.\nRunning simulations is recommended for further exploration.\n")
  )
  
  # Adding extra info
  output$extra <- list(estimator = estimator,
                      mplus.version = mplus.version,
                      final.message = final.message,
                      var.list = used.vars.list, # var.list.model, # some names could be shrunk
                      categorical.var.list  = categorical.var.list,
                      parameterization = parameterization,
                      type = type.analysis,
                      string = b.string)
  
  
  # Printing output
  if(!silent) {
    print(output$summary, row.names=T)
    cat("\n",final.message)
  }
  
  if(nice.tables && !silent) {
    nice.tab1 <- knitr::kable(output$non.invariant.pars, format = "html")
   trash <-  capture.output(kableExtra::kable_styling(nice.tab1, bootstrap_options=c("striped", "bordered"), position = "left", font_size = 12))
  }
  
  invisible(output)
}
