#' Extracts summaries from Mplus alignment output file
#' 
#' @param file filename of the Mplus alignment output file
#' @param nice.tables Logical. If tables should be send to RStudio Viewer pane using `kable` and `kableExtra` packages.
#' @param silent Logical. Used for debugging.
#' @examples
#' \dontrun{  
#'    align.summ <- extractAlignment("fixed.out") 
#'    }
#' @return A list of summary tables.
#' @seealso \code{\link[MIE]{runAlignment}}  and \code{\link[MIE]{extractAlignmentSim}} 
#' @export
extractAlignment <-  function(file = "fixed.out", nice.tables = FALSE, silent = FALSE ) {
  
  
  # Basic extraction function  
  extractBetween <- function(begin, end, string) {
    mapply(function(a, b) substr(string, a, b),
           gregexpr(begin, string)[[1]]+nchar(begin),
           gregexpr(end, string)[[1]]-1
    )  
  }
  
  
  
  # Read file
  b.string<-  paste(readLines(file), collapse="\n")
  
  
  # Extract pairwise comparisons ########
  
  # extract alignment part
  align.outp <- extractBetween("ALIGNMENT OUTPUT", "Average Invariance index", b.string)
  #separate intercepts/thresholds and loadings
  align.outp <-strsplit(align.outp, "Loadings\n")[[1]]
  # drop first header
  align.outp <- sub("\n\nINVARIANCE ANALYSIS\n\n Intercepts/Thresholds\n ", " ", align.outp)
  
  al.pw.i <- strsplit(align.outp[1], "Intercept for |Threshold ")[[1]]
  al.pw.l <- unlist(strsplit(align.outp[2:length(align.outp)], "Loadings for "))
  
  # this is reqired to name the fit contribution
  loading.names.by.factor <- {
    a <- strsplit(align.outp[2:length(align.outp)], "Loadings for ")
    lapply(a, function(y) {
      m <- sapply(y, function(b) substr( b, 1, regexpr("\n", b)-1))
      m <- m[!m==""]
      unname(m)
    })
  }
  
  
  al.pw.i <-al.pw.i[!al.pw.i %in% c(" ",NA,"")]
  al.pw.l <-al.pw.l[!al.pw.l %in% c(" ",NA,"")]
  
  al.pw.i.names1 <- sapply(1:length(al.pw.i), function(b) substr( al.pw.i[b], 1, regexpr("\n", al.pw.i[b])-1))
  al.pw.i.names2 <- sapply(al.pw.i.names1, function(nmz) ifelse( grepl("\\$", nmz), "Threshold", "Intercept"))
  
  al.pw <- c(paste(al.pw.i.names2, al.pw.i), paste("Loadings", al.pw.l))
  al.pw.names <- sapply(1:length(al.pw), function(b) substr( al.pw[b], 1, regexpr("\n", al.pw[b])-1))
  
  
  
  al.pw.list <-strsplit(al.pw, " Approximate Measurement Invariance Holds For Groups:")
  names(al.pw.list)<- al.pw.names
  
  align.outp <- lapply(al.pw.list, function(x) { 
    
    pairwise.tab <- read.table(text=x[1], stringsAsFactors = FALSE, skip=2)
    colnames(pairwise.tab) <- c("Group1", "Group2", "Est_in_G1", "Est_in_G2", "Difference", "SE", "P_value" )
    
    invariant.groups <- substr(x[2], 2, regexpr("Weighted Average Value Across Invariant Groups:",x[2])-1)
    invariant.groups <- unlist(strsplit(readLines(textConnection(invariant.groups)), " "))
    invariant.groups <- invariant.groups[!invariant.groups==""]
    
    all.groups <- unique(unlist(pairwise.tab[,1:2]))
    non.invariant.groups <- all.groups[!all.groups %in% invariant.groups]
    
    
    AlignedPar = as.numeric(sub(".*Weighted Average Value Across Invariant Groups: *(.*?) *\n.*", "\\1", x[2]))
    R2 = as.numeric(sub(".*R-square/Explained variance/Invariance index: *(.*?) *\n.*", "\\1", x[2]))
    
    inv.comparison <- substr(x[2], regexpr("Invariant Group Values, Difference to Average and Significance", x[2]), nchar(x[2]))
    inv.comparison <- read.table(text = inv.comparison, skip=2)
    colnames(inv.comparison) <- c("Group","Value", "Difference","SE","P.value")
    
    
    list(
      "Pairwise comparison" = pairwise.tab,
      "Aligned parameter" = AlignedPar,
      "R2" = R2,
      "Comparison of aligned" = inv.comparison,
      "Invariant groups" = invariant.groups,
      "Non-invariant groups" = non.invariant.groups
    )
    
  })
  output <- list()
  output[["alignment.output"]] <- align.outp
  
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
    t(sapply(align.outp, function(x)  c(AlignedParameter = x[["Aligned parameter"]],
                                        R2 = x[["R2"]], 
                                        invariant.gr = paste(x[["Invariant groups"]], collapse=" "), 
                                        non.invar.gr = paste(x[["Non-invariant groups"]], collapse = " "))))
  
  
  
  output[["summary"]] <- summ
  
  
  # Fit contribution ######
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
       contrib.l.tab <- contrib.l.tab[!duplicated(contrib.l.tab),]
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
    output$summary <- merge(summ, fit.contrib, 
          by = "row.names")
    
    output$summary <- output$summary #[order(output$summary$Factor),]
  }
  
  
  # Extract mean comparison ######
  
  mean.comparison <- extractBetween("FACTOR MEAN COMPARISON AT THE 5% SIGNIFICANCE LEVEL IN DESCENDING ORDER", "\n\n\n\n\n", b.string)
  mean.comparison<-mean.comparison[!mean.comparison==""]
  mean.comparison<- strsplit(mean.comparison,"(Results for Factor)")[[1]][-1]
  
  mean.comparison<- gsub("\n\n$", "", mean.comparison)
  
  names(mean.comparison) <- sapply(mean.comparison, function(x)  substr(x, 2, regexpr("\n", x)-1))
  
  mean.comp <- lapply(mean.comparison, function(x) {
    read.fwf(file=textConnection(x), skip=4, widths = c(7, 10, 10, 12, 1000),
             col.names = c("Ranking", "Latent class", "Group value", "Factor mean", "Groups With Significantly Smaller Factor Mean")) 
    
  })
  
  output[["mean.comparison"]] <- mean.comp
  
  
  
  
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
      message("Could not find ranking file.")
    }
    }
  
  
  if(!silent) print(output$summary, row.names=FALSE)
  if(nice.tables && !silent) {
    nice.tab1 <- knitr::kable(output$non.invariant.pars, format = "html")
   trash <-  capture.output(kableExtra::kable_styling(nice.tab1, bootstrap_options=c("striped", "bordered"), position = "left", font_size = 12))
  }
  invisible(output)
}
