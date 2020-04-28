#' Runs exact invariance tests in Mplus
#' 
#' @param model Character. Formula in Mplus format, e.g. "Factor1 BY item1 item2 item3 item4; item4 WITH item3;", see example.
#' @param group Character, name of the grouping variable.
#' @param dat Data frame containing data. 
#' @param categorical Character vector of variable names. Indicators that are binary or ordinal. Default is NULL.
#' @param filename Optional custom filename
#' @param Mplus_com  Sometimes you don't have a direct access to Mplus, so this argument specifies what to send to a system command line. Default value is "mplus".
#' 
#' @examples 
#' \dontrun{
#'    measurementInvarianceMplus(model = 
#'      "F1 BY F114 F115 F116 F117;
#'       F2 BY F118 F119  F120  F121 F122  F123;",
#'       group = "S003",
#'       data = wevs.subset)
#'       }
#' 
#' @export
measurementInvarianceMplus <- function(model, group, data, categorical=NULL, filename = NULL, Mplus_com = "mplus") {
  
  
  
  datfile = ifelse(is.null(filename), "mplus_temp.tab", paste0(filename, ".dat"))
  inpfile = ifelse(is.null(filename), "mplus_temp.inp", paste0(filename, ".inp"))
  outfile = ifelse(is.null(filename), "mplus_temp.out", paste0(filename, ".out"))
  
  
  var.list <- strsplit(model, "\n") [[1]]
  var.list <- gsub("!.*$", "",  var.list)
  var.list <- unlist(strsplit(var.list, ";|\n"))
  var.list <- gsub("\\s+", " ",  var.list)
  var.list <- var.list[!var.list==""]
  var.list <- var.list[!var.list==" "]
  var.list <- gsub("^.*((by)|(BY))", "", var.list)
  var.list <- unlist(strsplit(var.list, " "))
  var.list <- var.list[!var.list=="WITH"]
  var.list <- var.list[!var.list==""]
  var.list <- unique(var.list)
  
  if(any(!c(group, var.list) %in% colnames(data))) {
    print("The data doesn't contain the following variables:")
    print(c(group, var.list)[!c(group, var.list) %in% colnames(data)])
  }
  
  d <- data[c(group, var.list)]
  for(i in colnames(d)) d[,i] <- unclass(d[,i])
  rm(i)
  class(d)<-"data.frame"
  if(!is.numeric(d[,group])) {
    #d[,group] <- gsub(" ", "_", as.character( d[,group] )  )
    warning("The group variable must be numeric!")
    
    
  }
  
  write.table(d, datfile, quote=F, sep="\t", row.names=F, col.names=F, na=".")
  
  
  
  cat(paste("
DATA: file = '", datfile, "';
VARIABLE:
  NAMES = ", gsub("\\.", "_", group), paste(sub("\\.", "_", var.list), collapse = "\n          "), 
";\n  MISSING=.;\n",

ifelse(any(is.null(categorical)),
       "",
       paste("   categorical = ", paste(categorical, collapse = "\n"), ";\n")
),

"GROUPING IS",  gsub("\\.", "_", group), " (", paste(unique(d[,group]), collapse = "\n            "), ");\n",
"ANALYSIS:
  estimator = ", ifelse(any(is.null(categorical)),"MLR", "WLSMV"), ";
  model= configural ", ifelse(any(is.null(categorical)),"metric", ""), " scalar;

MODEL:\n", model), 
      
      file = inpfile)
  
  
  
  message("Run conventional MI test in Mplus.")
  trash <- system(paste(Mplus_com, inpfile))
  
  #outfile="/Users/maksimrudnev/Dropbox/STAT/political efficacy/mplus_temp.out"
  
# SUMMARIZE  
out <- paste(readLines(outfile), collapse = "\n")
  
blocks = sub(".*Invariance Testing *(.*?) *MODEL RESULTS FOR.*", "\\1", out)

#cat(blocks)

blocks = strsplit(blocks, "MODEL FIT INFORMATION FOR THE ")[[1]]

 blks <- lapply(blocks[2:length(blocks)], function(block2) {
    
      block2 = strsplit(block2, "\n")[[1]]
      
      rmsea = block2[which(block2=="RMSEA (Root Mean Square Error Of Approximation)")+2:4]
      rmsea = Reduce("rbind", strsplit(rmsea, "\\s\\s+"))[,-1]
      
      cfitli = block2[which(block2=="CFI/TLI")+2:3]
      cfitli = Reduce("rbind", strsplit(cfitli, "\\s\\s+"))[,c(2,3)]
      
      srmr = block2[which(block2=="SRMR (Standardized Root Mean Square Residual)")+2]
      srmr = Reduce("rbind", strsplit(srmr, "\\s\\s+"))[-1]
      
      chisq = block2[which(block2=="Chi-Square Test of Model Fit")+2:4]
      chisq = Reduce("rbind", strsplit(chisq, "\\s\\s+"))[,c(2,3)]
      
      list(rmsea=rmsea, cfitli=cfitli, srmr=srmr, chisq=chisq)
    
    })

 
tb.out = sapply(blks, function(x) 
            c(RMSEA=as.numeric(x$rmsea[x$rmsea[,1]=="Estimate",2]),
                 CFI=as.numeric(x$cfitli[x$cfitli[,1]=="CFI",2]),
                 SRMR=as.numeric(x$srmr[2])))
colnames(tb.out)<-unname(sapply(blocks[-1],   function(x) strsplit(x, "\n")[[1]][1]))
tb.out = cbind(tb.out,
  diff=as.matrix(apply(tb.out, 1, function(x) ifelse(length(x)==2, x[1]-x[2], c(x[1]-x[2], x[2]-x[3])))))

tb.out

}
