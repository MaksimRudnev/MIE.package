#' Extracts summaries of alignment simulations from Mplus output file
#' 
#' @param sim.outputs a character vector of file names containing simulation results of alignment ran in Mplus. 
#' @param silent Logical. Used for debugging.
#' @details Best used as part of \code{\link[MIE]{runAlignment}}
#' 
#' @examples 
#' \dontrun{  
#'      align.sim.summ <- extractAlignmentSim (c("sim500.out", "sim100.out", "sim1000.out")) 
#'      }
#' 
#' @return Invisibly returns a summary table.
#' @seealso \code{\link[MIE]{runAlignment}}  and \code{\link[MIE]{extractAlignment}} 
#' @export
extractAlignmentSim <- function(sim.outputs = c("sim500.out", "sim100.out", "sim1000.out"), silent = FALSE) {
  
  
  otp<- lapply(sim.outputs, function(x) {
    f <- paste(readLines(x), collapse="\n")
    mplus.version <- sub(".*Mplus VERSION *(.*?) *\n.*", "\\1", f)
    mplus.version <- trimws(gsub("(\\(.*\\))", "", mplus.version))
    
    if(!grepl("CORRELATIONS AND MEAN SQUARE ERROR OF POPULATION AND ESTIMATE VALUES", f)) {
      
      #stop("Correlations output is not found.")
      warning("Correlations output is not found. Calculating correlations manually")
      
      # extract parameters and correlate averages
      partable <- MIE:::extractParameters(string = f)
      partable.means <- partable[partable$L2 == "Means",c("Population", "parameter", "ESTIMATES_Average")]
      
      
   corr.of.averages =   sapply(setNames(nm=unique(partable.means$parameter)), function(latent.factor)
        cor(partable.means$ESTIMATES_Average[partable.means$parameter == latent.factor], 
            partable.means$Population[partable.means$parameter == latent.factor]))
   
   corr.of.averages.N = length(unique(partable$L1))
      
      
      # extract each replication mean, correlate to population value, and then average
      
      results.file <- sub(".*RESULTS SAVING INFORMATION(.*?) *Save file format.*", "\\1", f)
      results.file <- trimws(sub(".*Save file\n(.*?) *\n\n*", "\\1", results.file))
     if(file.exists(results.file)) {
       results.string <- readLines(results.file)
      } else {
        warning("Results file was not found.")
      }
      
      replication.id.index <- which(!grepl("\\s", results.string))
      parameter.vectors <- 
        lapply(1:length(replication.id.index), function(i) {
                 begin = replication.id.index[[i]] + 1
                 end = ifelse(i == length(replication.id.index), 
                              length(results.string), 
                              replication.id.index[[i + 1]]-1)
                 paste(results.string[begin:end], collapse = " ")
               })
        
    names(parameter.vectors) = results.string[replication.id.index]
  
     # find the location of ALPHA matrix in the technical output
    tech1.output <- substr(f, regexec("TECHNICAL OUTPUT", f)[[1]]+16, regexec("STARTING VALUES", f)[[1]])
    
    param.spec.output <- strsplit(tech1.output, "PARAMETER SPECIFICATION FOR")[[1]][-1]
    names(param.spec.output) <- unname(trimws(sapply(param.spec.output, function(x) strsplit(x, "\n\n\n")[[1]][[1]])))
    
   alpha.output <- lapply(param.spec.output, function(x) {
     alpha.string = sub(".*ALPHA(.*?)\n\n.*", "\\1", x)
     alpha.table = read.table(text = alpha.string, header = T)
     alpha.table[-1,]
   })
   
   alpha.index <- reshape2::melt(alpha.output, level = "group", id.vars = NULL) 
      
   means.for.each.replication = 
     lapply(parameter.vectors, function(x) {
       suppressMessages({
        param.vector = scan(text = parameter.vectors[[1]], what = numeric())
       })
        cbind(alpha.index, est = param.vector[as.numeric(alpha.index$value)])
     })
   
   pop.values = partable[partable$L2=="Means", c("L1", "parameter", "Population")]
   pop.values$L1 <- gsub("Group ", "", pop.values$L1)
   
   
   merged.dats = lapply(means.for.each.replication, function(x) {
     merge(x, pop.values, by.x = c("Lgroup", "variable"), by.y = c("L1", "parameter"))
   })
   
   cors.df = sapply(merged.dats, function(x) 
     sapply(setNames(nm=unique(x$variable)), function(latent.factor)
       cor(x$est[x$variable == latent.factor], x$Population[x$variable == latent.factor])))
  
    
    # output
    list(
      "Correlations and mean square error of population and estimate values" =
        list("Correlations Mean Average"=apply(cors.df, 1, mean),
                     "Correlations Variance Average"=NA, 
                     "Correlations Mean SD" = apply(cors.df, 1, sd), 
                     "Correlations Variance SD" = NA, 
                     "MSE Mean Average" = NA, 
                     "MSE Variance Average"= NA, 
                     "MSE Mean SD" = NA, 
                     "MSE Variance SD" = NA,
                     N = ncol(cors.df)),
        
      "Correlation and mean square error of the average estimates" = 
        list("Correlation of average means with true" = corr.of.averages, 
          "Correlation of average variances with true" = NA, 
          "MSE of average means with true", 
          "MSE of average variances with true",
          N = corr.of.averages.N
        )
         
        
      )
      
      
    } else {
    
    if(mplus.version == "8.10") { 
      cor.tab1 =  sub(".*CORRELATIONS AND MEAN SQUARE ERROR OF POPULATION AND ESTIMATE VALUES(.*?)\n\n\n.*", "\\1", f)
      cor.tab2 =  sub(".*CORRELATION AND MEAN SQUARE ERROR OF THE AVERAGE ESTIMATES(.*?)\n\n\n.*", "\\1", f)
      cor.tabs = list(cor.tab1, cor.tab2)
    
    } else {
    
      cor.tabs <- sub(".*CORRELATIONS AND MEAN SQUARE ERROR OF POPULATION AND ESTIMATE VALUES *(.*?)*\n\n\n\n.*", "\\1", f)
      cor.tabs <- strsplit(cor.tabs, "CORRELATION AND MEAN SQUARE ERROR OF THE AVERAGE ESTIMATES")[[1]]
      
      }
      
    cor.tabs1 <- strsplit(cor.tabs[[1]], "\n")[[1]][-c(1:4)]
    cor.tabs1 <- cor.tabs1[!cor.tabs1==""]
    
    dt1 <-sapply(seq(1, length(cor.tabs1), by=3), function(f.id) unlist(read.table(text=cor.tabs1[f.id+1:2])[-1]))
    colnames(dt1) = gsub(" ", "", cor.tabs1[seq(1, length(cor.tabs1), by=3)])
    row.names(dt1)<-c("Correlations Mean Average", "Correlations Variance Average", "Correlations Mean SD", "Correlations Variance SD",
                      "MSE Mean Average", "MSE Variance Average", "MSE Mean SD", "MSE Variance SD")
    
    
    cor.tabs2 <- strsplit(cor.tabs[[2]], "\n")[[1]]
    cor.tabs2 <- cor.tabs2[!cor.tabs2==""]
    dt2 <- read.table(text=gsub(".*(Mean|Variance) *", "", cor.tabs2))
    dt2 <- sapply(seq(1, nrow(dt2), by=2), function(y) unlist(dt2[y:(y+1),]) )
    colnames(dt2) <- gsub(" *", "", gsub("(Mean|Variance).* *", "",   cor.tabs2)) [seq(1, length(cor.tabs2), by=2)]
    rownames(dt2) <- c("Correlation of average means with true", "Correlation of average variances with true",
                       "MSE of average means with true","MSE of average variances with true")
    
    list('Correlations and mean square error of population and estimate values' = dt1,
         'Correlation and mean square error of the average estimates' = dt2)
    
    
  }})
  
  names(otp)<- sim.outputs
  
  
  if(!silent) {
    
    
    for(i in colnames(otp[[1]][[1]])  ) {
      cat("\n", "⎯⎯⎯⎯⎯⎯⎯⎯⎯ ", i," ", rep("⎯", getOption("width", 80)-nchar(i)-2),  "\n", sep="") 
      print(sapply(otp, function(x) x[[1]][,i] ), digits = 2)
      print(sapply(otp, function(x) x[[2]][,i] ), digits = 2)
    }
    }
  
  invisible(otp)
  
  
}


extractParameters <- function(file = NULL, string = NULL) {
  
  if(!is.null(file)) {
    if(file.exists(file))
        string <- paste(readLines(file), collapse="\n")
       else
         stop("File not found.")
    }

# Reading parameter table
all.pars.string <- sub(".*MODEL RESULTS *(.*?) *QUALITY OF NUMERICAL RESULTS.*", "\\1", string)
all.pars.split <- strsplit(all.pars.string, "\n")[[1]]
#all.pars.split <- gsub("\\s+", " ", all.pars.split)
all.pars.split <- all.pars.split[nchar(all.pars.split)!=0]

# identifying sections
full.string.length = max(nchar(all.pars.split))
headers.ind = which(nchar(all.pars.split)!=full.string.length)
headers1.ind =  headers.ind[c(headers.ind[-1] - headers.ind[-length(headers.ind)] == 1, F)]
headers2.ind = headers.ind[!headers.ind %in% headers1.ind]

# trimming space from sections names
all.pars.split[headers1.ind] <- trimws(gsub("\\s+", " ", all.pars.split[headers1.ind]))
all.pars.split[headers2.ind] <- trimws(gsub("\\s+", " ", all.pars.split[headers2.ind]))


# detecting the varnames (header of the partable)
varnames.string <- sub(".*\n\n(.*?) *\n\n.*", "\\1", all.pars.string)
if(grepl("\n", varnames.string)) {
  twoliner = read.fwf(textConnection(varnames.string), 
                      widths = c(17, 11, 11, 11, 11, 11, 6, 11))
  varnames.vector <- trimws(paste(twoliner[1,], twoliner[2,]))
  varnames.vector[[1]]<-"parameter"
} else {
  varnames.vector <- strsplit(varnames.string, split = "\\s+{2}")[[1]]
  varnames.vector[[1]]<-"parameter"
}

varnames.vector <- gsub("\\s+", "_", varnames.vector)

# making a list of parameters
all.pars.list <- list()
for(i in 1:length(headers1.ind)) {
  
  begin.section = headers1.ind[i]+1
  end.section = ifelse(i == length(headers1.ind),
                       length(all.pars.split),
                       headers1.ind[i+1]-1)
  
  section = all.pars.split[begin.section:end.section]
  
  headers2.in.section = headers2.ind[headers2.ind >= begin.section & 
                                       headers2.ind <= end.section]
  
  subsection <- list()
  for(k in 1:length(headers2.in.section)) {
    
    begin.subsection = headers2.in.section[k]+1
    end.subsection   = ifelse(k == length(headers2.in.section),
                              end.section,
                              headers2.in.section[k+1]-1)
    
    subsection[[ all.pars.split [[headers2.in.section[[k]]]] ]] <-
      read.table(text = all.pars.split[begin.subsection:end.subsection],
                 col.names = varnames.vector)
  }
  
  all.pars.list[[all.pars.split[[headers1.ind[[i]]]]]] <- subsection
  
}

# list to partable
partable <-  reshape2::melt(all.pars.list, id.vars = "parameter") %>%
  reshape2::dcast(L1 + L2 + parameter ~ variable, value.var = "value")

return(partable)
}
