#' Extracts summaries from Mplus Bayesian approximate invariance output file
#' 
#' @param file filename of the Mplus alignment output file
#' @examples
#' \dontrun{  
#'    align.summ <- extractDiffOutput("BSEM_prior01.out") 
#'    }
#' @return A summary table and a parsed difference outpout to data frame.
#' @seealso \code{\link[MIE]{extractAlignment}}  and \code{\link[MIE]{extractAlignmentSim}} 
#' @export
extractDiffOutput <- function(file) {
  
  str.bdiff <- paste(readLines(file),
                     collapse = "\n")
  
  
  diff.output <- sub(".*DIFFERENCE OUTPUT\n\n\n *(.*?) *\n\n\n\n.*", "\\1", str.bdiff)
  diff.output <- sub("Average   Std. Dev.     Deviations from the Mean\n\n", "", diff.output)
  
  # Extract PPP by group table to provide group labels
  ppp.by.group <- sub(".*Parameters \\(pD\\) From Each Group\\s*(.*?) *\\s*RMSEA.*", "\\1", str.bdiff)
  ppp.by.group <- gsub("\\n\\s+(Group)", "\nGroup", ppp.by.group)
  ppp.by.group = read.table(text=ppp.by.group, header = F)
  
  names(ppp.by.group)<- c("V1", "classes", "labels", "PPP", "Chi_LB", "Chi_UB", "DIC", "pD")
  
  ppp.by.group$labels <- as.numeric(gsub("\\(|\\)", "", ppp.by.group$labels))
  ppp.by.group$Chi_LB <- as.numeric(gsub("\\(|,", "", ppp.by.group$Chi_LB))
  ppp.by.group$Chi_UB <- as.numeric(gsub("\\)|,", "", ppp.by.group$Chi_UB))
  ppp.by.group <- ppp.by.group[,-1]
  
  
 # cat(diff.output)
  
#  strsplit(diff.output, "\n")[[1]]
  
  ff <- tempfile()
  cat(file = ff, diff.output, sep = "\n")
  #cat(diff.output, sep = "\n")
  average.pars <- read.fwf(file=ff, widths = c(11, 11, 11), flush=T)
  diff.pars <- read.fwf(file=ff, widths = c(11, 11, 11, 11, 12, 12, 12, 12), flush=T)
  unlink(ff)
  
  average.pars= na.omit(average.pars)
  names(average.pars) <- c("par", "average", "SD")
  
 # diff.pars= na.omit(diff.pars[-c(1:3)])
  diff.pars= diff.pars[-c(1:3)]
  diff.pars = diff.pars[!rowSums(is.na(diff.pars))==ncol(diff.pars),]
  
  diff.pars.names = diff.pars[which(1:nrow(diff.pars) %% 2 ==1),]
  diff.pars.values = diff.pars[which(1:nrow(diff.pars) %% 2 ==0),]
  
  diff.pars.values = trimws(as.vector(t(diff.pars.values)))
  diff.pars = data.frame(value=diff.pars.values,
                         labels = trimws(as.vector(t(diff.pars.names))))
  diff.pars <-na.omit(diff.pars)
  
  diff.pars$param <- gsub("^.+_",  "", diff.pars$labels)
  diff.pars$paramKind <- gsub("[^a-zA-Z]",  "", diff.pars$labels)
  diff.pars$group <- gsub("[a-zA-Z]+(.*?)_[0-9]", "\\1", diff.pars$labels)
  diff.pars$par <- gsub("[0-9]+_",  "", diff.pars$labels)
  average.pars$par <- unique(diff.pars$par)
  
  diff.pars$group.label <- ppp.by.group$labels[match(diff.pars$group, ppp.by.group$classes)]
  
  
  
  diff.pars <- merge(diff.pars,average.pars, by = "par", all.x= T )
  
  
  diff.tab <- reshape2::dcast(diff.pars, group + group.label ~ par, value.var = "value")
  average <- c("average", NA, reshape2::dcast(diff.pars, group + group.label ~ par, value.var = "average")[1,-c(1,2)])
  names(average) <- colnames(diff.tab)

  diff.tab <- rbind(diff.tab, average)
  
  list(diff.table = diff.pars, 
       ppp.by.group=ppp.by.group,
       summary = diff.tab)
  
}