#' Extracts summaries from Mplus Bayesian approximate invariance output file
#' 
#' @param file filename of the Mplus alignment output file
#' @examples
#' \dontrun{  
#'    align.summ <- extractBayesDiff("BSEM_prior01.out") 
#'    }
#' @return A summary table and a parsed difference outpout to data frame.
#' @seealso \code{\link[MIE]{extractAlignment}}  and \code{\link[MIE]{extractAlignmentSim}} 
#' @export
extractBayesDiff <- function(file) {
  
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


# Extract alignment within PSEM framework - TURNS OUT IT'S ALL ALREADY AVAILABLE  IN extractAlfDiff()
# 
# Extracts part of Mplus output created by `OUTPUT: ALIGNMENT;` when used in combination with `MODEL PRIORS:` and `~ALF(x, y)`. This output describes pairwise as well as omnibus (chi-square) tests of measurement invariance of parameters after applying the alignment procedure.
# 
# @param file filename of the Mplus output file.
# @export
# extractPSEMDiff <- function(file) {
#   
#   b.string <- paste(readLines(file),
#                      collapse = "\n")
#   
#   if(!any(grepl("ALIGNMENT OUTPUT", b.string)))
#     stop("No alignment output in the file. Check if OUTPUT: ALIGNMENT; was specified in the Mplus input.")
#   
#   # extract alignment part
#   align.outp <- MIE:::extractBetween("ALIGNMENT OUTPUT", "Beginning Time", b.string)
#   align.outp.l <- base::strsplit(align.outp, "(DIFF ANALYSIS FOR PARAMETERS)")[[1]]
# 
#   align.outp.l = align.outp.l[!sapply(align.outp.l, function(x) gsub("\n|\\s", "", x)=="")]
#   #length(align.outp.l)
#   
#   param.names = lapply(align.outp.l, function(x) {
#      pars = stringr::str_extract(gsub("\n", " ", x), "^(.*?)(?=Chi)")
#      pars = pars[!trimws(pars) %in% c("", "\n")]
#      pars = paste(unlist(pars), collapse = " ")
#      pars = gsub("\\s+", " ", trimws(pars))
#      strsplit(pars, " ")[[1]]
#   })
#   
#   parse_mplus_diff <- function(text_str) {
#     # Split text into trimmed lines
#     lines <- unlist(strsplit(text_str, "\r?\n"))
#     lines <- trimws(lines)
#     
#     ## 1. Global Fit Statistics Dataframe ####
#     chi_sq <- suppressWarnings(as.numeric(sub(".*Chi-square value\\s+", "", grep("Chi-square value", lines, value = TRUE)[1])))
#     df_val <- suppressWarnings(as.numeric(sub(".*Degrees of freedom\\s+", "", grep("Degrees of freedom", lines, value = TRUE)[1])))
#     
#     # Locate the global P-value specifically under Degrees of Freedom
#     df_idx <- grep("Degrees of freedom", lines)
#     p_val <- NA_real_
#     if (length(df_idx) > 0) {
#       sub_lines <- lines[df_idx[1]:min(df_idx[1] + 3, length(lines))]
#       p_line <- grep("^P-value", sub_lines, value = TRUE)
#       if (length(p_line) > 0) {
#         p_val <- suppressWarnings(as.numeric(sub("^P-value\\s+", "", p_line[1])))
#       }
#     }
#     
#     df_global <- data.frame(
#       ChiSquare = chi_sq,
#       DF = df_val,
#       PValue = p_val,
#       stringsAsFactors = FALSE
#     )
#     
# 
#     ## 2. Pairwise Comparisons Dataframe ####
# 
#     pw_start <- grep("^Param\\s+Param\\s+Value", lines)
#     pw_end   <- grep("Approximate Invariance Holds For:", lines)
#     df_pairwise <- data.frame()
#     
#     if (length(pw_start) > 0 && length(pw_end) > 0) {
#       pw_lines <- lines[(pw_start[1] + 1):(pw_end[1] - 1)]
#       
#       # Filter for lines that contain exactly 7 columns (Param1 Param2 Val1 Val2 Diff SE PVal)
#       tokens <- strsplit(pw_lines, "\\s+")
#       valid_pw <- sapply(tokens, length) == 7
#       
#       if (any(valid_pw)) {
#         df_pairwise <- read.table(
#           text = paste(pw_lines[valid_pw], collapse = "\n"),
#           header = FALSE,
#           col.names = c("Param1", "Param2", "Value1", "Value2", "Difference", "SE", "P_Value"),
#           stringsAsFactors = FALSE
#         )
#       }
#     }
#     
#     ## 3. Summary & Invariant Parameter Names ####
# 
#     avg_line <- grep("Average Value Across Invariant Parameters:", lines, value = TRUE)
#     avg_val <- if (length(avg_line) > 0) suppressWarnings(as.numeric(sub(".*Average Value Across Invariant Parameters:\\s+", "", avg_line[1]))) else NA_real_
#     
#     inv_start <- grep("Approximate Invariance Holds For:", lines)
#     inv_params <- character(0)
#     if (length(inv_start) > 0 && length(avg_line) > 0) {
#       avg_idx <- grep("Average Value Across Invariant Parameters:", lines)[1]
#       if (avg_idx > inv_start[1] + 1) {
#         param_lines <- lines[(inv_start[1] + 1):(avg_idx - 1)]
#         inv_params <- unlist(strsplit(paste(param_lines, collapse = " "), "\\s+"))
#         inv_params <- inv_params[inv_params != ""]
#       }
#     }
#     
#     df_summary <- data.frame(
#       Average_Value = avg_val,
#       Invariant_Parameters_Count = length(inv_params),
#       Invariant_Parameters = paste(inv_params, collapse = ", "),
#       stringsAsFactors = FALSE
#     )
#     
#     # 4. Invariant Values Dataframe
#     inv_tbl_start <- grep("^Param\\s+Value\\s+Difference", lines)
#     df_invariant <- data.frame()
#     
#     if (length(inv_tbl_start) > 0) {
#       inv_tbl_lines <- lines[(inv_tbl_start[1] + 1):length(lines)]
#       
#       # Filter for lines containing exactly 5 columns (Param Value Diff SE PVal)
#       tokens <- strsplit(inv_tbl_lines, "\\s+")
#       valid_inv <- sapply(tokens, length) == 5
#       
#       if (any(valid_inv)) {
#         df_invariant <- read.table(
#           text = paste(inv_tbl_lines[valid_inv], collapse = "\n"),
#           header = FALSE,
#           col.names = c("Param", "Value", "Difference", "SE", "P_Value"),
#           stringsAsFactors = FALSE
#         )
#       }
#     }
#     
#     
#     # Return Named List of Data Frames
#     return(list(
#       global_stats = df_global,
#       pairwise_comparisons = df_pairwise,
#       summary_stats = df_summary,
#       invariant_values = df_invariant
#     ))
#   }
#   
#   
#   align.outp.out = lapply(align.outp.l, function(x) parse_mplus_diff(x))
#   
#   # param.names.short = sapply(param.names, function(x) {
#   #   
#   #   chr.list = strsplit(x, character(0) )
#   #   shortest = min(lengths(chr.list))
#   #   ch.mat = sapply(1:shortest, function(i) sapply(chr.list, `[[`, i))
#   #   recurring.beginning = apply(ch.mat, 2, function(v) all(v == v[1]))
#   #   ch.mat[1, recurring.beginning]
#   #     
#   #   })
#   align.outp.out= lapply(1:length(align.outp.out), function(el) append(align.outp.out[[el]], list(parnames = param.names[[el]])))
#   names(align.outp.out) <- sapply(param.names, paste, collapse = "___")
#   
#   
#   align.outp.summary = sapply(unname(align.outp.out), function(x) { 
#     cbind(x$global_stats,
#           All_pars_N = length(x$parnames),
#           x$summary_stats,
#           all_pars = paste(x$parnames, collapse = ", "))
#   })
#   
#   return(list(summary = align.outp.summary,
#               details = align.outp.out))
#   
#   
# }
  
  
