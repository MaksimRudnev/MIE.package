#' Extracts summaries from Mplus alignment output file (version-robust parser)
#'
#' Handles thee output layouts of Mplus 8.x - 9.1.1, four estimators (ML, MLR, BAYES, WLSMV), both
#' multiple-group (GROUPING) and mixture (KNOWNCLASS) set-ups, and continuous, categorical
#' (binary / ordinal) or mixed indicators. 
#' 
#' 
#' @details 
#' Only alignment optimization (\code{ALIGNMENT = FIXED} or \code{FREE}) is handled here.
#' Outputs of ALF priors (\code{DIFF(...)~ALF(0,1)} in MODEL PRIORS) are read by
#' \code{\link[MIE]{extractAlfDiff2}}.
#' 
#' Variables, categorical variables, factors,
#' estimator, grouping variable and alignment settings are read from the SUMMARY OF ANALYSIS
#' section of the output, not from the echo of the input.
#' @section Differences between Mplus versions that the function accommodates:
#' \tabular{lll}{
#'   \strong{Part of the output} \tab \strong{Layouts found} \tab \strong{Handling} \cr
#'   Loadings section \tab Mplus 8.9 and later (ML, MLR, WLSMV): \code{Loadings for F1}, then \code{Loadings for X1} for each indicator.
#'      Mplus 8.4 and Bayes (all versions): a bare \code{Loadings}, then \code{Loadings for X1}, repeated for each factor \tab
#'      Factor names are taken from the summary of analysis; a bare header is counted as the next factor \cr
#'   Factor mean comparison title \tab \code{FACTOR MEAN COMPARISON} (8.4, 9.x), \code{FACTOR INTERCEPT COMPARISON} (8.9, 8.10; the column is then called Intercept),
#'      \code{FACTOR MEAN/INTERCEPT COMPARISON} (older versions) \tab Any of the three is accepted \cr
#'   Mean comparison columns \tab \code{Latent Class} (mixture) or \code{Group} (GROUPING); long lists of groups wrap over lines \tab
#'      Always returned as \code{Latent.class}; wrapped lists are joined \cr
#'   Bayes tables \tab 9 instead of 7 columns in the pairwise tables and 7 instead of 5 in the table of invariant values (two credibility limits);
#'      wording \code{Approximate Invariance (Noninvariance) Holds For Groups} \tab Columns are named \code{lowerCI}, \code{upperCI}; wording is ignored \cr
#'   Order of blocks \tab Thresholds, then intercepts, then loadings; categorical indicators come before continuous ones within a factor \tab
#'      Parameters are matched by name, never by position \cr
#'   TECH8 alignment results \tab Printed for TYPE = MIXTURE in all versions; for GROUPING only in Mplus 8.4 (not printed in 8.9 and later). Not tested with Bayes or WLSMV \tab
#'      \code{Fit.contribution} is added when present; otherwise a warning is given if TECH8 was requested \cr
#'   TECH8 tables \tab No row labels; rows repeat a variable once per threshold; a separate set of tables for each factor \tab
#'      Rows are assigned factor by factor; loadings + intercepts = total contribution holds for all assigned values \cr
#'   Estimator for continuous items \tab WLSMV with only continuous indicators is run by Mplus as ML \tab
#'      \code{extra$estimator} reports what Mplus used (\code{ML}) \cr
#'   Input echo \tab Keywords in any case, ranges in NAMES, no USEVARIABLES, comments \tab Not needed; summary of analysis is used
#' }
#'
#' @section What to expect in different conditions:
#' \tabular{lll}{
#'   \strong{Condition} \tab \strong{Rows in \code{$summary} / \code{$alignment.output}} \tab \strong{Other differences} \cr
#'   Continuous indicators \tab \code{Intercept X1} and \code{Loadings F1 by X1} for each indicator (2 rows) \tab
#'      \code{extra$categorical.var.names} is NULL; no \code{parameterization} or \code{link} \cr
#'   Binary / ordinal indicators, ML or MLR \tab \code{Threshold X1$1} (and \code{$2}, ... for more categories) and \code{Loadings}; no intercepts \tab
#'      Logit scale (\code{extra$link} is LOGIT); numerical integration is required by Mplus, hence TYPE = MIXTURE \cr
#'   Binary / ordinal indicators, WLSMV \tab Same as above \tab
#'      Probit scale; \code{extra$parameterization} is THETA or DELTA; GROUPING only; Mplus 8.9 or later \cr
#'   Mixed continuous and categorical \tab Thresholds of the categorical indicators, intercepts of the continuous ones, loadings of all; printed in this order \tab
#'      Both \code{extra$continuous.var.names} and \code{extra$categorical.var.names} are filled \cr
#'   Bayes \tab Categorical indicators (Mplus allows only binary ones with alignment) are reported as \code{Intercept X1}, not as thresholds \tab
#'      Pairwise tables and invariant-value tables have \code{lowerCI}, \code{upperCI}; probit link; Mplus accepts Bayes alignment only with TYPE = MIXTURE \cr
#'   ML / MLR versus WLSMV and Bayes \tab Same rows (see above) \tab
#'      \code{Fit.contribution} is read only for ML / MLR / MLF, and only if \code{OUTPUT: TECH8} was requested and Mplus printed alignment results (mixtures) \cr
#'   GROUPING (multiple-group) \tab Same \tab \code{extra$type} is MG; \code{extra$grouping.variable} is the grouping variable; tables list group values \cr
#'   KNOWNCLASS (mixture) \tab Same \tab \code{extra$type} is MIXTURE; \code{extra$grouping.variable} is the variable in KNOWNCLASS; tables list latent class numbers \cr
#'   ALIGNMENT = FIXED \tab Same \tab \code{extra$alignment} is FIXED; the factor mean of the reference group (\code{extra$reference.group}) is 0 in \code{$mean.comparison} \cr
#'   ALIGNMENT = FREE \tab Same \tab \code{extra$alignment} is FREE; the factor mean of the reference group is estimated \cr
#'   ALF priors (DIFF ... ~ALF) \tab Not an alignment optimization: no mean comparison, no aligned groups; one block per DIFF statement \tab
#'      Use \code{extractAlfDiff2}; this function stops with a message \cr
#'   SAVEDATA: RANKING \tab \code{$ranking.table}: one matrix of pairwise '>' / '<' signs per factor \tab
#'      Looked up in the folder of the output file; without it a warning is given \cr
#'   Several factors \tab Loadings of every factor, named \code{Loadings F2 by X2}; one table per factor in \code{$mean.comparison} and \code{$ranking.table} \tab
#'      \code{Fit.contribution} is matched factor by factor
#' }
#' @returns \describe{
#'   \item{mean.comparison}{(\code{"comparisons"}) A list with a data frame per factor: \code{Ranking}, \code{Latent.class}, \code{Group.value},
#'     \code{Factor.mean}, \code{Groups.With.Significantly.Smaller.Factor.Mean}, \code{NGroups.With.Smaller.Factor.Mean}.}
#'   \item{ranking.table}{(\code{"ranking"}) A list with a matrix per factor, from the ranking file saved with \code{SAVEDATA: RANKING}.}
#'   \item{alignment.output}{(\code{"summary"}) A list with an element per parameter (\code{Intercept X1}, \code{Threshold X1$1}, \code{Loadings F1 by X1}) with
#'     \code{Pairwise comparison}, \code{Aligned parameter}, \code{R2}, \code{Comparison of aligned pars to average}, \code{N_noninvariant}, \code{N_invariant},
#'     \code{Invariant groups}, \code{Non-invariant groups}.}
#'   \item{non.invariant.pars}{A parameters by groups matrix with "X" for non-invariant parameter-group combinations.}
#'   \item{summary}{A data frame with a row per parameter: \code{AlignedParameter}, \code{R2}, \code{N_invariant}, \code{N_noninvariant},
#'     \code{invariant.gr}, \code{non.invar.gr} and, if available, \code{Fit.contribution}.}
#'   \item{savedata}{(\code{"savedata"}) Data saved by the SAVEDATA command.}
#'   \item{extra}{\code{estimator}, \code{mplus.version}, \code{final.message}, \code{var.names}, \code{categorical.var.names}, \code{parameterization},
#'     \code{type}, \code{string}; and, new in this function, \code{alignment}, \code{reference.group}, \code{link}, \code{grouping.variable}, \code{n.groups},
#'     \code{factors}, \code{continuous.var.names}.}
#' }
#' Row names use one space (\code{Loadings F1 by X1}) in all versions; \code{extractAlignment} put two spaces after \code{Loadings} in some versions.
#'
#' @param file filename of the Mplus alignment output file
#' @param nice.tables Logical. If tables should be send to RStudio Viewer pane using `kable` and `kableExtra` packages.
#' @param silent Logical. Used for debugging.
#' @param what Character vector. What exactly to extract. Possible values: "summary", "ranking", "comparisons", "contributions", "savedata".
#' @examples
#' \dontrun{
#'    align.summ <- extractAlignment2("fixed.out")
#'    }
#' @return A list of summary tables, see the section "Returned list".
#' @seealso \code{\link[MIE]{extractAlignment}}, \code{\link[MIE]{extractAlfDiff2}}
#' @export
extractAlignment <- function(file = "fixed.out",
                              nice.tables = FALSE,
                              silent = FALSE,
                              what = c("summary", "ranking", "comparisons", "contributions")) {

  if (!file.exists(file)) stop("File not found: ", file)

  ## Read file ####
  lines <- sub("\r$", "", readLines(file, warn = FALSE))
  b.string <- paste(lines, collapse = "\n")
  out.dir  <- dirname(normalizePath(file))

  mplus.version <- trimws(gsub("(Mplus version)|(\\(.*\\))", "", lines[[1]], ignore.case = TRUE))

  if (!any(grepl("^ALIGNMENT OUTPUT", lines)))
    stop("No ALIGNMENT OUTPUT section in '", file, "'. Did Mplus terminate with an error, ",
         "or is 'OUTPUT: ALIGN;' missing from the input?")

  if (any(grepl("DIFF ANALYSIS FOR PARAMETERS", lines, fixed = TRUE)))
    stop("'", file, "' contains an ALF (DIFF priors) analysis, not an alignment optimization. Use extractAlfDiff2().")

  ## Input echo (only used for what the summary does not tell) ####
  input <- .ea2_input(lines)

  ## Summary of analysis ####
  sa <- .ea2_summary_of_analysis(lines)
  estimator <- sa$estimator
  type.analysis <- if (sa$knownclass || grepl("TYPE\\s*(IS|ARE|=)?\\s*MIXTURE", toupper(input))) "MIXTURE" else "MG"

  grp.var <- if (type.analysis == "MG") {
    sa$grouping
  } else {
    g <- regmatches(input, regexec("(?i)knownclass\\s*(?:is|are|=)?\\s*\\w+\\s*\\(\\s*(\\w+)\\s*=", input, perl = TRUE))[[1]]
    if (length(g) == 2) g[2] else NA_character_
  }

  # variables, in the order of the NAMES statement when it can be parsed
  used.vars.list <- c(sa$continuous, sa$categorical)
  nms <- .ea2_names(input)
  if (length(nms)) {
    ord <- match(toupper(used.vars.list), toupper(nms))
    if (!anyNA(ord)) used.vars.list <- used.vars.list[order(ord)]
  }
  categorical.var.names <- if (length(sa$categorical)) sa$categorical else NULL

  ## Factor mean comparison ####
  mean.comp <- .ea2_mean_comparison(lines)
  if (is.null(mean.comp)) {
    warning("Could not extract mean comparison table. Check if your Mplus input had a 'OUTPUT: align;' line.")
    mean.comp <- NA
  }

  output <- list()
  if ("comparisons" %in% what) output[["mean.comparison"]] <- mean.comp

  ## Ranking tables ####
  if ("ranking" %in% what && grepl("Factor Mean Ranking Tables", b.string, fixed = TRUE)) {
    ranking.tab <- .ea2_ranking(lines, out.dir)
    if (is.null(ranking.tab)) warning("Could not find ranking file.") else output[["ranking.table"]] <- ranking.tab
  }

  ## Alignment output: pairwise comparisons ####
  factors <- sa$factors
  if (!length(factors) && is.list(mean.comp)) factors <- names(mean.comp)
  blocks <- .ea2_alignment_blocks(lines, factors)
  if (!length(blocks)) stop("The ALIGNMENT OUTPUT section of '", file, "' contains no parameter tables.")
  align.outp <- lapply(blocks, `[[`, "result")
  meta <- do.call(rbind, lapply(blocks, `[[`, "meta"))
  rownames(meta) <- NULL

  if ("summary" %in% what) output[["alignment.output"]] <- align.outp

  # non-invariant parameters by group
  all.groups <- unique(unlist(lapply(align.outp, function(x) c(x[["Invariant groups"]], x[["Non-invariant groups"]]))))
  all.groups <- all.groups[order(suppressWarnings(as.numeric(all.groups)), all.groups)]
  non.invariant.groups.pars.tab <- sapply(align.outp, function(x) {
    temp <- rep("", length(all.groups))
    temp[!all.groups %in% x[["Invariant groups"]]] <- "X"
    temp
  })
  non.invariant.groups.pars.tab <- matrix(non.invariant.groups.pars.tab, nrow = length(all.groups),
                                          dimnames = list(all.groups, names(align.outp)))
  output[["non.invariant.pars"]] <- t(non.invariant.groups.pars.tab)

  # Combine summaries: R2, aligned parameters, list of invariant and non-invariant groups
  summ <- lapply(align.outp, function(x) data.frame(
    AlignedParameter = x[["Aligned parameter"]],
    R2 = x[["R2"]],
    N_invariant = x[["N_invariant"]],
    N_noninvariant = x[["N_noninvariant"]],
    invariant.gr = paste(x[["Invariant groups"]], collapse = " "),
    non.invar.gr = paste(x[["Non-invariant groups"]], collapse = " "),
    stringsAsFactors = FALSE))
  summ1 <- Reduce("rbind", summ)
  rownames(summ1) <- names(summ)
  output[["summary"]] <- summ1

  ## Fit contribution (TECH8) ####
  if (estimator %in% c("MLR", "ML", "MLF") && "contributions" %in% what) {
    fc <- .ea2_fit_contrib(lines, meta)
    if (is.null(fc)) {
      if (grepl("TECH8", toupper(input)))
        warning("Extracting fit contributions failed: the output has no TECHNICAL 8 alignment results ",
                "(Mplus prints them for TYPE = MIXTURE, and for GROUPING only in some versions).")
    } else {
      output$summary[["Fit.contribution"]] <- fc[rownames(output$summary)]
    }
  }

  ## Savedata ####
  if ("savedata" %in% what) {
    if (!grepl("SAVEDATA INFORMATION", b.string)) {
      warning("Couldn't find the Savedata information.")
    } else {
      saved <- tryCatch(.ea2_savedata(lines, out.dir), error = function(e) {
        warning("Reading the saved data failed: ", conditionMessage(e)); NULL })
      if (!is.null(saved)) output$savedata <- saved
    }
  }

  ## Final message ####
  sum.noninv <- sum(output$summary[, "N_noninvariant"])
  n.params.total <- sum(output$summary[, c("N_invariant", "N_noninvariant")])
  final.message <- paste("There are", sum.noninv,
        "non-invariant parameters out of",
        n.params.total,
        "which is",
        scales::percent(sum.noninv / n.params.total),
        ".\n This is", ifelse(sum.noninv / n.params.total < .25,
                              "smaller than the recommended cutoff of 25%,\n thereby the results support approximate invariance. \nRunning simuations is still recommended.\n",
                              "greater than the recommended cutoff of 25%,\nthereby the results do NOT support approximate invariance.\nRunning simulations is recommended for further exploration.\n"))

  output$extra <- list(estimator = estimator,
                       mplus.version = mplus.version,
                       final.message = final.message,
                       var.names = used.vars.list,
                       categorical.var.names = categorical.var.names,
                       parameterization = if (is.na(sa$parameterization)) NA else sa$parameterization,
                       type = type.analysis,
                       string = b.string,
                       # additional information, not available in extractAlignment()
                       alignment = sa$alignment,
                       reference.group = sa$reference.group,
                       link = sa$link,
                       grouping.variable = grp.var,
                       n.groups = length(all.groups),
                       factors = factors,
                       continuous.var.names = if (length(sa$continuous)) sa$continuous else NULL)

  if (!silent) {
    print(output$summary, row.names = TRUE)
    cat("\n", final.message)
  }

  if (nice.tables && !silent) {
    nice.tab1 <- knitr::kable(output$non.invariant.pars, format = "html")
    trash <- capture.output(kableExtra::kable_styling(nice.tab1, bootstrap_options = c("striped", "bordered"),
                                                      position = "left", font_size = 12))
  }

  invisible(output)
}


# ---- helpers: input echo and summary of analysis ------------------------------------------------

.ea2_num <- function(x) suppressWarnings(as.numeric(x))

.ea2_tokens <- function(x) {
  x <- trimws(x)
  if (!length(x)) return(character())
  out <- unlist(strsplit(x, "\\s+"))
  out[nzchar(out)]
}

# input echo as a single squished string (comments and Mplus' own messages removed)
.ea2_input <- function(lines) {
  i0 <- grep("^INPUT INSTRUCTIONS", lines)[1]
  i1 <- grep("^SUMMARY OF ANALYSIS", lines)[1]
  if (is.na(i0) || is.na(i1)) return("")
  x <- lines[(i0 + 1):(i1 - 1)]
  x <- gsub("!.*", "", x)
  x <- x[!grepl("^\\s*\\*\\*\\*", x)]
  gsub("\\s+", " ", paste(x, collapse = " "))
}

# NAMES statement of the input, ranges such as x1_1-x1_4 expanded; NULL if it cannot be parsed
.ea2_names <- function(input) {
  m <- regmatches(input, regexec("(?i)\\bnames\\s*(?:are|is|=)?\\s*([^;]*);", input, perl = TRUE))[[1]]
  if (length(m) != 2) return(NULL)
  tok <- .ea2_tokens(m[2])
  tok <- tok[!toupper(tok) %in% c("ARE", "IS", "=")]
  out <- unlist(lapply(tok, function(t) {
    if (!grepl("-", t, fixed = TRUE)) return(t)
    ends <- strsplit(t, "-", fixed = TRUE)[[1]]
    a <- regmatches(ends[1], regexec("^(.*?)(\\d+)$", ends[1]))[[1]]
    b <- regmatches(ends[2], regexec("^(.*?)(\\d+)$", ends[2]))[[1]]
    if (length(a) == 3 && length(b) == 3 && a[2] == b[2])
      paste0(a[2], seq(as.integer(a[3]), as.integer(b[3])))
    else NA_character_
  }))
  if (anyNA(out)) NULL else out
}

.ea2_summary_of_analysis <- function(lines) {
  i0 <- grep("^SUMMARY OF ANALYSIS", lines)[1]
  i1 <- grep("^Input data format", lines)
  i1 <- if (length(i1)) i1[i1 > i0][1] else NA
  if (is.na(i0)) stop("No SUMMARY OF ANALYSIS section found.")
  s  <- lines[i0:(if (is.na(i1)) length(lines) else i1)]
  get1 <- function(pat) {
    x <- grep(pat, s, value = TRUE)[1]
    if (is.na(x)) NA_character_ else trimws(sub(pat, "", x))
  }

  # indented lists below a top-level heading
  block <- function(heading) {
    h <- grep(paste0("^", heading), s)[1]
    if (is.na(h)) return(character())
    rest <- s[-seq_len(h)]
    e <- grep("^[^ ]", rest)[1]
    if (!is.na(e)) rest <- rest[seq_len(e - 1)]
    rest
  }

  odv <- block("Observed dependent variables")
  cont <- cat <- character()
  kind <- NULL
  for (l in odv) {
    if (grepl("^  [^ ]", l)) {                       # two spaces: type of variables
      kind <- if (grepl("^  Continuous", l)) "cont"
              else if (grepl("^  Binary and ordered categorical", l)) "cat"
              else "other"
    } else if (grepl("^   +[^ ]", l)) {
      v <- .ea2_tokens(l)
      if (identical(kind, "cont")) cont <- c(cont, v)
      if (identical(kind, "cat")) cat <- c(cat, v)
    }
  }
  fac <- .ea2_tokens(block("Continuous latent variables"))

  list(
    estimator = {e <- get1("^Estimator\\s+"); if (is.na(e)) NA_character_ else e},
    parameterization = get1("^Parameterization\\s+"),
    link = get1("^Link\\s+"),
    grouping = get1("^\\s*Grouping variable\\s+"),
    knownclass = any(grepl("^\\s*Knownclass\\s", s)),
    alignment = get1("^\\s*Factor mean for reference group\\s+"),
    reference.group = get1("^\\s*Reference group\\s+"),
    continuous = cont, categorical = cat, factors = fac)
}


# ---- helpers: factor mean comparison and ranking -------------------------------------------------

.ea2_mean_comparison <- function(lines) {
  h <- grep("^FACTOR (MEAN|INTERCEPT|MEAN/INTERCEPT) COMPARISON AT THE", lines)[1]
  if (is.na(h)) return(NULL)
  rest <- lines[-seq_len(h)]
  e <- grep("^[A-Za-z]", rest)[1]
  sec <- if (is.na(e)) rest else rest[seq_len(e - 1)]
  fi <- grep("^\\s*Results for Factor", sec)
  if (!length(fi)) return(NULL)
  fi_end <- c(fi[-1] - 1, length(sec))

  res <- lapply(seq_along(fi), function(k) {
    blk <- sec[fi[k]:fi_end[k]]
    nm <- trimws(sub("^\\s*Results for Factor\\s*", "", blk[1]))
    rows <- list(); last <- 0
    for (l in blk[-1]) {
      m <- regmatches(l, regexec("^\\s*(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(-?\\d*\\.\\d+)(.*)$", l))[[1]]
      if (length(m) == 6) {
        rows[[length(rows) + 1]] <- list(m[2:5], .ea2_tokens(m[6]))
        last <- length(rows)
      } else if (last > 0 && nzchar(trimws(l)) && !grepl("Ranking|Factor|Group|Class", l)) {
        rows[[last]][[2]] <- c(rows[[last]][[2]], .ea2_tokens(l))     # wrapped list of groups
      }
    }
    smaller <- lapply(rows, `[[`, 2)
    tab <- data.frame(
      Ranking = .ea2_num(sapply(rows, function(r) r[[1]][1])),
      Latent.class = .ea2_num(sapply(rows, function(r) r[[1]][2])),
      Group.value = .ea2_num(sapply(rows, function(r) r[[1]][3])),
      Factor.mean = .ea2_num(sapply(rows, function(r) r[[1]][4])),
      Groups.With.Significantly.Smaller.Factor.Mean = sapply(smaller, paste, collapse = ", "),
      NGroups.With.Smaller.Factor.Mean = sapply(smaller, length),
      stringsAsFactors = FALSE)
    list(name = nm, tab = tab)
  })
  setNames(lapply(res, `[[`, "tab"), sapply(res, `[[`, "name"))
}

.ea2_ranking <- function(lines, out.dir) {
  i <- grep("Factor Mean Ranking Tables", lines)[1]
  j <- grep("Save file", lines)
  j <- j[j > i][1]
  if (is.na(j)) return(NULL)
  nxt <- lines[(j + 1):length(lines)]
  fn <- trimws(nxt[nzchar(trimws(nxt))][1])
  path <- file.path(out.dir, fn)
  if (!file.exists(path)) path <- fn
  if (!file.exists(path)) return(NULL)
  ranking.tabs <- paste(readLines(path, warn = FALSE), collapse = "\n")
  ranking.tabs <- strsplit(ranking.tabs, "Ranking table for ")[[1]][-1]
  names(ranking.tabs) <- sapply(ranking.tabs, function(x) substr(x, 1, regexpr("\n", x) - 1))
  lapply(ranking.tabs, function(x) {
    rt <- utils::read.csv(text = x, skip = 2, row.names = 1)
    colnames(rt) <- rownames(rt)
    rt[, -ncol(rt)]
  })
}


# ---- helpers: pairwise alignment tables ----------------------------------------------------------

# numeric rows (all tokens numeric) of a character vector of lines, as a list of numeric vectors
.ea2_numrows <- function(x) {
  tk <- lapply(x, .ea2_tokens)
  tk <- tk[vapply(tk, function(t) length(t) >= 2 && !anyNA(.ea2_num(t)), logical(1))]
  lapply(tk, .ea2_num)
}

.ea2_rows_to_df <- function(rows, nms_by_ncol) {
  if (!length(rows)) return(NULL)
  nc <- max(vapply(rows, length, integer(1)))
  rows <- rows[vapply(rows, length, integer(1)) == nc]
  m <- do.call(rbind, rows)
  df <- as.data.frame(m)
  nm <- nms_by_ncol[[as.character(nc)]]
  colnames(df) <- if (is.null(nm)) paste0("V", seq_len(nc)) else nm
  df
}

.ea2_alignment_blocks <- function(lines, factors) {
  i0 <- grep("^ALIGNMENT OUTPUT", lines)[1]
  ends <- grep("Average Invariance index|^SAVEDATA INFORMATION|^TECHNICAL 8 OUTPUT|^Beginning Time", lines)
  ends <- ends[ends > i0]
  i1 <- if (length(ends)) ends[1] - 1 else length(lines)
  reg <- lines[i0:i1]

  re.int  <- "^\\s*Intercept for (\\S+)\\s*$"
  re.thr  <- "^\\s*Threshold (\\S+)\\s*$"
  re.load <- "^\\s*Loadings for (\\S+)\\s*$"
  re.bare <- "^\\s*Loadings\\s*$"
  is.hdr <- grepl(re.int, reg) | grepl(re.thr, reg) | grepl(re.load, reg) | grepl(re.bare, reg)
  hdr <- which(is.hdr)

  next_nonblank <- function(i) { j <- i + 1; while (j <= length(reg) && !nzchar(trimws(reg[j]))) j <- j + 1; j }

  cur.factor <- NA_character_; n.bare <- 0
  blocks <- list()
  for (k in seq_along(hdr)) {
    i <- hdr[k]; l <- reg[i]
    if (grepl(re.bare, l)) {                       # older versions / Bayes: factor is not named
      n.bare <- n.bare + 1
      cur.factor <- if (n.bare <= length(factors)) factors[n.bare] else paste0("F", n.bare)
      next
    }
    if (grepl(re.load, l)) {
      j <- next_nonblank(i)
      if (j <= length(reg) && grepl(re.load, reg[j])) {      # "Loadings for F1" followed by "Loadings for X1"
        cur.factor <- sub(re.load, "\\1", l)
        next
      }
    }
    end <- if (k < length(hdr)) hdr[k + 1] - 1 else length(reg)
    body <- reg[(i + 1):max(i + 1, end)]
    if (end < i + 1) body <- character()

    if (grepl(re.int, l)) {
      type <- "Intercept"; var <- sub(re.int, "\\1", l); nm <- paste("Intercept", var)
    } else if (grepl(re.thr, l)) {
      type <- "Threshold"; var <- sub(re.thr, "\\1", l); nm <- paste("Threshold", var)
      var <- sub("\\$.*$", "", var)
    } else {
      type <- "Loadings"; var <- sub(re.load, "\\1", l)
      fac <- if (is.na(cur.factor)) "F?" else cur.factor
      nm <- paste("Loadings", fac, "by", var)
    }
    blocks[[nm]] <- list(
      result = .ea2_parse_param_block(body),
      meta = data.frame(name = nm, type = type, var = var,
                        factor = if (type == "Loadings") fac else NA_character_,
                        stringsAsFactors = FALSE))
  }
  blocks
}

.ea2_parse_param_block <- function(b) {
  hold <- grep("Holds For", b)[1]
  pre  <- if (is.na(hold)) b else b[seq_len(hold - 1)]
  pw <- .ea2_rows_to_df(.ea2_numrows(pre), list(
    "7" = c("Group1", "Group2", "Est_in_G1", "Est_in_G2", "Difference", "SE", "P_value"),
    "9" = c("Group1", "Group2", "Est_in_G1", "Est_in_G2", "Difference", "SE", "P_value", "lowerCI", "upperCI")))
  if (is.null(pw)) pw <- data.frame()

  wa <- grep("Weighted Average Value Across Invariant Groups", b)[1]
  inv.groups <- character()
  if (!is.na(hold)) {
    stop.at <- if (is.na(wa)) length(b) + 1 else wa
    if (stop.at > hold + 1) inv.groups <- .ea2_tokens(b[(hold + 1):(stop.at - 1)])
  }
  grab <- function(pat) {
    x <- grep(pat, b, value = TRUE)[1]
    if (is.na(x)) NA_real_ else .ea2_num(trimws(sub(paste0(".*", pat, "\\s*"), "", x)))
  }
  aligned <- grab("Weighted Average Value Across Invariant Groups:")
  r2 <- grab("R-square/Explained variance/Invariance index:")

  ic <- grep("Invariant Group Values, Difference to Average and Significance", b)[1]
  inv.comp <- NULL
  if (!is.na(ic)) {
    rows <- list()
    for (l in b[-seq_len(ic)]) {
      t <- .ea2_tokens(l)
      if (length(t) >= 2 && !anyNA(.ea2_num(t))) rows[[length(rows) + 1]] <- .ea2_num(t)
      else if (length(rows)) break
    }
    inv.comp <- .ea2_rows_to_df(rows, list(
      "5" = c("Group", "Value", "Difference", "SE", "P.value"),
      "7" = c("Group", "Value", "Difference", "SE", "P.value", "lowerCI", "upperCI")))
  }
  if (is.null(inv.comp)) {
    ncl <- if (ncol(pw) == 9) 7 else 5
    inv.comp <- as.data.frame(t(rep(NA_real_, ncl)))
    colnames(inv.comp) <- c("Group", "Value", "Difference", "SE", "P.value", "lowerCI", "upperCI")[seq_len(ncl)]
  }

  all.groups <- unique(as.character(c(unlist(pw[, c("Group1", "Group2")[c("Group1", "Group2") %in% names(pw)]]), inv.groups)))
  non.invariant.groups <- all.groups[!all.groups %in% inv.groups]
  ord <- function(g) g[order(.ea2_num(g), g)]

  list(
    "Pairwise comparison" = pw,
    "Aligned parameter" = aligned,
    "R2" = r2,
    "Comparison of aligned pars to average" = inv.comp,
    "N_noninvariant" = length(non.invariant.groups),
    "N_invariant" = length(inv.groups),
    "Invariant groups" = inv.groups,
    "Non-invariant groups" = ord(non.invariant.groups)
  )
}


# ---- helpers: fit contributions (TECHNICAL 8) ----------------------------------------------------

# Returns a named numeric vector (names = parameter names in `meta`), or NULL if TECH8 has no alignment results.
# The tables in the output have no row labels. Within "ALIGNMENT RESULTS FOR <factor>" the rows follow the
# order in which the thresholds / intercepts of that factor's indicators are printed in the alignment output
# (the loadings table repeats a variable once per threshold).
.ea2_fit_contrib <- function(lines, meta) {
  i0 <- grep("^TECHNICAL 8 OUTPUT", lines)[1]
  if (is.na(i0)) return(NULL)
  ends <- grep("^SAVEDATA INFORMATION|^Beginning Time", lines)
  ends <- ends[ends > i0]
  i1 <- if (length(ends)) ends[1] - 1 else length(lines)
  t8 <- lines[i0:i1]
  fi <- grep("ALIGNMENT RESULTS FOR", t8)
  if (!length(fi)) return(NULL)
  fi_end <- c(fi[-1] - 1, length(t8))

  one_col <- function(sec, title) {
    j <- grep(paste0("Fit Function ", title, " Contribution By Variable"), sec)[1]
    if (is.na(j)) return(NULL)
    v <- numeric()
    for (l in sec[-seq_len(j)]) {
      t <- .ea2_tokens(l)
      if (length(t) == 1 && !is.na(.ea2_num(t))) v <- c(v, .ea2_num(t)) else break
    }
    v
  }

  contrib <- setNames(rep(NA_real_, nrow(meta)), meta$name)
  par.meta <- meta[meta$type %in% c("Intercept", "Threshold"), ]
  for (k in seq_along(fi)) {
    sec <- t8[fi[k]:fi_end[k]]
    fac <- trimws(sub(".*ALIGNMENT RESULTS FOR\\s*", "", sec[1]))
    lm <- meta[meta$type == "Loadings" & meta$factor == fac, ]
    pm <- par.meta[par.meta$var %in% lm$var, ]
    ci <- one_col(sec, "Intercepts"); cl <- one_col(sec, "Loadings")
    if (is.null(ci) || is.null(cl)) next
    if (length(ci) != nrow(pm) || length(cl) != nrow(pm)) {
      warning("Fit contributions for ", fac, ": ", length(ci), " rows in TECH8 but ", nrow(pm),
              " intercepts/thresholds in the alignment output; skipped.")
      next
    }
    contrib[pm$name] <- ci
    contrib[lm$name] <- cl[match(lm$var, pm$var)]
  }
  if (all(is.na(contrib))) NULL else contrib
}


# ---- helper: saved data --------------------------------------------------------------------------

.ea2_savedata <- function(lines, out.dir) {
  save.dat.str <- lines[grep("SAVEDATA INFORMATION", lines):grep("Save missing symbol", lines)]
  filename <- trimws(save.dat.str[grep("Save file$", save.dat.str) + 1])
  if (!file.exists(filename) && file.exists(file.path(out.dir, filename))) filename <- file.path(out.dir, filename)
  fileformat <- save.dat.str[grep("Order and format of variables", save.dat.str):grep("Save file format", save.dat.str)]
  fileformat <- fileformat[3:(length(fileformat) - 2)]
  fileformat <- utils::read.table(text = paste(fileformat, sep = "\n"), col.names = c("varname", "type"))
  utils::read.fwf(filename,
                  widths = as.numeric(gsub("F|I|\\.\\d", "", fileformat$type)),
                  col.names = fileformat$varname,
                  strip.white = TRUE,
                  na.strings = trimws(gsub("Save missing symbol", "", save.dat.str[grepl("Save missing symbol", save.dat.str)])),
                  colClasses = "numeric")
}


#' Extracts DIFF alignment (ALF priors) results from Mplus output (version-robust parser)
#'
#' A line-based re-implementation of \code{\link[MIE]{extractAlfDiff}}. It copes with the
#' typo in the Mplus 8.10 output (\code{Invariant Parametetrs}), with unequal wrapping of the
#' parameter lists, and stops with an informative error if the output has no DIFF analysis
#' (Mplus 8.9 ignores \code{OUTPUT: ALIGN} without \code{ALIGNMENT=}, and Mplus 8.10 prints
#' no DIFF analysis for TYPE = MIXTURE).
#' @section What to expect in different conditions:
#' \tabular{lll}{
#'   \strong{Condition} \tab \strong{Output} \tab \strong{Notes} \cr
#'   Any indicators, GROUPING, ML / MLR / WLSMV \tab A DIFF ANALYSIS block per DIFF statement: chi-square test of equality, pairwise differences between the labelled parameters, the parameters for which approximate invariance holds, their average, and their differences from the average \tab
#'      Mplus 8.10 or later (8.9 ignores OUTPUT: ALIGN here). WLSMV with only continuous indicators is run as ML \cr
#'   TYPE = MIXTURE (KNOWNCLASS) \tab Same \tab Mplus 9 or later; Mplus 8.10 runs the model but prints no DIFF analysis, then an error is raised. With MLR, Mplus 9.x rejected some models with categorical indicators ("MIXED TYPE PRIOR") \cr
#'   BAYES \tab Not available: Mplus rejects ALF priors with ESTIMATOR = BAYES \tab  \cr
#'   Mplus 8.10 \tab Same, but \code{Average Value Across Invariant Parametetrs} is misspelled in the output \tab Handled
#' }
#' The tables list the labelled parameters (for instance \code{L1_1}, \code{A1_1}, \code{T4A_1}) and not the groups; non-invariance is
#' therefore reported per parameter set. Unlike \code{extractAlignment2} there are no factor mean comparisons.
#' @param input  Mplus output either as a single string or a filename
#' @param file Logical. If TRUE, `input` is treated as a filename.
#' @details  Mplus input file has to include the `OUTPUT: align;` in the OUTPUT section and `DIFF()` statements in the `Model priors:` section.
#' @returns A list with \code{outputs} (one element per DIFF statement: \code{param.list},
#'  \code{chisq}, \code{comparison.table}, \code{Average.invariant.value}, \code{invariant.param.list},
#'  \code{compared.to.invar.value}), \code{summaries} (one row per DIFF statement) and \code{message}.
#'  The message reports non-invariant parameters out of all parameters (extractAlfDiff divides
#'  by the number of invariant parameters).
#' @examples
#' # alf.out = extractAlfDiff2("alf_MLR_cont_1f.out")
#' @export
extractAlfDiff <- function(input, file = TRUE) {

  if (file) {
    if (!file.exists(input)) stop("File not found: ", input)
    lines <- sub("\r$", "", readLines(input, warn = FALSE))
  } else {
    lines <- sub("\r$", "", unlist(strsplit(input, "\n", fixed = TRUE)))
  }
  lines <- gsub("!.*", "", lines)

  i0 <- grep("^ALIGNMENT OUTPUT", lines)[1]
  if (is.na(i0)) stop("No ALIGNMENT OUTPUT section found.")
  starts <- grep("DIFF ANALYSIS FOR PARAMETERS", lines)
  starts <- starts[starts > i0]
  if (!length(starts))
    stop("The output has no DIFF ANALYSIS section. Check that the input has DIFF() priors in MODEL PRIORS and ",
         "OUTPUT: ALIGN; and that the Mplus version supports it (ALF priors need Mplus 8.10 or later; ",
         "for TYPE = MIXTURE Mplus 9).")
  stop.at <- grep("^Beginning Time|^SAVEDATA INFORMATION|^TECHNICAL", lines)
  stop.at <- stop.at[stop.at > max(starts)]
  last <- if (length(stop.at)) stop.at[1] - 1 else length(lines)
  ends <- c(starts[-1] - 1, last)

  outputs <- lapply(seq_along(starts), function(k) {
    x <- lines[(starts[k] + 1):ends[k]]

    # parameter list: lines up to the first blank line after the first non-blank one
    nz <- which(nzchar(trimws(x)))
    first <- nz[1]
    blank.after <- which(!nzchar(trimws(x)) & seq_along(x) > first)[1]
    param.list <- .ea2_tokens(gsub(":", " ", x[first:(blank.after - 1)]))

    grabn <- function(pat) {
      v <- grep(pat, x, value = TRUE)[1]
      if (is.na(v)) NA_real_ else .ea2_num(trimws(sub(paste0(".*", pat), "", v)))
    }
    chisq <- grabn("Chi-square value")
    chisq.df <- grabn("Degrees of freedom")
    chisq.p <- grabn("P-value")

    # pairwise comparison table: rows with 7 numeric columns after the "Param Param ..." header
    h1 <- grep("Param\\s+Param\\s+Value\\s+Value", x)[1]
    comp.tab1 <- NULL
    if (!is.na(h1)) {
      rows <- list()
      for (l in x[-seq_len(h1)]) {
        t <- .ea2_tokens(l)
        if (length(t) == 7 && !anyNA(.ea2_num(t[3:7]))) rows[[length(rows) + 1]] <- t
        else if (length(rows)) break
      }
      if (length(rows)) {
        comp.tab1 <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
        colnames(comp.tab1) <- c("Param", "Param.1", "Value", "Value.1", "Difference", "SE", "P.value")
        comp.tab1[3:7] <- lapply(comp.tab1[3:7], .ea2_num)
      }
    }

    ih <- grep("Approximate Invariance Holds For", x)[1]
    ia <- grep("Average Value Across Invariant Param", x)[1]
    invar.param.list <- if (!is.na(ih) && !is.na(ia) && ia > ih + 1) .ea2_tokens(gsub(":", " ", x[(ih + 1):(ia - 1)])) else character()
    AVAIP <- if (is.na(ia)) NA_real_ else .ea2_num(trimws(sub(".*:", "", x[ia])))

    h2 <- grep("Invariant Values, Difference to Average and Significance", x)[1]
    comp.tab2 <- NULL
    if (!is.na(h2)) {
      rows <- list()
      for (l in x[-seq_len(h2 + 1)]) {          # skip the "Param Value Difference SE P-value" header
        t <- .ea2_tokens(l)
        if (length(t) == 5 && !anyNA(.ea2_num(t[2:5]))) rows[[length(rows) + 1]] <- t
        else if (length(rows)) break
      }
      if (length(rows)) {
        comp.tab2 <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
        colnames(comp.tab2) <- c("Param", "Value", "Difference", "SE", "P.value")
        comp.tab2[2:5] <- lapply(comp.tab2[2:5], .ea2_num)
      }
    }

    list(
      param.list = param.list,
      chisq = c(statistic = chisq, df = chisq.df, p.value = chisq.p),
      comparison.table = comp.tab1,
      Average.invariant.value = AVAIP,
      invariant.param.list = invar.param.list,
      compared.to.invar.value = comp.tab2
    )
  })

  summaries <- lapply(outputs, function(x) data.frame(
    al.pars = paste(substr(paste(x$param.list, collapse = ", "), 1, 20), ">>>"),
    N_invariant = length(x$invariant.param.list),
    N_noninvariant = length(x$param.list) - length(x$invariant.param.list),
    Noninvariant.params = paste(setdiff(x$param.list, x$invariant.param.list), collapse = ", "),
    Chi.square = x$chisq[["statistic"]],
    df = x$chisq[["df"]],
    p.value = x$chisq[["p.value"]],
    Average.invariant.value = x$Average.invariant.value,
    stringsAsFactors = FALSE))
  summaries <- Reduce(rbind, summaries)

  n.non <- sum(summaries$N_noninvariant)
  n.tot <- sum(summaries$N_noninvariant + summaries$N_invariant)
  msg <- paste("ALF DIFF alignment estimated", n.non, "non-invariant parameters out of", n.tot,
               "parameters, which is ", scales::percent(n.non / n.tot))

  list(outputs = outputs,
       summaries = summaries,
       message = msg)
}
