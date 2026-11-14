#' Stratified measurement invariance tests with two grouping variables
#'
#' Tests measurement invariance when groups are defined by two variables,
#' either nested (`mode = "hierarchical"`, e.g. classrooms within countries) or
#' crossed (`mode = "intersectional"`, e.g. gender by ethnicity). Besides the
#' usual global tests, it fits *stratified* models, in which the parameters are
#' constrained to equality only among the groups that share a level of one of
#' the grouping variables. Everything is run for the metric stage (loadings)
#' and then for the scalar stage (loadings and intercepts/thresholds).
#'
#' @param model character, lavaan syntax model.
#' @param data data frame.
#' @param group1 character, name of a grouping variable. In the hierarchical
#'   mode it is the **higher-level** variable (e.g. country); in the
#'   intersectional mode, the first of the two crossed variables.
#' @param group2 character, name of a grouping variable. In the hierarchical
#'   mode it is the **lower-level** variable nested in `group1` (e.g.
#'   classroom); in the intersectional mode, the second crossed variable.
#' @param mode `"hierarchical"` or `"intersectional"`, see Details.
#' @param stages character vector, `"metric"`, `"scalar"`, or both (default).
#' @param print logical, whether the result should be printed.
#' @param ... further arguments passed to \code{\link[lavaan]{cfa}}, e.g.
#'   `estimator`, `ordered`, `missing`, `std.lv`, `group.partial`. The arguments
#'   `group` and `group.equal` are set by the function and cannot be used.
#'
#' @details
#' In both modes the finest grouping is the combination of `group1` and
#' `group2` (the *cells*). A stratified model uses the cells as groups and
#' applies equality constraints within strata, that is, among the cells that
#' share a level of a grouping variable. Each stratum gets its own reference
#' group for the scale of the latent variables (latent means at the scalar
#' stage, latent variances if `std.lv = TRUE`), so the latent means and
#' variances stay free across the cells.
#'
#' **Hierarchical.** Three sets of models are fitted at each stage:
#' \enumerate{
#'   \item `lower.*`: MGCFA across all lower-level groups (cells);
#'   \item `higher.*`: MGCFA across the higher-level groups, the lower-level
#'     groups being pooled;
#'   \item `within.*`: lower-level groups as groups, parameters equal across
#'     lower-level groups within each higher-level group, and free across the
#'     higher-level groups.
#' }
#' `lower.configural`, `within.metric`, and `lower.metric` are nested, so the
#' overall test of metric invariance splits into a within-part (across
#' lower-level groups within higher-level groups) and a between-part (across
#' higher-level groups, given invariance within them). The `higher.*` models
#' are fitted to differently grouped data: their fit cannot be compared with
#' the fit of the other models, and they take invariance, equal latent
#' means, and equal (co)variances of the lower-level groups for granted.
#'
#' **Intersectional.** All models use the cells as groups:
#' \enumerate{
#'   \item `configural`, `metric`, `scalar`: the usual global models across all
#'     cells;
#'   \item `*.within.<group1>`: parameters equal among the cells sharing a level
#'     of `group1`, i.e. invariance across `group2` within every level of
#'     `group1`; parameters differ by `group1` only;
#'   \item `*.within.<group2>`: the same for `group2`.
#' }
#' The two stratified models are not nested in each other; both are nested
#' between the less and the more constrained global models.
#'
#' At the scalar stage the stratified models constrain both loadings and
#' intercepts within strata, and are compared with the stratified metric
#' models.
#'
#' Comparisons report the (scaled, if applicable) chi-square difference test,
#' differences in CFI and RMSEA, and `rmsea.d`, the RMSEA of the difference
#' test (Savalei, Brace, & Fouladi, 2023), which, unlike the difference in
#' CFI, does not depend on the fit of the baseline model and is expressed
#' per degree of freedom.
#'
#' @returns Invisibly, a list of class `stratifiedMI2` with elements
#'   `models` (fitted lavaan objects), `fit` (data frame of fit measures),
#'   `comparisons` (data frame of model comparisons), `cells` (data frame
#'   mapping the cells to the levels of the grouping variables), `mode`, and
#'   `group` (names of the grouping variables).
#'
#' @examples
#' \dontrun{
#' # classrooms nested in countries
#' stratifiedMI2("F =~ v1 + v2 + v3 + v4", data = Dat1,
#'               group1 = "country", group2 = "classroom",
#'               mode = "hierarchical")
#'
#' # gender crossed with ethnicity
#' stratifiedMI2("F =~ v1 + v2 + v3 + v4", data = Dat1,
#'               group1 = "gender", group2 = "ethnicity",
#'               mode = "intersectional")
#' }
#' @seealso \code{\link{stratifiedMI}}, \code{\link{globalMI}}
#' @export
stratifiedMI2 <- function(model, data, group1, group2,
                          mode = c("hierarchical", "intersectional"),
                          stages = c("metric", "scalar"),
                          print = TRUE, ...) {

  mode   <- match.arg(mode)
  stages <- match.arg(stages, several.ok = TRUE)
  dots   <- list(...)

  if (any(c("group", "group.equal") %in% names(dots)))
    stop("Arguments 'group' and 'group.equal' are set by stratifiedMI2() and cannot be passed.")
  dots$do.fit <- NULL
  if (!is.character(group1) || !is.character(group2) ||
      length(group1) != 1 || length(group2) != 1 || group1 == group2)
    stop("'group1' and 'group2' should be names of two different variables.")
  data <- as.data.frame(data)
  if (!all(c(group1, group2) %in% names(data)))
    stop("Grouping variables are not found in the data.")

  # Cells: combinations of the two grouping variables --------------------
  miss <- is.na(data[[group1]]) | is.na(data[[group2]])
  if (any(miss)) {
    message(sum(miss), " rows with missing values on the grouping variables were dropped.")
    data <- data[!miss, , drop = FALSE]
  }
  cell <- make.unique(c(names(data), paste0(group1, "_", group2)))[ncol(data) + 1]
  data[[cell]] <- paste(data[[group1]], data[[group2]], sep = "_")
  cells <- unique(data[, c(cell, group1, group2)])
  rownames(cells) <- NULL
  names(cells)[1] <- "cell"
  cells$n <- as.vector(table(data[[cell]])[cells$cell])

  if (mode == "hierarchical") {
    if (any(duplicated(unique(data[, c(group1, group2)])[[group2]])))
      message("Some labels of '", group2, "' occur in more than one level of '",
              group1, "'; they are treated as different groups.")
    if (length(unique(data[[group1]])) < 2)
      stop("'", group1, "' should have at least two levels.")
  } else {
    full.cross <- length(unique(data[[group1]])) * length(unique(data[[group2]]))
    if (nrow(cells) < full.cross)
      message(full.cross - nrow(cells), " of ", full.cross, " combinations of '",
              group1, "' and '", group2, "' are empty.")
  }

  # Fitting -------------------------------------------------------------
  equal <- list(configural = NULL,
                metric     = "loadings",
                scalar     = c("loadings", "intercepts", "thresholds"))

  fit <- function(level, group = cell, strata.by = NULL) {
    args <- c(list(model = model, data = data, group = group), dots)
    args$group.equal <- equal[[level]]
    tryCatch({
      if (!is.null(strata.by)) {
        full   <- do.call(lavaan::cfa, c(args, list(do.fit = FALSE)))
        args$group.equal <- args$group.partial <- NULL
        config <- do.call(lavaan::cfa, c(args, list(do.fit = FALSE)))
        labs   <- lavaan::lavInspect(full, "group.label")
        args$model <- stratify_partable(
          full, config, stratum = cells[match(labs, cells$cell), strata.by])
      }
      do.call(lavaan::cfa, args)
    },
    error = function(e) {
      warning("Model could not be fitted: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  }

  st <- c("metric", "scalar")[c("metric", "scalar") %in% stages]
  models <- list()
  comps  <- list()
  add.comp <- function(stage, from, to, test) {
    comps[[length(comps) + 1]] <<- data.frame(stage = stage, from = from, to = to,
                                              test = test, stringsAsFactors = FALSE)
  }

  if (mode == "hierarchical") {

    models$lower.configural  <- fit("configural")
    models$higher.configural <- fit("configural", group = group1)
    if ("metric" %in% st || "scalar" %in% st) {
      models$lower.metric  <- fit("metric")
      models$higher.metric <- fit("metric", group = group1)
      models$within.metric <- fit("metric", strata.by = group1)
    }
    if ("metric" %in% st) {
      add.comp("metric", "lower.configural", "within.metric",
               paste0("loadings across ", group2, " within ", group1))
      add.comp("metric", "within.metric", "lower.metric",
               paste0("loadings across ", group1, ", given the above"))
      add.comp("metric", "lower.configural", "lower.metric",
               paste0("loadings across all ", group2, " (overall)"))
      add.comp("metric", "higher.configural", "higher.metric",
               paste0("loadings across ", group1, ", ", group2, " pooled"))
    }
    if ("scalar" %in% st) {
      models$lower.scalar  <- fit("scalar")
      models$higher.scalar <- fit("scalar", group = group1)
      models$within.scalar <- fit("scalar", strata.by = group1)
      add.comp("scalar", "within.metric", "within.scalar",
               paste0("intercepts across ", group2, " within ", group1))
      add.comp("scalar", "within.scalar", "lower.scalar",
               paste0("loadings & intercepts across ", group1, ", given the above"))
      add.comp("scalar", "lower.metric", "lower.scalar",
               paste0("intercepts across all ", group2, " (overall)"))
      add.comp("scalar", "higher.metric", "higher.scalar",
               paste0("intercepts across ", group1, ", ", group2, " pooled"))
    }

  } else {

    w1 <- paste0("within.", group1)
    w2 <- paste0("within.", group2)
    models$configural <- fit("configural")
    models$metric <- fit("metric")
    models[[paste0("metric.", w1)]] <- fit("metric", strata.by = group1)
    models[[paste0("metric.", w2)]] <- fit("metric", strata.by = group2)

    if ("metric" %in% st) {
      add.comp("metric", "configural", "metric", "loadings across all cells")
      add.comp("metric", "configural", paste0("metric.", w1),
               paste0("loadings across ", group2, " within ", group1))
      add.comp("metric", paste0("metric.", w1), "metric",
               paste0("loadings across ", group1, ", given the above"))
      add.comp("metric", "configural", paste0("metric.", w2),
               paste0("loadings across ", group1, " within ", group2))
      add.comp("metric", paste0("metric.", w2), "metric",
               paste0("loadings across ", group2, ", given the above"))
    }
    if ("scalar" %in% st) {
      models$scalar <- fit("scalar")
      models[[paste0("scalar.", w1)]] <- fit("scalar", strata.by = group1)
      models[[paste0("scalar.", w2)]] <- fit("scalar", strata.by = group2)
      add.comp("scalar", "metric", "scalar", "intercepts across all cells")
      add.comp("scalar", paste0("metric.", w1), paste0("scalar.", w1),
               paste0("intercepts across ", group2, " within ", group1))
      add.comp("scalar", paste0("scalar.", w1), "scalar",
               paste0("loadings & intercepts across ", group1, ", given the above"))
      add.comp("scalar", paste0("metric.", w2), paste0("scalar.", w2),
               paste0("intercepts across ", group1, " within ", group2))
      add.comp("scalar", paste0("scalar.", w2), "scalar",
               paste0("loadings & intercepts across ", group2, ", given the above"))
    }
  }

  # Summaries -----------------------------------------------------------
  fits <- as.data.frame(t(sapply(models, smi_fit)))
  comps <- do.call(rbind, comps)
  comps <- cbind(comps, t(mapply(
    function(from, to) smi_compare(models[[from]], models[[to]], fits[from, ], fits[to, ]),
    comps$from, comps$to)))
  rownames(comps) <- NULL

  out <- structure(list(models = models, fit = fits, comparisons = comps,
                        cells = cells, mode = mode, group = c(group1, group2)),
                   class = "stratifiedMI2")
  if (print) print(out)
  invisible(out)
}


#' Turn a globally constrained parameter table into a stratified one
#'
#' @param full lavaan object (may be unfitted) with `group.equal` constraints
#'   across all groups.
#' @param config lavaan object, the same model without `group.equal`.
#' @param stratum vector with the stratum of every group, in the order of the
#'   groups in the lavaan objects.
#' @returns A parameter table in which the `group.equal` constraints hold
#'   within strata only, and the first group of every stratum identifies the
#'   scale of the latent variables in the same way the first group does it in
#'   `full`.
#' @noRd
stratify_partable <- function(full, config, stratum) {

  pt <- as.data.frame(lavaan::parTable(full),   stringsAsFactors = FALSE)
  pc <- as.data.frame(lavaan::parTable(config), stringsAsFactors = FALSE)
  pt$start <- pt$est <- pt$se <- NULL

  ref  <- match(stratum, stratum)   # reference (first) group of each group's stratum
  pos  <- stats::ave(seq_len(nrow(pt)), pt$group, FUN = seq_along)
  auto <- function(x) x == "" | grepl("^\\.p[0-9]+\\.$", x)

  # 1. Re-point the equality constraints to the stratum's reference group
  eq   <- which(pt$op == "==" & pt$user == 2)
  l    <- match(pt$lhs[eq], pt$plabel)
  r    <- match(pt$rhs[eq], pt$plabel)
  ok   <- !is.na(l) & !is.na(r)
  eq <- eq[ok]; l <- l[ok]; r <- r[ok]
  target <- ref[pt$group[r]]
  new.l  <- match(paste(target, pos[l]), paste(pt$group, pos))
  pt$lhs[eq] <- pt$plabel[new.l]
  relabel <- auto(pt$label[r])
  pt$label[r[relabel]] <- pt$plabel[new.l[relabel]]
  # the reference groups are constrained to nothing but themselves
  drop <- eq[new.l == r]

  # 2. Reference groups of the strata set the scale as in the configural model
  for (g in setdiff(unique(ref), 1)) {
    it <- which(pt$group == g)
    ic <- which(pc$group == g)
    if (length(it) != length(ic))
      stop("Constrained and unconstrained models have different parameters.")
    fix <- pt$free[it] > 0 & pc$free[ic] == 0
    pt$free[it[fix]]   <- 0
    pt$ustart[it[fix]] <- pc$ustart[ic[fix]]
  }

  if (length(drop) > 0) pt <- pt[-drop, ]
  free <- pt$free > 0
  pt$free[free] <- seq_len(sum(free))
  pt$id <- seq_len(nrow(pt))
  rownames(pt) <- NULL
  pt
}


#' Fit measures of one model, scaled versions if available
#' @noRd
smi_fit <- function(x) {
  nms <- c("chisq", "df", "cfi", "rmsea", "srmr", "bic")
  out <- setNames(rep(NA_real_, length(nms)), nms)
  if (is.null(x) || !lavaan::lavInspect(x, "converged")) return(out)
  fm <- tryCatch(lavaan::fitMeasures(x), error = function(e) NULL)
  if (is.null(fm)) return(out)
  sc <- if ("chisq.scaled" %in% names(fm) && !is.na(fm[["chisq.scaled"]])) ".scaled" else ""
  get <- function(nm) if (nm %in% names(fm)) fm[[nm]] else NA_real_
  out[] <- c(get(paste0("chisq", sc)), get(paste0("df", sc)), get(paste0("cfi", sc)),
             get(paste0("rmsea", sc)), get("srmr"), get("bic"))
  out
}


#' Compare two nested models
#' @param a,b lavaan objects, `b` is nested in `a`.
#' @param fa,fb one-row data frames of their fit measures.
#' @noRd
smi_compare <- function(a, b, fa, fb) {
  out <- c(chisq.diff = NA_real_, df.diff = NA_real_, p = NA_real_,
           rmsea.d = NA_real_, cfi.diff = NA_real_, rmsea.diff = NA_real_)
  if (is.na(fa$chisq) || is.na(fb$chisq)) return(out)
  lrt <- tryCatch(suppressWarnings(lavaan::lavTestLRT(a, b)), error = function(e) NULL)
  if (!is.null(lrt)) {
    out[c("chisq.diff", "df.diff", "p")] <-
      unlist(lrt[2, c("Chisq diff", "Df diff", "Pr(>Chisq)")])
  } else {
    out[c("chisq.diff", "df.diff")] <- c(fb$chisq - fa$chisq, fb$df - fa$df)
  }
  d <- out[["chisq.diff"]]; ddf <- out[["df.diff"]]
  if (!is.na(d) && !is.na(ddf) && ddf > 0) {
    n <- lavaan::lavInspect(b, "ntotal")
    g <- lavaan::lavInspect(b, "ngroups")
    out["rmsea.d"] <- sqrt(max(0, (d - ddf) / (ddf * n))) * sqrt(g)
  }
  out["cfi.diff"]   <- fb$cfi - fa$cfi
  out["rmsea.diff"] <- fb$rmsea - fa$rmsea
  out
}


#' @export
print.stratifiedMI2 <- function(x, digits = 3, ...) {
  g <- x$group
  cat("\nStratified measurement invariance,", x$mode, "mode\n")
  if (x$mode == "hierarchical") {
    cat("  lower.*  : groups are the", nrow(x$cells), "levels of", g[2], "\n")
    cat("  higher.* : groups are the", length(unique(x$cells[[g[1]]])), "levels of",
        g[1], "(fit is not comparable with the other models)\n")
    cat("  within.* : groups as in lower.*, constraints across", g[2], "within", g[1], "\n")
  } else {
    cat("  Groups are the", nrow(x$cells), "combinations of", g[1], "and", g[2], "\n")
    cat("  *.within.", g[1], " : constraints across ", g[2], " within ", g[1],
        " (parameters differ by ", g[1], " only)\n", sep = "")
    cat("  *.within.", g[2], " : constraints across ", g[1], " within ", g[2],
        " (parameters differ by ", g[2], " only)\n", sep = "")
  }
  cat("  Smallest group: n =", min(x$cells$n), "\n")

  cat("\nModel fit\n")
  print(round(x$fit, digits), na.print = "")

  for (s in unique(x$comparisons$stage)) {
    cat("\nComparisons,", s, "stage\n")
    d <- x$comparisons[x$comparisons$stage == s, -1]
    num <- sapply(d, is.numeric)
    d[num] <- round(d[num], digits)
    print(d, row.names = FALSE, right = FALSE, na.print = "")
  }
  cat("\n")
  invisible(x)
}
