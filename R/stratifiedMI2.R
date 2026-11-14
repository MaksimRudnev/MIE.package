#' Stratified measurement invariance tests with two grouping variables
#'
#' Tests measurement invariance when groups are defined by two variables,
#' either nested (`mode = "hierarchical"`, e.g. classrooms within countries) or
#' crossed (`mode = "intersectional"`, e.g. gender by ethnicity). Besides the
#' usual global tests, it fits *stratified* models, in which parameters are
#' constrained to equality only among the groups that share a level of one of
#' the grouping variables. Everything is run for the metric stage (loadings)
#' and then for the scalar stage (loadings and intercepts).
#'
#' @param model character, lavaan syntax model. Factors must be identified by
#'   marker indicators (lavaan's default); first-order factors and continuous
#'   indicators only.
#' @param data data frame.
#' @param group1 character, name of a grouping variable. In the hierarchical
#'   mode it is the **higher-level** variable (e.g. country); in the
#'   intersectional mode, the first of the two crossed variables.
#' @param group2 character, name of a grouping variable. In the hierarchical
#'   mode it is the **lower-level** variable nested in `group1` (e.g.
#'   classroom); in the intersectional mode, the second crossed variable.
#' @param mode `"hierarchical"` or `"intersectional"`, see Details.
#' @param stages character vector, `"metric"`, `"scalar"`, or both (default).
#' @param additive logical, intersectional mode only: should the models without
#'   interaction of the two grouping variables be fitted? They use nonlinear
#'   constraints and require all combinations of the two variables to be
#'   present.
#' @param pooled logical, hierarchical mode only: should the conventional MGCFA
#'   across the higher-level groups (lower-level groups pooled) be fitted?
#' @param min.n the smallest acceptable number of observations in a cell (a
#'   combination of `group1` and `group2`). Smaller cells are dropped with a
#'   warning. By default, the number of free parameters of the model in one
#'   group.
#' @param group.partial character vector of parameters that are never
#'   constrained, e.g. `c("F=~v2", "v3~1")`.
#' @param print logical, whether the result should be printed.
#' @param ... further arguments passed to \code{\link[lavaan]{cfa}}, e.g.
#'   `estimator` or `missing`. The arguments `group`, `group.equal`, `std.lv`,
#'   and `ordered` cannot be used.
#'
#' @details
#' In both modes the finest grouping is the combination of `group1` and
#' `group2` (the *cells*), and all the main models use the cells as groups.
#' A stratified model applies equality constraints within strata, that is,
#' among the cells sharing a level of a grouping variable.
#'
#' **Identification.** In every cell and every model, the loading of the
#' marker indicator (the first indicator of a factor) is fixed to 1 and its
#' intercept to 0; latent means, variances, and covariances are free. The
#' identification is thus the same in all the models and does not depend on
#' the strata. The global models are equivalent to the ones lavaan fits by
#' default (same chi-square and df).
#'
#' **Hierarchical mode.** Loadings / intercepts are:
#' \tabular{lll}{
#'   `configural` \tab free \tab free \cr
#'   `metric.within` \tab equal within `group1` \tab free \cr
#'   `metric.all` \tab equal across all cells \tab free \cr
#'   `scalar.within` \tab equal within `group1` \tab equal within `group1` \cr
#'   `scalar.within.metric.all` \tab equal across all cells \tab equal within `group1` \cr
#'   `scalar.all` \tab equal across all cells \tab equal across all cells \cr
#' }
#' The chain `configural` - `metric.within` - `metric.all` splits the overall
#' test of metric invariance into a part within the higher-level groups and a
#' part between them; `scalar.within.metric.all` does the same for the
#' intercepts. In addition, `pooled.configural`, `pooled.metric`, and
#' `pooled.scalar` are the conventional MGCFA across `group1` with the
#' lower-level groups pooled. Their test statistics are corrected for
#' clustering in the lower-level groups if the estimator allows it and every
#' higher-level group has at least as many lower-level groups as the model
#' has free parameters.
#'
#' **Intersectional mode.** For each stage:
#' \tabular{ll}{
#'   `configural` \tab parameters differ across all cells \cr
#'   `*.additive` \tab parameters differ by `group1` and by `group2`, without their interaction \cr
#'   `*.within.<group1>` \tab equal among cells sharing a level of `group1`, i.e. differ by `group1` only \cr
#'   `*.within.<group2>` \tab equal among cells sharing a level of `group2`, i.e. differ by `group2` only \cr
#'   `*.all` \tab equal across all cells \cr
#' }
#' Note that misfit of `*.within.<group1>` points to non-invariance across
#' `group2`, and the other way round. The two `within` models are not nested
#' in each other, but both are nested in the additive model, so that
#' `configural` vs. `additive` tests the interaction (non-invariance specific
#' to combinations of the two variables), and `additive` vs. `within` tests
#' the main effect of the other variable. No interaction is defined
#' multiplicatively for loadings
#' (\eqn{\lambda_{ab}\lambda_{11} = \lambda_{a1}\lambda_{1b}}), which does not
#' depend on the choice of marker, and additively for intercepts
#' (\eqn{\nu_{ab} + \nu_{11} = \nu_{a1} + \nu_{1b}}), which is independent of
#' the marker only when loadings are invariant.
#'
#' At the scalar stage, the stratified and additive models constrain both
#' loadings and intercepts, and are compared with their metric counterparts.
#'
#' **Output.** Comparisons include the (scaled, if applicable) chi-square
#' difference test, differences in CFI and RMSEA, and `rmsea.d`, the RMSEA
#' of the difference test (Savalei, Brace, & Fouladi, 2023). The tests of
#' invariance within strata are also broken down by stratum (`strata`), which
#' is exact because the strata have no parameters in common. `cell.discrepancy`
#' gives every cell's contribution to the chi-square of each model divided by
#' the cell's sample size; unlike the global fit, its change between models
#' is not dominated by large cells.
#'
#' **Comparability of fit indices.** All models that use cells as groups are
#' fitted to the same data, with the same saturated and baseline models, so
#' their chi-squares, CFI, RMSEA, SRMR, and BIC are comparable, and nested
#' ones can be tested against each other. The `pooled.*` models are fitted
#' to differently grouped data: their fit measures are comparable only among
#' themselves and should never be set beside those of the other models.
#' Differences in CFI and RMSEA depend on the number of constraints tested,
#' which varies a lot between the within and between steps; `rmsea.d` and
#' `chisq.diff/df.diff` are per degree of freedom and are better suited for
#' comparing steps. Per-stratum tests use the unscaled chi-square.
#'
#' **Limitations.** (1) Every cell is a group with its own parameters, so
#' cells must be large enough to estimate the model; with many small
#' lower-level units (e.g. actual classrooms) multilevel CFA is preferable.
#' (2) Conventional cutoffs for differences in CFI and RMSEA were
#' established for a few groups and global constraints and are not
#' calibrated for stratified models. (3) Tests between strata are
#' interpretable only if invariance holds within strata.
#'
#' @returns Invisibly, a list of class `stratifiedMI2` with elements
#'   `models` (fitted lavaan objects), `fit` (fit measures), `comparisons`
#'   (model comparisons), `strata` (within-strata tests by stratum),
#'   `cell.discrepancy` (per-cell misfit, see Details), `cells` (the cells used
#'   and their sizes), `dropped` (the cells dropped), `mode`, and `group`.
#'
#' @references Savalei, V., Brace, J. C., & Fouladi, R. T. (2023). We need to
#'   change how we compute RMSEA for nested model comparisons in structural
#'   equation modeling. *Psychological Methods*.
#'
#' @examples
#' \dontrun{
#' # schools nested in countries
#' stratifiedMI2("F =~ v1 + v2 + v3 + v4", data = Dat1,
#'               group1 = "country", group2 = "school",
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
                          additive = TRUE, pooled = TRUE,
                          min.n = NULL, group.partial = NULL,
                          print = TRUE, ...) {

  mode   <- match.arg(mode)
  stages <- match.arg(stages, several.ok = TRUE)
  dots   <- list(...)

  banned <- c("group", "group.equal", "std.lv", "ordered", "cluster")
  if (any(banned %in% names(dots)))
    stop("Arguments ", paste(sQuote(intersect(banned, names(dots))), collapse = ", "),
         " cannot be used with stratifiedMI2().")
  dots$do.fit <- dots$meanstructure <- NULL
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
  names(cells)[1] <- "cell"
  cells <- cells[order(cells[[group1]], cells[[group2]]), ]
  cells$n <- as.vector(table(data[[cell]])[cells$cell])
  rownames(cells) <- NULL

  if (mode == "hierarchical" && any(duplicated(cells[[group2]])))
    message("Some labels of '", group2, "' occur in more than one level of '",
            group1, "'; they are treated as different groups.")
  if (mode == "intersectional") {
    n.all <- length(unique(cells[[group1]])) * length(unique(cells[[group2]]))
    if (nrow(cells) < n.all)
      warning(n.all - nrow(cells), " of ", n.all, " combinations of '", group1,
              "' and '", group2, "' are empty.", call. = FALSE)
  }

  # Drop the cells that are too small to estimate the model
  one.group <- lavaan::lavaanify(model, meanstructure = TRUE, int.ov.free = TRUE,
                                 int.lv.free = FALSE, auto.fix.first = TRUE,
                                 auto.fix.single = TRUE, auto.var = TRUE,
                                 auto.cov.lv.x = TRUE, auto.cov.y = TRUE)
  npar <- sum(one.group$free > 0)
  if (is.null(min.n)) min.n <- npar
  small   <- cells$n < min.n
  dropped <- cells[small, ]
  if (any(small)) {
    warning(sum(small), " cell(s) with fewer than ", min.n, " observations were dropped: ",
            paste0(dropped$cell, " (n = ", dropped$n, ")", collapse = ", "), call. = FALSE)
    cells <- cells[!small, ]
    rownames(cells) <- NULL
    data  <- data[data[[cell]] %in% cells$cell, , drop = FALSE]
  }
  if (length(unique(cells[[group1]])) < 2 ||
      (mode == "intersectional" && length(unique(cells[[group2]])) < 2))
    stop("Each grouping variable should have at least two levels.")
  if (nrow(cells) <= length(unique(cells[[group1]])))
    stop("There should be more cells than levels of '", group1, "'.")

  crossed <- nrow(cells) ==
    length(unique(cells[[group1]])) * length(unique(cells[[group2]]))
  additive <- additive && mode == "intersectional"
  if (additive && !crossed) {
    warning("Additive models are skipped because some combinations of '", group1,
            "' and '", group2, "' are missing.", call. = FALSE)
    additive <- FALSE
  }

  # Baseline parameter table with marker identification -----------------
  base <- smi_cfa(c(list(model = model, data = data, group = cell,
                                      meanstructure = TRUE, do.fit = FALSE), dots))
  if (length(lavaan::lavNames(base, "ov.ord")) > 0)
    stop("Categorical indicators are not supported.")
  labs  <- lavaan::lavInspect(base, "group.label")
  cells <- cells[match(labs, cells$cell), ]
  rownames(cells) <- NULL
  pt0 <- smi_marker_partable(base)

  key  <- paste0(pt0$lhs, pt0$op, pt0$rhs)
  keep <- !key %in% gsub("\\s", "", group.partial)
  ov   <- lavaan::lavNames(base, "ov")
  keys <- list(
    loadings   = unique(key[pt0$op == "=~" & pt0$free > 0 & keep]),
    intercepts = unique(key[pt0$op == "~1" & pt0$lhs %in% ov & pt0$free > 0 & keep]))

  # A model is defined by the pattern of loadings and of intercepts:
  # "free", "all", "additive", or the name of a stratifying variable
  fit <- function(loadings = "free", intercepts = "free") {
    pt <- pt0
    spec <- list(loadings = loadings, intercepts = intercepts)
    for (what in names(spec)) {
      s <- spec[[what]]
      if (s == "free" || length(keys[[what]]) == 0) next
      pt <- rbind(pt, if (s == "additive") {
        smi_additive_rows(pt0, key, keys[[what]], cells[[group1]], cells[[group2]],
                          sign = if (what == "loadings") "*" else "+")
      } else {
        smi_equality_rows(pt0, key, keys[[what]],
                          stratum = if (s == "all") rep(1, nrow(cells)) else cells[[s]])
      })
    }
    pt$id <- seq_len(nrow(pt))
    smi_try(smi_cfa(c(list(model = pt, data = data, group = cell), dots)))
  }

  models <- list()
  comps  <- list()
  compare <- function(stage, from, to, test) {
    comps[[length(comps) + 1]] <<- data.frame(stage = stage, from = from, to = to,
                                              test = test, stringsAsFactors = FALSE)
  }
  scalar <- "scalar" %in% stages
  breakdown <- list()   # within-strata tests to be split by stratum

  # Models ----------------------------------------------------------------
  if (mode == "hierarchical") {

    models$configural    <- fit()
    models$metric.within <- fit(group1)
    models$metric.all    <- fit("all")
    compare("metric", "configural", "metric.within",
            paste0("loadings across ", group2, " within ", group1))
    compare("metric", "metric.within", "metric.all",
            paste0("loadings across ", group1))
    compare("metric", "configural", "metric.all", "loadings across all cells (overall)")
    breakdown$metric <- list("configural", "metric.within", group1, "loadings")

    if (scalar) {
      models$scalar.within            <- fit(group1, group1)
      models$scalar.within.metric.all <- fit("all", group1)
      models$scalar.all               <- fit("all", "all")
      compare("scalar", "metric.within", "scalar.within",
              paste0("intercepts across ", group2, " within ", group1))
      compare("scalar", "scalar.within", "scalar.within.metric.all",
              paste0("loadings across ", group1, ", given scalar within"))
      compare("scalar", "metric.all", "scalar.within.metric.all",
              paste0("intercepts across ", group2, " within ", group1, ", given metric overall"))
      compare("scalar", "scalar.within.metric.all", "scalar.all",
              paste0("intercepts across ", group1))
      compare("scalar", "metric.all", "scalar.all", "intercepts across all cells (overall)")
      breakdown$scalar <- list("metric.within", "scalar.within", group1, "intercepts")
    }

    if (pooled) {
      # correction for clustering needs enough clusters in every group
      robust <- min(table(cells[[group1]])) >= npar
      if (!robust)
        message("Pooled models are not corrected for clustering in '", group2,
                "': fewer than ", npar, " (number of free parameters) levels of '",
                group2, "' in some levels of '", group1, "'.")
      pool <- function(equal = NULL) {
        args <- c(list(model = model, data = data, group = group1,
                       group.partial = group.partial), dots)
        args$group.equal <- equal
        fit.robust <- if (robust) suppressWarnings(tryCatch(
          smi_cfa(c(args, list(cluster = cell))), error = function(e) NULL))
        if (!is.null(fit.robust)) fit.robust else smi_try(smi_cfa(args))
      }
      models$pooled.configural <- pool()
      models$pooled.metric     <- pool("loadings")
      compare("metric", "pooled.configural", "pooled.metric",
              paste0("pooled: loadings across ", group1))
      if (scalar) {
        models$pooled.scalar <- pool(c("loadings", "intercepts"))
        compare("scalar", "pooled.metric", "pooled.scalar",
                paste0("pooled: intercepts across ", group1))
      }
    }

  } else {

    w <- paste0("within.", c(group1, group2))
    g <- c(group1, group2)
    m <- function(stage, suffix) paste0(stage, ".", suffix)

    models$configural <- fit()
    if (additive) models$metric.additive <- fit("additive")
    models[[m("metric", w[1])]] <- fit(g[1])
    models[[m("metric", w[2])]] <- fit(g[2])
    models$metric.all <- fit("all")
    if (scalar) {
      if (additive) models$scalar.additive <- fit("additive", "additive")
      models[[m("scalar", w[1])]] <- fit(g[1], g[1])
      models[[m("scalar", w[2])]] <- fit(g[2], g[2])
      models$scalar.all <- fit("all", "all")
    }

    for (stage in c("metric", if (scalar) "scalar")) {
      par   <- if (stage == "metric") "loadings" else "intercepts"
      par2  <- if (stage == "metric") "loadings" else "loadings & intercepts"
      start <- function(suffix) if (stage == "metric") "configural" else m("metric", suffix)
      if (additive) {
        compare(stage, start("additive"), m(stage, "additive"),
                paste0(par, ": interaction of ", g[1], " and ", g[2]))
        for (i in 1:2)
          compare(stage, m(stage, "additive"), m(stage, w[i]),
                  paste0(par2, ": main effect of ", g[3 - i]))
      }
      for (i in 1:2) {
        compare(stage, start(w[i]), m(stage, w[i]),
                paste0(par, " across ", g[3 - i], " within ", g[i]))
        compare(stage, m(stage, w[i]), m(stage, "all"),
                paste0(par2, " across ", g[i], ", given no effect of ", g[3 - i]))
        breakdown[[paste0(stage, i)]] <- list(start(w[i]), m(stage, w[i]), g[i], par)
      }
      compare(stage, start("all"), m(stage, "all"),
              paste0(par, " across all cells (overall)"))
    }
  }

  # Summaries -------------------------------------------------------------
  fits  <- as.data.frame(t(sapply(models, smi_fit)))
  comps <- do.call(rbind, comps)
  comps <- cbind(comps, t(mapply(
    function(from, to) smi_compare(models[[from]], models[[to]], fits[from, ], fits[to, ]),
    comps$from, comps$to)))
  rownames(comps) <- NULL

  cell.stat <- sapply(models[!grepl("^pooled", names(models))], smi_group_stat,
                      n.groups = nrow(cells))
  rownames(cell.stat) <- cells$cell

  strata <- do.call(rbind, lapply(names(breakdown), function(nm) {
    b <- breakdown[[nm]]
    d <- cell.stat[, b[[2]]] - cell.stat[, b[[1]]]
    s <- cells[[b[[3]]]]
    k <- tapply(d, s, length)
    out <- data.frame(stage = sub("[12]$", "", nm), by = b[[3]], stratum = names(k),
                      cells = as.vector(k), n = as.vector(tapply(cells$n, s, sum)),
                      chisq.diff = as.vector(tapply(d, s, sum)),
                      df.diff = as.vector(k - 1) * length(keys[[b[[4]]]]),
                      stringsAsFactors = FALSE)
    out <- out[out$df.diff > 0, ]
    out$p <- stats::pchisq(out$chisq.diff, out$df.diff, lower.tail = FALSE)
    out$rmsea.d <- sqrt(pmax(0, (out$chisq.diff - out$df.diff) / (out$df.diff * out$n))) *
      sqrt(out$cells)
    out
  }))
  rownames(strata) <- NULL

  out <- structure(list(models = models, fit = fits, comparisons = comps,
                        strata = strata, cell.discrepancy = cell.stat / cells$n,
                        cells = cells, dropped = dropped, mode = mode,
                        group = c(group1, group2)),
                   class = "stratifiedMI2")
  if (print) print(out)
  invisible(out)
}


#' Fit a model, turning an error into a warning
#' @noRd
smi_try <- function(expr) {
  tryCatch(expr, error = function(e) {
    warning("A model could not be fitted: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}


#' Parameter table with marker identification of loadings and intercepts
#'
#' @param base lavaan object (may be unfitted) of a multiple group model with
#'   mean structure and without equality constraints.
#' @returns Parameter table where in each group the intercept of every factor's
#'   marker indicator is fixed to 0 and the factor's mean is free.
#' @noRd
smi_marker_partable <- function(base) {

  pt <- as.data.frame(lavaan::parTable(base), stringsAsFactors = FALSE)
  pt$start <- pt$est <- pt$se <- NULL
  lv <- lavaan::lavNames(base, "lv")
  if (any(pt$op == "=~" & pt$rhs %in% lv))
    stop("Higher-order factors are not supported.")

  load   <- pt$op == "=~" & pt$group == pt$group[1]
  fixed  <- load & pt$free == 0 & !is.na(pt$ustart) & pt$ustart != 0
  marker <- sapply(lv, function(f) pt$rhs[fixed & pt$lhs == f][1])
  if (anyNA(marker))
    stop("Every factor should have a marker indicator, i.e. a loading fixed to a non-zero value: ",
         paste(lv[is.na(marker)], collapse = ", "))
  if (any(duplicated(marker)))
    stop("Factors should have different marker indicators.")

  int <- pt$op == "~1"
  pt$free[int & pt$lhs %in% marker]   <- 0
  pt$ustart[int & pt$lhs %in% marker] <- 0
  pt$free[int & pt$lhs %in% lv]       <- 1
  pt$ustart[int & pt$lhs %in% lv]     <- NA

  free <- pt$free > 0
  pt$free[free] <- seq_len(sum(free))
  pt
}


#' Rows of a parameter table with constraints
#' @param pt parameter table, `key` its parameter names (lhs, op, rhs pasted).
#' @param keys names of the parameters to constrain.
#' @noRd
smi_constraint_rows <- function(pt, lhs, rhs) {
  new <- pt[rep(1, length(lhs)), ]
  new[] <- lapply(new, function(x) if (is.character(x)) "" else x * 0L)
  new$lhs <- lhs
  new$rhs <- rhs
  new$op  <- "=="
  new$user   <- 1L
  new$ustart <- NA
  new
}

#' Equality of parameters across groups within strata
#' @param stratum stratum of every group, in the order of the groups.
#' @noRd
smi_equality_rows <- function(pt, key, keys, stratum) {
  ref <- match(stratum, stratum)   # first group of each group's stratum
  g   <- which(ref != seq_along(ref))
  d   <- expand.grid(g = g, key = keys, stringsAsFactors = FALSE)
  plabel <- function(group, k) pt$plabel[match(paste(group, k), paste(pt$group, key))]
  smi_constraint_rows(pt, plabel(ref[d$g], d$key), plabel(d$g, d$key))
}

#' No interaction of two crossed grouping variables in parameters
#' @param a,b levels of the two grouping variables for every group.
#' @param sign "*" for multiplicative (loadings), "+" for additive (intercepts).
#' @noRd
smi_additive_rows <- function(pt, key, keys, a, b, sign) {
  a1 <- a == a[1]
  b1 <- b == b[1]
  g  <- which(!a1 & !b1)
  d  <- expand.grid(g = g, key = keys, stringsAsFactors = FALSE)
  g11 <- which(a1 & b1)
  ga1 <- match(paste(a[d$g], b[1]), paste(a, b))
  g1b <- match(paste(a[1], b[d$g]), paste(a, b))
  plabel <- function(group, k) pt$plabel[match(paste(group, k), paste(pt$group, key))]
  smi_constraint_rows(pt,
                      paste0(plabel(d$g, d$key), sign, plabel(g11, d$key)),
                      paste0(plabel(ga1, d$key), sign, plabel(g1b, d$key)))
}


#' Fit measures of one model, scaled versions if available
#' @noRd
smi_fit <- function(x) {
  nms <- c("chisq", "df", "cfi", "rmsea", "srmr", "bic")
  out <- setNames(rep(NA_real_, length(nms)), nms)
  if (is.null(x) || !lavaan::lavInspect(x, "converged")) return(out)
  fm <- tryCatch(lavaan::fitMeasures(x), error = function(e) NULL)
  if (is.null(fm)) return(out)
  # scaled version of a measure if it is available, the standard one otherwise
  get <- function(nm) {
    val <- fm[intersect(c(paste0(nm, ".scaled"), nm), names(fm))]
    val <- val[!is.na(val)]
    if (length(val) == 0) NA_real_ else val[[1]]
  }
  out[] <- sapply(nms, get)
  out
}

#' Contributions of the groups to the (unscaled) chi-square
#' @noRd
smi_group_stat <- function(x, n.groups) {
  out <- rep(NA_real_, n.groups)
  if (is.null(x) || !lavaan::lavInspect(x, "converged")) return(out)
  stat <- tryCatch(lavaan::lavInspect(x, "test")[[1]]$stat.group, error = function(e) NULL)
  if (length(stat) == n.groups) out[] <- stat
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
  op <- options(width = max(getOption("width"), 200))
  on.exit(options(op))
  show <- function(d) {
    num <- sapply(d, is.numeric)
    d[num] <- round(d[num], digits)
    print(d, row.names = FALSE, right = FALSE, na.print = "")
  }
  cat("\nStratified measurement invariance,", x$mode, "mode\n")
  cat("  Groups: ", nrow(x$cells), " combinations of ", g[1], " (",
      length(unique(x$cells[[g[1]]])), ") and ", g[2], ", n from ", min(x$cells$n),
      " to ", max(x$cells$n), "\n", sep = "")
  if (nrow(x$dropped) > 0)
    cat("  Dropped as too small:", paste(x$dropped$cell, collapse = ", "), "\n")
  if (x$mode == "hierarchical") {
    cat("  *.within : equal across ", g[2], " within ", g[1], ", free across ", g[1], "\n", sep = "")
    if (any(grepl("^pooled", rownames(x$fit))))
      cat("  pooled.* : groups are levels of ", g[1],
          "; fit is NOT comparable with the other models\n", sep = "")
  } else {
    for (i in 1:2)
      cat("  *.within.", g[i], " : equal across ", g[3 - i], " within ", g[i],
          " (differ by ", g[i], " only)\n", sep = "")
    if ("metric.additive" %in% rownames(x$fit))
      cat("  *.additive : differ by ", g[1], " and by ", g[2], ", no interaction\n", sep = "")
  }

  cat("\nModel fit\n")
  print(round(x$fit, digits), na.print = "")

  for (s in unique(x$comparisons$stage)) {
    cat("\nComparisons,", s, "stage\n")
    show(x$comparisons[x$comparisons$stage == s, -1])
    for (b in unique(x$strata$by[x$strata$stage == s])) {
      cat("\n  Invariance within", b, "by stratum,", s, "stage\n")
      show(x$strata[x$strata$stage == s & x$strata$by == b, -(1:2)])
    }
  }
  cat("\n")
  invisible(x)
}


#' Call lavaan::cfa with a list of arguments
#' (cfa() inspects its own call, so it has to be called by name)
#' @noRd
smi_cfa <- function(args) {
  do.call("cfa", args, envir = asNamespace("lavaan"))
}
