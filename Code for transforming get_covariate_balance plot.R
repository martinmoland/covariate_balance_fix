#####Open trace window for function and allow for edits#######
trace("get_covariate_balance", edit = TRUE)



###Copy the following into the trace######
function (matched.sets, data, covariates, use.equal.weights = FALSE, 
          verbose = TRUE, plot = FALSE, reference.line = TRUE, legend = TRUE, 
          ylab = "SD", ...) 
{
  calculate.network.proportion.balance = FALSE
  calculate.network.count.balance = FALSE
  adjacency.matrix = NULL
  neighborhood.degree = NULL
  continuous.treatment = FALSE
  if (is.null(covariates)) {
    stop("please specify the covariates for which you would like to check the balance")
  }
  if (!all(covariates %in% colnames(data))) {
    stop("Some of the specified covariates are not columns in the data set.")
  }
  if (!inherits(matched.sets, "matched.set")) 
    stop("Please pass a matched.set object")
  unit.id <- attr(matched.sets, "id.var")
  time.id <- attr(matched.sets, "t.var")
  lag <- attr(matched.sets, "lag")
  treatment <- attr(matched.sets, "treatment.var")
  if (!inherits(data[, unit.id], "integer") && !inherits(data[, 
                                                              unit.id], "numeric")) 
    stop("please convert unit id column to integer or numeric")
  if (!inherits(data[, time.id], "integer")) 
    stop("please convert time id to consecutive integers")
  if (any(table(data[, unit.id]) != max(table(data[, unit.id])))) {
    testmat <- data.table::dcast(data.table::as.data.table(data), 
                                 formula = paste0(unit.id, "~", time.id), value.var = treatment)
    d <- data.table::melt(data.table(testmat), id = unit.id, 
                          variable = time.id, value = treatment, variable.factor = FALSE, 
                          value.name = treatment)
    d <- data.frame(d)[, c(1, 2)]
    class(d[, 2]) <- "integer"
    data <- merge(data.table(d), data.table(data), all.x = TRUE, 
                  by = c(unit.id, time.id))
    data <- as.data.frame(data)
  }
  ordered.data <- data[order(data[, unit.id], data[, time.id]), 
  ]
  if (calculate.network.proportion.balance || calculate.network.count.balance) {
    if (is.null(adjacency.matrix)) {
      stop("Please provide adjacency matrix")
    }
    ordered.data <- calculate_neighbor_treatment(data = ordered.data, 
                                                 edge.matrix = adjacency.matrix, n.degree = neighborhood.degree, 
                                                 unit.id = unit.id, time.id = time.id, treatment.variable = treatment)
    if (!is.null(calculate.network.proportion.balance)) {
      covariates <- c(covariates, make.names(paste0("neighborhood_t_prop", 
                                                    ".", 1:neighborhood.degree)))
    }
    if (!is.null(calculate.network.count.balance)) {
      covariates <- c(covariates, make.names(paste0("neighborhood_t_count", 
                                                    ".", 1:neighborhood.degree)))
    }
  }
  matched.sets <- matched.sets[sapply(matched.sets, length) > 
                                 0]
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% 
                                        c(time.id, unit.id, treatment)]
  othercols <- othercols[othercols %in% covariates]
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, 
                                   othercols), drop = FALSE]
  ordered.data <- ordered.data[, unique(c(unit.id, time.id, 
                                          treatment, covariates)), drop = FALSE]
  if (is.null(attr(matched.sets[[1]], "weights")) | use.equal.weights) {
    for (i in 1:length(matched.sets)) {
      attr(matched.sets[[i]], "weights") <- rep(1/length(matched.sets[[i]]), 
                                                length(matched.sets[[i]]))
      names(attr(matched.sets[[i]], "weights")) <- matched.sets[[i]]
    }
  }
  treated.ts <- as.integer(sub(".*\\.", "", names(matched.sets)))
  treated.ids <- as.integer(sub("\\..*", "", names(matched.sets)))
  tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
  idxlist <- get_yearly_dmats(as.matrix(ordered.data), treated.ids, 
                              tlist, matched_sets = matched.sets, lag)
  balance_mats <- build_balance_mats(ordered_expanded_data = ordered.data, 
                                     idx = idxlist, msets = matched.sets)
  unlistedmats <- unlist(balance_mats, recursive = F)
  plotpoints <- list()
  for (k in 1:(lag + 1)) {
    var.points <- list()
    for (i in 1:length(covariates)) {
      variable <- covariates[i]
      sd.val <- sd(sapply(unlistedmats[seq(from = k, to = (length(matched.sets) * 
                                                             (lag + 1)), by = lag + 1)], function(x) {
                                                               x[nrow(x), variable]
                                                             }), na.rm = T)
      if (isTRUE(all.equal(sd.val, 0))) {
        sd.val <- NA
      }
      tprd <- unlistedmats[seq(from = k, to = (length(matched.sets) * 
                                                 (lag + 1)), by = lag + 1)]
      get_mean_difs <- function(x, variable) {
        return(x[nrow(x), variable] - sum(x[1:(nrow(x) - 
                                                 1), "weights"] * x[1:(nrow(x) - 1), variable], 
                                          na.rm = T))
      }
      diffs <- sapply(tprd, get_mean_difs, variable = variable)
      var.points[[i]] <- mean(diffs/sd.val, na.rm = T)
    }
    names(var.points) <- covariates
    plotpoints[[k]] <- var.points
  }
  names(plotpoints) <- paste0("t_", lag:0)
  pointmatrix <- apply((as.matrix(do.call(rbind, plotpoints))), 
                       2, function(x) {
                         (as.numeric(x))
                       }, simplify = TRUE)
  rownames(pointmatrix) <- names(plotpoints)
  remove.vars.idx <- apply(apply(pointmatrix, 2, is.nan), 2, 
                           any)
  if (sum(remove.vars.idx) > 0) {
    removed.vars <- names(which(apply(apply(pointmatrix, 
                                            2, is.nan), 2, any)))
    pointmatrix <- pointmatrix[, !remove.vars.idx]
    warning(paste0("Some variables were removed due to low variation, inadequate data needed for calculation: ", 
                   removed.vars))
  }
  pointmatrix <- pointmatrix[-nrow(pointmatrix), , drop = FALSE]
  if (!plot) 
    return(pointmatrix)
  if (plot) { 
    treated.included <- treatment %in% colnames(pointmatrix)
    if (!continuous.treatment) {
      if (treated.included) {
        treated.data <- pointmatrix[, which(colnames(pointmatrix) == 
                                              treatment)]
        pointmatrix <- pointmatrix[, -which(colnames(pointmatrix) == 
                                              treatment)]
        graphics::matplot(pointmatrix, type = "l", col = "black", #This transforms all lines to black
                          lty = 1:nrow(pointmatrix)  #This changes linetypes for n of the pointmatrix, 
                          ylab = ylab, xaxt = "n", 
                          ...)
        graphics::lines(x = 1:nrow(pointmatrix), y = as.numeric(treated.data), 
                        type = "l", lty = 2, lwd = 3)
        graphics::axis(side = 1, labels = paste0("t-", 
                                                 (nrow(pointmatrix)):1), at = 1:nrow(pointmatrix))
      }
      else {
        graphics::matplot(pointmatrix, type = "l", col = "black", #This transforms all lines to black 
                          lty = 1:nrow(pointmatrix), #This changes linetypes for n of the pointmatrix 
                          ylab = ylab, xaxt = "n", 
                          ...)
        graphics::axis(side = 1, labels = paste0("t-", 
                                                 (nrow(pointmatrix)):1), at = 1:nrow(pointmatrix))
      }
    }
    else {
      if (treated.included) {
        treated.data <- pointmatrix[, which(colnames(pointmatrix) == 
                                              treatment)]
        pointmatrix <- pointmatrix[, -which(colnames(pointmatrix) == 
                                              treatment)]
        graphics::matplot(pointmatrix, type = "l", col = "black", #Same as above 
                          lty = 1:ncol(pointmatrix), #Same as above 
                          ylab = ylab, xaxt = "n", 
                          ...)
        graphics::lines(x = 1:(nrow(pointmatrix) - 1), 
                        y = as.numeric(treated.data)[1:(nrow(pointmatrix) - 
                                                          1)], type = "l", lty = 2, lwd = 3)
        graphics::axis(side = 1, labels = paste0("t-", 
                                                 (nrow(pointmatrix)):1), at = 1:nrow(pointmatrix))
      }
      else {
        graphics::matplot(pointmatrix, type = "l", col = 1:ncol(pointmatrix), 
                          lty = 1, ylab = ylab, xaxt = "n", ...)
        graphics::axis(side = 1, labels = paste0("t-", 
                                                 (nrow(pointmatrix)):1), at = 1:nrow(pointmatrix))
      }
    }
    if (legend) {
      if (treated.included) {
        legend("topleft", legend = c(colnames(pointmatrix), 
                                     treatment), col = "black", lty = c(1:ncol(pointmatrix), #This varies legend solely by linetype
                                                                        2))
      }
      else {
        legend("topleft", legend = colnames(pointmatrix), 
               col = "black", lty = 1:(ncol(pointmatrix))) #Legend varies solely by linetype
      }
    }
    if (reference.line) 
      graphics::abline(h = 0, lty = "longdash")
  }
}
