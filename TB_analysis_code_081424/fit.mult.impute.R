fit.mult.impute =
  function (formula, data, weights, subset, na.action, init, control, 
          ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
          robust, model = FALSE, x = FALSE, y = TRUE, tt, method = ties, 
          id, cluster, istate, statedata, nocenter = c(-1, 0, 1), ...) 
{
  ties <- match.arg(ties)
  Call <- match.call()
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxph.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(indx == 0L)) 
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx == 
                                                                  0L]), domain = NA)
  }
  if (missing(control)) 
    control <- coxph.control(...)
  if (missing(formula)) 
    stop("a formula argument is required")
  ss <- c("cluster", "offset")
  if (is.list(formula)) 
    Terms <- if (missing(data)) 
      terms(formula[[1]], specials = ss)
  else terms(formula[[1]], specials = ss, data = data)
  else Terms <- if (missing(data)) 
    terms(formula, specials = ss)
  else terms(formula, specials = ss, data = data)
  tcl <- attr(Terms, "specials")$cluster
  if (length(tcl) > 1) 
    stop("a formula cannot have multiple cluster terms")
  if (length(tcl) > 0) {
    factors <- attr(Terms, "factors")
    if (any(factors[tcl, ] > 1)) 
      stop("cluster() cannot be in an interaction")
    if (attr(Terms, "response") == 0) 
      stop("formula must have a Surv response")
    temp <- attr(Terms, "term.labels")
    oo <- attr(Terms, "specials")$offset
    if (!is.null(oo)) {
      ooterm <- rownames(factors)[oo]
      if (oo < tcl) 
        temp <- c(ooterm, temp)
      else temp <- c(temp, ooterm)
    }
    if (is.null(Call$cluster)) 
      Call$cluster <- attr(Terms, "variables")[[1 + tcl]][[2]]
    else warning("cluster appears both in a formula and as an argument, formula term ignored")
    if (is.list(formula)) 
      formula[[1]][[3]] <- reformulate(temp[1 - tcl])[[2]]
    else formula[[3]] <- reformulate(temp[1 - tcl])[[2]]
    Call$formula <- formula
  }
  indx <- match(c("formula", "data", "weights", "subset", "na.action", 
                  "cluster", "id", "istate"), names(Call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A formula argument is required")
  tform <- Call[c(1, indx)]
  tform[[1L]] <- quote(stats::model.frame)
  if (is.list(formula)) {
    multiform <- TRUE
    dformula <- formula[[1]]
    if (missing(statedata)) 
      covlist <- parsecovar1(formula[-1])
    else {
      if (!inherits(statedata, "data.frame")) 
        stop("statedata must be a data frame")
      if (is.null(statedata$state)) 
        stop("statedata data frame must contain a 'state' variable")
      covlist <- parsecovar1(formula[-1], names(statedata))
    }
    tlab <- unlist(lapply(covlist$rhs, function(x) attr(terms.formula(x$formula), 
                                                        "term.labels")))
    tlab <- c(attr(terms.formula(dformula), "term.labels"), 
              tlab)
    newform <- reformulate(tlab, dformula[[2]])
    environment(newform) <- environment(dformula)
    formula <- newform
    tform$na.action <- na.pass
  }
  else {
    multiform <- FALSE
    covlist <- NULL
    dformula <- formula
  }
  special <- c("strata", "tt", "frailty", "ridge", "pspline")
  tform$formula <- if (missing(data)) 
    terms(formula, special)
  else terms(formula, special, data = data)
  if (!is.null(attr(tform$formula, "specials")$tt)) {
    coxenv <- new.env(parent = environment(formula))
    assign("tt", function(x) x, envir = coxenv)
    environment(tform$formula) <- coxenv
  }
  mf <- eval(tform, parent.frame())
  Terms <- terms(mf)
  n <- nrow(mf)
  Y <- model.response(mf)
  isSurv2 <- inherits(Y, "Surv2")
  if (isSurv2) {
    if (length(attr(mf, "na.action"))) {
      tform$na.action <- na.pass
      mf <- eval.parent(tform)
    }
    if (!is.null(attr(Terms, "specials")$cluster)) 
      stop("cluster() cannot appear in the model statement")
    new <- surv2data(mf)
    mf <- new$mf
    istate <- new$istate
    id <- new$id
    Y <- new$y
    n <- nrow(mf)
  }
  else {
    if (!is.Surv(Y)) 
      stop("Response must be a survival object")
    id <- model.extract(mf, "id")
    istate <- model.extract(mf, "istate")
  }
  if (n == 0) 
    stop("No (non-missing) observations")
  type <- attr(Y, "type")
  multi <- FALSE
  if (type == "mright" || type == "mcounting") 
    multi <- TRUE
  else if (type != "right" && type != "counting") 
    stop(paste("Cox model doesn't support \"", type, "\" survival data", 
               sep = ""))
  data.n <- nrow(Y)
  if (!multi && multiform) 
    stop("formula is a list but the response is not multi-state")
  if (multi && length(attr(Terms, "specials")$frailty) > 0) 
    stop("multi-state models do not currently support frailty terms")
  if (multi && length(attr(Terms, "specials")$pspline) > 0) 
    stop("multi-state models do not currently support pspline terms")
  if (multi && length(attr(Terms, "specials")$ridge) > 0) 
    stop("multi-state models do not currently support ridge penalties")
  if (control$timefix) 
    Y <- aeqSurv(Y)
  if (length(attr(Terms, "variables")) > 2) {
    ytemp <- terms.inner(formula[1:2])
    xtemp <- terms.inner(formula[-2])
    if (any(!is.na(match(xtemp, ytemp)))) 
      warning("a variable appears on both the left and right sides of the formula")
  }
  strats <- attr(Terms, "specials")$strata
  hasinteractions <- FALSE
  dropterms <- NULL
  if (length(strats)) {
    stemp <- untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) == 1) 
      strata.keep <- mf[[stemp$vars]]
    else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
    istrat <- as.integer(strata.keep)
    for (i in stemp$vars) {
      if (any(attr(Terms, "order")[attr(Terms, "factors")[i, 
                                                          ] > 0] > 1)) 
        hasinteractions <- TRUE
    }
    if (!hasinteractions) 
      dropterms <- stemp$terms
  }
  else istrat <- NULL
  if (hasinteractions && multi) 
    stop("multi-state coxph does not support strata*covariate interactions")
  timetrans <- attr(Terms, "specials")$tt
  if (missing(tt)) 
    tt <- NULL
  if (length(timetrans)) {
    if (multi || isSurv2) 
      stop("the tt() transform is not implemented for multi-state or Surv2 models")
    timetrans <- untangle.specials(Terms, "tt")
    ntrans <- length(timetrans$terms)
    if (is.null(tt)) {
      tt <- function(x, time, riskset, weights) {
        obrien <- function(x) {
          r <- rank(x)
          (r - 0.5)/(0.5 + length(r) - r)
        }
        unlist(tapply(x, riskset, obrien))
      }
    }
    if (is.function(tt)) 
      tt <- list(tt)
    if (is.list(tt)) {
      if (any(!sapply(tt, is.function))) 
        stop("The tt argument must contain function or list of functions")
      if (length(tt) != ntrans) {
        if (length(tt) == 1) {
          temp <- vector("list", ntrans)
          for (i in 1:ntrans) temp[[i]] <- tt[[1]]
          tt <- temp
        }
        else stop("Wrong length for tt argument")
      }
    }
    else stop("The tt argument must contain a function or list of functions")
    if (ncol(Y) == 2) {
      if (length(strats) == 0) {
        sorted <- order(-Y[, 1], Y[, 2])
        newstrat <- rep.int(0L, nrow(Y))
        newstrat[1] <- 1L
      }
      else {
        sorted <- order(istrat, -Y[, 1], Y[, 2])
        newstrat <- as.integer(c(1, 1 * (diff(istrat[sorted]) != 
                                           0)))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount1, Y[sorted, ], as.integer(newstrat))
      tindex <- sorted[counts$index]
    }
    else {
      if (length(strats) == 0) {
        sort.end <- order(-Y[, 2], Y[, 3])
        sort.start <- order(-Y[, 1])
        newstrat <- c(1L, rep(0, nrow(Y) - 1))
      }
      else {
        sort.end <- order(istrat, -Y[, 2], Y[, 3])
        sort.start <- order(istrat, -Y[, 1])
        newstrat <- c(1L, as.integer(diff(istrat[sort.end]) != 
                                       0))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount2, Y, as.integer(sort.start - 
                                                  1L), as.integer(sort.end - 1L), as.integer(newstrat))
      tindex <- counts$index
    }
    Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
    type <- "right"
    mf <- mf[tindex, ]
    istrat <- rep(1:length(counts$nrisk), counts$nrisk)
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights))) 
      stop("weights must be finite")
    tcall <- attr(Terms, "variables")[timetrans$terms + 2]
    pvars <- attr(Terms, "predvars")
    pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
    for (i in 1:ntrans) {
      newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[, 1], 
                         istrat, weights)
      mf[[timetrans$var[i]]] <- newtt
      nclass <- class(newtt)
      if (any(nclass %in% pmethod)) {
        dummy <- as.call(list(as.name(class(newtt)[1]), 
                              tcall[[i]][[2]]))
        ptemp <- makepredictcall(newtt, dummy)
        pvars[[timetrans$terms[i] + 2]] <- ptemp
      }
    }
    attr(Terms, "predvars") <- pvars
  }
  xlevels <- .getXlevels(Terms, mf)
  cluster <- model.extract(mf, "cluster")
  weights <- model.weights(mf)
  has.cluster <- !(missing(cluster) || length(cluster) == 0)
  has.id <- !(missing(id) || length(id) == 0)
  has.rwt <- (!is.null(weights) && any(weights != floor(weights)))
  has.robust <- (!missing(robust) && !is.null(robust))
  if (has.id) 
    id <- as.factor(id)
  if (missing(robust) || is.null(robust)) {
    if (has.cluster || has.rwt || (has.id && anyDuplicated(id[Y[, 
                                                                ncol(Y)] == 1]))) 
      robust <- TRUE
    else robust <- FALSE
  }
  if (!is.logical(robust)) 
    stop("robust must be TRUE/FALSE")
  if (has.cluster) {
    if (!robust) {
      warning("cluster specified with robust=FALSE, cluster ignored")
      ncluster <- 0
      clname <- NULL
    }
    else {
      if (is.factor(cluster)) {
        clname <- levels(cluster)
        cluster <- as.integer(cluster)
      }
      else {
        clname <- sort(unique(cluster))
        cluster <- match(cluster, clname)
      }
      ncluster <- length(clname)
    }
  }
  else {
    if (robust && has.id) {
      clname <- levels(id)
      cluster <- as.integer(id)
      ncluster <- length(clname)
    }
    else {
      ncluster <- 0
    }
  }
  if (robust && is.null(cluster)) {
    if (ncol(Y) == 2 || !has.robust) 
      cluster <- seq.int(1, nrow(mf))
    else stop("one of cluster or id is needed")
  }
  contrast.arg <- NULL
  attr(Terms, "intercept") <- 1
  if (multi) {
    if (length(id) == 0) 
      stop("an id statement is required for multi-state models")
    mcheck <- survcheck2(Y, id, istate)
    if (mcheck$flag["overlap"] > 0) 
      stop("data set has overlapping intervals for one or more subjects")
    transitions <- mcheck$transitions
    istate <- mcheck$istate
    states <- mcheck$states
    if (missing(statedata)) 
      covlist2 <- parsecovar2(covlist, NULL, dformula = dformula, 
                              Terms, transitions, states)
    else covlist2 <- parsecovar2(covlist, statedata, dformula = dformula, 
                                 Terms, transitions, states)
    tmap <- covlist2$tmap
    if (!is.null(covlist)) {
      good.tran <- bad.tran <- rep(FALSE, nrow(Y))
      termname <- rownames(attr(Terms, "factors"))
      trow <- (!is.na(match(rownames(tmap), termname)))
      termiss <- matrix(0L, nrow(mf), ncol(mf))
      for (i in 1:ncol(mf)) {
        xx <- is.na(mf[[i]])
        if (is.matrix(xx)) 
          termiss[, i] <- apply(xx, 1, any)
        else termiss[, i] <- xx
      }
      for (i in levels(istate)) {
        rindex <- which(istate == i)
        j <- which(covlist2$mapid[, 1] == match(i, states))
        for (jcol in j) {
          k <- which(trow & tmap[, jcol] > 0)
          bad.tran[rindex] <- (bad.tran[rindex] | apply(termiss[rindex, 
                                                                k, drop = FALSE], 1, any))
          good.tran[rindex] <- (good.tran[rindex] | apply(!termiss[rindex, 
                                                                   k, drop = FALSE], 1, all))
        }
      }
      n.partially.used <- sum(good.tran & bad.tran & !is.na(Y))
      omit <- (!good.tran & bad.tran) | is.na(Y)
      if (all(omit)) 
        stop("all observations deleted due to missing values")
      temp <- setNames(seq(omit)[omit], attr(mf, "row.names")[omit])
      attr(temp, "class") <- "omit"
      mf <- mf[!omit, , drop = FALSE]
      attr(mf, "na.action") <- temp
      Y <- Y[!omit]
      id <- id[!omit]
      if (length(istate)) 
        istate <- istate[!omit]
    }
  }
  if (length(dropterms)) {
    Terms2 <- Terms[-dropterms]
    X <- model.matrix(Terms2, mf, constrasts.arg = contrast.arg)
    temp <- attr(X, "assign")
    shift <- sort(dropterms)
    for (i in seq(along = shift)) temp <- temp + 1 * (shift[i] <= 
                                                        temp)
    attr(X, "assign") <- temp
  }
  else X <- model.matrix(Terms, mf, contrasts.arg = contrast.arg)
  Xatt <- attributes(X)
  if (hasinteractions) 
    adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  else adrop <- 0
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset == 0)) 
    offset <- rep(0, nrow(mf))
  else if (any(!is.finite(offset) | !is.finite(exp(offset)))) 
    stop("offsets must lead to a finite risk score")
  weights <- model.weights(mf)
  if (!is.null(weights) && any(!is.finite(weights))) 
    stop("weights must be finite")
  assign <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  if (sum(Y[, ncol(Y)]) == 0) {
    ncoef <- ncol(X)
    ctemp <- rep(NA, ncoef)
    names(ctemp) <- colnames(X)
    concordance = c(concordant = 0, discordant = 0, tied.x = 0, 
                    tied.y = 0, tied.xy = 0, concordance = NA, std = NA, 
                    timefix = FALSE)
    rval <- list(coefficients = ctemp, var = matrix(0, ncoef, 
                                                    ncoef), loglik = c(0, 0), score = 0, iter = 0, linear.predictors = offset, 
                 residuals = rep(0, data.n), means = colMeans(X), 
                 method = method, n = data.n, nevent = 0, terms = Terms, 
                 assign = assign, concordance = concordance, wald.test = 0, 
                 y = Y, call = Call)
    class(rval) <- "coxph"
    return(rval)
  }
  if (multi) {
    if (length(strats) > 0) {
      stratum_map <- tmap[c(1L, strats), ]
      stratum_map[-1, ] <- ifelse(stratum_map[-1, ] > 0, 
                                  1L, 0L)
      if (nrow(stratum_map) > 2) {
        temp <- stratum_map[-1, ]
        if (!all(apply(temp, 2, function(x) all(x == 
                                                0) || all(x == 1)))) {
          strata.keep <- mf[, strats]
          istrat <- sapply(strata.keep, as.numeric)
        }
      }
    }
    else stratum_map <- tmap[1, , drop = FALSE]
    cmap <- parsecovar3(tmap, colnames(X), attr(X, "assign"), 
                        covlist2$phbaseline)
    xstack <- stacker(cmap, stratum_map, as.integer(istate), 
                      X, Y, strata = istrat, states = states)
    rkeep <- unique(xstack$rindex)
    transitions <- survcheck2(Y[rkeep, ], id[rkeep], istate[rkeep])$transitions
    X <- xstack$X
    Y <- xstack$Y
    istrat <- xstack$strata
    if (length(offset)) 
      offset <- offset[xstack$rindex]
    if (length(weights)) 
      weights <- weights[xstack$rindex]
    if (length(cluster)) 
      cluster <- cluster[xstack$rindex]
    t2 <- tmap[-c(1, strats), , drop = FALSE]
    r2 <- row(t2)[!duplicated(as.vector(t2)) & t2 != 0]
    c2 <- col(t2)[!duplicated(as.vector(t2)) & t2 != 0]
    a2 <- lapply(seq(along = r2), function(i) {
      cmap[assign[[r2[i]]], c2[i]]
    })
    tab <- table(r2)
    count <- tab[r2]
    names(a2) <- ifelse(count == 1, row.names(t2)[r2], paste(row.names(t2)[r2], 
                                                             colnames(cmap)[c2], sep = "_"))
    assign <- a2
  }
  if (!all(is.finite(X))) 
    stop("data contains an infinite predictor")
  if (missing(init)) 
    init <- NULL
  else {
    if (length(init) != ncol(X)) 
      stop("wrong length for init argument")
    temp <- X %*% init - sum(colMeans(X) * init) + offset
    if (any(exp(temp) > .Machine$double.xmax) || all(exp(temp) == 
                                                     0)) 
      stop("initial values lead to overflow or underflow of the exp function")
  }
  pterms <- sapply(mf, inherits, "coxph.penalty")
  if (any(pterms)) {
    pattr <- lapply(mf[pterms], attributes)
    pname <- names(pterms)[pterms]
    ord <- attr(Terms, "order")[match(pname, attr(Terms, 
                                                  "term.labels"))]
    if (any(ord > 1)) 
      stop("Penalty terms cannot be in an interaction")
    pcols <- assign[match(pname, names(assign))]
    fit <- coxpenal.fit(X, Y, istrat, offset, init = init, 
                        control, weights = weights, method = method, row.names(mf), 
                        pcols, pattr, assign, nocenter = nocenter)
  }
  else {
    rname <- row.names(mf)
    if (multi) 
      rname <- rname[xstack$rindex]
    if (method == "breslow" || method == "efron") {
      if (grepl("right", type)) 
        fit <- coxph.fit(X, Y, istrat, offset, init, 
                         control, weights = weights, method = method, 
                         rname, nocenter = nocenter)
      else fit <- agreg.fit(X, Y, istrat, offset, init, 
                            control, weights = weights, method = method, 
                            rname, nocenter = nocenter)
    }
    else if (method == "exact") {
      if (type == "right") 
        fit <- coxexact.fit(X, Y, istrat, offset, init, 
                            control, weights = weights, method = method, 
                            rname, nocenter = nocenter)
      else fit <- agexact.fit(X, Y, istrat, offset, init, 
                              control, weights = weights, method = method, 
                              rname, nocenter = nocenter)
    }
    else stop(paste("Unknown method", method))
  }
  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "coxph"
  }
  else {
    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable", 
                   paste(vars, collapse = " "))
      if (!singular.ok) 
        stop(msg)
    }
    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)])
    fit$terms <- Terms
    fit$assign <- assign
    class(fit) <- fit$class
    fit$class <- NULL
    if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {
      fit$naive.var <- fit$var
      fit2 <- c(fit, list(x = X, y = Y, weights = weights))
      if (length(istrat)) 
        fit2$strata <- istrat
      if (length(cluster)) {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
                                collapse = cluster, weighted = TRUE)
        if (is.null(init)) 
          fit2$linear.predictors <- 0 * fit$linear.predictors
        else fit2$linear.predictors <- c(X %*% init)
        temp0 <- residuals.coxph(fit2, type = "score", 
                                 collapse = cluster, weighted = TRUE)
      }
      else {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
                                weighted = TRUE)
        fit2$linear.predictors <- 0 * fit$linear.predictors
        temp0 <- residuals.coxph(fit2, type = "score", 
                                 weighted = TRUE)
      }
      fit$var <- t(temp) %*% temp
      u <- apply(as.matrix(temp0), 2, sum)
      fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u, 
                                control$toler.chol)$test
    }
    if (length(fit$coefficients) && is.null(fit$wald.test)) {
      nabeta <- !is.na(fit$coefficients)
      if (is.null(init)) 
        temp <- fit$coefficients[nabeta]
      else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
      fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta], 
                                   temp, control$toler.chol)$test
    }
    if (length(cluster)) 
      temp <- concordancefit(Y, fit$linear.predictors, 
                             istrat, weights, cluster = cluster, reverse = TRUE, 
                             timefix = FALSE)
    else temp <- concordancefit(Y, fit$linear.predictors, 
                                istrat, weights, reverse = TRUE, timefix = FALSE)
    if (is.matrix(temp$count)) 
      fit$concordance <- c(colSums(temp$count), concordance = temp$concordance, 
                           std = sqrt(temp$var))
    else fit$concordance <- c(temp$count, concordance = temp$concordance, 
                              std = sqrt(temp$var))
    na.action <- attr(mf, "na.action")
    if (length(na.action)) 
      fit$na.action <- na.action
    if (model) {
      if (length(timetrans)) {
        stop("'model=TRUE' not supported for models with tt terms")
      }
      fit$model <- mf
    }
    if (x) {
      fit$x <- X
      if (length(timetrans)) 
        fit$strata <- istrat
      else if (length(strats)) 
        fit$strata <- strata.keep
    }
    if (y) 
      fit$y <- Y
    fit$timefix <- control$timefix
  }
  if (!is.null(weights) && any(weights != 1)) 
    fit$weights <- weights
  if (multi) {
    fit$transitions <- transitions
    fit$states <- states
    fit$cmap <- cmap
    fit$stratum_map <- stratum_map
    fit$resid <- rowsum(fit$resid, xstack$rindex)
    single <- apply(cmap, 1, function(x) all(x %in% c(0, 
                                                      max(x))))
    cindx <- col(cmap)[match(1:length(fit$coefficients), 
                             cmap)]
    rindx <- row(cmap)[match(1:length(fit$coefficients), 
                             cmap)]
    suffix <- ifelse(single[rindx], "", paste0("_", colnames(cmap)[cindx]))
    names(fit$coefficients) <- paste0(names(fit$coefficients), 
                                      suffix)
    if (x) 
      fit$strata <- istrat
    class(fit) <- c("coxphms", class(fit))
  }
  names(fit$means) <- names(fit$coefficients)
  fit$formula <- formula(Terms)
  if (length(xlevels) > 0) 
    fit$xlevels <- xlevels
  fit$contrasts <- contr.save
  if (any(offset != 0)) 
    fit$offset <- offset
  fit$call <- Call
  fit
}




###########################
fit.mult.impute =
  function (formula,
            fitter,
            xtrans,
            data,
            n.impute = xtrans$n.impute,
            fit.reps = FALSE,
            dtrans,
            derived,
            vcovOpts = NULL,
            pr = F,
            subset,
            ...)
  {
    call <- match.call()
    if (deparse(substitute(fitter)) == "lm")
      warning(
        "If you use print, summary, or anova on the result, lm methods use the\nsum of squared residuals rather than the Rubin formula for computing\nresidual variance and standard errors.  It is suggested to use ols\ninstead of lm."
      )
    using.Design <- FALSE
    fits <- if (fit.reps)
      vector("list", n.impute)
    used.mice <- any(class(xtrans) == "mids")
    if (used.mice && missing(n.impute))
      n.impute <- xtrans$m
    stats.ok2average <- c(
      "linear.predictors",
      "fitted.values",
      "stats",
      "means",
      "icoef",
      "scale",
      "center",
      "y.imputed"
    )
    for (i in 1:n.impute) {
      if (used.mice) {
        completed.data <- mice::complete(xtrans, i)
        for (impvar in names(completed.data))
          if (length(attr(completed.data[[impvar]],
                          "contrasts")))
            attr(completed.data[[impvar]], "contrasts") <- NULL
      }
      else {
        completed.data <- data
        imputed.data <- impute.transcan(
          xtrans,
          imputation = i,
          data = data,
          list.out = TRUE,
          pr = FALSE,
          check = FALSE
        )
        completed.data[names(imputed.data)] <- imputed.data
      }
      if (!missing(dtrans))
        completed.data <- dtrans(completed.data)
      if (!missing(derived)) {
        stop("derived variables in fit.mult.imputed not yet implemented")
        eval(derived, completed.data)
      }
      if (using.Design)
        options(Design.attr = da)
      f <- if (missing(subset))
        fitter(formula, data = completed.data, ...)
      else
        fitter(formula, data = completed.data, subset = subset,
               ...)
      if (fit.reps)
        fits[[i]] <- f
      cof <- f$coef
      v <- do.call("vcov", c(list(
        object = f, intercepts = "all"
      ),
      vcovOpts))
      if (i == 1) {
        if (inherits(f, "orm") && length(f$na.action) &&
            length(f$na.action$nmiss) && f$na.action$nmiss[1] >
            0)
          warning(
            "When using fit.mult.impute with orm, there should not be any missing\nY values because different imputations will result in differing numbers\nof intercepts"
          )
        assign <- f$assign
        ns <- num.intercepts(f)
        ik <- coef.intercepts <- NULL
        if (ns > 0) {
          ik <- attr(v, "intercepts")
          if (length(ik)) {
            if (ik == "all")
              ik <- 1:ns
            else if (ik == "none")
              ik <- 0
            lenik <- length(ik)
            if (length(ik) == 1 && ik == 0)
              lenik <- 0
            if (lenik != ns) {
              for (j in 1:length(assign))
                assign[[j]] <- assign[[j]] -
                  (ns - lenik)
              coef.intercepts <- ik
            }
          }
        }
      }
      if (length(ik))
        cof <- c(cof[ik], cof[-(1:ns)])
      nvar0 <- length(cof)
      nvar <- nrow(v)
      if (nvar > nvar0) {
        cof <- c(cof, log(f$scale))
        names(cof) <- c(names(f$coef), if ((nvar - nvar0) ==
                                           1)
          "Log(scale)"
          else
            names(f$scale))
      }
      if (i == 1) {
        vavg <- 0 * v
        p <- length(cof)
        bar <- rep(0, p)
        vname <- names(cof)
        cov <- matrix(
          0,
          nrow = p,
          ncol = p,
          dimnames = list(vname,
                          vname)
        )
        astats <- NULL
        fitcomp <- names(f)[names(f) %in% stats.ok2average]
        if (length(fitcomp))
          for (ncomp in fitcomp)
            astats[[ncomp]] <- f[[ncomp]]
        if (inherits(f, "Design") | inherits(f, "rms")) {
          using.Design <- TRUE
          da <- f$Design
        }
      }
      vavg <- vavg + v
      bar <- bar + cof
      cof <- as.matrix(cof)
      cov <- cov + cof %*% t(cof)
      if (i > 1 && length(fitcomp))
        for (ncomp in fitcomp)
          astats[[ncomp]] <- astats[[ncomp]] +
        f[[ncomp]]
    }
    vavg <- vavg / n.impute
    bar <- bar / n.impute
    bar <- as.matrix(bar)
    cov <- (cov - n.impute * bar %*% t(bar)) / (n.impute - 1)
    U <- diag(vavg)
    B <- diag(cov)
    cov <- vavg + (n.impute + 1) / n.impute * cov
    r <- diag(cov) / diag(vavg)
    names(r) <- vname
    tau <- (1 + 1 / n.impute) * B / U
    missingInfo <- tau / (1 + tau)
    dfmi <- (n.impute - 1) * ((1 + 1 / tau) ^ 2)
    if (length(fitcomp))
      for (ncomp in fitcomp)
        f[[ncomp]] <- astats[[ncomp]] / n.impute
    if (pr) {
      cat("\nVariance Inflation Factors Due to Imputation:\n\n")
      print(round(r, 2))
      cat("\nRate of Missing Information:\n\n")
      print(round(missingInfo, 2))
      cat("\nd.f. for t-distribution for Tests of Single Coefficients:\n\n")
      print(round(dfmi, 2))
      if (length(fitcomp)) {
        cat(
          "\nThe following fit components were averaged over the",
          n.impute,
          "model fits:\n\n"
        )
        cat(" ", fitcomp, "\n\n")
      }
    }
    f$coefficients <- drop(bar)
    if (length(coef.intercepts))
      attr(f$coefficients, "intercepts") <- coef.intercepts
    attr(cov, "intercepts") <- ik
    f$var <- cov
    f$variance.inflation.impute <- r
    f$missingInfo <- missingInfo
    f$dfmi <- dfmi
    f$fits <- fits
    f$formula <- formula
    f$assign <- assign
    f$call <- call
    if (using.Design)
      options(Design.attr = NULL)
    class(f) <- c("fit.mult.impute", class(f))
    f
  }