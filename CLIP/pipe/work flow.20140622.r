
# R version 3.0.0 (2013-04-03) -- "Masked Marvel"

install.packages('VGAM')

setwd('E:/Yyx/13学年第二学期/CBI/魏丽萍/20140621_ref_From_HouMei/ZTNB')

### 主要目的就是对len和read进行回归

suppressMessages(require(MASS));
suppressMessages(require(VGAM));

args <- c("2.filter.rehead.merge", "0.05")
# args = commandArgs(TRUE);
data = read.table(args[1], sep = "\t");
len = data[,3] - data[,2];
read = data[,5];
cut = as.numeric(as.character(args[2]))
# tabularize the data
wts = table(read,len);
temp.coln = as.integer(names(wts[1,]));
temp.rown = as.integer(names(wts[,1]));

drow = length(temp.rown);
dcol = length(temp.coln);

rown = matrix(rep(temp.rown,dcol),drow,dcol);
coln = t(matrix(rep(temp.coln,drow),dcol,drow));

wdx = which(wts>0);
tread = rown[wdx];
tlen = coln[wdx];
twts = wts[wdx];

# local polynomial regression: read ~ len
# smoothing parameter: 
alp = 0.95
# polynomial degree:
pd = 1
# 
lregfit = loess(tread ~ tlen, weights = twts, span = alp, family ="symmetric", degree = pd);
mu = lregfit$fitted;
mu[mu<=0] = min(exp(-4),mu[mu>0]); # bounded mean function
logmu = log(mu);

length(tlen)   # 2221
plot(tread ~ tlen, pch=16, cex=0.5, col='gray')
xx <- seq(from=min(tlen),to=max(tlen),by=1)
lines(xx, predict(lregfit, newdata=xx), col='red', lwd=2)
plot(tread ~ tlen, pch=16, cex=0.5, col='gray', ylim=c(0,5000))
lines(xx, predict(lregfit, newdata=xx), col='red', lwd=2)

# compute p-values using N(0,1) for (tread[-cdx] - mu[-cdx]) / sqrt(mu[-cdx])
nb.pvalue=numeric(length(tread));
zscore = (tread-mu)/sqrt(mu); 
cdx = which(zscore < 4);
nb.pvalue[-cdx] = pnorm(zscore[-cdx],lower.tail=FALSE);

my_col_vec <- rep('gray', length(tlen))
my_col_vec[-cdx] <- 'orange'
plot(tread ~ tlen, pch=16, cex=0.5, col=my_col_vec)
lines(xx, predict(lregfit, newdata=xx), col='red', lwd=2)
plot(tread ~ tlen, pch=16, cex=0.5, col=my_col_vec, ylim=c(0,5000))
lines(xx, predict(lregfit, newdata=xx), col='red', lwd=2)
plot(tread ~ tlen, pch=16, cex=0.5, col=my_col_vec, ylim=c(0,1000))
lines(xx, predict(lregfit, newdata=xx), col='red', lwd=2)
range(nb.pvalue[-cdx])   # 0.000000e+00 3.067879e-05

tread_fit = tread[cdx];
tlen_fit = tlen[cdx];
twts_fit = twts[cdx];

lregfit = loess(tread_fit ~ tlen_fit, weights = twts_fit, span = alp, family ="symmetric", degree = pd);
mu = lregfit$fitted;
mu[mu<=0] = min(exp(-4),mu[mu>0]); # bounded mean function
logmu = log(mu);

length(tlen_fit)   # 630
plot(tread_fit ~ tlen_fit, pch=16, cex=0.5, col='gray')
lines(xx, predict(lregfit, newdata=xx), col='red', lwd=2)

# negative binomail regression with the known predictor log(mu)
nb = vglm(tread_fit ~ 1, posnegbinomial(), weights = twts_fit, 
          maxit = 200, trace = FALSE, offset = logmu,epsilon=10^(-3))
# model parameters
(khat <- Coef(nb)["size"]);
(muhat <- Coef(nb)["munb"]);
nbsize = khat * mu;
nbmu = muhat * mu;

# pvalues
nb.pvalue[cdx] = pnbinom(tread_fit - 1, size = nbsize, mu = nbmu, lower.tail = FALSE) / (1-dnbinom(0, size = nbsize, mu = nbmu))
nb.fdr = p.adjust(nb.pvalue,method='BH')

hist(nb.fdr, breaks=20)

# output
nbdx = nb.fdr <= cut;
#sum(twts[nbdx]) # the number of clusters whose FDR <= cut
# read counts and cluster length of clusters whose FDR <= cut
nb.out = as.data.frame(matrix(c(tread[nbdx],tlen[nbdx],nb.pvalue[nbdx],nb.fdr[nbdx]),sum(nbdx),4))
names(nb.out) = c("read","length","p","fdr")
outname = paste(args[1],".ztnb",sep="")
write.table(nb.out,outname,sep="\t",quote=F,row.names=F);

nb.all = as.data.frame(matrix(c(tread,tlen,nb.pvalue,nb.fdr), ncol=4))
names(nb.all) = c("read","length","p","fdr")
outname = paste(args[1],".ztnb_all",sep="")
write.table(nb.all,outname,sep="\t",quote=F,row.names=F);



stop('end of the program')

> vglm
function (formula, family, data = list(), weights = NULL, subset = NULL, 
    na.action = na.fail, etastart = NULL, mustart = NULL, coefstart = NULL, 
    control = vglm.control(...), offset = NULL, method = "vglm.fit", 
    model = FALSE, x.arg = TRUE, y.arg = TRUE, contrasts = NULL, 
    constraints = NULL, extra = list(), form2 = NULL, qr.arg = TRUE, 
    smart = TRUE, ...) 
{
    dataname <- as.character(substitute(data))
    function.name <- "vglm"
    ocall <- match.call()
    if (smart) 
        setup.smart("write")
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), vglm.fit = 1, stop("invalid 'method': ", 
        method))
    mt <- attr(mf, "terms")
    xlev <- .getXlevels(mt, mf)
    y <- model.response(mf, "any")
    x <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(y), 0)
    attr(x, "assign") <- attrassigndefault(x, mt)
    if (!is.null(form2)) {
        if (!is.null(subset)) 
            stop("argument 'subset' cannot be used when ", "argument 'form2' is used")
        retlist <- shadowvglm(formula = form2, family = family, 
            data = data, na.action = na.action, control = vglm.control(...), 
            method = method, model = model, x.arg = x.arg, y.arg = y.arg, 
            contrasts = contrasts, constraints = constraints, 
            extra = extra, qr.arg = qr.arg)
        Ym2 <- retlist$Ym2
        Xm2 <- retlist$Xm2
        if (length(Ym2)) {
            if (nrow(as.matrix(Ym2)) != nrow(as.matrix(y))) 
                stop("number of rows of 'y' and 'Ym2' are unequal")
        }
        if (length(Xm2)) {
            if (nrow(as.matrix(Xm2)) != nrow(as.matrix(x))) 
                stop("number of rows of 'x' and 'Xm2' are unequal")
        }
    }
    else {
        Xm2 <- Ym2 <- NULL
    }
    offset <- model.offset(mf)
    if (is.null(offset)) 
        offset <- 0
    w <- model.weights(mf)
    if (!length(w)) {
        w <- rep(1, nrow(mf))
    }
    else if (ncol(as.matrix(w)) == 1 && any(w < 0)) 
        stop("negative weights not allowed")
    if (is.character(family)) 
        family <- get(family)
    if (is.function(family)) 
        family <- family()
    if (!inherits(family, "vglmff")) {
        stop("'family = ", family, "' is not a VGAM family function")
    }
    eval(vcontrol.expression)
    if (length(slot(family, "first"))) 
        eval(slot(family, "first"))
    vglm.fitter <- get(method)
    fit <- vglm.fitter(x = x, y = y, w = w, offset = offset, 
        Xm2 = Xm2, Ym2 = Ym2, etastart = etastart, mustart = mustart, 
        coefstart = coefstart, family = family, control = control, 
        constraints = constraints, criterion = control$criterion, 
        extra = extra, qr.arg = qr.arg, Terms = mt, function.name = function.name, 
        ...)
    fit$misc$dataname <- dataname
    if (smart) {
        fit$smart.prediction <- get.smart.prediction()
        wrapup.smart()
    }
    answer <- new(Class = "vglm", assign = attr(x, "assign"), 
        call = ocall, coefficients = fit$coefficients, constraints = fit$constraints, 
        criterion = fit$crit.list, df.residual = fit$df.residual, 
        df.total = fit$df.total, dispersion = 1, effects = fit$effects, 
        family = fit$family, misc = fit$misc, model = if (model) 
            mf
        else data.frame(), R = fit$R, rank = fit$rank, residuals = as.matrix(fit$residuals), 
        res.ss = fit$res.ss, smart.prediction = as.list(fit$smart.prediction), 
        terms = list(terms = mt))
    if (!smart) 
        answer@smart.prediction <- list(smart.arg = FALSE)
    if (qr.arg) {
        class(fit$qr) <- "list"
        slot(answer, "qr") <- fit$qr
    }
    if (length(attr(x, "contrasts"))) 
        slot(answer, "contrasts") <- attr(x, "contrasts")
    if (length(fit$fitted.values)) 
        slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
    slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action"))) 
        list(aaa)
    else list()
    if (length(offset)) 
        slot(answer, "offset") <- as.matrix(offset)
    if (length(fit$weights)) 
        slot(answer, "weights") <- as.matrix(fit$weights)
    if (x.arg) 
        slot(answer, "x") <- fit$x
    if (x.arg && length(Xm2)) 
        slot(answer, "Xm2") <- Xm2
    if (y.arg && length(Ym2)) 
        slot(answer, "Ym2") <- as.matrix(Ym2)
    if (!is.null(form2)) 
        slot(answer, "callXm2") <- retlist$call
    answer@misc$formula <- formula
    answer@misc$form2 <- form2
    if (length(xlev)) 
        slot(answer, "xlevels") <- xlev
    if (y.arg) 
        slot(answer, "y") <- as.matrix(fit$y)
    slot(answer, "control") <- fit$control
    slot(answer, "extra") <- if (length(fit$extra)) {
        if (is.list(fit$extra)) 
            fit$extra
        else {
            warning("'extra' is not a list, therefore placing ", 
                "'extra' into a list")
            list(fit$extra)
        }
    }
    else list()
    slot(answer, "iter") <- fit$iter
    slot(answer, "post") <- fit$post
    fit$predictors <- as.matrix(fit$predictors)
    if (length(fit$misc$predictors.names) == ncol(fit$predictors)) 
        dimnames(fit$predictors) <- list(dimnames(fit$predictors)[[1]], 
            fit$misc$predictors.names)
    slot(answer, "predictors") <- fit$predictors
    if (length(fit$prior.weights)) 
        slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)
    answer
}
<environment: namespace:VGAM>
attr(,"smart")
[1] TRUE


> vglm.fit
function (x, y, w = rep(1, length(x[, 1])), X.vlm.arg = NULL, 
    Xm2 = NULL, Ym2 = NULL, etastart = NULL, mustart = NULL, 
    coefstart = NULL, offset = 0, family, control = vglm.control(), 
    criterion = "coefficients", qr.arg = FALSE, constraints = NULL, 
    extra = NULL, Terms = Terms, function.name = "vglm", ...) 
{
    eff.n <- nrow(x)
    specialCM <- NULL
    post <- list()
    check.rank <- TRUE
    check.rank <- control$Check.rank
    nonparametric <- FALSE
    epsilon <- control$epsilon
    maxit <- control$maxit
    save.weight <- control$save.weight
    trace <- control$trace
    orig.stepsize <- control$stepsize
    minimize.criterion <- control$min.criterion
    fv <- NULL
    n <- dim(x)[1]
    copy.X.vlm <- FALSE
    stepsize <- orig.stepsize
    old.coeffs <- coefstart
    intercept.only <- ncol(x) == 1 && dimnames(x)[[2]] == "(Intercept)"
    y.names <- predictors.names <- NULL
    n.save <- n
    if (length(slot(family, "initialize"))) 
        eval(slot(family, "initialize"))
    if (length(etastart)) {
        eta <- etastart
        mu <- if (length(mustart)) 
            mustart
        else if (length(body(slot(family, "linkinv")))) 
            slot(family, "linkinv")(eta, extra)
        else warning("argument 'etastart' assigned a value ", 
            "but there is no 'linkinv' slot to use it")
    }
    if (length(mustart)) {
        mu <- mustart
        if (length(body(slot(family, "linkfun")))) {
            eta <- slot(family, "linkfun")(mu, extra)
        }
        else {
            warning("argument 'mustart' assigned a value ", "but there is no 'link' slot to use it")
        }
    }
    M <- if (is.matrix(eta)) 
        ncol(eta)
    else 1
    if (length(slot(family, "constraints"))) 
        eval(slot(family, "constraints"))
    Hlist <- process.constraints(constraints, x, M, specialCM = specialCM)
    ncolHlist <- unlist(lapply(Hlist, ncol))
    dimB <- sum(ncolHlist)
    X.vlm.save <- if (length(X.vlm.arg)) {
        X.vlm.arg
    }
    else {
        lm2vlm.model.matrix(x, Hlist, xij = control$xij, Xm2 = Xm2)
    }
    if (length(coefstart)) {
        eta <- if (ncol(X.vlm.save) > 1) {
            matrix(X.vlm.save %*% coefstart, n, M, byrow = TRUE) + 
                offset
        }
        else {
            matrix(X.vlm.save * coefstart, n, M, byrow = TRUE) + 
                offset
        }
        if (M == 1) 
            eta <- c(eta)
        mu <- slot(family, "linkinv")(eta, extra)
    }
    if (criterion != "coefficients") {
        tfun <- slot(family, criterion)
    }
    iter <- 1
    new.crit <- switch(criterion, coefficients = 1, tfun(mu = mu, 
        y = y, w = w, res = FALSE, eta = eta, extra))
    old.crit <- ifelse(minimize.criterion, 10 * new.crit + 10, 
        -10 * new.crit - 10)
    deriv.mu <- eval(slot(family, "deriv"))
    wz <- eval(slot(family, "weight"))
    if (control$checkwz) 
        wz <- checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon)
    U <- vchol(wz, M = M, n = n, silent = !trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
    z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset
    c.list <- list(z = as.double(z), fit = as.double(t(eta)), 
        one.more = TRUE, coeff = as.double(rep(1, ncol(X.vlm.save))), 
        U = as.double(U), copy.X.vlm = copy.X.vlm, X.vlm = if (copy.X.vlm) as.double(X.vlm.save) else double(3))
    dX.vlm <- as.integer(dim(X.vlm.save))
    nrow.X.vlm <- dX.vlm[[1]]
    ncol.X.vlm <- dX.vlm[[2]]
    if (nrow.X.vlm < ncol.X.vlm) 
        stop(ncol.X.vlm, " parameters but only ", nrow.X.vlm, 
            " observations")
    while (c.list$one.more) {
        tfit <- vlm.wfit(xmat = X.vlm.save, z, Hlist = NULL, 
            U = U, matrix.out = FALSE, is.vlmX = TRUE, qr = qr.arg, 
            xij = NULL)
        c.list$coeff <- tfit$coefficients
        tfit$predictors <- tfit$fitted.values
        c.list$fit <- tfit$fitted.values
        if (!c.list$one.more) {
            break
        }
        fv <- c.list$fit
        new.coeffs <- c.list$coeff
        if (length(slot(family, "middle"))) 
            eval(slot(family, "middle"))
        eta <- fv + offset
        mu <- slot(family, "linkinv")(eta, extra)
        if (length(slot(family, "middle2"))) 
            eval(slot(family, "middle2"))
        old.crit <- new.crit
        new.crit <- switch(criterion, coefficients = new.coeffs, 
            tfun(mu = mu, y = y, w = w, res = FALSE, eta = eta, 
                extra))
        if (trace && orig.stepsize == 1) {
            cat("VGLM    linear loop ", iter, ": ", criterion, 
                "= ")
            UUUU <- switch(criterion, coefficients = format(new.crit, 
                dig = round(1 - log10(epsilon))), format(new.crit, 
                dig = max(4, round(-0 - log10(epsilon) + log10(sqrt(eff.n))))))
            switch(criterion, coefficients = {
                if (length(new.crit) > 2) cat("\n")
                cat(UUUU, fill = TRUE, sep = ", ")
            }, cat(UUUU, fill = TRUE, sep = ", "))
        }
        take.half.step <- (control$half.stepsizing && length(old.coeffs)) && 
            ((orig.stepsize != 1) || (criterion != "coefficients" && 
                (if (minimize.criterion) 
                  new.crit > old.crit
                else new.crit < old.crit)))
        if (!is.logical(take.half.step)) 
            take.half.step <- TRUE
        if (take.half.step) {
            stepsize <- 2 * min(orig.stepsize, 2 * stepsize)
            new.coeffs.save <- new.coeffs
            if (trace) 
                cat("Taking a modified step")
            repeat {
                if (trace) {
                  cat(".")
                  flush.console()
                }
                stepsize <- stepsize/2
                if (too.small <- stepsize < 0.001) 
                  break
                new.coeffs <- (1 - stepsize) * old.coeffs + stepsize * 
                  new.coeffs.save
                if (length(slot(family, "middle"))) 
                  eval(slot(family, "middle"))
                fv <- X.vlm.save %*% new.coeffs
                if (M > 1) 
                  fv <- matrix(fv, n, M, byrow = TRUE)
                eta <- fv + offset
                mu <- slot(family, "linkinv")(eta, extra)
                if (length(slot(family, "middle2"))) 
                  eval(slot(family, "middle2"))
                new.crit <- switch(criterion, coefficients = new.coeffs, 
                  tfun(mu = mu, y = y, w = w, res = FALSE, eta = eta, 
                    extra))
                if ((criterion == "coefficients") || (minimize.criterion && 
                  new.crit < old.crit) || (!minimize.criterion && 
                  new.crit > old.crit)) 
                  break
            }
            if (trace) 
                cat("\n")
            if (too.small) {
                warning("iterations terminated because ", "half-step sizes are very small")
                one.more <- FALSE
            }
            else {
                if (trace) {
                  cat("VGLM    linear loop ", iter, ": ", criterion, 
                    "= ")
                  UUUU <- switch(criterion, coefficients = format(new.crit, 
                    dig = round(1 - log10(epsilon))), format(new.crit, 
                    dig = max(4, round(-0 - log10(epsilon) + 
                      log10(sqrt(eff.n))))))
                  switch(criterion, coefficients = {
                    if (length(new.crit) > 2) cat("\n")
                    cat(UUUU, fill = TRUE, sep = ", ")
                  }, cat(UUUU, fill = TRUE, sep = ", "))
                }
                one.more <- eval(control$convergence)
            }
        }
        else {
            one.more <- eval(control$convergence)
        }
        flush.console()
        if (!is.logical(one.more)) 
            one.more <- FALSE
        if (one.more) {
            iter <- iter + 1
            deriv.mu <- eval(slot(family, "deriv"))
            wz <- eval(slot(family, "weight"))
            if (control$checkwz) 
                wz <- checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon)
            U <- vchol(wz, M = M, n = n, silent = !trace)
            tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
            z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset
            c.list$z <- z
            c.list$U <- U
            if (copy.X.vlm) 
                c.list$X.vlm <- X.vlm.save
        }
        c.list$one.more <- one.more
        c.list$coeff <- runif(length(new.coeffs))
        old.coeffs <- new.coeffs
    }
    if (maxit > 1 && iter >= maxit && !control$noWarning) 
        warning("convergence not obtained in ", maxit, " iterations")
    dnrow.X.vlm <- labels(X.vlm.save)
    xnrow.X.vlm <- dnrow.X.vlm[[2]]
    ynrow.X.vlm <- dnrow.X.vlm[[1]]
    if (length(slot(family, "fini"))) 
        eval(slot(family, "fini"))
    if (M > 1) 
        tfit$predictors <- matrix(tfit$predictors, n, M)
    coefs <- tfit$coefficients
    asgn <- attr(X.vlm.save, "assign")
    names(coefs) <- xnrow.X.vlm
    rank <- tfit$rank
    cnames <- xnrow.X.vlm
    if (check.rank && rank < ncol.X.vlm) 
        stop("vglm only handles full-rank models (currently)")
    R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    attributes(R) <- list(dim = c(ncol.X.vlm, ncol.X.vlm), dimnames = list(cnames, 
        cnames), rank = rank)
    effects <- tfit$effects
    neff <- rep("", nrow.X.vlm)
    neff[seq(ncol.X.vlm)] <- cnames
    names(effects) <- neff
    dim(tfit$predictors) <- c(n, M)
    dn <- labels(x)
    yn <- dn[[1]]
    xn <- dn[[2]]
    residuals <- z - tfit$predictors
    if (M == 1) {
        tfit$predictors <- as.vector(tfit$predictors)
        residuals <- as.vector(residuals)
        names(residuals) <- names(tfit$predictors) <- yn
    }
    else {
        dimnames(residuals) <- dimnames(tfit$predictors) <- list(yn, 
            predictors.names)
    }
    if (is.matrix(mu)) {
        if (length(dimnames(y)[[2]])) {
            y.names <- dimnames(y)[[2]]
        }
        if (length(dimnames(mu)[[2]])) {
            y.names <- dimnames(mu)[[2]]
        }
        dimnames(mu) <- list(yn, y.names)
    }
    else {
        names(mu) <- names(fv)
    }
    df.residual <- nrow.X.vlm - rank
    fit <- list(assign = asgn, coefficients = coefs, constraints = Hlist, 
        df.residual = df.residual, df.total = n * M, effects = effects, 
        fitted.values = mu, offset = offset, rank = rank, residuals = residuals, 
        R = R, terms = Terms)
    if (qr.arg) {
        fit$qr <- tfit$qr
        dimnames(fit$qr$qr) <- dnrow.X.vlm
    }
    if (M == 1) {
        wz <- as.vector(wz)
    }
    fit$weights <- if (save.weight) 
        wz
    else NULL
    misc <- list(colnames.x = xn, colnames.X.vlm = xnrow.X.vlm, 
        criterion = criterion, function.name = function.name, 
        intercept.only = intercept.only, predictors.names = predictors.names, 
        M = M, n = n, nonparametric = nonparametric, nrow.X.vlm = nrow.X.vlm, 
        orig.assign = attr(x, "assign"), p = ncol(x), ncol.X.vlm = ncol.X.vlm, 
        ynames = dimnames(y)[[2]])
    crit.list <- list()
    if (criterion != "coefficients") 
        crit.list[[criterion]] <- fit[[criterion]] <- new.crit
    for (ii in names(.min.criterion.VGAM)) {
        if (ii != criterion && any(slotNames(family) == ii) && 
            length(body(slot(family, ii)))) {
            fit[[ii]] <- crit.list[[ii]] <- (slot(family, ii))(mu = mu, 
                y = y, w = w, res = FALSE, eta = eta, extra)
        }
    }
    if (w[1] != 1 || any(w != w[1])) 
        fit$prior.weights <- w
    if (length(slot(family, "last"))) 
        eval(slot(family, "last"))
    structure(c(fit, list(predictors = tfit$predictors, contrasts = attr(x, 
        "contrasts"), control = control, crit.list = crit.list, 
        extra = extra, family = family, iter = iter, misc = misc, 
        post = post, res.ss = tfit$res.ss, x = x, y = y)), vclass = slot(family, 
        "vfamily"))
}
<environment: namespace:VGAM>
