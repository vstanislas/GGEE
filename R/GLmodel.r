
#' Fonction to create interaction variable. 
#'
#' The fonction create interaction variable from a genetic data set structured by group. 
#'
#' @param X a genotype matrix where columns are genetic markers and rows are samples. 
#' @param Y phenotype values can be logical (1/0) or numericGLmodel.
#' @param method method to consider for interactions variable  construction :
#' \itemize{
#' \item{"GGEE":}{ Compute the gene-gene interaction variables for each gene pair by finding its Eigen-epistasis Component defined as the linear combination of Gene SNPs having the highest correlation with the phenotype.} 
#' \item{"PCA":}{ Compute the gene-gene interaction variables for each gene pair with Principal Component Analysis.}
#' \item{"PLS":}{ Compute the gene-gene interaction variables for each gene pair with Partial Least Square.}
#' \item{"CCA":}{ Compute the gene-gene interaction variables for a given gene pair with Canonical Correlation Analysis.}
#' }
#' @param listGenesSNP list containing the names of genetic marker for each group.
#' @param nbcomp number of components to consider to built interaction variables for "PCA", "PLS" and "CCA" methods (2 by default).
# @param seuilVarPCA an optional threshold of percentage of variance for the "PCA" method, if noticed the interaction variables are built with the PCA component that explain at least this threshold.
#' @return Returns a list including:
#' \item{XBet}{a matrix where columns are interaction variables and rows are samples.}
#' \item{interLength}{a vector that indicate the number of interaction variable for each gene couple.}
#' @export
BuiltEpiVar <- function(X, Y, method, listGenesSNP, nbcomp =2) {
 
 	genesLength <- sapply(listGenesSNP, function(x) length(x))
    genes <- names(listGenesSNP)
    groupsGenes <- unlist(lapply(names(genesLength), function(x) rep(x, genesLength[x])))
    
    
    if (method == "CCA") {
        bet <- between.mat(X = X, G = groupsGenes, nbcomp = nbcomp)
        XBet <- as.matrix(bet$XBet)
        interLength <- bet$interLength
    } else if (method == "PCA") {
        PCA.L <- PCAGenes(X = X, listGenesSNP = listGenesSNP, nbcomp = nbcomp)
        XBet <- PCA.L$I
        interLength <- PCA.L$interLength
    } else if (method == "PLS") {
        PLS.L <- PLSGenes(X = X, Y, listGenesSNP = listGenesSNP, nbcomp = nbcomp)
        XBet <- PLS.L$Int
        interLength <- PLS.L$interLength
    } else if (method == "GGEE") {
        GGEE.L <- GGEE(X=X, Y, listGenesSNP=listGenesSNP)
        XBet <- GGEE.L$Int
        interLength <- GGEE.L$interLength
    } else {
        warning("Interaction variables construction method unknown")
    }
    
  list(XBet = XBet, interLength=interLength)

}


#' Fonction to add interaction variable and to fit a group lasso model.
#'
#' The GLmodel fonction first add interaction variable to a genetic data set structured by group and then fit a group lasso model with an adaptive ridge cleaning approach.
#'
#' @param X a genotype matrix where columns are genetic markers and rows are samples. 
#' @param Y phenotype values can be logical (1/0) or numericGLmodel.
#' @param XBet a matrix where columns are interaction variables and rows are samples.
#' @param interLength a vector that indicate the number of interaction variable for each gene couple.
#' @param listGenesSNP list containing the names of genetic marker for each group.
#' @param nlambda length of the grid of lambda values (100 by default).
#' @param limitLambda number of lambda values among the grid to consider for the cross validation (30 by default).
#' @param lambda.cri criteria for lambda selection are "min" or "oneSE" ("min" by default).
#' @return Returns a list including:
#' \item{res_GL.min}{a matrix nx1 where row are single marker or interaction variables and column the coefficient values obtain with the lambda criteria chosen.}
#' \item{pval.adj}{a matrix where row are single marker or interaction variables for which res_GL.min is no equal to zero and column the pvalues obtained after the cleaning procedure and adjusted with Benjamini & Hochberg correction.}
#' \item{vc}{a list with result of the cross validation.}
#' @export
GLmodel <- function(X, Y, XBet, interLength, listGenesSNP, nlambda=100, limitLambda=30, lambda.cri="min") {
    
   genesLength <- sapply(listGenesSNP, function(x) length(x))
   genes <- names(listGenesSNP)
    
           
    XwBet <- cbind(scale(X, center = TRUE, scale = TRUE), scale(XBet, center = TRUE, 
        scale = TRUE))
    nbInt <- length(genes) * (length(genes) - 1)/2
    groups.wBet <- c(rep(seq(genes), sapply(genes, function(x) genesLength[x])), 
        rep((length(genes) + 1):(nbInt + length(genes)), interLength))
    namesgroups.wBet <- c(rep(genes, sapply(genes, function(x) genesLength[x])), 
        rep(names(interLength), interLength))
    
    
    subs <- split(1:dim(XwBet)[1], 1:2)
    vc <- valid.crois(X = XwBet[subs[[1]], ], y = Y[subs[[1]]], groups = groups.wBet, 
        K = 10, nlambda = nlambda, limitLambda = limitLambda)
    
    if (lambda.cri =="min"){
    lambda <- vc$lambda.min
    }else if (lambda.cri =="oneSE")
    {
    lambda <- vc$lambda.oneSE
    }
    
    if (all(Y == 1 | Y == 0)) {
        coefs.min <- grplassoLog(X = XwBet[subs[[1]], ], y = Y[subs[[1]]], lambda = lambda, 
            groups = groups.wBet)
    } else {
        coefs.min <- grplassoLin(X = XwBet[subs[[1]], ], y = Y[subs[[1]]], lambda = lambda, 
            groups = groups.wBet)
    }
    
    
     res_GL.min <- as.matrix(coefs.min$coefs)
     rownames(res_GL.min) <- namesgroups.wBet
    colnames(res_GL.min) <- "Coefs"
    
    if (sum(res_GL.min) == 0) {
        clean.grp <- NULL
        pval <- c()
    } else {
        clean.grp <- ridgeAdap::cleaning.ridge(XwBet[subs[[2]], ], Y[subs[[2]]], lambda1 = lambda, 
            beta = coefs.min$coefs, center = FALSE, group = groups.wBet, penalty = "grplasso")
        pval <- clean.grp$pval
    }
    
    namesGroups <- c(genes, names(interLength))
    pval.adj <- as.matrix(p.adjust(pval, "BH"))
    dimnames(pval.adj) <- list(namesGroups[as.numeric(names(pval))], "pval.adj")
    
    
    results <- list(res_GL.min = res_GL.min, pval.adj = pval.adj, vc = vc)
    class(results)<-"GLmodel"
	return(results)
}


.seqLog <- function(from, to, length) {
  10^(seq(log10(from), log10(to), length=length))
}



valid.crois <- function(X, y, groups, K, nlambda, center = TRUE, limitLambda = nlambda) {
    
    if (all(y == 1 | y == 0)) {
        print("LogReg()")
        theta.fit <- function(x, y, lambda, groups) {
            grplasso::grplasso(x, y, index = groups, lambda = lambda, model = grplasso::LogReg(), control = grplasso::grpl.control(update.hess = "lambda", 
                trace = 0), center = center, standardize = FALSE)
        }
        theta.predict <- function(fit, x) {
            res <- x %*% as.numeric(fit$coef)
            exp(res)/(1 + exp(res))
        }
        groups <- c(NA, (groups + 1))
        X <- cbind(1, X)
        lambdamax <- grplasso::lambdamax(X, y, index = groups, penscale = sqrt, model = grplasso::LogReg(), 
            center = center, standardize = FALSE)
        lambda <- .seqLog(lambdamax, lambdamax * 1e-10, length = nlambda)[1:limitLambda]
    } else {
        cat("\nLinReg()")
        theta.fit <- function(x, y, lambda, groups) {
            grplasso::grplasso(x, y, index = groups, lambda = lambda, model = grplasso::LinReg(), control = grplasso::grpl.control(update.hess = "lambda", 
                trace = 0), center = center, standardize = FALSE)
        }
        theta.predict <- function(fit, x) {
            x %*% as.numeric(fit$coef)
        }
        groups <- c(NA, (groups + 1))
        X <- cbind(1, X)
        lambdamax <- grplasso::lambdamax(X, y, index = groups, penscale = sqrt, model = grplasso::LinReg(), 
            center = center, standardize = FALSE)
        lambda <- .seqLog(lambdamax, lambdamax * 1e-10, length = nlambda)[1:limitLambda]
    }
    
    cv.error <- c()
    cat("\nLambda =")
    for (l in lambda) {
        cat("", l)
        r <- bootstrap::crossval(X, y, theta.fit, theta.predict, lambda = l, groups = groups, 
            ngroup = K)
        err <- (y - r$cv.fit)^2
        sd.err <- sd(err)/sqrt(length(err))
        cv.error <- rbind(cv.error, c(mean(err) - sd.err, mean(err), mean(err) + 
            sd.err))
    }
    
    id.lambda.min <- which.min(cv.error[, 2])
    lambda.min <- lambda[id.lambda.min]
    lambda.oneSE <- lambda[min(which(cv.error[, 2] <= cv.error[id.lambda.min, 3]))]
    
    return(list(cv.error = cv.error, lambda = lambda, lambda.min = lambda.min, lambda.oneSE = lambda.oneSE, 
        id.lambda.min = id.lambda.min))
}






grplassoLog <- function(X, y, lambda, groups, ...) {
    ## Take care of the intercept
    groups <- c(NA, (groups + 1))
    X <- cbind(1, X)
    res <- grplasso::grplasso(X, y, index = groups, lambda = lambda, model = grplasso::LogReg(), control = grplasso::grpl.control(update.hess = "lambda", 
        trace = 0), standardize = FALSE)
    coefs <- as.matrix(res$coefficients[-1, ])
    groups <- groups[-1]  ## delete the intercept !
    groups <- groups - 1
    list(coefs = coefs, groups = groups)
}




grplassoLin <- function(X, y, lambda, groups, ...) {
    ## Take care of the intercept
    groups <- c(NA, (groups + 1))
    X <- cbind(1, X)
    res <- grplasso::grplasso(X, y, index = groups, lambda = lambda, model = grplasso::LinReg(), control = grplasso::grpl.control(update.hess = "lambda", 
        trace = 0), standardize = FALSE)
    coefs <- as.matrix(res$coefficients[-1, ])
    groups <- groups[-1]  ## delete the intercept !
    groups <- groups - 1
    list(coefs = coefs, groups = groups)
}


#' Plot GLmodel cross-validation results
#'
#' Plot results from the GLmodel cross-validation with indication of minimal and oneSE lamba values 
#'
#' @param x an object inheriting from class "GLmodel" that contains cross-validation results.
#' @param zoom if zoom=TRUE, only a reduced number of lambda values to plot. The number depending of dotNb.
#' @param dotNb number of lambda values to plot if zoom=TRUE (10 by default)
#' @export
plotGLmodel <- function(x, zoom = FALSE, dotNb = 10) {
	is.GLmodel<-function(x){if (class(x)=="GLmodel") TRUE else FALSE}
	if (!is.GLmodel(x)) stop("Not a GLmodel object")
	vc <- x$vc

    i1 <- which.min(vc$cv.error[, 2])
    i2 <- min(which(vc$cv.error[, 2] <= vc$cv.error[i1, 3]))
    if (zoom == FALSE) {
        xlim = range(vc$lambda)
        ylim = range(c(vc$cv.error[, 1], vc$cv.error[, 3]))
    } else {
        limitLambda = length(vc$lambda)
        xlim = c(0, vc$lambda[limitLambda - (dotNb-1)])
        ylim = range(c(vc$cv.error[c((limitLambda - (dotNb-1)):limitLambda), 1], vc$cv.error[c((limitLambda - 
            (dotNb-1)):limitLambda), 3]))
    }
    plot(vc$lambda, vc$cv.error[, 2], type = "l", ylim = ylim, xlab = expression(lambda), 
        xlim = xlim, ylab = "CV error")
    points(vc$lambda, vc$cv.error[, 2], pch = 20, cex = 0.7)
    lines(vc$lambda, vc$cv.error[, 1], lty = 3)
    lines(vc$lambda, vc$cv.error[, 3], lty = 3)
    abline(h = vc$cv.error[i1, 3], col = "green")
    abline(v = vc$lambda[i1], col = "red", lty = 2)
    abline(v = vc$lambda[i2], col = "blue", lty = 2)
    exprs <- list(bquote(hat(lambda)[min] == .(round(vc$lambda[i1], 2))), bquote(hat(lambda)[oneSE] == 
        .(round(vc$lambda[i2], 2))))
    exprs <- sapply(exprs, as.expression)
    legend("bottomright", legend = exprs, col = c("red", "blue"), lty = 2, cex = 1.2)

}
 