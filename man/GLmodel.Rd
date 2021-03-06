% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLmodel.r
\name{GLmodel}
\alias{GLmodel}
\title{Fonction to add interaction variable and to fit a group lasso model.}
\usage{
GLmodel(X, Y, XIntTrain, XIntTest, idSubs, interLength, listGenesSNP,
  nlambda = 100, limitLambda = 30, lambda.cri = "min")
}
\arguments{
\item{X}{a genotype matrix where columns are genetic markers and rows are samples.}

\item{Y}{phenotype values can be logical (1/0) or numericGLmodel.}

\item{interLength}{a vector that indicate the number of interaction variable for each gene couple.}

\item{listGenesSNP}{list containing the names of genetic marker for each group.}

\item{nlambda}{length of the grid of lambda values (100 by default).}

\item{limitLambda}{number of lambda values among the grid to consider for the cross validation (30 by default).}

\item{lambda.cri}{criteria for lambda selection are "min" or "oneSE" ("min" by default).}

\item{XBet}{a matrix where columns are interaction variables and rows are samples.}
}
\value{
Returns a list including:
\item{res_GL.min}{a matrix nx1 where row are single marker or interaction variables and column the coefficient values obtain with the lambda criteria chosen.}
\item{pval.adj}{a matrix where row are single marker or interaction variables for which res_GL.min is no equal to zero and column the pvalues obtained after the cleaning procedure and adjusted with Benjamini & Hochberg correction.}
\item{vc}{a list with result of the cross validation.}
}
\description{
The GLmodel fonction first add interaction variable to a genetic data set structured by group and then fit a group lasso model with an adaptive ridge cleaning approach.
}

