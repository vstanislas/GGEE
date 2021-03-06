% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{simGeno}
\alias{simGeno}
\title{Fonction to simulate genotype data}
\usage{
simGeno(N = 600, corr = 0.8, sizeGenesMain, sizeGenesPair, sizeGenesRemain,
  SameMainPair, MAFcausalSNP = NULL, causalSNPnb, causalSNPportion,
  AllMainPair = TRUE, minMAF = 0.05, maxMAF = 0.5)
}
\arguments{
\item{N}{integer; number of individuals.}

\item{corr}{numeric value between 0 and 1. Correlation level among the SNPs within each gene.}

\item{sizeGenesMain}{a vector of integers specifying the sizes of the genes having main effects.}

\item{sizeGenesPair}{a vector of integers specifying the sizes of the genes having interaction effects. The size of the vector has to be an even number, pairs being defined with the gene taken two by two along the vector.}

\item{sizeGenesRemain}{a vector of integers specifying the sizes of the genes having none effects.}

\item{SameMainPair}{logical. If TRUE, genes with interation effects will also have main effects}

\item{MAFcausalSNP}{MAF value for causal SNPs}

\item{causalSNPnb}{integer. number of SNP in each gene to considered as causal}

\item{causalSNPportion}{value between 0 and 1. Portion of SNP in each gene to considered as causal (if causalSNPnb is NULL).}

\item{AllMainPair}{logical. If FALSE, only the first gene of each interaction pair will have a main effect. SameMainPair has to be TRUE.}

\item{minMAF}{minimum value of the Minor Allele Frequency. Default \eqn{0.05}.}

\item{maxMAF}{maximum value of the Minor Allele Frequency. Default \eqn{0.5}.}
}
\value{
Returns a list including:
\item{X}{a matrix where columns are SNPs and rows are individuals.}
\item{listGenesSNP}{a list that indicate the names of the SNPs composing each gene.}
\item{MainEff}{a vector containing the names of genes having main effects.}
\item{GenePair}{a vector containing the names of genes having interaction effects.}
\item{MAF}{a vector that give the minor allele frequency observed for each simulated SNP.}
}
\description{
The fonction create a SNP data set structured by genes with genotypes coded by 1, 2 and 3.
}

