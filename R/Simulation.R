
#' Fonction to simulate genotype data 
#'
#' The fonction create a SNP data set structured by genes with genotypes coded by 1, 2 and 3.
#' @param N  integer; number of individuals.
#' @param corr numeric value between 0 and 1. Correlation level among the SNPs within each gene.
#' @param sizeGenesMain a vector of integers specifying the sizes of the genes having main effects.
#' @param sizeGenesPair  a vector of integers specifying the sizes of the genes having interaction effects. The size of the vector has to be an even number, pairs being defined with the gene taken two by two along the vector.
#' @param sizeGenesRemain a vector of integers specifying the sizes of the genes having none effects.
#' @param SameMainPair logical. If TRUE, genes with interation effects will also have main effects 
#' @param MAFcausalSNP MAF value for causal SNPs
#' @param causalSNPnb integer. number of SNP in each gene to considered as causal
#' @param causalSNPportion value between 0 and 1. Portion of SNP in each gene to considered as causal (if causalSNPnb is NULL). 
#' @param AllMainPair logical. If FALSE, only the first gene of each interaction pair will have a main effect. SameMainPair has to be TRUE.
#' @param minMAF minimum value of the Minor Allele Frequency. Default \eqn{0.05}.
#' @param maxMAF maximum value of the Minor Allele Frequency. Default \eqn{0.5}.
#' @return Returns a list including:
#' \item{X}{a matrix where columns are SNPs and rows are individuals.}
#' \item{listGenesSNP}{a list that indicate the names of the SNPs composing each gene.}
#' \item{MainEff}{a vector containing the names of genes having main effects.}
#' \item{GenePair}{a vector containing the names of genes having interaction effects.}
#' \item{MAF}{a vector that give the minor allele frequency observed for each simulated SNP.}
#' @export
simGeno <- function(N = 600, corr = 0.8, sizeGenesMain, sizeGenesPair, sizeGenesRemain, 
                    SameMainPair, MAFcausalSNP = NULL, causalSNPnb, causalSNPportion, AllMainPair = TRUE, 
                    minMAF = 0.05, maxMAF = 0.5) {
  genesSize <- c(sizeGenesMain, sizeGenesPair, sizeGenesRemain)
  nbSNP <- sum(genesSize)
  nbGenes <- length(genesSize)
  
  Z <- matrix(0, N, nbSNP)
  count <- 0
  matList <- lapply(genesSize, FUN = function(x) {
    m <- matrix(corr, x, x)
    diag(m) <- 1
    return(m)
  })
  for (bb in matList) {
    l.bb <- ncol(bb)
    Z.bb <- matrix(rnorm(l.bb * N), N, l.bb) %*% chol(bb)
    Z[, (count + 1):(count + l.bb)] <- Z.bb
    count <- count + l.bb
  }
  
  
  MAF <- runif(nbSNP, min = minMAF, max = maxMAF)
  if (is.null(MAFcausalSNP)) {
    
  } else {
    MAFlist <- split(MAF, rep(cumsum(genesSize), genesSize))
    if (is.null(causalSNPnb)) {
      MAFlist <- lapply(MAFlist, function(x) replace(x, 1:(causalSNPportion * 
                                                             length(x)), MAFcausalSNP))
      
    } else {
      MAFlist <- lapply(MAFlist, function(x) replace(x, 1:causalSNPnb, MAFcausalSNP))
    }
    MAF <- as.vector(unlist(MAFlist))
  }
  
  
  X <- sapply(1:nbSNP, FUN = function(ii) {
    maf <- MAF[ii]
    c1 <- qnorm(maf^2)
    c2 <- qnorm(2 * maf - maf^2)
    vect <- 1 * (Z[, ii] > c2) + 3 * (Z[, ii] < c1) + 2 * ((Z[, ii] < c2) & (Z[,ii] > c1))
    return(vect)
  })
  
  
  
  listGenesSNP <- vector("list", nbGenes)
  namesGenes <- c()
  groupsGenes <- c()
  for (i in seq(nbGenes)) {
    namesGenes = c(namesGenes, paste("Genes", i, sep = ""))
    groupsGenes <- c(groupsGenes, rep(namesGenes[i], genesSize[i]))
    namesSNP <- c()
    for (j in seq(genesSize[i])) {
      namesSNP <- c(namesSNP, paste("Gene", i, "SNP", j, sep = "."))
    }
    listGenesSNP[[i]] <- namesSNP
  }
  names(listGenesSNP) <- namesGenes
  colnames(X) <- as.vector(unlist(listGenesSNP))
  
  MAFobs <- apply(X, 2, function(x) (table(x)[2] + 2 * table(x)[3])/(N * 2))
  names(MAFobs) <- as.vector(unlist(listGenesSNP))
  
  
  
  
  if (is.null(sizeGenesMain)) {
    MainEff <- NULL
    idMainEff <- NULL
  } else {
    idMainEff <- seq(length(sizeGenesMain))
    MainEff <- names(listGenesSNP[idMainEff])
  }
  
  if (SameMainPair == TRUE) {
    if (is.null(sizeGenesPair)) {
      warning("sizeGenesPair is NULL, working genes selection impossible")
    } else {
      idGenePair <- (length(sizeGenesMain) + 1):(length(sizeGenesMain) + length(sizeGenesPair))
      GenePair <- names(listGenesSNP[idGenePair])
      
      if (AllMainPair == TRUE) {
        MainEff <- c(MainEff, GenePair)
        idMainEff <- c(idMainEff, idGenePair)
        sizeGenesMain <- c(sizeGenesMain, sizeGenesPair)
        
      } else {
        MainEff <- c(MainEff, GenePair[seq(1, length(GenePair), 2)])
        idMainEff <- c(idMainEff, idGenePair[seq(1, length(idGenePair), 2)])
        sizeGenesMain <- c(sizeGenesMain, sizeGenesPair[seq(1, length(GenePair), 2)])
        
      }
      
      
      
      
    }
  } else {
    if (is.null(sizeGenesPair)) {
      GenePair <- NULL
      idGenePair <- NULL
    } else {
      idGenePair <- (length(sizeGenesMain) + 1):(length(sizeGenesMain) + length(sizeGenesPair))
      GenePair <- names(listGenesSNP[idGenePair])
    }
  }
  
  if (is.null(sizeGenesRemain)) {
    GeneRemain <- NULL
    idGeneRemain <- NULL
  } else {
    if (is.null(sizeGenesMain)&is.null(sizeGenesPair)){
      GeneRemain <- names(listGenesSNP)
    }else{
      GeneRemain <- names(listGenesSNP[-c(idMainEff, idGenePair)])
    }
  }
  
  names(sizeGenesMain) <- MainEff 
  names(sizeGenesPair) <- GenePair
  names(sizeGenesRemain) <- GeneRemain
  
  sizeGenes <- list(sizeGenesMain=sizeGenesMain,sizeGenesPair=sizeGenesPair,sizeGenesRemain=sizeGenesRemain)
  
  
  list(X = X, listGenesSNP = listGenesSNP, MainEff = MainEff, GenePair = GenePair, GeneRemain=GeneRemain, sizeGenes=sizeGenes,
       MAF = MAFobs)
}



#' Fonction to simulate genotype data from a real data set
#'
#' @param data a matrix of real data with SNP in columns and individuals in rows.
#' @param listGenes a list that indicate the names of the SNPs composing each gene.
#' @param sizeGenesMain a vector of integers specifying the sizes of the genes having main effects.
#' @param sizeGenesPair  a vector of integers specifying the sizes of the genes having interaction effects. The size of the vector has to be an even number, pairs being defined with the gene taken two by two along the vector.
#' @param sizeGenesRemain a vector of integers specifying the sizes of the genes having none effects.
#' @param SameMainPair logical. If TRUE, genes with interation effects will also have main effects 
#' @param AllMainPair logical. If FALSE, only the first gene of each interaction pair will have a main effect. SameMainPair has to be TRUE.
#' @param sampling character. Can either be "gene" or "size". If "gene" is chosen, the genes are randomly selected in the initial data set whereas if it's "size", the genes are selected in the data set depending on their size, defined randomly.  
#' @return Returns a list including:
#' \item{X}{a matrix where columns are SNPs and rows are individuals.}
#' \item{listGenesSNP}{a list that indicate the names of the SNPs composing each gene.}
#' \item{MainEff}{a vector containing the names of genes having main effects.}
#' \item{GenePair}{a vector containing the names of genes having interaction effects.}
#' @export
simGenoDataR <- function(data, listGenesData, GenesMainNb, GenesPairNb, GenesRemainNb, sizeGenesMain, sizeGenesPair, 
                         sizeGenesRemain, SameMainPair, AllMainPair=TRUE, sampling="gene", GenesLengthMin=NULL, GenesActifLengthMin=NULL){
  
  genesNames <- names(listGenesData)
  genesLength <- sapply(listGenesData, function(x) length(x))
  
  # Genes actif selection
  # With sizeGenesMain and sizeGenesPair 
  if(is.null(GenesMainNb)| is.null(GenesPairNb) | is.null(GenesRemainNb)){
    
    #Genes Main
    if(is.null(sizeGenesMain)){
      MainEff <- NULL
      idMainEff <- NULL
      genesLength2 <- genesLength
      genesNames2 <- genesNames
      
    }else{
      selection <- selectGenesBySize(sizeGenesMain, genesNames, genesLength)
      MainEff <- selection$SelectedGenes
      idMainEff <- selection$idSelectedGenes
      sizeGenesMain <- genesLength[idMainEff]
      genesLength2 <- genesLength[- idMainEff]
      genesNames2 <- genesNames[- idMainEff]
    }
    
    #Genes Pair
    if(is.null(sizeGenesPair)){
      GenePair <- NULL
      idGenePair <- NULL 
    }else{
      selection <- selectGenesBySize(sizeGenesPair, genesNames, genesLength)
      GenePair <- selection$SelectedGenes
      idGenePair <- selection$idSelectedGenes
      sizeGenesPair <- genesLength[idGenePair]
    }
    
    
  }else{
    
    if(is.null(GenesActifLengthMin)){
      GenesActifLengthMin <- GenesLengthMin
    } 
    
    # With GenesMainNb and GenesPairNb  
    #Genes Main
    if(GenesMainNb == 0){
      idMainEff <- NULL
      MainEff <- NULL
      sizeGenesMain <- NULL
      genesLength2 <- genesLength
      genesNames2 <- genesNames
    }else{
      if(sampling=="gene"){
        idGenes <- seq(length(genesNames))
        if(!is.null(GenesActifLengthMin)){
          idSmallGenes <- which(genesLength< GenesActifLengthMin)
          idGenes <- idGenes[-idSmallGenes]
        }
        idMainEff <- sample(idGenes,GenesMainNb)
      }else if(sampling=="size"){
        sizeGenes <- unique(genesLength)
        if(!is.null(GenesActifLengthMin)){
          sizeGenes <- sizeGenes[which(sizeGenes>= GenesActifLengthMin)]
        }
        sizeGenesMain <- sample(sizeGenes, GenesMainNb, replace = TRUE)
        selection <- selectGenesBySize(sizeGenesMain, genesNames, genesLength)
        idMainEff <- selection$idSelectedGenes
      }
      MainEff <- genesNames[idMainEff]
      sizeGenesMain <- genesLength[idMainEff]
      genesLength2 <- genesLength[- idMainEff]
      genesNames2 <- genesNames[- idMainEff]
    }
    
    #Genes Pair 
    if(GenesPairNb == 0){
      idGenePair <- NULL
      GenePair <- NULL
      sizeGenesPair <- NULL
    }else{
      if(sampling=="gene"){
        idGenes <- seq(length(genesNames2))
        if(!is.null(GenesActifLengthMin)){
          idSmallGenes <- which(genesLength2< GenesActifLengthMin)
          idGenes <- idGenes[-idSmallGenes]
        }
        idGenePair <- sample(idGenes,GenesPairNb)
      }else if(sampling=="size"){
        sizeGenes <- unique(genesLength2)
        if(!is.null(GenesActifLengthMin)){
          sizeGenes <- sizeGenes[which(sizeGenes>= GenesActifLengthMin)]
        }
        sizeGenesPair <- sample(sizeGenes, GenesPairNb, replace = TRUE)
        selection <- selectGenesBySize(sizeGenesPair, genesNames2, genesLength2)
        idGenePair <- selection$idSelectedGenes
      }
      GenePair <- genesNames2[idGenePair]
      sizeGenesPair <- genesLength2[idGenePair]
    }
    
    
  }
  
  genesActifs <- names(genesLength[c(MainEff,GenePair)])
  
  if(SameMainPair == FALSE){	
  }else{
    if(AllMainPair == TRUE){
      MainEff <- c(MainEff, GenePair)
      sizeGenesMain <- c(sizeGenesMain, sizeGenesPair)
      
    }else{
      MainEff <- c(MainEff, GenePair[seq(1, length(GenePair), 2)])
      sizeGenesMain <- c(sizeGenesMain, sizeGenesPair[seq(1, length(GenePair), 2)])
    }
  }
  
  
  # GenesRemain selection
  idGenesActifs<- unlist(sapply(genesActifs, function(x) which(names(genesLength) == x)))
  if(is.null(idGenesActifs)){
    genesLength3 <- genesLength
    genesNames3 <- genesNames
  }else{
    genesLength3 <- genesLength[- idGenesActifs]
    genesNames3 <- genesNames[- idGenesActifs]
  }
  
  # With sizeGenesRemain   
  if(is.null(GenesMainNb)| is.null(GenesPairNb) | is.null(GenesPairNb)){ 
    
    if(is.null(sizeGenesRemain)){
      GeneRemain <- NULL
      idGeneRemain <- NULL
    }else{
      selection <- selectGenesBySize(sizeGenesRemain, genesNames, genesLength)
      GeneRemain <- selection$SelectedGenes
      idGeneRemain <- selection$idSelectedGenes
      sizeGenesRemain <- genesLength[idGeneRemain]
    }
    
  }else{
    
    # With GenesRemainNb 
    if(GenesRemainNb == 0){
      idGeneRemain <- NULL
      GeneRemain <- NULL
      sizeGenesRemain <- NULL
    }else{
      if(sampling=="gene"){
        idGenes <- seq(length(genesNames3))
        if(!is.null(GenesLengthMin)){
          idSmallGenes <- which(genesLength3< GenesLengthMin)
          idGenes <- idGenes[-idSmallGenes]
        }
        idGeneRemain <- sample(idGenes,GenesRemainNb)
      }else if(sampling=="size"){
        sizeGenes <- unique(genesLength3)
        if(!is.null(GenesLengthMin)){
          sizeGenes <- sizeGenes[which(sizeGenes>= GenesLengthMin)]
        }
        sizeGenesRemain <- sample(sizeGenes, GenesRemainNb, replace = TRUE)
        selection <- selectGenesBySize(sizeGenesRemain, genesNames3, genesLength3)
        idGeneRemain <- selection$idSelectedGenes
        
      }
      GeneRemain <- genesNames3[idGeneRemain]
      sizeGenesRemain <- genesLength3[idGeneRemain]
    }
    
  }
  genesRemain <- names(genesLength[GeneRemain])
  
  
  # Definition of the reduced data set  
  genes <- c(genesActifs, genesRemain)
  red <- choixGenes(genes, data, listGenesData)
  X <- red$X
  listGenesSNP <- red$listGenesRed
  
  sizeGenes <- list(sizeGenesMain=sizeGenesMain,sizeGenesPair=sizeGenesPair,sizeGenesRemain=sizeGenesRemain)
  
  MAFobs <- apply(X, 2, function(x) (table(x)[2] + 2 * table(x)[3])/(dim(X)[1] * 2))
  names(MAFobs) <- as.vector(unlist(listGenesSNP))
  
  
  
  list(X=X, listGenesSNP=listGenesSNP, MainEff=MainEff, GenePair=GenePair, GeneRemain=GeneRemain, sizeGenes=sizeGenes, 
       MAF = MAFobs)	
  
}

selectGenesBySize <- function(sizeGenes, genesNames, genesLength){
  sizeGenesUN <- unique(sizeGenes)
  occsizeGenes <- sapply(sizeGenesUN, function(x) length(which(sizeGenes==x)))
  if(length(sizeGenesUN) ==1){
    genesPossByLength <- genesNames[which(genesLength == sizeGenesUN)]
    SelectedGenes <- sample(genesPossByLength, occsizeGenes)
    idSelectedGenes <- unlist(sapply(SelectedGenes, function(x) which(names(genesLength) == x)))
  }else{
    genesPossByLength <- sapply(sizeGenesUN, function(x) genesNames[which(genesLength == x)])
    SelectedGenes <- unlist(sapply(seq(length(genesPossByLength)), function(x) sample(genesPossByLength[[x]], occsizeGenes[x])))
    idSelectedGenes <- unlist(sapply(SelectedGenes, function(x) which(names(genesLength) == x)))
  }
  list(SelectedGenes=SelectedGenes, idSelectedGenes=idSelectedGenes)
}


#' Fonction to simulate phenotype values
#'
#' The fonction simulate a quantitative phenotype assuming a gene-structure among the SNPs. 
#' @param X a matrix where columns are SNPs and rows are individuals.
#' @param listGenes a list that indicate the names of the SNPs composing each gene.
#' @param MainEff a vector containing the names of genes having main effects.
#' @param GenePair a vector containing the names of genes having interaction effects. The size of the vector GenePair has to be an even number, pairs being defined with the gene taken two by two along the vector.
#' @param model the model to consider to simulate interaction effects, either "SNPproduct" or "PCproduct".
#' @param pvBeta numerical vector of possible values for main effects regression coefficients.
#' @param pvGamma  numerical vector of possible values for interaction effects regression coefficients. 
#' @param r2 numeric value between 0 and 1. Coefficient of determination.
#' @param causalSNPnb numerical value. number of SNP in each gene to considered as causal.
#' @param causalSNPportion value between 0 and 1. Portion of SNP in each gene to considered as causal (if causalSNPnb is NULL). 
#' @param beta0 numeric value; intercept coefficient. Default 0.
#' @return Returns a list including:
#' \item{y}{vector of simulated phenotype continuous values}
#' \item{G}{matrix of the simulated main effects.}
#' \item{GG}{matrix of the simulated interaction effects.}
#' \item{R2T}{numerical value for the coefficient of determination R2 when considering the model containing simulated main and interaction effects.} 
#' \item{R2I}{numerical value for the coefficient of determination R2 when considering the model containing only simulated interaction effects.} 
#' \item{R2S}{numerical value for the coefficient of determination R2 when considering the model containing only simulated main effects.} 
#' \item{caract}{a list with the caracteristic elements of the simulation.}
#' @export
simPheno <- function(X, listGenes, MainEff, GenePair, model = "SNPproduct", pvBeta = c(2,2),
                     pvGamma = c(2, 2), r2 = 0.2, causalSNPnb, causalSNPportion = NULL, beta0 = 0, logical=FALSE) {
  
  if(is.null(MainEff) & is.null(GenePair)){
    stop("No actifs genes")
  }else{
    
    if (length(GenePair)%%2 == 1) {
      print("The number of gene in GenePair has to be an even number")
    }
    
    nbMainEff <- length(MainEff)
    nbSNPparMainEff <- sapply(MainEff, function(x) length(listGenes[[x]]))
    nbGenePair <- length(GenePair)/2
    nbGeneInt <- length(GenePair)
    nbSNPparGeneEnInter <- sapply(GenePair, function(x) length(listGenes[[x]]))
    genes <- names(listGenes)
    GeneRemain <- genes[-which(genes%in%c(MainEff,GenePair))]
    nbSNPparGeneRemain <- sapply(GeneRemain, function(x) length(listGenes[[x]]))
    
    # Main Gene effect
    if (is.null(MainEff)) {
      NamesCausalSNP_M <- NULL
      G <- matrix(0, nrow = dim(X)[1])
      beta <- c(0)
      
    } else {
      NamesCausalSNP_M <- sapply(MainEff, function(x) SNPinGene(x, portionSNP = causalSNPportion, 
                                                                causalSNPnb, listGenes))
      if (is.matrix(NamesCausalSNP_M)) {
        isMat <- TRUE
        m <- NamesCausalSNP_M 
        NamesCausalSNP_M = as.list(data.frame(NamesCausalSNP_M))
      } else {
        isMat <- FALSE
      }
      XCausalSNP_M <- X[, as.character(unlist(NamesCausalSNP_M))]
      causalSNPnbbyGene <- sapply(NamesCausalSNP_M, function(x) length(x))
      XCausalSNPagg <- aggregate(t(XCausalSNP_M), by = list(gene = rep(MainEff, causalSNPnbbyGene)), FUN = sum, na.rm = TRUE)
      nomsGene <- XCausalSNPagg[, 1]
      XCausalSNPagg <- XCausalSNPagg[, -1]
      G <- t(XCausalSNPagg)
      colnames(G) <- nomsGene
      G <- scale(G, center = TRUE, scale = TRUE)
      beta <- replicate(nbMainEff, sample(pvBeta, 1))
      names(beta) <- nomsGene
      
      if(isMat==FALSE){
        caM<-NamesCausalSNP_M 
        max <- max(sapply(caM, function(x) length(x)))
        m <- matrix(ncol=length(caM),nrow=max)
        colnames(m) <- names(caM)
        for(i in seq(length(caM))){
          m[1:length(caM[[i]]),i]<- unlist(caM[i])
        }
      }
      NamesCausalSNP_M <- m
    }
    
    
    if (is.null(GenePair)) {
      NamesCausalSNP_I <- NULL
      GG <- matrix(0, nrow = dim(X)[1])
      gamma <- c(0)
      
    } else {
      # Gene Pair effect
      GenePair2 <- matrix(GenePair, ncol = 2, byrow = TRUE)
      GenePair2 <- apply(GenePair2, 1, function(x) paste("X", x[1], x[2], sep = "."))
      
      
      if (model == "SNPproduct") {
        NamesCausalSNP_I <- sapply(GenePair, function(x) SNPinGene(x, causalSNPportion, 
                                                                   causalSNPnb, listGenes))
        if (is.matrix(NamesCausalSNP_I)) {
          isMat <- TRUE
          m <- NamesCausalSNP_I
          NamesCausalSNP_I = as.list(data.frame(NamesCausalSNP_I))
        } else {
          isMat <- FALSE
        }
        
        causalSNPnbbyGene <- sapply(NamesCausalSNP_I, function(x) length(x))
        GG <- matrix(, nrow = dim(X)[1])
        for (i in seq(1:(nbGeneInt/2)) * 2 - 1) {
          X1 <- as.matrix(X[, as.character(NamesCausalSNP_I[[i]])])
          X2 <- as.matrix(X[, as.character(NamesCausalSNP_I[[i + 1]])])
          X3 <- matrix(0)
          for (j in seq(dim(X1)[2])) {
            for (k in seq(dim(X2)[2])) {
              prod <- data.frame(X1[, j] * X2[, k])
              colnames(prod) <- paste("pair", i, "snp", j, k, sep = ".")
              X3 <- cbind(X3, prod)
            }
          }
          gg <- apply(X3, 1, sum)
          gg <- data.frame(gg)
          colnames(gg) <- paste("pair", i, i + 1, sep = ".")
          GG <- cbind(GG, gg)
        }
        GG <- as.matrix(GG[, -1])
        
        
        
      } else if (model == "PCproduct") {
       
        red <- choixGenes(genes = GenePair, X, listGenesSNP = listGenes)
        Xred <- red$Xred
        
        listGenesred <- listGenes[sapply(GenePair, function(x) which(names(listGenes) == 
                                                                       x))]
        
        PCA.L <- PCAGenes(X = Xred, listGenesSNP = listGenesred, nbcomp = 1)
        XBet <- PCA.L$I
        interLength <- PCA.L$interLength
        
        
        GG <- XBet[, sapply(GenePair2, function(x) which(names(interLength) == 
                                                           x))]
        
        NamesCausalSNP_I <- sapply(GenePair, function(x) SNPinGene(x, 1, NULL, 
                                                                   listGenes))
        if (is.matrix(NamesCausalSNP_I)) {
          isMat <- TRUE
          m <- NamesCausalSNP_I
          NamesCausalSNP_I = as.list(data.frame(NamesCausalSNP_I))
        } else {
          isMat <- FALSE
        }
        
        
      }
      
      GG <- scale(GG, center = TRUE, scale = TRUE)
      gamma <- replicate(nbGenePair, sample(pvGamma, 1))
      names(gamma) <- GenePair2
      colnames(GG) <- GenePair2
      
      if(isMat==FALSE){
        caI<-NamesCausalSNP_I 
        max <- max(sapply(caI, function(x) length(x)))
        m <- matrix(ncol=length(caI),nrow=max)
        colnames(m) <- names(caI)
        for(i in seq(length(caI))){
          m[1:length(caI[[i]]),i]<- unlist(caI[i])
        }
      }
      NamesCausalSNP_I <- m
      
    }
    
 
    XComp <- cbind(G, GG)
    
    if(logical==TRUE){
      #y binaire
      y.ex <- beta0 + G %*% beta + GG %*% gamma 
      prob<- exp(y.ex)/(1+exp(y.ex))
      y <- sapply(prob, function(x) rbinom(1,1,x)) 
      
      nullmod <- glm(y~1, family="binomial")
      modT <- glm(y ~ XComp - 1, family=binomial(link='logit'))
      modS <- glm(y ~ G - 1, family=binomial(link='logit'))
      modI <- glm(y ~ GG - 1, family=binomial(link='logit'))
      
      #McFadden Pseudo-R-squared
      R2T <- 1-logLik(modT)/logLik(nullmod)
      R2S <- 1-logLik(modS)/logLik(nullmod) 
      R2I <- 1-logLik(modI)/logLik(nullmod) 
      
    }else{
      # Terme d'erreur
      Coef <- c(beta, gamma)
      ybar <- mean(XComp %*% Coef)
      sigma <- sum((XComp %*% Coef - ybar)^2) * (r2 - 1)/((2 - dim(X)[1]) * r2)
      epsilon <- rnorm(dim(X)[1], 0, sqrt(sigma))
      
      # y continu
      y <- beta0 + G %*% beta + GG %*% gamma + epsilon
      
      regT <- lm(y ~ XComp - 1)
      regS <- lm(y ~ G - 1)
      regI <- lm(y ~ GG - 1)
      
      R2T <- summary(regT)$r.squared
      R2S <- summary(regS)$r.squared
      R2I <- summary(regI)$r.squared
    }
    
    
    # R2I/R2T
    PartR2I <- R2I/R2T * 100
    PartR2S <- R2S/R2T * 100
    
  }
   
  
  list(y = as.matrix(y), G = G, GG = GG, R2T = R2T, R2I = R2I, R2S = R2S, caract = list(MainEff = MainEff, 
                                                                                        nbSNPbyMainEff = nbSNPparMainEff, Coef_MainEff = beta, causalSNPMainEff = NamesCausalSNP_M, 
                                                                                        GenePair = GenePair, nbSNPbyInterGene = nbSNPparGeneEnInter, Coef_GenePair = gamma, 
                                                                                        causalSNPInter = NamesCausalSNP_I, GeneRemain=GeneRemain, nbSNPparGeneRemain=nbSNPparGeneRemain,
                                                                                        beta0 = beta0, r2 = r2, causalSNPportion = causalSNPportion, 
                                                                                        causalSNPnb = causalSNPnb, R2T = R2T, PartR2I = PartR2I, PartR2S = PartR2S))
}


#' Fonction to simulate phenotype values without any genetic effects
simPhenoWE <- function(data, listGenesData, listGenesRemain, r2 = 0.2, causalSNPportion = NULL, beta0 = 0, logical=FALSE) {
  
  genes <- names(listGenesData)
  genesRemain <- names(listGenesRemain)
  nbSNPparGeneRemain <- sapply(genesRemain, function(x) length(listGenesRemain[[x]]))
  idGenesRemain <- sapply(genesRemain, function(x) which(genes==x))
  SelectableGenes <- genes[-idGenesRemain]
  listSelectableGenes <- listGenesData[-idGenesRemain]
  
  idSelectableGenes <- 1:length(SelectableGenes)
  idMainEff <- sample(idSelectableGenes,2)
  MainEff <- SelectableGenes[idMainEff]
  idGenePair <- sample(idSelectableGenes[-idMainEff],2)
  GenePair <- SelectableGenes[idGenePair]
  
  # Main Eff
  NamesCausalSNP_M <- sapply(MainEff, function(x) SNPinGene(x, portionSNP = causalSNPportion, 
                                                            nbSNPcausaux=2, listSelectableGenes))
  if (is.matrix(NamesCausalSNP_M)) {
    isMat <- TRUE
    m <- NamesCausalSNP_M 
    NamesCausalSNP_M = as.list(data.frame(NamesCausalSNP_M))
  } else {
    isMat <- FALSE
  }
  XCausalSNP_M <- data[, as.character(unlist(NamesCausalSNP_M))]
  causalSNPnbbyGene <- sapply(NamesCausalSNP_M, function(x) length(x))
  XCausalSNPagg <- aggregate(t(XCausalSNP_M), by = list(gene = rep(MainEff, causalSNPnbbyGene)), FUN = sum, na.rm = TRUE)
  nomsGene <- XCausalSNPagg[, 1]
  XCausalSNPagg <- XCausalSNPagg[, -1]
  G <- t(XCausalSNPagg)
  colnames(G) <- nomsGene
  G <- scale(G, center = TRUE, scale = TRUE)
  beta <- c(2, 2)
  names(beta) <- nomsGene
  
  if(isMat==FALSE){
    caM<-NamesCausalSNP_M 
    max <- max(sapply(caM, function(x) length(x)))
    m <- matrix(ncol=length(caM),nrow=max)
    colnames(m) <- names(caM)
    for(i in seq(length(caM))){
      m[1:length(caM[[i]]),i]<- unlist(caM[i])
    }
  }
  NamesCausalSNP_M <- m
  
  
  # Gene Pair (SNP product)
  GenePair2 <- matrix(GenePair, ncol = 2, byrow = TRUE)
  GenePair2 <- apply(GenePair2, 1, function(x) paste("X", x[1], x[2], sep = "."))
  NamesCausalSNP_I <- sapply(GenePair, function(x) SNPinGene(x, causalSNPportion, 
                                                               nbSNPcausaux = 2, listSelectableGenes))
  if (is.matrix(NamesCausalSNP_I)) {
    isMat <- TRUE
    m <- NamesCausalSNP_I
    NamesCausalSNP_I = as.list(data.frame(NamesCausalSNP_I))
  } else {
    isMat <- FALSE
  }
    
  causalSNPnbbyGene <- sapply(NamesCausalSNP_I, function(x) length(x))
  X1 <- as.matrix(data[, as.character(NamesCausalSNP_I[[1]])])
  X2 <- as.matrix(data[, as.character(NamesCausalSNP_I[[2]])])
  X3 <- matrix(0)
  for (j in seq(dim(X1)[2])) {
    for (k in seq(dim(X2)[2])) {
      prod <- data.frame(X1[, j] * X2[, k])
      colnames(prod) <- paste("snp", j, k, sep = ".")
      X3 <- cbind(X3, prod)
      }
  }
  GG <- apply(X3, 1, sum)
  GG <- data.frame(GG)
  colnames(GG) <- "pair.1.2"
  GG <- scale(GG, center = TRUE, scale = TRUE)
  gamma <- 2
  names(gamma) <- GenePair2
  colnames(GG) <- GenePair2
  
  if(isMat==FALSE){
    caI<-NamesCausalSNP_I 
    max <- max(sapply(caI, function(x) length(x)))
    m <- matrix(ncol=length(caI),nrow=max)
    colnames(m) <- names(caI)
    for(i in seq(length(caI))){
      m[1:length(caI[[i]]),i]<- unlist(caI[i])
    }
  }
  NamesCausalSNP_I <- m
  
  # Permutation
  XComp1 <- cbind(G, GG)
  XComp <- XComp1[sample(nrow(XComp1)),]
  
  # Terme d'erreur
  Coef <- c(beta, gamma)
  ybar <- mean(XComp %*% Coef)
  sigma <- sum((XComp %*% Coef - ybar)^2) * (r2 - 1)/((2 - dim(X)[1]) * r2)
  epsilon <- rnorm(dim(X)[1], 0, sqrt(sigma))
  

  if(logical==TRUE){
    #y binaire
    y.ex <- beta0 + XComp %*% Coef 
    prob<- exp(y.ex)/(1+exp(y.ex))
    y <- sapply(prob, function(x) rbinom(1,1,x)) 
    
  }else{
    # y continu
    y <- beta0 + XComp %*% Coef + epsilon
  }
  
  
  
  # R2I/R2T
  regT <- lm(y ~ XComp - 1)
  R2T <- summary(regT)$r.squared
  
  R2S <- 0
  R2I <- 0
  PartR2I <- 0
  PartR2S <- 0
  
  caract <- list(MainEff = NULL, nbSNPbyMainEff = 0, Coef_MainEff = NULL, causalSNPMainEff = NULL, 
                 GenePair = NULL, nbSNPbyInterGene = 0, Coef_GenePair = NULL, causalSNPInter = NULL, 
                 GeneRemain=genesRemain, nbSNPparGeneRemain=nbSNPparGeneRemain,
                 beta0 = beta0, r2 = r2, causalSNPportion = causalSNPportion, 
                 causalSNPnb = NULL, R2T = R2T, PartR2I = PartR2I, PartR2S = PartR2S, 
                 namesG = colnames(G), namesGG=colnames(GG))
  
  
  list(y = as.matrix(y), G = G, GG = GG, R2T = R2T, R2I = R2I, R2S = R2S, caract = caract)
  }




#' Fonction to simulate phenotype values as in GGI package
simPhenoGGI <- function(X, listGenes, GenePair, PairCausalSNPnb, beta0 = 0, tau=2) {
  
  if (length(GenePair)%%2 == 1) {
    print("The number of gene in GenePair has to be an even number")
  }
  
  nbGenePair <- length(GenePair)/2
  nbGeneInt <- length(GenePair)
  nbSNPparGeneEnInter <- sapply(GenePair, function(x) length(listGenes[[x]]))
  genes <- names(listGenes)
  GeneRemain <- genes[-which(genes%in%GenePair)]
  nbSNPparGeneRemain <- sapply(GeneRemain, function(x) length(listGenes[[x]]))
  N <- dim(X)[1]
  
  if (is.null(GenePair)) {
    NamesCausalSNP_I <- NULL
    GG <- matrix(0, nrow = N)
    gamma <- c(0)
    
  } else {
    # Gene Pair effect
    GenePair2 <- matrix(GenePair, ncol = 2, byrow = TRUE)
    GenePair2 <- apply(GenePair2, 1, function(x) paste("X", x[1], x[2], sep = "."))
    
    
    
    NamesCausalSNP_I <- sapply(GenePair, function(x) SNPinGene(x, nbSNPcausaux=PairCausalSNPnb, 
                                                               listGenes=listGenes))
    if (is.matrix(NamesCausalSNP_I)) {
      isMat <- TRUE
      m <- NamesCausalSNP_I
      NamesCausalSNP_I = as.list(data.frame(NamesCausalSNP_I))
    } else {
      isMat <- FALSE
    }
    
    GG <- matrix(, nrow = N)
    gamma <- c()
    for (i in seq(1:(nbGeneInt/2)) * 2 - 1) {
      X1 <- as.matrix(X[, as.character(NamesCausalSNP_I[[i]])])
      X2 <- as.matrix(X[, as.character(NamesCausalSNP_I[[i + 1]])])
      X3 <- matrix(0)
      for (j in seq(dim(X1)[2])) {
        SNP1 <- X1[, j]
        SNP2 <- X2[, j]
        prod <- data.frame(SNP1 * SNP2)
        colnames(prod) <- paste("pair", i, i+1, "snp", j, j, sep = ".")
        X3 <- cbind(X3, prod)
        
        tSNP1 <- table(SNP1)
        tSNP2 <- table(SNP2)
        
        if(length(unique(SNP1))==2){
          namestSNP1 <- names(tSNP1)
          
          alelePres <- unique(SNP1)
          aleleMiss <-  which(c(1,2,3) %in% alelePres==FALSE)
          
          t2SNP1 <- matrix(nrow = 1, ncol = 3)
          t2SNP1[,as.numeric(namestSNP1[1])] <-  tSNP1[namestSNP1[1]]
          t2SNP1[,as.numeric(namestSNP1[2])] <-  tSNP1[namestSNP1[2]]
          t2SNP1[,aleleMiss] <- 0
          colnames(t2SNP1) <- c(1,2,3)
          
          MAF1 <- (t2SNP1[2] + 2 * t2SNP1[3])/(N * 2) 
          
          }else{
          MAF1 <- (tSNP1[2] + 2 * tSNP1[3])/(N * 2) 
          }
        

        if(length(unique(SNP2))==2){
          namestSNP2 <- names(tSNP2)
          
          alelePres <- unique(SNP2)
          aleleMiss <-  which(c(1,2,3) %in% alelePres==FALSE)
          
          t2SNP2 <- matrix(nrow = 1, ncol = 3)
          t2SNP2[,as.numeric(namestSNP2[1])] <-  tSNP2[namestSNP2[1]]
          t2SNP2[,as.numeric(namestSNP2[2])] <-  tSNP2[namestSNP2[2]]
          t2SNP2[,aleleMiss] <- 0
          colnames(t2SNP2) <- c(1,2,3)
          
          MAF2 <- (t2SNP2[2] + 2 * t2SNP2[3])/(N * 2) 
          
        }else{
          MAF2 <- (tSNP2[2] + 2 * tSNP2[3])/(N * 2) 
        }
        
        MAF1 <- min(MAF1, 1 - MAF1)
        MAF2 <- min(MAF2, 1 - MAF2)
        
        gamma_jj <- log(tau)/16 * abs(log10(MAF1 * MAF2 ))
        names(gamma_jj) <- paste("pair", i, i+1, "snp", j, j, sep = ".") 
        gamma <- c(gamma, gamma_jj)
        
      }
      GG <- cbind(GG, X3[,-1])
    }
    GG <- as.matrix(GG[, -1])
    
     
    
    #GG <- scale(GG, center = TRUE, scale = TRUE)
    #gamma <- replicate(nbGenePair, sample(pvGamma, 1))
    
    #names(gamma) <- GenePair2
    #colnames(GG) <- GenePair2
    
    if(isMat==FALSE){
      caI<-NamesCausalSNP_I 
      max <- max(sapply(caI, function(x) length(x)))
      m <- matrix(ncol=length(caI),nrow=max)
      colnames(m) <- names(caI)
      for(i in seq(length(caI))){
        m[1:length(caI[[i]]),i]<- unlist(caI[i])
      }
    }
    NamesCausalSNP_I <- m
    
  }
  
  
  #y binaire
  y.ex <- beta0 + GG %*% gamma 
  prob<- exp(y.ex)/(1+exp(y.ex))
  y <- sapply(prob, function(x) rbinom(1,1,x)) 
  
  nullmod <- glm(y~1, family="binomial")
  modI <- glm(y ~ GG - 1, family=binomial(link='logit'))
  
  #McFadden Pseudo-R-squared
  R2I <- 1-logLik(modI)/logLik(nullmod) 
  
  
  
  list(y = as.matrix(y), GG = GG, R2I = R2I, prob=prob,
       caract = list(GenePair = GenePair, nbSNPbyInterGene = nbSNPparGeneEnInter, 
                     Coef_GenePair = gamma, causalSNPInter = NamesCausalSNP_I, GeneRemain=GeneRemain, 
                     nbSNPparGeneRemain=nbSNPparGeneRemain, beta0 = beta0, 
                     PairCausalSNPnb=PairCausalSNPnb, tau=tau))
}

