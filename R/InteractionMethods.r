#' @useDynLib GGEE
#' @importFrom Rcpp sourceCpp



GGEE <- function(X, Y, listGenesSNP) {
    
    genes <- names(listGenesSNP)
    
    NamesSNP <- sapply(genes, function(x) SNPinGene(x, 1, nbSNPcausaux = NULL, listGenes = listGenesSNP))
    
    if (is.matrix(NamesSNP)) {
        NamesSNP = as.list(data.frame(NamesSNP))
    }
    
    
    nbSNPbyGene <- sapply(NamesSNP, function(x) length(x))
    Z <- matrix(, nrow = dim(X)[1])
    for (i in 1:(length(genes) - 1)) {
        X1 <- as.matrix(X[, as.character(NamesSNP[[i]])])
        for (k in (i + 1):(length(genes))) {
            X2 <- as.matrix(X[, as.character(NamesSNP[[k]])])
			print(paste("X", genes[i], genes[k], sep="."))
  		
			resProd <- IntProd(X1,X2, i, k)
			W <- as.matrix(resProd$W)
			colnames(W) <- resProd$names
            W <- scale(W, center = TRUE, scale = TRUE)
            
            
            A <- t(W) %*% Y
            A <-as.vector(A)
			u <- A/sqrt(A%*%A) 
			#u <- A/sqrt(sum(A^2)) 

            z <- W %*% u
            colnames(z) <- paste("X", genes[i], genes[k], sep = ".")
            Z <- cbind(Z, z)
 
            
        }
    }
    Z <- as.matrix(Z[, -1])
    Z <- scale(Z, center = TRUE, scale = TRUE)
    
    
    interLength <- rep(1, dim(Z)[2])
    names(interLength) <- colnames(Z)
    return(list(Int = Z, interLength = interLength))
    
}


SNPinGene <- function(gene, portionSNP, nbSNPcausaux, listGenes) {
   	nbSNP <- length(listGenes[[gene]])
    if(nbSNP == 1){
		listGenes[[gene]]
    }else{
    	if(is.null(nbSNPcausaux)){
  	 	nbCausalSNP <- nbSNP*portionSNP
  	 	}
		else{
		nbCausalSNP <- nbSNPcausaux
		}
		listGenes[[gene]][1:nbCausalSNP]
     }
}
  


between <- function(X1, X2, nbcomp, nameX1, nameX2) {
    nbcomp1 = nbcomp
    nbcomp2 = nbcomp
    nbrow = length(X1)
    if (is.vector(X1) == FALSE & is.vector(X2) == FALSE) {
        nbrow = dim(X1)[1]
        colnames(X1) <- c(1:dim(X1)[2])
        colnames(X2) <- c(1:dim(X2)[2])
        cca <- cancor(X1, X2)
        a1 <- cca$xcoef[, 1]
        a2 <- cca$ycoef[, 1]
        if (length(a1) < dim(X1)[2]) {
            n1 <- which( !(1:dim(X1)[2] %in% as.numeric(names(a1))))
            X1 <- X1[, -n1]
            # suppression des variables qui n'ont pas de coef dans a1
        }
        
        if (length(a2) < dim(X2)[2]) {
            n2 <- which(!(1:dim(X2)[2] %in% as.numeric(names(a2))))
            X2 <- X2[, -n2]
        }
    }
    
    
    if (is.vector(X1) == TRUE & is.vector(X2) == FALSE) {
        nbrow = length(X1)
        colnames(X2) <- c(1:dim(X2)[2])
        cca <- cancor(X1, X2)
        a1 <- cca$xcoef[, 1]
        a2 <- cca$ycoef[, 1]
        
        if (length(a2) < dim(X2)[2]) {
            n2 <- which(!(1:dim(X2)[2] %in% as.numeric(names(a2))))
            X2 <- X2[, -n2]
        }
    }
    
    if (is.vector(X1) == FALSE & is.vector(X2) == TRUE) {
        nbrow = dim(X1)[1]
        colnames(X1) <- c(1:dim(X1)[2])
        cca <- cancor(X1, X2)
        a1 <- cca$xcoef[, 1]
        a2 <- cca$ycoef[, 1]
        if (length(a1) < dim(X1)[2]) {
            n1 <- which(!(1:dim(X1)[2] %in% as.numeric(names(a1))))
            X1 <- X1[, -n1]
        }
    }
    
    cca2 <- cancor(X1, X2)
    
    
    if (dim(cca2$xcoef)[2] < nbcomp1) {
        nbcomp1 = dim(cca2$xcoef)[2]
    }
    if (dim(cca2$ycoef)[2] < nbcomp2) {
        nbcomp2 = dim(cca2$ycoef)[2]
    }
    
    A1 <- cca2$xcoef[, 1:nbcomp1]
    A2 <- cca2$ycoef[, 1:nbcomp2]
    
    I <- data.frame(matrix(1, nrow = nbrow))
    
    for (i in seq(min(nbcomp1, nbcomp2))) {
        a1 <- as.matrix(A1)[, i]
        a2 <- as.matrix(A2)[, i]
        Z1 <- t(a1 %*% t(X1))
        Z2 <- t(a2 %*% t(X2))
        int <- Z1 * Z2
        int <- data.frame(int)
        names(int) <- paste(nameX1, nameX2, i, i, sep = ".")
        I <- cbind(I, int)
    }
    
    I <- I[, -1]
    if (is.vector(I)) {
        I <- as.matrix(I)
        colnames(I) <- names(int)
    }
   # nbVar <- nbcomp
    nbVar <- dim(I)[2]
    names(nbVar) <- paste("X", nameX1, nameX2, sep = ".")
    
    return(list(Int = I, nbVar = nbVar))
}




between.mat <- function(X, G, nbcomp) {
    if (is.null(nbcomp)) {
        nbcomp = 1
    }
    genes <- unique(G)
    nbGroup <- length(levels(as.factor(G)))
    B <- data.frame(matrix(1, nrow = dim(X)[1]))
    d <- data.frame(matrix(1, ncol = 2))
    interLength <- c()
    for (i in 1:(nbGroup - 1)) {
        f <- i + 1
        for (k in f:nbGroup) {
            X1 <- X[, G == genes[i]]  #1er groupe de variables
            X2 <- X[, G == genes[k]]  #2em groupe de variables
            betw <- between(X1, X2, nbcomp, nameX1 = genes[i], nameX2 = genes[k])
            newVar <- betw$Int
            # newVar<-data.frame(newVar) colnames(newVar) <-paste('X', genes[i], genes[k],
            # sep='.')
            B <- data.frame(B, newVar)
            d <- rbind(d, c(i, k))
            nbVar <- betw$nbVar
            # names(nbVar) <- paste(genes[i], genes[k], sep='.')
            interLength <- c(interLength, nbVar)
            
        }
    }
    nam <- colnames(B)[-1]
    B <- as.data.frame(B[, -1])
    colnames(B) <- nam
    d <- d[-1, ]
    return(list(XBet = B, ident = d, interLength = interLength, nbVar = nbVar))
}




PCAGenes <- function(X, listGenesSNP, nbcomp) {
    nbGenes <- length(listGenesSNP)
    namesGenes <- names(listGenesSNP)
    
    Xacp <- c()
    allnbComp <- c()
    ResACPvar <- list()
    for (i in seq(namesGenes)) {
        red <- choixGenes(genes = namesGenes[i], X, listGenesSNP = listGenesSNP)
        Xred <- red$Xred
        
        if(dim(Xred)[2]== 1){
			newVar <- Xred
			colnames(newVar) <- paste(namesGenes[i], 1, sep =".")
			allnbComp <-c(allnbComp,1)
		}else{
		
			ResACP <- FactoMineR::PCA(Xred, scale.unit = TRUE, ncp = dim(Xred)[2], graph = F)
        	eig <- list(ResACP$eig)
       		names(eig) <- namesGenes[i]
       	 	ResACPvar <- c(ResACPvar, eig)
        	
        	if (dim(Xred)[2]< nbcomp){
				nbComp <- dim(Xred)[2]
			}else{
				nbComp <- nbcomp
			}	
        	allnbComp <- c(allnbComp, nbComp)
        
       	 	newVar <- as.matrix(ResACP$ind$coord[, c(1:nbComp)])
        	namesVar <- c()
        	for (j in seq(nbComp)) {
            	namesVar <- c(namesVar, paste(namesGenes[i], j, sep = "."))
        	}
        	colnames(newVar) <- namesVar
        }
        
        Xacp <- cbind(Xacp, newVar)
    }
    
    G <- rep(seq(nbGenes), allnbComp)
    I <- data.frame(matrix(1, nrow = dim(Xacp)[1]))
    interLength <- c()
    d <- data.frame(matrix(1, ncol = 2))
    c <- c()
    for (g1 in 1:(nbGenes - 1)) {
        for (g2 in (g1 + 1):nbGenes) {
            nbVar1 <- allnbComp[g1]
            nbVar2 <- allnbComp[g2]
            
            for (i in seq(nbVar1)) {
                for (k in seq(nbVar2)) {
                  X1 <- as.matrix(Xacp[, G == g1])[, i]
                  X2 <- as.matrix(Xacp[, G == g2])[, k]
                  Y <- X1 * X2
                  Y <- data.frame(Y)
                  colnames(Y) <- paste(namesGenes[g1], namesGenes[g2], i, k, sep = ".")
                  I <- cbind(I, Y)
                }
            }
            nbVar <- nbVar1 * nbVar2
            names(nbVar) <- paste("X", namesGenes[g1], namesGenes[g2], sep = ".")
            interLength <- c(interLength, nbVar)
            d <- rbind(d, c(g1, g2))
        }
    }
    I <- I[, -1]
    d <- d[-1, ]
    
    if (is.vector(I)) {
        I <- as.matrix(I)
        colnames(I) <- colnames(Y)
    }
    
    XacpwInt <- as.matrix(cbind(Xacp, I))
    
    return(list(XacpwInt = XacpwInt, Int = I, interLength = interLength, ResACP = ResACPvar))
}


# Allow to select a reduce number of gene among a dataset 

choixGenes <- function(genes, X, listGenesSNP){
	genesLength <- sapply(listGenesSNP, function(x) length(x))
	nbSNPNames <- sum(genesLength)
	nbSNP <- dim(X)[2]
	if(nbSNPNames < nbSNP){print("Warning, SNPs from matrix X are not referenced in the gene list")}
	
	SNPaGarder <-c()
	groupsGenes <- c()
	listGenesRed <- c()
	
	for(i in seq(genes)){
		SNPaGarder <- c(SNPaGarder,unlist(listGenesSNP[genes[[i]]]))
		groupsGenes <- c(groupsGenes,rep(genes[[i]], genesLength[genes[[i]]]))
		listGenesRed <- c(listGenesRed, listGenesSNP[genes[[i]]])
	}
	Xred <- as.matrix(X[,SNPaGarder])
	colnames(Xred)<- SNPaGarder
	list(Xred=Xred, groupsGenes=groupsGenes, listGenesRed =listGenesRed)
}



PLSGenes <- function(X, Y, listGenesSNP, nbcomp = NULL) {
    nbGenes <- length(listGenesSNP)
    namesGenes <- names(listGenesSNP)
    
    I <- data.frame(matrix(1, nrow = dim(X)[1]))
    d <- data.frame(matrix(1, ncol = 2))
    interLength <- c()
    
    for (i in 1:(nbGenes - 1)) {
        for (j in (i + 1):nbGenes) {
            
            red1 <- choixGenes(genes = namesGenes[i], X, listGenesSNP = listGenesSNP)
            g1 <- red1$Xred
            red2 <- choixGenes(genes = namesGenes[j], X, listGenesSNP = listGenesSNP)
            g2 <- red2$Xred
            
             
            if (is.null(nbcomp)) {
                pls = pls::plsr(as.matrix(cbind(Y, g1)) ~ g2, validation = "LOO")
            } else {
                if (dim(cbind(Y, g1))[2] < nbcomp | dim(g2)[2] < nbcomp) {
                    nbcomp = min(dim(cbind(Y, g1))[2], dim(g2))
                  }
                  pls = pls::plsr(as.matrix(cbind(Y, g1)) ~ g2, ncomp = nbcomp, validation = "LOO")
            }
            scPLS <- pls::scores(pls)
            nbCompPLS <- dim(scPLS)[2]
            scPLS <- data.frame(scPLS[1:dim(X)[1], 1:nbCompPLS])
                
             
            noms <- c()
            for (k in seq(nbCompPLS)) {
                noms <- c(noms, paste(namesGenes[i], namesGenes[j], k, sep = "."))
            }
            colnames(scPLS) <- noms
            names(nbCompPLS) <- paste("X", namesGenes[i], namesGenes[j], sep = ".")
            interLength <- c(interLength, nbCompPLS)
            I = cbind(I, scPLS)
            
            d <- rbind(d, c(i, j))
        }
    }
    I <- I[, -1]
    d <- d[-1, ]
    
    return(list(Int = I, interLength = interLength))
}
 
