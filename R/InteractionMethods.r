#' @useDynLib GGEE
#' @importFrom Rcpp sourceCpp



GGEE <- function(X, Y, listGenesSNP, idSubs) {
  
  genes <- names(listGenesSNP)
  
  NamesSNP <- sapply(genes, function(x) SNPinGene(x, 1, nbSNPcausaux = NULL, listGenes = listGenesSNP))
  
  if (is.matrix(NamesSNP)) {
    NamesSNP = as.list(data.frame(NamesSNP))
  }
  
  
  nbSNPbyGene <- sapply(NamesSNP, function(x) length(x))
  nb_int <- 0
  Ztrain <- data.frame(matrix(, nrow = length(idSubs$idTrain)), row.names = rownames(X)[idSubs$idTrain])
  Ztest <- data.frame(matrix(, nrow = length(idSubs$idTest)), row.names = rownames(X)[idSubs$idTest])
  
  
  for (i in 1:(length(genes) - 1)) {
    X1 <- as.matrix(X[, as.character(NamesSNP[[i]])])
    for (k in (i + 1):(length(genes))) {
      X2 <- as.matrix(X[, as.character(NamesSNP[[k]])])
      print(paste(genes[i], genes[k], sep="."))
      
      resProd <- IntProd(X1,X2, i, k)
      W <- as.matrix(resProd$W)
      colnames(W) <- resProd$names
      W <- scale(W, center = TRUE, scale = TRUE)
      
      W_tr <- W[idSubs$idTrain,]
      W_te <- W[idSubs$idTest,]
      Y_tr <- Y[idSubs$idTrain]
      
      A <- t(W_tr) %*% Y_tr
      A <-as.vector(A)
      u <- A/sqrt(A%*%A) 
      #u <- A/sqrt(sum(A^2)) 
      
      z_tr <- W_tr %*% u
      z_te <- W_te %*% u
      colnames(z_tr) <- paste("X", genes[i], genes[k], sep = ".")
      colnames(z_te) <- paste("X", genes[i], genes[k], sep = ".")
      Ztrain <- cbind(Ztrain, z_tr)
      Ztest <- cbind(Ztest, z_te)
      
      nb_int <- nb_int + dim(W)[2]
    }
  }
  namesINT <- colnames(Ztrain)[-1]
  Ztrain <- data.frame(as.matrix(Ztrain[, -1]))
  Ztest <- data.frame(as.matrix(Ztest[, -1]))
  #Z <- scale(Z, center = TRUE, scale = TRUE)
  #scale realized in GLmodel
  
  
  interLength <- rep(1, dim(Ztrain)[2])
  names(interLength) <- namesINT
  return(list(IntTrain = Ztrain, IntTest = Ztest, interLength = interLength, nb_int=nb_int))
  
}




RFW <- function(X, Y, listGenesSNP, idSubs) {
  
  genes <- names(listGenesSNP)
  NamesSNP <- sapply(genes, function(x) SNPinGene(x, 1, nbSNPcausaux = NULL, listGenes = listGenesSNP))
  if (is.matrix(NamesSNP)) {
    NamesSNP = as.list(data.frame(NamesSNP))
  }
  
  
  nbSNPbyGene <- sapply(NamesSNP, function(x) length(x))
  nb_int <- 0
  Ztrain <- data.frame(matrix(, nrow = length(idSubs$idTrain)), row.names = rownames(X)[idSubs$idTrain])
  Ztest <- data.frame(matrix(, nrow = length(idSubs$idTest)), row.names = rownames(X)[idSubs$idTest])
  
  for (i in 1:(length(genes) - 1)) {
    X1 <- as.matrix(X[, as.character(NamesSNP[[i]])])
    for (k in (i + 1):(length(genes))) {
      X2 <- as.matrix(X[, as.character(NamesSNP[[k]])])
      print(paste(genes[i], genes[k], sep="."))
      
      resProd <- IntProd(X1,X2, i, k)
      W <- as.matrix(resProd$W)
      colnames(W) <- resProd$names
      W <- scale(W, center = TRUE, scale = TRUE)
      dat <- cbind(Y, W)
      
      #Interaction on training dataset
      rf.couple <- randomForest(Y~.,data=dat, subset=idSubs$idTrain, ntree=500)
      ztrain <- as.matrix(rf.couple$predicted)
      
      #Interaction on test dataset
      ztest <- as.matrix(predict(rf.couple, newdata=dat[idSubs$idTest,]))
      
      colnames(ztrain) <- paste("X", genes[i], genes[k], sep = ".")
      colnames(ztest) <- paste("X", genes[i], genes[k], sep = ".")
      
      Ztrain <- cbind(Ztrain, ztrain)
      Ztest <- cbind(Ztest, ztest)
    }
  }
  namesINT <- colnames(Ztrain)[-1]
  Ztrain <- data.frame(as.matrix(Ztrain[, -1]))
  Ztest <- data.frame(as.matrix(Ztest[, -1]))
  
  interLength <- rep(1, dim(Ztrain)[2])
  names(interLength) <- namesINT
  return(list(IntTrain = Ztrain, IntTest = Ztest, interLength = interLength, nb_int=nb_int))
}


RF <- function(X, Y, listGenesSNP, idSubs) {
  
  genes <- names(listGenesSNP)
  NamesSNP <- sapply(genes, function(x) SNPinGene(x, 1, nbSNPcausaux = NULL, listGenes = listGenesSNP))
  if (is.matrix(NamesSNP)) {
    NamesSNP = as.list(data.frame(NamesSNP))
  }
  
  
  nbSNPbyGene <- sapply(NamesSNP, function(x) length(x))
  nb_int <- 0
  Ztrain <- data.frame(matrix(, nrow = length(idSubs$idTrain)), row.names = rownames(X)[idSubs$idTrain])
  Ztest <- data.frame(matrix(, nrow = length(idSubs$idTest)), row.names = rownames(X)[idSubs$idTest])
  
  for (i in 1:(length(genes) - 1)) {
    X1 <- as.matrix(X[, as.character(NamesSNP[[i]])])
    for (k in (i + 1):(length(genes))) {
      X2 <- as.matrix(X[, as.character(NamesSNP[[k]])])
      print(paste(genes[i], genes[k], sep="."))
      
      dat <- data.frame(Y,X1,X2)
      
      #Interaction on training dataset
      rf.couple <- randomForest(Y~.,data=dat, subset=idSubs$idTrain, ntree=500)
      ztrain <- as.matrix(rf.couple$predicted)
      
      #Interaction on test dataset
      ztest <- as.matrix(predict(rf.couple, newdata=dat[idSubs$idTest,]))
      
      colnames(ztrain) <- paste("X", genes[i], genes[k], sep = ".")
      colnames(ztest) <- paste("X", genes[i], genes[k], sep = ".")
      
      Ztrain <- cbind(Ztrain, ztrain)
      Ztest <- cbind(Ztest, ztest)
    }
  }
  namesINT <- colnames(Ztrain)[-1]
  Ztrain <- data.frame(as.matrix(Ztrain[, -1]))
  Ztest <- data.frame(as.matrix(Ztest[, -1]))
  
  interLength <- rep(1, dim(Ztrain)[2])
  names(interLength) <- namesINT
  return(list(IntTrain = Ztrain, IntTest = Ztest, interLength = interLength, nb_int=nb_int))
}




SVMmet <- function(X, Y, listGenesSNP, kernel, degree=3, idSubs) {
  
  genes <- names(listGenesSNP)
  
  NamesSNP <- sapply(genes, function(x) SNPinGene(x, 1, nbSNPcausaux = NULL, listGenes = listGenesSNP))
  
  if (is.matrix(NamesSNP)) {
    NamesSNP = as.list(data.frame(NamesSNP))
  }
  
  
  nbSNPbyGene <- sapply(NamesSNP, function(x) length(x))
  nb_int <- 0
  Ztrain <- data.frame(matrix(, nrow = length(idSubs$idTrain)), row.names = rownames(X)[idSubs$idTrain])
  Ztest <- data.frame(matrix(, nrow = length(idSubs$idTest)), row.names = rownames(X)[idSubs$idTest])
  
  for (i in 1:(length(genes) - 1)) {
    X1 <- as.matrix(X[, as.character(NamesSNP[[i]])])
    for (k in (i + 1):(length(genes))) {
      X2 <- as.matrix(X[, as.character(NamesSNP[[k]])])
      print(paste(genes[i], genes[k], sep="."))
      
      dat <- data.frame(Y,X1,X2)
      
      #Interaction on training dataset
      #kernel="polynomial"
      #kernel="radial"
      #kernel="sigmoid"
      svm.couple <- svm(Y~.,data=dat, subset=idSubs$idTrain, kernel=kernel, degree=degree)
      #tune.out=tune(svm, Y~., data=dat[idSubs$idTrain,], kernel=kernel,
      #             ranges=list(cost=c(0.01, 0.1,1,10,100),
      #                          gamma=c(0.01, 0.05, 0.1, 0.5,1) ))
      #summary(tune.out)
      #svm.couple <- svm(Y~.,data=dat, subset=idSubs$idTrain, kernel=kernel, gamma=tune.out$best.model$gamma, cost=tune.out$best.model$cost)     
      
      ztrain <- as.matrix(svm.couple$fitted)
      #dec <- attributes(predict(svm.couple,dat[idSubs$idTrain,], decision.values=TRUE))$decision.values
      #sum((Y[idSubs$idTrain]-ztrain)^2)
      #sum((Y[idSubs$idTrain]-dec)^2)
      
      #Interaction on test dataset
      ztest <- as.matrix(predict(svm.couple, newdata=dat[idSubs$idTest,]))
      #dect <-attributes(predict(svm.couple,dat[idSubs$idTest,], decision.values=TRUE))$decision.values
      #sum((Y[idSubs$idTest]-ztest)^2)
      #sum((Y[idSubs$idTest]-dect)^2)
      
      #plot(Y[idSubs$idTest], ztest)
      #plot(Y[idSubs$idTest],dect)
      
      #plot(Y[idSubs$idTrain], ztrain)
      #plot(Y[idSubs$idTrain],dec)
      
      
      colnames(ztrain) <- paste("X", genes[i], genes[k], sep = ".")
      colnames(ztest) <- paste("X", genes[i], genes[k], sep = ".")
      
      Ztrain <- cbind(Ztrain, ztrain)
      Ztest <- cbind(Ztest, ztest)
      
    }
  }
  
  namesINT <- colnames(Ztrain)[-1]
  Ztrain <- data.frame(as.matrix(Ztrain[, -1]))
  Ztest <- data.frame(as.matrix(Ztest[, -1]))
  
  interLength <- rep(1, dim(Ztrain)[2])
  names(interLength) <- namesINT
  return(list(IntTrain = Ztrain, IntTest = Ztest, interLength = interLength, nb_int=nb_int))
  
}





BOOST <- function(X, Y, listGenesSNP, idSubs) {
  
  genes <- names(listGenesSNP)
  
  NamesSNP <- sapply(genes, function(x) SNPinGene(x, 1, nbSNPcausaux = NULL, listGenes = listGenesSNP))
  
  if (is.matrix(NamesSNP)) {
    NamesSNP = as.list(data.frame(NamesSNP))
  }
  
  
  nbSNPbyGene <- sapply(NamesSNP, function(x) length(x))
  nb_int <- 0
  Ztrain <- data.frame(matrix(, nrow = length(idSubs$idTrain)), row.names = rownames(X)[idSubs$idTrain])
  Ztest <- data.frame(matrix(, nrow = length(idSubs$idTest)), row.names = rownames(X)[idSubs$idTest])
  
   for (i in 1:(length(genes) - 1)) {
    X1 <- as.matrix(X[, as.character(NamesSNP[[i]])])
    for (k in (i + 1):(length(genes))) {
      X2 <- as.matrix(X[, as.character(NamesSNP[[k]])])
    print(paste(genes[i], genes[k], sep="."))
      
      dat <- data.frame(Y,X1,X2)
      
      boost.couple <- gbm(Y~.,data=dat[idSubs$idTrain,], distribution="gaussian", n.trees=2000, interaction.depth = 4)
      ztrain <- as.matrix(boost.couple$fit)
      #sum((Y[idSubs$idTrain]-ztrain)^2)
      
      #Interaction on test dataset
      ztest <- as.matrix(predict(boost.couple, newdata=dat[idSubs$idTest,], n.trees=2000))
      #sum((Y[idSubs$idTest]-ztest)^2)
      
      #par(mfrow=c(1,2)) 
      #plot(Y[idSubs$idTest], ztest)
      #plot(Y[idSubs$idTrain], ztrain)

      colnames(ztrain) <- paste("X", genes[i], genes[k], sep = ".")
      colnames(ztest) <- paste("X", genes[i], genes[k], sep = ".")
      
      Ztrain <- cbind(Ztrain, ztrain)
      Ztest <- cbind(Ztest, ztest)
      
      }
   }
  
  namesINT <- colnames(Ztrain)[-1]
  Ztrain <- data.frame(as.matrix(Ztrain[, -1]))
  Ztest <- data.frame(as.matrix(Ztest[, -1]))
  
  interLength <- rep(1, dim(Ztrain)[2])
  names(interLength) <- namesINT
  return(list(IntTrain = Ztrain, IntTest = Ztest, interLength = interLength, nb_int=nb_int))

}




NN <- function(X, Y, listGenesSNP, idSubs) {
  
  genes <- names(listGenesSNP)
  
  NamesSNP <- sapply(genes, function(x) SNPinGene(x, 1, nbSNPcausaux = NULL, listGenes = listGenesSNP))
  
  if (is.matrix(NamesSNP)) {
    NamesSNP = as.list(data.frame(NamesSNP))
  }
  
  
  nbSNPbyGene <- sapply(NamesSNP, function(x) length(x))
  nb_int <- 0
  Ztrain <- data.frame(matrix(, nrow = length(idSubs$idTrain)), row.names = rownames(X)[idSubs$idTrain])
  Ztest <- data.frame(matrix(, nrow = length(idSubs$idTest)), row.names = rownames(X)[idSubs$idTest])
  
  for (i in 1:(length(genes) - 1)) {
    X1 <- as.matrix(X[, as.character(NamesSNP[[i]])])
    for (k in (i + 1):(length(genes))) {
      X2 <- as.matrix(X[, as.character(NamesSNP[[k]])])
      print(paste(genes[i], genes[k], sep="."))
      
      
      dat <- data.frame(Y,X1,X2)
      dat <- as.data.frame(scale(dat, center = TRUE, scale = TRUE))
      n <- names(dat)
      f <- as.formula(paste("Y~", paste(n[!n %in% "Y"], collapse = " + ")))
      nn <- neuralnet(f,data=dat[idSubs$idTrain,],hidden=length(listGenesSNP),linear.output=T)
      ztrain <- as.matrix(nn$net.result[[1]])
      ztest <- as.matrix(compute(nn, dat[idSubs$idTest,-1])$net.result)
      colnames(ztrain) <- paste("X", genes[i], genes[k], sep = ".")
      colnames(ztest) <- paste("X", genes[i], genes[k], sep = ".")
      
      Ztrain <- cbind(Ztrain, ztrain)
      Ztest <- cbind(Ztest, ztest)
      
    }
  }
  
  namesINT <- colnames(Ztrain)[-1]
  Ztrain <- data.frame(as.matrix(Ztrain[, -1]))
  Ztest <- data.frame(as.matrix(Ztest[, -1]))
  
  interLength <- rep(1, dim(Ztrain)[2])
  names(interLength) <- namesINT
  return(list(IntTrain = Ztrain, IntTest = Ztest, interLength = interLength, nb_int=nb_int))
  
}




SNPinGene <- function(gene, portionSNP, nbSNPcausaux, listGenes) {
  nbSNP <- length(listGenes[[gene]])
  if(is.null(nbSNPcausaux)){
    nbCausalSNP <- nbSNP*portionSNP
    listGenes[[gene]][1:nbCausalSNP]
  }
  else{
    if(nbSNP <= nbSNPcausaux){
      listGenes[[gene]]
    }else{
      nbCausalSNP <- nbSNPcausaux
      listGenes[[gene]][1:nbCausalSNP]
    }
  }
}




IntVarCCA <- function(X1, X2, nbcomp, nameX1, nameX2, idSubs) {
  nbcomp1 = nbcomp
  nbcomp2 = nbcomp
  
  
  X1_tr <- as.matrix(X1[idSubs$idTrain,])
  X2_tr <- as.matrix(X2[idSubs$idTrain,])
  
  X1_te <- as.matrix(X1[idSubs$idTest,])
  X2_te <- as.matrix(X2[idSubs$idTest,])
  
  if (dim(X1_tr)[2] > 1 & dim(X2_tr)[2] > 1) {
    colnames(X1_tr) <- c(1:dim(X1_tr)[2])
    colnames(X2_tr) <- c(1:dim(X2_tr)[2])
    cca <- cancor(X1_tr, X2_tr)
    a1 <- cca$xcoef[, 1]
    a2 <- cca$ycoef[, 1]
    if (length(a1) < dim(X1_tr)[2]) {
      n1 <- which( !(1:dim(X1_tr)[2] %in% as.numeric(names(a1))))
      X1_tr <- X1_tr[, -n1]
      X1_te <- X1_te[, -n1]
      # suppression des variables qui n'ont pas de coef dans a1
    }
    
    if (length(a2) < dim(X2_tr)[2]) {
      n2 <- which(!(1:dim(X2_tr)[2] %in% as.numeric(names(a2))))
      X2_tr <- X2_tr[, -n2]
      X2_te <- X2_te[, -n2]
    }
  }
  
  
  if (dim(X1_tr)[2] == 1 & dim(X2_tr)[2] > 1) {
    colnames(X2_tr) <- c(1:dim(X2_tr)[2])
    cca <- cancor(X1_tr, X2_tr)
    a1 <- cca$xcoef[, 1]
    a2 <- cca$ycoef[, 1]
    
    if (length(a2) < dim(X2_tr)[2]) {
      n2 <- which(!(1:dim(X2_tr)[2] %in% as.numeric(names(a2))))
      X2_tr <- X2_tr[, -n2]
      X2_te <- X2_te[, -n2]
    }
  }
  
  if (dim(X1_tr)[2] > 1 & dim(X2_tr)[2] == 1) {
    colnames(X1_tr) <- c(1:dim(X1_tr)[2])
    cca <- cancor(X1_tr, X2_tr)
    a1 <- cca$xcoef[, 1]
    a2 <- cca$ycoef[, 1]
    if (length(a1) < dim(X1_tr)[2]) {
      n1 <- which(!(1:dim(X1_tr)[2] %in% as.numeric(names(a1))))
      X1_tr <- X1_tr[, -n1]
      X1_te <- X1_te[, -n1]
    }
  }
  
  cca2 <- cancor(X1_tr, X2_tr)
  
  
  if (dim(cca2$xcoef)[2] < nbcomp1) {
    nbcomp1 = dim(cca2$xcoef)[2]
  }
  if (dim(cca2$ycoef)[2] < nbcomp2) {
    nbcomp2 = dim(cca2$ycoef)[2]
  }
  
  A1 <- cca2$xcoef[, 1:nbcomp1]
  A2 <- cca2$ycoef[, 1:nbcomp2]
  
  I_tr <- data.frame(matrix(1, nrow = length(idSubs$idTrain)))
  I_te <- data.frame(matrix(1, nrow = length(idSubs$idTest)))
  
  for (i in seq(min(nbcomp1, nbcomp2))) {
    a1 <- as.matrix(A1)[, i]
    a2 <- as.matrix(A2)[, i]
    Z1_tr <- t(a1 %*% t(X1_tr))
    Z2_tr <- t(a2 %*% t(X2_tr))
    int_tr <- Z1_tr * Z2_tr
    int_tr <- data.frame(int_tr)
    names(int_tr) <- paste(nameX1, nameX2, i, i, sep = ".")
    I_tr <- cbind(I_tr, int_tr)
    
    Z1_te <- t(a1 %*% t(X1_te))
    Z2_te <- t(a2 %*% t(X2_te))
    int_te <- Z1_te * Z2_te
    int_te <- data.frame(int_te)
    names(int_te) <- paste(nameX1, nameX2, i, i, sep = ".")
    I_te <- cbind(I_te, int_te)
    
  }
  
  I_tr <- I_tr[, -1]
  I_te <- I_te[, -1]
  if (is.vector(I_tr)) {
    I_tr <- as.matrix(I_tr)
    colnames(I_tr) <- names(int_tr)
    rownames(I_tr) <- rownames(int_tr)
    I_te <- as.matrix(I_te)
    colnames(I_te) <- names(int_te)
    rownames(I_te) <- rownames(int_te)
  }

  nbVar <- dim(I_tr)[2]
  names(nbVar) <- paste("X", nameX1, nameX2, sep = ".")
  
  return(list(I_tr = I_tr, I_te = I_te, nbVar = nbVar))
}




CCAGenes <- function(X, G, nbcomp, idSubs) {
  if (is.null(nbcomp)) {
    nbcomp = 1
  }
  genes <- unique(G)
  nbGroup <- length(levels(as.factor(G)))
  Ztrain <- data.frame(matrix(1, nrow = length(idSubs$idTrain)))
  Ztest <- data.frame(matrix(1, nrow = length(idSubs$idTest)))
  
  d <- data.frame(matrix(1, ncol = 2))
  interLength <- c()
  for (i in 1:(nbGroup - 1)) {
    for (k in (i + 1):nbGroup) {
      X1 <- as.matrix(X[, G == genes[i]])  #1er groupe de variables
      X2 <- as.matrix(X[, G == genes[k]])  #2em groupe de variables
      print(paste(genes[i], genes[k], sep="."))
      
      resInt <- IntVarCCA(X1, X2, nbcomp, nameX1 = genes[i], nameX2 = genes[k], idSubs)
      z_tr <- resInt$I_tr
      z_te <- resInt$I_te

      Ztrain <- data.frame(Ztrain, z_tr)
      Ztest <- data.frame(Ztest, z_te)
      
      d <- rbind(d, c(i, k))
      nbVar <- resInt$nbVar
      
      interLength <- c(interLength, nbVar)
      
    }
  }
  nam <- colnames(Ztrain)[-1]
  Ztrain <- as.data.frame(Ztrain[, -1])
  Ztest <- as.data.frame(Ztest[, -1])
  colnames(Ztrain) <- nam
  colnames(Ztest) <- nam
  d <- d[-1, ]
  return(list(IntTrain = Ztrain, IntTest = Ztest, ident = d, interLength = interLength, nbVar = nbVar))
}




PCAGenes <- function(X, listGenesSNP, nbcomp, idSubs) {
  nbGenes <- length(listGenesSNP)
  namesGenes <- names(listGenesSNP)
  
  ## PC construction for each gene
  XacpTrain <- c()
  XacpTest <- c()
  allnbComp <- c()
  ResACPvar <- list()
  for (i in seq(namesGenes)) {
    red <- choixGenes(genes = namesGenes[i], X, listGenesSNP = listGenesSNP)
    Xred <- red$Xred
    
    if(dim(Xred)[2]== 1){
      CP_train <- as.matrix(Xred[idSubs$idTrain,])
      CP_test <- as.matrix(Xred[idSubs$idTest,])
      colnames(CP_train) <- paste(namesGenes[i], 1, sep =".")
      colnames(CP_test) <- paste(namesGenes[i], 1, sep =".")
      allnbComp <-c(allnbComp,1)
    }else{
      
      ResACP <- FactoMineR::PCA(Xred, scale.unit = TRUE, ncp = dim(Xred)[2], graph = F, ind.sup =idSubs$idTest)
      eig <- list(ResACP$eig)
      names(eig) <- namesGenes[i]
      ResACPvar <- c(ResACPvar, eig)
      
      if (dim(Xred)[2]< nbcomp){
        nbComp <- dim(Xred)[2]
      }else{
        nbComp <- nbcomp
      }	
      allnbComp <- c(allnbComp, nbComp)
      
      CP_train <- as.matrix(ResACP$ind$coord[, c(1:nbComp)])
      CP_test <- as.matrix(ResACP$ind.sup$coord[, c(1:nbComp)])
      
      namesVar <- c()
      for (j in seq(nbComp)) {
        namesVar <- c(namesVar, paste(namesGenes[i], j, sep = "."))
      }
      colnames(CP_train) <- namesVar
      colnames(CP_test) <- namesVar
    }
    
    XacpTrain <- cbind(XacpTrain, CP_train)
    XacpTest <- cbind(XacpTest, CP_test)
  }
  
  
  ## Interaction variables construction for each couple
  G <- rep(seq(nbGenes), allnbComp)
  IntTrain <- data.frame(matrix(1, nrow = dim(XacpTrain)[1]))
  IntTest <- data.frame(matrix(1, nrow = dim(XacpTest)[1]))
 
  interLength <- c()
  d <- data.frame(matrix(1, ncol = 2))
  c <- c()
  for (g1 in 1:(nbGenes - 1)) {
    for (g2 in (g1 + 1):nbGenes) {
      print(paste(namesGenes[g1], namesGenes[g2], sep="."))
      nbVar1 <- allnbComp[g1]
      nbVar2 <- allnbComp[g2]
      
      for (i in seq(nbVar1)) {
        for (k in seq(nbVar2)) {
          CPi_tr <- as.matrix(XacpTrain[, G == g1])[, i]
          CPk_tr <- as.matrix(XacpTrain[, G == g2])[, k]
          zTrain <- CPi_tr * CPk_tr
          zTrain <- data.frame(zTrain)
          
          CPi_te <- as.matrix(XacpTest[, G == g1])[, i]
          CPk_te <- as.matrix(XacpTest[, G == g2])[, k]
          zTest <- CPi_te * CPk_te
          zTest <- data.frame(zTest)
          
          colnames(zTrain) <- paste(namesGenes[g1], namesGenes[g2], i, k, sep = ".")
          colnames(zTest) <- paste(namesGenes[g1], namesGenes[g2], i, k, sep = ".")
          
          IntTrain <- cbind(IntTrain, zTrain)
          IntTest <- cbind(IntTest, zTest)
        }
      }
      nbVar <- nbVar1 * nbVar2
      names(nbVar) <- paste("X", namesGenes[g1], namesGenes[g2], sep = ".")
      interLength <- c(interLength, nbVar)
      d <- rbind(d, c(g1, g2))
    }
  }
  IntTrain <- IntTrain[, -1]
  IntTest <- IntTest[, -1]
  d <- d[-1, ]
  
  if (is.vector(IntTrain)) {
    IntTrain <- as.matrix(IntTrain)
    colnames(IntTrain) <- colnames(zTrain)
    
    IntTest <- as.matrix(IntTest)
    colnames(IntTest) <- colnames(zTest)
  }
  
  XCPwInt_tr <- as.matrix(cbind(XacpTrain, IntTrain))
  XCPwInt_te <- as.matrix(cbind(XacpTest, IntTest))
  
  return(list(XCPwInt_tr = XCPwInt_tr, IntTrain = IntTrain, IntTest = IntTest, interLength = interLength, ResACP = ResACPvar))
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



PLSGenes <- function(X, Y, listGenesSNP, nbcomp = NULL, idSubs) {
  nbGenes <- length(listGenesSNP)
  namesGenes <- names(listGenesSNP)
  
  IntTrain <- data.frame(matrix(1, nrow = length(idSubs$idTrain)))
  IntTest <- data.frame(matrix(1, nrow = length(idSubs$idTest)))
  
  d <- data.frame(matrix(1, ncol = 2))
  interLength <- c()
  
  for (i in 1:(nbGenes - 1)) {
    for (j in (i + 1):nbGenes) {
      print(paste(namesGenes[i], namesGenes[j], sep="."))
      
      red1 <- choixGenes(genes = namesGenes[i], X, listGenesSNP = listGenesSNP)
      g1 <- red1$Xred
      
      red2 <- choixGenes(genes = namesGenes[j], X, listGenesSNP = listGenesSNP)
      g2 <- red2$Xred
        
      dat<-data.frame(Yg1=I(cbind(Y, g1)), g2=I(g2))
      
      if (is.null(nbcomp)) {
        pls = pls::plsr(Yg1 ~ g2, data=dat[idSubs$idTrain,], validation = "LOO")
        
      } else {
        if (dim(cbind(Y, g1))[2] < nbcomp | dim(g2)[2] < nbcomp) {
          nbcomp2 = min(dim(cbind(Y, g1))[2], dim(g2))
        }else{
          nbcomp2=nbcomp
        }
        pls = pls::plsr(Yg1 ~ g2, ncomp = nbcomp2, data=dat[idSubs$idTrain,], validation = "LOO")

      }
      scPLS_tr <- pls::scores(pls)
      nbCompPLS <- dim(scPLS_tr)[2]
      scPLS_tr <- data.frame(scPLS_tr[1:length(idSubs$idTrain), 1:nbCompPLS])
    
      scPLS_te <- predict(pls, newdata=dat[idSubs$idTest,], type="scores")
      #scPLS_te <- predict(pls, newdata=g2_te, comps=1:nbcomp2, type="scores")
      #meme resultats, seul les nouveaux individus du gene2 sont pris en compte
      scPLS_te <- data.frame(scPLS_te)
      
      
      noms <- c()
      for (k in seq(nbCompPLS)) {
        noms <- c(noms, paste(namesGenes[i], namesGenes[j], k, sep = "."))
      }
      colnames(scPLS_tr) <- noms
      colnames(scPLS_te) <- noms
      names(nbCompPLS) <- paste("X", namesGenes[i], namesGenes[j], sep = ".")
      interLength <- c(interLength, nbCompPLS)
      
      IntTrain = cbind(IntTrain, scPLS_tr)
      IntTest = cbind(IntTest, scPLS_te)
      
      d <- rbind(d, c(i, j))
    }
  }
  IntTrain <- IntTrain[, -1]
  IntTest <- IntTest[, -1]
  d <- d[-1, ]
  
  return(list(IntTrain = IntTrain, IntTest = IntTest, interLength = interLength))
}

