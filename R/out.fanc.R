out <- function(x,rho, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	 if(scores==TRUE && is.null(x$x)==TRUE) stop("Data matrix is needed for computing the factor score in fitting procedure by fanc")
	 
	 gamma_vec <- x$gamma
	 gamma_length <- length(gamma_vec)
	 if(gamma==Inf) gamma_index <- 1
	 if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
	 rhovec <- x$rho[,gamma_index]
	 rho_index <- which.min(abs(rho-rhovec))	
	 Lambda <- x$loadings[[gamma_index]][[rho_index]]
	 diagPsi <- x$uniquenesses[rho_index,,gamma_index]
	 if(x$cor.factor==TRUE){
		Phi <- x$Phi[,,rho_index,gamma_index]
		Phi <- as.matrix(Phi)
	 }
	 if(scores==TRUE){
		Lambda_mat <- as.matrix(Lambda)
		 diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
		 diagPsiinvLambda <- diagPsiinvrep * Lambda_mat
		 M0 <- crossprod(Lambda_mat,diagPsiinvLambda)
		 if(x$cor.factor==TRUE) M <- M0 + solve(Phi)
		 if(x$cor.factor==FALSE) M <- M0 + diag(x$factors)
		 solveM <- solve(M)
		 PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
		 ans_scores <- x$x %*% PsiinvLambdaMinv
	 }
	 df <- x$df[rho_index,gamma_index]
	 if(is.null(x$AIC)==FALSE){
		AIC <- x$AIC[rho_index,gamma_index]
		 BIC <- x$BIC[rho_index,gamma_index]
		 CAIC <- x$CAIC[rho_index,gamma_index]
		 criteria <- c(AIC,BIC,CAIC)
		 names(criteria) <- c("AIC","BIC","CAIC")
	 }
	 
	 gamma0 <- gamma_vec[gamma_index]
	 rho0 <- rhovec[rho_index]
	 
	 
	if(is.null(x$GFI)==FALSE){
		 GFI <- x$GFI[rho_index,gamma_index];
		 AGFI <- x$AGFI[rho_index,gamma_index];
		 GOF <- c(GFI,AGFI)
		 names(GOF) <- c("GFI","AGFI")
	 }
	 
	 ans <- list(loadings=Lambda, uniquenesses=diagPsi)
	 if(x$cor.factor==TRUE) ans <- append(ans,list(Phi=Phi))
	 if(scores==TRUE) ans <- append(ans,list(scores=ans_scores)) 
	 ans <- append(ans,list(df=df)) 
	 if(is.null(x$AIC)==FALSE) ans <- append(ans,list(criteria=criteria))
	 if(is.null(x$GFI)==FALSE) ans <- append(ans,list(goodness.of.fit=GOF))
	 ans <- append(ans,list(rho=rho0, gamma=gamma0)) 
	 ans
 }
 
 
 out.aic <- function(x, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(gamma<=1) stop("gamma must be greater than 1")
	if(scores==TRUE && is.null(x$x)==TRUE) stop("Data matrix is needed for computing the factor score in fitting procedure by fanc")
	if(is.null(x$AIC)==TRUE) stop("AIC was not able to be calculated. Data matrix or the number of observations is needed in fitting procedure by fanc.")
	gamma_vec <- x$gamma
	gamma_length <- length(gamma_vec)
	if(gamma==Inf) gamma_index <- 1
	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
	 
	if(gamma_length == 1) AIC_vec=c(x$AIC)
	else AIC_vec=x$AIC[,gamma_index]
	
	rho_index <- which.min(AIC_vec)
	Lambda <- x$loadings[[gamma_index]][[rho_index]]
	diagPsi <- x$uniquenesses[rho_index,,gamma_index]
	 if(x$cor.factor==TRUE){
		Phi <- x$Phi[,,rho_index,gamma_index]
		Phi <- as.matrix(Phi)
	 }
	rho0 <- x$rho[rho_index,gamma_index]
	gamma0 <- gamma_vec[gamma_index]
	AIC <- min(AIC_vec)
	df <- x$df[rho_index,gamma_index]
	
	 if(scores==TRUE){
		Lambda_mat <- as.matrix(Lambda)
		 diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
		 diagPsiinvLambda <- diagPsiinvrep * Lambda_mat
		 M0 <- crossprod(Lambda_mat,diagPsiinvLambda)
		 if(x$cor.factor==TRUE) M <- M0 + solve(Phi)
		 if(x$cor.factor==FALSE) M <- M0 + diag(x$factors)
		 solveM <- solve(M)
		 PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
		 ans_scores <- x$x %*% PsiinvLambdaMinv
	 }
	 
	if(is.null(x$GFI)==FALSE){
		 GFI <- x$GFI[rho_index,gamma_index];
		 AGFI <- x$AGFI[rho_index,gamma_index];
		 GOF <- c(GFI,AGFI)
		 names(GOF) <- c("GFI","AGFI")
	 }
	 
	 	 
	 ans <- list(loadings=Lambda, uniquenesses=diagPsi)
	 if(x$cor.factor==TRUE) ans <- append(ans,list(Phi=Phi))
	 if(scores==TRUE) ans <- append(ans,list(scores=ans_scores)) 
	 ans <- append(ans,list(df=df, AIC=AIC))
	 if(is.null(x$GFI)==FALSE) ans <- append(ans,list(goodness.of.fit=GOF))
	 ans <- append(ans,list(rho=rho0, gamma=gamma0)) 
	 

	 ans
 }
 
 
 out.bic <- function(x, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(gamma<=1) stop("gamma must be greater than 1")
	if(scores==TRUE && is.null(x$x)==TRUE) stop("Data matrix is needed for computing the factor score in fitting procedure by fanc")
	if(is.null(x$BIC)==TRUE) stop("BIC was not able to be calculated. Data matrix or the number of observations is needed in fitting procedure by fanc.")
	gamma_vec <- x$gamma
	gamma_length <- length(gamma_vec)
	if(gamma==Inf) gamma_index <- 1
	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
	 
	if(gamma_length == 1) BIC_vec=c(x$BIC)
	else BIC_vec=x$BIC[,gamma_index]
	
	rho_index <- which.min(BIC_vec)
	Lambda <- x$loadings[[gamma_index]][[rho_index]]
	diagPsi <- x$uniquenesses[rho_index,,gamma_index]
	 if(x$cor.factor==TRUE){
		Phi <- x$Phi[,,rho_index,gamma_index]
		Phi <- as.matrix(Phi)
	 }
	rho0 <- x$rho[rho_index,gamma_index]
	gamma0 <- gamma_vec[gamma_index]
	BIC <- min(BIC_vec)
	df <- x$df[rho_index,gamma_index]
	
	 if(scores==TRUE){
		Lambda_mat <- as.matrix(Lambda)
		 diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
		 diagPsiinvLambda <- diagPsiinvrep * Lambda_mat
		 M0 <- crossprod(Lambda_mat,diagPsiinvLambda)
		 if(x$cor.factor==TRUE) M <- M0 + solve(Phi)
		 if(x$cor.factor==FALSE) M <- M0 + diag(x$factors)
		 solveM <- solve(M)
		 PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
		 ans_scores <- x$x %*% PsiinvLambdaMinv
	 }
	 
	if(is.null(x$GFI)==FALSE){
		 GFI <- x$GFI[rho_index,gamma_index];
		 AGFI <- x$AGFI[rho_index,gamma_index];
		 GOF <- c(GFI,AGFI)
		 names(GOF) <- c("GFI","AGFI")
	 }
	 
	 	 
	 ans <- list(loadings=Lambda, uniquenesses=diagPsi)
	 if(x$cor.factor==TRUE) ans <- append(ans,list(Phi=Phi))
	 if(scores==TRUE) ans <- append(ans,list(scores=ans_scores)) 
	 ans <- append(ans,list(df=df, BIC=BIC))
	 if(is.null(x$GFI)==FALSE) ans <- append(ans,list(goodness.of.fit=GOF))
	 ans <- append(ans,list(rho=rho0, gamma=gamma0)) 
	 
	 ans
}
 
 
 out.caic <- function(x, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(gamma<=1) stop("gamma must be greater than 1")
	if(scores==TRUE && is.null(x$x)==TRUE) stop("Data matrix is needed for computing the factor score in fitting procedure by fanc")
	if(is.null(x$CAIC)==TRUE) stop("CAIC was not able to be calculated. Data matrix or the number of observations is needed in fitting procedure by fanc.")
	gamma_vec <- x$gamma
	gamma_length <- length(gamma_vec)
	if(gamma==Inf) gamma_index <- 1
	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
	 
	if(gamma_length == 1) CAIC_vec=c(x$CAIC)
	else CAIC_vec=x$CAIC[,gamma_index]
	
	rho_index <- which.min(CAIC_vec)
	Lambda <- x$loadings[[gamma_index]][[rho_index]]
	diagPsi <- x$uniquenesses[rho_index,,gamma_index]
	 if(x$cor.factor==TRUE){
		Phi <- x$Phi[,,rho_index,gamma_index]
		Phi <- as.matrix(Phi)
	 }
	rho0 <- x$rho[rho_index,gamma_index]
	gamma0 <- gamma_vec[gamma_index]
	CAIC <- min(CAIC_vec)
	df <- x$df[rho_index,gamma_index]
	
	 if(scores==TRUE){
		Lambda_mat <- as.matrix(Lambda)
		 diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
		 diagPsiinvLambda <- diagPsiinvrep * Lambda_mat
		 M0 <- crossprod(Lambda_mat,diagPsiinvLambda)
		 if(x$cor.factor==TRUE) M <- M0 + solve(Phi)
		 if(x$cor.factor==FALSE) M <- M0 + diag(x$factors)
		 solveM <- solve(M)
		 PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
		 ans_scores <- x$x %*% PsiinvLambdaMinv
	 }
	 
	if(is.null(x$GFI)==FALSE){
		 GFI <- x$GFI[rho_index,gamma_index];
		 AGFI <- x$AGFI[rho_index,gamma_index];
		 GOF <- c(GFI,AGFI)
		 names(GOF) <- c("GFI","AGFI")
	 }
	 
	 	 
	 ans <- list(loadings=Lambda, uniquenesses=diagPsi)
	 if(x$cor.factor==TRUE) ans <- append(ans,list(Phi=Phi))
	 if(scores==TRUE) ans <- append(ans,list(scores=ans_scores)) 
	 ans <- append(ans,list(df=df, CAIC=CAIC))
	 if(is.null(x$GFI)==FALSE) ans <- append(ans,list(goodness.of.fit=GOF))
	 ans <- append(ans,list(rho=rho0, gamma=gamma0)) 
	 
	 ans
 }


 out.ebic <- function(x, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(gamma<=1) stop("gamma must be greater than 1")
	if(scores==TRUE && is.null(x$x)==TRUE) stop("Data matrix is needed for computing the factor score in fitting procedure by fanc")
	if(is.null(x$EBIC)==TRUE) stop("Extended BIC was not able to be calculated. Data matrix or the number of observations is needed in fitting procedure by fanc.")
	gamma_vec <- x$gamma
	gamma_length <- length(gamma_vec)
	if(gamma==Inf) gamma_index <- 1
	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
	 
	if(gamma_length == 1) EBIC_vec=c(x$EBIC)
	else EBIC_vec=x$EBIC[,gamma_index]
	
	rho_index <- which.min(EBIC_vec)
	Lambda <- x$loadings[[gamma_index]][[rho_index]]
	diagPsi <- x$uniquenesses[rho_index,,gamma_index]
	 if(x$cor.factor==TRUE){
		Phi <- x$Phi[,,rho_index,gamma_index]
		Phi <- as.matrix(Phi)
	 }
	rho0 <- x$rho[rho_index,gamma_index]
	gamma0 <- gamma_vec[gamma_index]
	EBIC <- min(EBIC_vec)
	df <- x$df[rho_index,gamma_index]
	
	 if(scores==TRUE){
		Lambda_mat <- as.matrix(Lambda)
		 diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
		 diagPsiinvLambda <- diagPsiinvrep * Lambda_mat
		 M0 <- crossprod(Lambda_mat,diagPsiinvLambda)
		 if(x$cor.factor==TRUE) M <- M0 + solve(Phi)
		 if(x$cor.factor==FALSE) M <- M0 + diag(x$factors)
		 solveM <- solve(M)
		 PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
		 ans_scores <- x$x %*% PsiinvLambdaMinv
	 }
	 
	if(is.null(x$GFI)==FALSE){
		 GFI <- x$GFI[rho_index,gamma_index];
		 AGFI <- x$AGFI[rho_index,gamma_index];
		 GOF <- c(GFI,AGFI)
		 names(GOF) <- c("GFI","AGFI")
	 }
	 
	 	 
	 ans <- list(loadings=Lambda, uniquenesses=diagPsi)
	 if(x$cor.factor==TRUE) ans <- append(ans,list(Phi=Phi))
	 if(scores==TRUE) ans <- append(ans,list(scores=ans_scores)) 
	 ans <- append(ans,list(df=df, EBIC=EBIC))
	 if(is.null(x$GFI)==FALSE) ans <- append(ans,list(goodness.of.fit=GOF))
	 ans <- append(ans,list(rho=rho0, gamma=gamma0)) 
	 
	 ans
}


