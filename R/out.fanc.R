out <- function(x,rho, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
 	gamma_vec <- x$gamma
 	gamma_length <- length(gamma_vec)
 	if(gamma==Inf) gamma_index <- 1
 	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
 	rhovec <- x$rho[,gamma_index]
 	rho_index <- which.min(abs(rho-rhovec))	
 	Lambda <- x$loadings[,,rho_index,gamma_index]
	Lambda <- as.matrix(Lambda)
 	diagPsi <- x$uniquenesses[rho_index,,gamma_index]
 	if(scores==TRUE && is.null(x$x)==TRUE) stop("data matrix is needed for computing the factor score when fitting the data by fanc")
 	if(scores==TRUE){
 		diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
 		diagPsiinvLambda <- diagPsiinvrep * Lambda
 		M0 <- crossprod(Lambda,diagPsiinvLambda)
 		M <- M0 + diag(ncol(Lambda))
 		solveM <- solve(M)
 		PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
 		ans_scores <- x$x %*% PsiinvLambdaMinv
 	}
 	df <- x$df[rho_index,gamma_index]
 	AIC <- x$AIC[rho_index,gamma_index]
 	BIC <- x$BIC[rho_index,gamma_index]
 	CAIC <- x$CAIC[rho_index,gamma_index]
 	criteria <- c(AIC,BIC,CAIC)
 	names(criteria) <- c("AIC","BIC","CAIC")
 	
 	gamma0 <- gamma_vec[gamma_index]
 	rho0 <- rhovec[rho_index]
 	
 	class(Lambda) <- "loadings"
 	#calculate GFI
 	# if(GFI==TRUE){
 		# Sigma <- Lambda %*% t(Lambda) + diag(diagPsi)
 		# V <- solve(Sigma)
 		# S <- x$covmat
 		# GFI <- 1 - sum(diag( (V %*% (S-Sigma)) %*% (V %*% (S-Sigma)) )) / sum(diag(V %*% S %*% V %*% S))
 		# p <- nrow(Lambda)
 		# m <- ncol(Lambda)
 		# AGFI <- 1 - (p*(p+1)*(1-GFI)) / (p*(p+1)-2*df)
 		# if(scores==TRUE){
 			# ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, criteria=criteria, rho=rho0, gamma=gamma0, GFI=GFI, AGFI=AGFI)
 		# }else{
 			# ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, criteria=criteria, rho=rho0, gamma=gamma0, GFI=GFI, AGFI=AGFI)
 		# }
	p0 <- dim(x$loadings)[1]
  	if(p0<=500){
	 	GFI <- x$GFI[rho_index,gamma_index];
	 	AGFI <- x$AGFI[rho_index,gamma_index];
	 	GOF <- c(GFI,AGFI)
	 	names(GOF) <- c("GFI","AGFI")
		if(scores==TRUE){
 			ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, criteria=criteria, goodness.of.fit=GOF, rho=rho0, gamma=gamma0)
 		}else{
 			ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, criteria=criteria,  goodness.of.fit=GOF,rho=rho0, gamma=gamma0)
 		}
 	}else{
		if(scores==TRUE){
 			ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, criteria=criteria, rho=rho0, gamma=gamma0)
 		}else{
 			ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, criteria=criteria, rho=rho0, gamma=gamma0)
 		}
 	}
 	ans
}


out.aic <- function(x, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(gamma<=1) stop("gamma must be greater than 1")
 	gamma_vec <- x$gamma
 	gamma_length <- length(gamma_vec)
 	if(gamma==Inf) gamma_index <- 1
 	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
 	Lambda <- x$outAIC$loadings[,,gamma_index]
	Lambda <- as.matrix(Lambda)
 	diagPsi <- x$outAIC$uniquenesses[gamma_index,]
 	if(scores==TRUE && is.null(x$x)==TRUE) stop("data matrix is needed for computing the factor score when fitting the data by fanc")
 	if(scores==TRUE){
 		diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
 		diagPsiinvLambda <- diagPsiinvrep * Lambda
 		M0 <- crossprod(Lambda,diagPsiinvLambda)
 		M <- M0 + diag(ncol(Lambda))
 		solveM <- solve(M)
 		PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
 		ans_scores <- x$x %*% PsiinvLambdaMinv
 	}
 	df <- x$outAIC$df[gamma_index]
 	AIC <- x$outAIC$criteria[gamma_index]
 	gamma0 <- x$outAIC$gamma[gamma_index]
 	rho0 <- x$outAIC$rho[gamma_index]

 	class(Lambda) <- "loadings"
 	if(scores==TRUE){
 		ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, AIC=AIC, rho=rho0, gamma=gamma0)
 		#ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, AIC=AIC, rho=rho0)
 	}else{
 		ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, AIC=AIC, rho=rho0, gamma=gamma0)
 		#ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, AIC=AIC, rho=rho0)
 	}
 	ans
}


out.bic <- function(x, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(gamma<=1) stop("gamma must be greater than 1")
 	gamma_vec <- x$gamma
 	gamma_length <- length(gamma_vec)
 	if(gamma==Inf) gamma_index <- 1
 	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
 	Lambda <- x$outBIC$loadings[,,gamma_index]
	Lambda <- as.matrix(Lambda)
 	diagPsi <- x$outBIC$uniquenesses[gamma_index,]
 	if(scores==TRUE && is.null(x$x)==TRUE) stop("data matrix is needed for computing the factor score when fitting the data by fanc")
 	if(scores==TRUE){
 		diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
 		diagPsiinvLambda <- diagPsiinvrep * Lambda
 		M0 <- crossprod(Lambda,diagPsiinvLambda)
 		M <- M0 + diag(ncol(Lambda))
 		solveM <- solve(M)
 		PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
 		ans_scores <- x$x %*% PsiinvLambdaMinv
 	}
 	df <- x$outBIC$df[gamma_index]
 	BIC <- x$outBIC$criteria[gamma_index]
 	gamma0 <- x$outBIC$gamma[gamma_index]
 	rho0 <- x$outBIC$rho[gamma_index]

 	class(Lambda) <- "loadings"
 	if(scores==TRUE){
 		ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, BIC=BIC, rho=rho0, gamma=gamma0)
 		#ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, BIC=BIC, rho=rho0)
 	}else{
 		ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, BIC=BIC, rho=rho0, gamma=gamma0)
 		#ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, BIC=BIC, rho=rho0)
 	}
 	ans
}


out.caic <- function(x, gamma, scores=FALSE){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(gamma<=1) stop("gamma must be greater than 1")
 	gamma_vec <- x$gamma
 	gamma_length <- length(gamma_vec)
 	if(gamma==Inf) gamma_index <- 1
 	if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
 	Lambda <- x$outCAIC$loadings[,,gamma_index]
	Lambda <- as.matrix(Lambda)
 	diagPsi <- x$outCAIC$uniquenesses[gamma_index,]
 	if(scores==TRUE && is.null(x$x)==TRUE) stop("data matrix is needed for computing the factor score when fitting the data by fanc")
 	if(scores==TRUE){
 		diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
 		diagPsiinvLambda <- diagPsiinvrep * Lambda
 		M0 <- crossprod(Lambda,diagPsiinvLambda)
 		M <- M0 + diag(ncol(Lambda))
 		solveM <- solve(M)
 		PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
 		ans_scores <- x$x %*% PsiinvLambdaMinv
 	}
 	df <- x$outCAIC$df[gamma_index]
 	CAIC <- x$outCAIC$criteria[gamma_index]
 	gamma0 <- x$outCAIC$gamma[gamma_index]
 	rho0 <- x$outCAIC$rho[gamma_index]

 	class(Lambda) <- "loadings"
 	if(scores==TRUE){
 		ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, CAIC=CAIC, rho=rho0, gamma=gamma0)
 		#ans <- list(loadings=Lambda, uniquenesses=diagPsi, scores=ans_scores, df=df, CAIC=CAIC, rho=rho0)
 	}else{
 		ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, CAIC=CAIC, rho=rho0, gamma=gamma0)
 		#ans <- list(loadings=Lambda, uniquenesses=diagPsi, df=df, CAIC=CAIC, rho=rho0)
 	}
 	ans
}
