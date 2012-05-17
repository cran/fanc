#X:data matrix
#m:number of factors
#length.rho=50: number of candidates of tuning parameters
#gamma=10000: tuning parameter in MC+
#Lambdainitial:initial value of Lambda(maybe delete??)
#eta=0.0:hyper-parameter for preventing the occurrence of improper solutions
#max.count.em=2000:Maximum number of iteration in EM algorithm
#max.count.cd=2000:Maximum number of iteration in coordinate descent algorithm 
#max.count.initial=100
#Delta=0.001: min_rho=max_rho*Delta
#min.uniquevar=0.005: if(Psi<min.uniquevar) Psi=min.uniquevar
#tol1=1e-5: tol in EM
#tol2=1e-5: tol in Coordinate descent
#scale.x = TRUE
#initial.iter=5: the number of candidates of initial values (initial values are determined randomly)
fanc <- function(x, factors, covmat, n.obs, length.rho=20, length.gamma=5, max.gamma=20, min.gamma=1.1, eta=0.0, max.count.em=2000, max.count.cd=2000, Delta=0.001, min.uniquevar=0.005, tol1=1e-5, tol2=1e-5, scale.x = TRUE, min.rho.zero=FALSE){
	#START CHECK ERRORS
	if(missing(x)==FALSE){
		if (!is.matrix(x)) stop('"x" must be a matrix.')
    	if (mode(x) != "numeric") stop('"x" must be a matrix.')
		if(factors<1 || factors>=ncol(x) ) stop('"factors"  must be a positive integer less than the number of variables.')
	}
	if(missing(covmat)==FALSE){
		if (!is.matrix(covmat)) stop('"covmat" must be a covariance matrix.')
		if(factors<1 || factors>=ncol(covmat) ) stop('"factors"  must be a positive integer less than the number of variables.')
	}
	if(length.gamma<1) stop("length.gamma must be a positive integer.")
	if(missing(n.obs)==FALSE){
		if (n.obs<2) stop('"n.obs" must be an integer greater than 2.')
	}
	if(length.rho<1) stop('"length.rho"  must be a positive integer.')
	if(length.gamma<1) stop('"length.gamma"  must be a positive integer.')
	if(min.gamma<=1) stop('"min.gamma"  must be a real value greater than 1.')
	if(eta<0) stop('"eta"  must be a non-negative real vaule.')
	if(max.count.em<1) stop('"max.count.em"  must be a positive integer.')
	if(max.count.cd<1) stop('"max.count.cd"  must be a positive integer.')
	if(Delta<=0) stop('"Delta"  must be a positive real value.')
	if(min.uniquevar<=0) stop('"min.uniquevar"  must be a positive real value.')
	if(tol1<=0) stop('"tol1"  must be a positive real value.')
	if(tol2<=0) stop('"tol2"  must be a positive real value.')
	candidate_TF <- c(TRUE, FALSE)
	if(sum(candidate_TF == scale.x) != 1) stop('"scale.x"  must be TRUE or FALSE.')
	if(sum(candidate_TF == min.rho.zero) != 1) stop('"min.rho.zero"  must be TRUE or FALSE.')
	#END CHEDKING ERROR
	
	
	#START!!
	m <- factors
	xmissing <- 0
	if(missing(x)){
		if(missing(covmat)) stop("input data matrix or covariance matrix is needed.")
		p <- ncol(covmat)
		diagS <- diag(covmat)
		if(missing(n.obs)) N <- p+1
		else N <- n.obs
		x <- x0 <- 1
		xmissing <- 1
		#x0 <- 1
		Npflag <- as.integer(1)
	}else{
		p <- ncol(x)
		N <- nrow(x) 
		if(scale.x==TRUE){
			x <- scale(x, scale=TRUE)[1:N,1:p] /sqrt(N-1)*sqrt(N)
			#x <- x %*% solve(diag(sqrt(apply(x,2,function(x) var(x)*(N-1)/N ))))
		}else{
			x <- scale(x, scale=FALSE)[1:N,1:p] 
		}
		if(N>=p){
			covmat <- cov(x) * (N-1) / N
			Npflag <- as.integer(1)
			diagS <- diag(covmat)
		}
		if(N<p) {
			covmat <- 1
			Npflag <- as.integer(0)
		diagS <- apply(x,2,function(x) var(x) * (N-1) / N)
		}
		x0 <- x[1:N,1:p]/sqrt(N)
	}
	diagm <- diag(m)
	m0 <- m0new <- 0
	AIC <-rep(0, length.rho)
	CAIC <-rep(0, length.rho)
	BIC <-rep(0, length.rho)
	Lambda_index_allfree0 <- matrix(0,p,m)
	Lambda_index_allfree <- as.integer(Lambda_index_allfree0)
	max.count.em <- as.integer(max.count.em)
	max.count.cd <- as.integer(max.count.cd)
	dimS <- as.integer(length(covmat))
	dimx0 <- as.integer(length(x0))


	#set initial values
	rhos_length=5000
	max.count.initial=100
	max.count.initial <- as.integer(max.count.initial)
	initial.iter=10
	

	#gammmamax, Lambdaの初期値決め（New version）
	Lambdainitialij_for_rhoinitial <- seq(0.3,2,length=10)
	rho_cand <- rep(0,length(Lambdainitialij_for_rhoinitial))
	penalized_loglikelihood_rho_cand <- rep(0,length(Lambdainitialij_for_rhoinitial))
	penalized_loglikelihood_Lambdainitial <- 1e+100
	rho_max <-0
	rho <- 0
	Lambdainitial00 <- matrix(0,p,m)
	diagPsiinitial00 <- rep(1,p)
	Lambdainitial_rand <- matrix(runif(p* initial.iter),p, initial.iter)
	diagPsiinitial_rand <- matrix(runif(p* initial.iter),p, initial.iter)

#	dyn.load("/Users/hirosekei/Documents/阪大（2011年度から）/研究/Factoranalysis_LASSO/EMCD/fancC_v0.3.so")
 	parameterupdate=.Call("LamdPsiupdate_initialver1_R", Lambdainitial00, diagPsiinitial00, as.integer(N), tol1, tol2, max.count.initial ,eta, min.uniquevar,rho,Inf,covmat,diagS, x0,as.integer(0), Lambda_index_allfree, Lambdainitial_rand, diagPsiinitial_rand, Lambdainitialij_for_rhoinitial, Npflag, dimx0, dimS)
 #	dyn.unload("/Users/hirosekei/Documents/阪大（2011年度から）/研究/Factoranalysis_LASSO/EMCD/fancC_v0.3.so")

	Lambda <- Lambdainitial <- parameterupdate[[2]]
	diagPsi <- diagPsiinitial <- parameterupdate[[4]]
	rho_max <- parameterupdate[[5]]
	
	
	#Psiの初期値決め2
	rho_min <- rho_max * Delta
	rho_vec <- exp(seq(log(rho_max), log(rho_min), length= length.rho))
	if(min.rho.zero==TRUE) rho_vec[length.rho] = 0.0
	
	
	
	
	
	#Compute appropriate value of rho
	equation_rho <- function(rhos, rho, gamma){
		pnorm(rhos* gamma)-gamma*pnorm(rhos)+(gamma-1)*pnorm(rho)
	}

	equation_rho1 <- function(rhos, gamma){
		pnorm(rhos* gamma)
	}
	equation_rho2 <- function(rhos){
		pnorm(rhos)
	}
	equation_rho3 <- function(rho){
		pnorm(rho)
	}


	
	determine_rho <- function(rhos,rho_vec,gamma_vec){
		ans <- matrix(NA,length(rho_vec),length(gamma_vec))
		result02 <- equation_rho2(rhos)
		result03 <- equation_rho3(rho_vec)
		for(j in 1:length(gamma_vec)){
			gamma <- gamma_vec[j]
			result01 <- equation_rho1(rhos, gamma)
			for(i in 1:length(rho_vec)){
				rho <- rho_vec[i]
				#result0 <- equation_rho(rhos, rho, gamma)
				result0 <- result01 - gamma*result02 + (gamma-1)*result03[i]
				ans[i,j] <- rhos[which.min(abs(result0))]
			}
		}
		ans
	}
	if(length.gamma<2){
		rhomatrix <- as.matrix(rho_vec)
		gamma_vec=Inf
	}else if(length.gamma<3){
		rhos <- exp(seq(log(10*rho_max), log(rho_min/2), length= rhos_length))
		gammmamatrix0 <- determine_rho(rhos,rho_vec,min.gamma)	
		rhomatrix <- cbind(rho_vec, gammmamatrix0)
		gamma_vec=c(Inf, min.gamma)
	}else{
		gamma_vec <- exp(seq(from=(log(max.gamma)),to=(log(min.gamma)),length.out=(length.gamma-1)))
		rhos <- exp(seq(log(10*rho_max), log(rho_min/2), length= rhos_length))
		gammmamatrix0 <- determine_rho(rhos,rho_vec,gamma_vec)	
		rhomatrix <- cbind(rho_vec, gammmamatrix0)
		gamma_vec <- c(Inf,gamma_vec)
	}
	
	colnames(rhomatrix) <- NULL

	
	
	#####main#####
#	dyn.load("/Users/hirosekei/Documents/阪大（2011年度から）/研究/Factoranalysis_LASSO/EMCD/fancC_v0.3.so")
 	result =.Call("fancmain", Lambda, diagPsi, as.integer(N), tol1, tol2, max.count.em, max.count.cd , max.count.initial, eta, min.uniquevar, rhomatrix, as.matrix(gamma_vec), x0, diagS, covmat, Lambdainitial_rand, Npflag, dimx0, dimS)
#	dyn.unload("/Users/hirosekei/Documents/阪大（2011年度から）/研究/Factoranalysis_LASSO/EMCD/fancC_v0.3.so")



	#結果整理

     LambdaAIC_rhofixed <- array(0,dim=c(p,m,ncol(rhomatrix)))
	 LambdaBIC_rhofixed <- array(0,dim=c(p,m,ncol(rhomatrix)))
	 LambdaCAIC_rhofixed <- array(0,dim=c(p,m,ncol(rhomatrix)))
	 diagPsiAIC_rhofixed <- array(0,dim=c(ncol(rhomatrix),p))
	 diagPsiBIC_rhofixed <- array(0,dim=c(ncol(rhomatrix),p))
	 diagPsiCAIC_rhofixed <- array(0,dim=c(ncol(rhomatrix),p))
	 dfAIC_rhofixed <- dfBIC_rhofixed <- dfCAIC_rhofixed <- rep(0,ncol(rhomatrix))
	 AICvalue_rhofixed <- BICvalue_rhofixed <- CAICvalue_rhofixed <- rep(0,ncol(rhomatrix))
	 rhoAIC_rhofixed <- rhoBIC_rhofixed <- rhoCAIC_rhofixed <- rep(0,ncol(rhomatrix))
	 gammaAIC_rhofixed <- gammaBIC_rhofixed <- gammaCAIC_rhofixed <- rep(0,ncol(rhomatrix))
	
	Lambda_all <- array(result[[1]],dim=c(p,m,nrow(rhomatrix),ncol(rhomatrix)))
	diagPsi_all <- array(t(result[[2]]),dim=c(nrow(rhomatrix),p, ncol(rhomatrix)))
	diagPsi_all0 <- t(result[[2]])
	for(i in 1:ncol(rhomatrix)){
		diagPsi_all[,,i] <- diagPsi_all0[(i-1)*nrow(rhomatrix) + (1:nrow(rhomatrix)),]
	}
	AIC_all <- matrix(result[[4]], nrow(rhomatrix), ncol(rhomatrix))
	BIC_all <- matrix(result[[5]], nrow(rhomatrix), ncol(rhomatrix))
	CAIC_all <- matrix(result[[6]], nrow(rhomatrix), ncol(rhomatrix))
	dfmatrix <- matrix(result[[7]],nrow(rhomatrix), ncol(rhomatrix))
	GFI <- matrix(result[[8]],nrow(rhomatrix), ncol(rhomatrix))
	AGFI <- AGFI <- 1 - (p*(p+1)*(1-GFI)) / (p*(p+1)-2*dfmatrix)
	

	if(length.gamma == 1){
		AICindexr <- which.min(AIC_all); BICindexr <- which.min(BIC_all); CAICindexr <- which.min(CAIC_all); 
		LambdaAIC <- Lambda_all[,,AICindexr,1]; LambdaBIC <- Lambda_all[,,BICindexr,1]; LambdaCAIC <- Lambda_all[,,CAICindexr,1]; 
		diagPsiAIC <- diagPsi_all[AICindexr,,1]; diagPsiBIC <- diagPsi_all[BICindexr,,1]; diagPsiCAIC <- diagPsi_all[CAICindexr,,1]; 
		rhoAIC <- rhomatrix[AICindexr,1]; rhoBIC <- rhomatrix[BICindexr,1]; rhoCAIC <- rhomatrix[CAICindexr,1]; 
		gammaAIC <- gamma_vec[1]; gammaBIC <- gamma_vec[1]; gammaCAIC <- gamma_vec[1]; 
		AICvalue <- min(AIC_all); BICvalue <- min(BIC_all); CAICvalue <- min(CAIC_all) 
		dfAIC <- dfmatrix[AICindexr,1]; dfBIC <- dfmatrix[BICindexr,1]; dfCAIC <- dfmatrix[CAICindexr,1]; 
	}else{
		AICindexr0 <- apply(AIC_all,2,which.min); AICindexc <- which.min(diag(AIC_all[AICindexr0,1: length.gamma])); AICindexr <- AICindexr0[AICindexc]
		BICindexr0 <- apply(BIC_all,2,which.min); BICindexc <- which.min(diag(BIC_all[BICindexr0,1: length.gamma])); BICindexr <- BICindexr0[BICindexc]
		CAICindexr0 <- apply(CAIC_all,2,which.min); CAICindexc <- which.min(diag(CAIC_all[CAICindexr0,1: length.gamma])); CAICindexr <- CAICindexr0[CAICindexc]
		LambdaAIC <- Lambda_all[,,AICindexr,AICindexc]; LambdaBIC <- Lambda_all[,,BICindexr,BICindexc]; LambdaCAIC <- Lambda_all[,,CAICindexr,CAICindexc]; 
		diagPsiAIC <- diagPsi_all[AICindexr,,AICindexc]; diagPsiBIC <- diagPsi_all[BICindexr,,BICindexc]; diagPsiCAIC <- diagPsi_all[CAICindexr,,CAICindexc];
		rhoAIC <- rhomatrix[AICindexr,AICindexc]; rhoBIC <- rhomatrix[BICindexr,BICindexc]; rhoCAIC <- rhomatrix[CAICindexr,CAICindexc]; 
		gammaAIC <- gamma_vec[AICindexc]; gammaBIC <- gamma_vec[BICindexc]; gammaCAIC <- gamma_vec[CAICindexc]; 
		AICvalue <- min(AIC_all); BICvalue <- min(BIC_all); CAICvalue <- min(CAIC_all) 
		dfAIC <- dfmatrix[AICindexr,AICindexc]; dfBIC <- dfmatrix[BICindexr,BICindexc]; dfCAIC <- dfmatrix[CAICindexr,CAICindexc]; 
		 for(i in 1:length.gamma){
			 LambdaAIC_rhofixed[,,i] <- Lambda_all[,, AICindexr0[i],i];
			 LambdaBIC_rhofixed[,,i] <- Lambda_all[,, BICindexr0[i],i];
			 LambdaCAIC_rhofixed[,,i] <- Lambda_all[,, CAICindexr0[i],i];
			 diagPsiAIC_rhofixed[i,] <- diagPsi_all[AICindexr0[i],,i];
			 diagPsiBIC_rhofixed[i,] <- diagPsi_all[BICindexr0[i],,i];
		     diagPsiCAIC_rhofixed[i,] <- diagPsi_all[CAICindexr0[i],,i];
			dfAIC_rhofixed[i] <- dfmatrix[AICindexr0[i],i]
			dfBIC_rhofixed[i] <- dfmatrix[BICindexr0[i],i]
			dfCAIC_rhofixed[i] <- dfmatrix[CAICindexr0[i],i]
			AICvalue_rhofixed[i] <- AIC_all[AICindexr0[i],i]
			BICvalue_rhofixed[i] <- BIC_all[BICindexr0[i],i]
			CAICvalue_rhofixed[i] <- CAIC_all[CAICindexr0[i],i]
			rhoAIC_rhofixed[i] <- rhomatrix[AICindexr0[i],i]
			rhoBIC_rhofixed[i] <- rhomatrix[BICindexr0[i],i]
			rhoCAIC_rhofixed[i] <- rhomatrix[CAICindexr0[i],i]
			gammaAIC_rhofixed[i] <- gamma_vec[i]
			gammaBIC_rhofixed[i] <- gamma_vec[i]
			gammaCAIC_rhofixed[i] <- gamma_vec[i]
		 }
	}
	
	outAIC <- list(loadings=LambdaAIC, uniquenesses=diagPsiAIC, df=dfAIC, criteria=AICvalue, rho=rhoAIC, gamma=gammaAIC)
	outBIC <- list(loadings=LambdaBIC, uniquenesses=diagPsiBIC, df=dfBIC, criteria=BICvalue, rho=rhoBIC, gamma=gammaBIC)
	outCAIC <- list(loadings=LambdaCAIC, uniquenesses=diagPsiCAIC, df=dfCAIC, criteria=CAICvalue, rho=rhoCAIC, gamma=gammaCAIC)

	outAIC_rhofixed <- list(loadings=LambdaAIC_rhofixed, uniquenesses=diagPsiAIC_rhofixed, df=dfAIC_rhofixed, criteria=AICvalue_rhofixed, rho=rhoAIC_rhofixed, gamma=gammaAIC_rhofixed)
	outBIC_rhofixed <- list(loadings=LambdaBIC_rhofixed, uniquenesses=diagPsiBIC_rhofixed, df=dfBIC_rhofixed, criteria=BICvalue_rhofixed, rho=rhoBIC_rhofixed, gamma=gammaBIC_rhofixed)
	outCAIC_rhofixed <- list(loadings=LambdaCAIC_rhofixed, uniquenesses=diagPsiCAIC_rhofixed, df=dfCAIC_rhofixed, criteria=CAICvalue_rhofixed, rho=rhoCAIC_rhofixed, gamma=gammaCAIC_rhofixed)
	
	#gamma_mat <- matrix(gamma_vec,nrow(rhomatrix),ncol(rhomatrix),byrow=T)
		

	#if(xmissing==1) ans <-list(loadings=Lambda_all,uniquenesses=diagPsi_all,rho=rhomatrix, AIC=AIC_all, BIC=BIC_all, CAIC=CAIC_all, outAIC=outAIC, outBIC=outBIC, outCAIC=outCAIC, dfmatrix= dfmatrix, GFI=GFI, AGFI=AGFI, gamma=gamma_vec, factors=factors, Npflag=Npflag, call=match.call())

	if(xmissing==1) ans <-list(loadings=Lambda_all,uniquenesses=diagPsi_all,rho=rhomatrix, AIC=AIC_all, BIC=BIC_all, CAIC=CAIC_all, outAIC=outAIC_rhofixed, outBIC=outBIC_rhofixed, outCAIC=outCAIC_rhofixed, dfmatrix= dfmatrix, GFI=GFI, AGFI=AGFI, gamma=gamma_vec, factors=factors, Npflag=Npflag, call=match.call())



	#if(xmissing==0) ans <-list(loadings=Lambda_all,uniquenesses=diagPsi_all,rho=rhomatrix, AIC=AIC_all, BIC=BIC_all, CAIC=CAIC_all, outAIC=outAIC, outBIC=outBIC, outCAIC=outCAIC, dfmatrix= dfmatrix, GFI=GFI, AGFI=AGFI, gamma=gamma_vec, factors=factors, Npflag=Npflag, call=match.call(), x=x)

	if(xmissing==0) ans <-list(loadings=Lambda_all,uniquenesses=diagPsi_all,rho=rhomatrix, AIC=AIC_all, BIC=BIC_all, CAIC=CAIC_all, outAIC=outAIC_rhofixed, outBIC=outBIC_rhofixed, outCAIC=outCAIC_rhofixed, dfmatrix= dfmatrix, GFI=GFI, AGFI=AGFI, gamma=gamma_vec, factors=factors, Npflag=Npflag, call=match.call(), x=x)

	class(ans) <- "fanc"
	ans


}
