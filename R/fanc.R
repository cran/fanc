fanc <- function(x, factors, n.obs, cor.factor=FALSE, normalize = TRUE, rho.max, covmat, control=list()){
	#ORGANAIZE "CONTROL"
    con <- list(length.rho=30, length.gamma=9, max.gamma=100, min.gamma=1.01, eta=0.0, ncand.initial=10,  maxit.initial=500,  maxit.em=10000, maxit.cd=500, maxit.bfgs=500, start="warm", Delta=1e-5, min.uniquevar=0.005, tol.em=1e-5, tol.cd=1e-5, tol.bfgs=1e-5, min.rhozero=FALSE, zita=0.0, progress=FALSE, openmp=FALSE, num.threads=0, gamma.ebic=1)
    con[names(control)] <- control
    
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
	if(con$length.gamma<1) stop('"length.gamma" in control must be a positive integer.')
	if(missing(n.obs)==FALSE){
		if (n.obs<2) stop('"n.obs" must be an integer greater than 2.')
	}
	if(class(cor.factor) != "logical") stop('"cor.factor"  must be logical.')
	if(con$length.rho<1) stop('"length.rho" in control must be a positive integer.')
	if(missing(rho.max)==FALSE){
		if(rho.max<=0) stop('"rho.max"  must be a positive real value.')
		if(rho.max>8.25) warning('"rho.max" is greater than 8.25. In such cases, the reparametrization of the penalty funcion may be failed')
	}
	if(con$min.gamma<=1) stop('"min.gamma" in control must be a real value greater than 1.')
	if(con$max.gamma<=con$min.gamma) stop('"max.gamma"  must be a real value greater than min.gamma.')
	if(con$eta<0) stop('"eta"  must be a non-negative real vaule.')
	if(con$ncand.initial<1) stop('"ncand.initial" in control must be a positive integer.')
	if(con$maxit.em<1) stop('"maxit.em" in control must be a positive integer.')
	if(con$maxit.cd<1) stop('"maxit.cd" in control must be a positive integer.')
	if(con$maxit.bfgs<1) stop('"maxit.bfgs" in control must be a positive integer.')
    cand <- c("warm", "cold")
	if(sum(con$start==cand) != 1) stop('"start" in control must be "warm" or "cold".')
	if(con$Delta<=0) stop('"Delta" in control must be a positive real value.')
	if(con$min.uniquevar<=0) stop('"min.uniquevar" in control must be a positive real value.')
	if(con$tol.em<=0) stop('"tol.em" in control must be a positive real value.')
	if(con$tol.cd<=0) stop('"tol.cd" in control must be a positive real value.')
	if(con$tol.bfgs<=0) stop('"tol.bfgs" in control must be a positive real value.')
	if(con$zita<0) stop('"zita" in control must be a non-negative real vaule.')
	if(class(normalize) != "logical") stop('"normalize" must be logical.')
	if(class(con$min.rhozero) != "logical") stop('"min.rhozero" in control must be logical.')
	if(class(con$progress) != "logical") stop('"progress" in control must be logical.')
	if(class(con$openmp) != "logical") stop('"openmp" in control must be logical.')
	#END CHEDKING ERROR
	
	
	#START!!
	if(con$progress==TRUE){
		pb <- txtProgressBar(min = 0, max = (con$length.rho*con$length.gamma+1), style = 3)
		setTxtProgressBar(pb, 0)
	}else{
		pb <- NULL
	}
	pmax_for_S <- as.integer(500) #maximum number of vairables which calculates the GFI and AGFI
	m <- factors #number of factors
	xmissing <- 0 #missing indicator for data matrix
	likelihood.availale <-  TRUE #if both x and n.obs are missing, likelihood is not available.
	
	if(missing(x)){
		if(missing(covmat)) stop("input data matrix or covariance matrix is needed.")
		p <- ncol(covmat) #number of variables
		diagS <- diag(covmat) #diagonal elements of sample covariance matrix
		if(missing(n.obs)){
			N <- p+1 #dummy
			likelihood.availale <- FALSE
		}
		else N <- n.obs
		x <- x0 <- 1 #dummy
		xmissing <- 1 #xmissing is an indicator variable
		Npflag <- as.integer(1) #dummy
	}else{
		p <- ncol(x) #number of variables
		N <- nrow(x)   #number of observations
		if(normalize==TRUE){
			x <- scale(x, scale=TRUE)[1:N,1:p] /sqrt(N-1)*sqrt(N) #scaling
		}else{
			x <- scale(x, scale=FALSE)[1:N,1:p] #centering
		}
		if(N>=p || p <= pmax_for_S){
			covmat <- cov(x) * (N-1) / N
			if(N>=p) Npflag <- as.integer(1)
			if(N<p) Npflag <- as.integer(0)
			diagS <- diag(covmat)
		}else{
			covmat <- 1 #dummy
			Npflag <- as.integer(0)
            diagS <- apply(x,2,function(x) var(x) * (N-1) / N) #diagonal elements of sample covariance matrix
        }
		x0 <- x[1:N,1:p]/sqrt(N) #data matrix
	}
    
    
	dimS <- as.integer(length(covmat)) #p*p
	dimx0 <- as.integer(length(x0)) #p*m
	
	
	#set initial values
	rhos_length=5000
	
	
	#set gammmamax and initial value of  Lambda
	Lambdainitialij_for_rhoinitial <- seq(0.3,2,length=10)  #constant factor for initial value of Lambda
	rho_cand <- rep(0,length(Lambdainitialij_for_rhoinitial)) #candidates of initial value of Lambda
	
	#SET INITIAL VALUES
    parameterupdate=.Call("LamdPsiupdate_forinitial_C", as.integer(p), as.integer(N), as.double(con$tol.em), as.double(con$tol.cd), as.integer(con$ncand.initial), as.integer(con$maxit.initial), as.integer(con$maxit.cd), as.integer(con$maxit.bfgs), as.double(con$eta),  as.double(con$min.uniquevar),Inf,covmat,diagS, x0, as.integer(0), Lambdainitialij_for_rhoinitial, Npflag, dimx0, dimS, as.integer(con$openmp), as.integer(con$num.threads))
	if(con$progress==TRUE) setTxtProgressBar(pb, 1)
	
	Lambda <- matrix(0,p,m)
	Lambda[,1] <- parameterupdate[[1]]
	diagPsi <- parameterupdate[[2]]
	if(missing(rho.max)){
		rho.max <- parameterupdate[[3]]
		if(rho.max>8.25) warning('"rho.max" is greater than 8.25. In such cases, the reparametrization of the penalty funcion may be failed')
	}
	
	#Set initival value of Psi
	rho.min <- rho.max * con$Delta
	rho_vec <- exp(seq(log(rho.max), log(rho.min), length= con$length.rho))
	if(con$min.rhozero==TRUE) rho_vec[con$length.rho] = 0.0
	
	#Compute appropriate value of rho
	if(con$length.gamma<2) gamma_vec=Inf
	else if(con$length.gamma<3) gamma_vec=c(Inf, con$min.gamma)
	else gamma_vec <- c(Inf,  exp(seq(from=(log(con$max.gamma)),to=(log(con$min.gamma)),length.out=(con$length.gamma-1))) )
	rhomatrix <- matrix(NA,con$length.rho,con$length.gamma)
	rhomatrix[,1] <- rho_vec #rho vector for lasso
	if(con$length.gamma!=1){
		for(i in 1:con$length.rho){
			for(j in 2:con$length.gamma){
				reparametrization.f <- function(rhos) pnorm(gamma_vec[j]*rhos) - gamma_vec[j]*pnorm(rhos) + (gamma_vec[j]-1)*pnorm(rho_vec[i])
				rhomatrix[i,j] <- uniroot(reparametrization.f,c(rho_vec[i], max(rho.max,8.3)),tol=.Machine$double.eps)$root
			}
		}
	}
	rhomatrix[,1] <- rho_vec #rho vector for lasso
	
	if(con$min.rhozero==TRUE) rhomatrix[con$length.rho,] = 0.0
	colnames(rhomatrix) <- NULL
	
	if(con$start=="warm") allrandom.initial <- as.integer(0)
	if(con$start=="cold") allrandom.initial <- as.integer(1)
	
	#####main#####
    result =.Call("fancmain", Lambda, diagPsi, as.integer(N), as.double(con$tol.em), as.double(con$tol.cd), as.double(con$tol.bfgs), as.integer(con$ncand.initial), as.integer(con$maxit.em), as.integer(con$maxit.cd), as.integer(con$maxit.bfgs), as.integer(con$maxit.initial), as.double(con$eta), as.double(con$min.uniquevar), rhomatrix, as.matrix(gamma_vec), x0, diagS, covmat, Npflag, dimx0, dimS, pmax_for_S, as.integer(cor.factor),as.double(con$zita), as.integer(allrandom.initial), as.integer(con$progress), as.integer(con$openmp), as.integer(con$num.threads), pb)
	if(con$progress==TRUE) close(pb)
	#name
	name_factor <- paste("Factor",1:factors,sep="")
	if(is.null(colnames(x))==FALSE) name_variables <- colnames(x)
	if(is.null(colnames(x))==TRUE) name_variables <- paste("V",1:p,sep="")
	if(xmissing==1&&is.null(colnames(covmat))==FALSE) name_variables <- colnames(covmat)
	
	#make sparse matrix
	index_resultLambda0all <- arrayInd((result[[1]]+1), .dim=c(p,con$length.rho*con$length.gamma*m))
	result_Lambda0all <- sparseMatrix(i=index_resultLambda0all[,1],j=index_resultLambda0all[,2], x=result[[2]], dims=c(p,con$length.rho*con$length.gamma*m))
	Lambda_all <- vector("list", con$length.gamma)
	names(Lambda_all) <- paste("gamma",1:con$length.gamma,sep="")
	for(i in 1:con$length.gamma){
		Lambda_all[[i]] <- vector("list", con$length.rho)
		names(Lambda_all[[i]]) <- paste("rho",1:con$length.rho,sep="")
		for(j in 1:con$length.rho){
            i0 <- (i-1)*con$length.rho*m+(j-1)*m + 1
			if(m==1) Lambda_all[[i]][[j]] <- as(as.matrix(result_Lambda0all[,i0:(i0+m-1)]),"dgCMatrix")
			else Lambda_all[[i]][[j]] <- result_Lambda0all[,i0:(i0+m-1)]
			rownames(Lambda_all[[i]][[j]]) <- name_variables
			colnames(Lambda_all[[i]][[j]]) <- name_factor
		}
	}
	
	#result of Psi
	diagPsi_all <- array(0,dim=c(nrow(rhomatrix),p, ncol(rhomatrix)))
	diagPsi_all0 <- t(result[[3]])
	for(i in 1:ncol(rhomatrix)){
		diagPsi_all[,,i] <- diagPsi_all0[(i-1)*nrow(rhomatrix) + (1:nrow(rhomatrix)),]
	}
	colnames(diagPsi_all) <- name_variables
    
	#result of Phi
	Phi_all <- array(result[[10]],dim=c(m,m,nrow(rhomatrix),ncol(rhomatrix)))
    colnames(Phi_all) <- name_factor
    rownames(Phi_all) <- name_factor
    
	#result of model selection criteria, goodness-of-fit-index
	AIC_all <- matrix(result[[5]], nrow(rhomatrix), ncol(rhomatrix))
	BIC_all <- matrix(result[[6]], nrow(rhomatrix), ncol(rhomatrix))
	CAIC_all <- matrix(result[[7]], nrow(rhomatrix), ncol(rhomatrix))
	dfmatrix <- matrix(result[[8]],nrow(rhomatrix), ncol(rhomatrix))
	nonzero.loadings <- matrix(result[[12]],nrow(rhomatrix), ncol(rhomatrix))
	GFI <- matrix(result[[9]],nrow(rhomatrix), ncol(rhomatrix))
	AGFI <- AGFI <- 1 - (p*(p+1)*(1-GFI)) / (p*(p+1)-2*dfmatrix)


    
	#result of likelihood:  log(det(Sigma)) + tr(Sigmainv*S) + penalty
	likelihood <- array(0,dim=c(nrow(rhomatrix),3, ncol(rhomatrix)))
	likelihood0 <- t(result[[4]])
	likelihood0[,2] <- likelihood0[,2] * 2
	for(i in 1:ncol(rhomatrix)){
		likelihood[,,i] <- likelihood0[(i-1)*nrow(rhomatrix) + (1:nrow(rhomatrix)),]
	}
	colnames(likelihood) <- c("logF", "penalty", "logF+penalty")
	

	#extended BIC
	EBIC_all <- BIC_all
#	if(cor.factor==FALSE) EBIC_all <- EBIC_all + 2*con$gamma.ebic*lchoose(p*factors,dfmatrix-p)
#	else EBIC_all <- EBIC_all + 2*con$gamma.ebic*lchoose(p*factors,dfmatrix-p-factors*(factors-1)/2)
	if(cor.factor==FALSE) EBIC_all <- EBIC_all + 2*con$gamma.ebic*(dfmatrix-p)*log(p*factors)
	else EBIC_all <- EBIC_all + 2*con$gamma.ebic*(dfmatrix-p-factors*(factors-1)/2)*log(p*factors)
	
	#BIC for high dimensional data
	BICH_all <- matrix(NA, nrow(rhomatrix), ncol(rhomatrix))
	if(cor.factor==FALSE) for(i in 1:ncol(rhomatrix))  BICH_all[,i] <- N*likelihood[,1,i] + dfmatrix[,i]*log(N)*log(log(p*factors + p))
	else  for(i in 1:ncol(rhomatrix)) BICH_all[,i] <- N*likelihood[,1,i] + dfmatrix[,i]*log(N)*log(log(p*factors + p + factors*(factors-1)/2 ))


	#convergence
	convergence <- result[[11]]
	names(convergence) <- c("EM", "Coordinate descent", "BFGS")
    
	#answer
	ans <-list(loadings=Lambda_all,uniquenesses=diagPsi_all)
	if(cor.factor==TRUE) ans <- append(ans,list(Phi=Phi_all))
	ans <- append(ans,list(rho=rhomatrix, gamma=gamma_vec, df= dfmatrix, nonzero.loadings=nonzero.loadings, convergence=convergence, Npflag=Npflag))
	if(likelihood.availale==TRUE){
		ans <- append(ans, list(likelihood=likelihood, AIC=AIC_all, BIC=BIC_all, BIC.H=BICH_all, EBIC=EBIC_all, CAIC=CAIC_all))
		if(sum(abs(GFI)) != 0) ans <- append(ans, list(GFI=GFI, AGFI=AGFI))
	}
	ans <- append(ans,list(factors=factors, cor.factor=cor.factor))
	if(xmissing==0) ans <- append(ans,list(x=x))
	ans <- append(ans,list(call=match.call()))
    
	class(ans) <- "fanc"
	ans
	
	
}
