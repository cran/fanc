//fanc C
//#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Parse.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <omp.h>
#include "profiling.h"

//////////////////////////////////////////////////
double alphaOne = 1.0;
double alphaminusOne = -1.0;
double betaOne = 1.0;
double betaZero = 0.0;
double Zero = 0.0;
char TRANSN = 'N';
char TRANST = 'T';
char UPPER = 'U';
char jobzN = 'N';

/////////////


//ex; m,zita,A,

/* complete data log-likelihood function with respect to Phi used for vmmin */
double vmmin_Phi_f(int n, double *x, void *ex){
	double ans = 0;
	double *ex_double = (double *) ex;
	
	int m=(int) ex_double[0];
	int mm=m*m;
	int blocksize=64;
	int mblocksize=m*blocksize;
	int one=1;
	int i,j;
	int info=0;
	int ipiv_Phi[m];
	int count1=0;
	
	double zita=ex_double[1];
	double eigenPhivalues[m];
	double work1[mblocksize];
	double Phi[mm];
	double Phi_temp2[mm];
	double A[mm];
	double Phiinv[mm];
	double Phiinv_A[mm];
	
	//definition of log(det(Phi)), tr(inv(Phi))*A
	double logdetPhi = 0.0;
	double trPhiinvA = 0.0;
	
	//initialize
	for( i=0; i<mm; i++){
		Phi[i] = 0.0;
		Phi_temp2[i] = 0.0;
		Phiinv[i] = 0.0;
		Phiinv_A[i] = 0.0;
	}
	for( i=0; i<m; i++){
		eigenPhivalues[i]=0.0;
		Phi[i+i*m] = 1.0;
		Phi_temp2[i+i*m] = 1.0;
	}
	for( i=0; i<mblocksize; i++){
		work1[i]=0.0;
	}
	for( i=0; i<m; i++){
		ipiv_Phi[i]=0;
	}
	
	count1 =0;
	for( i=0; i<m-1; i++){
		for(j=i+1; j<m; j++){
			Phi[i+j*m]=x[count1];
			Phi[j+i*m]=x[count1];
			count1=count1+1;
		}
	}
	for( i=0; i<mm; i++){
		A[i]=ex_double[i+2];
	}
	
	
	//calculation of log(det(Phi))
	//eigen(Phi)
	F77_CALL(dcopy)(&mm, Phi, &one, Phi_temp2, &one);
	F77_CALL(dsyev)(&jobzN, &UPPER, &m, Phi_temp2, &m, eigenPhivalues, work1, &mblocksize, &info);
	logdetPhi = 0.0;
	for( i=0; i<m; i++){
		logdetPhi = logdetPhi + log(eigenPhivalues[i]);
	}
	//Phi^(-1)
	F77_CALL(dcopy)(&mm, Phi, &one, Phiinv, &one);
	F77_CALL(dgetrf)(&m, &m, Phiinv, &m, ipiv_Phi,  &info);
	F77_CALL(dgetri)(&m, Phiinv, &m, ipiv_Phi, work1, &mblocksize, &info);
	//Phi^(-1)A
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &m, &alphaOne, Phiinv, &m, A, &m, &betaZero, Phiinv_A, &m);
	
	for( i=0; i<m; i++){
		trPhiinvA = trPhiinvA + Phiinv_A[i+i*m];
	}
	
	//answer
	ans = (1.0 - zita) * logdetPhi + trPhiinvA;
	return ans;
}


/* first derivative of complete data log-likelihood function with respect to Phi used for vmmin */
void vmmin_Phi_gr(int n, double *x, double *gr, void *ex){
	double *ex_double = (double *) ex;
	
	int m=(int) ex_double[0];
	int mm=m*m;
	int blocksize=64;
	int mblocksize=m*blocksize;
	int one=1;
	int i,j;
	int info=0;
	int ipiv_Phi[m];
	int count1=0;
	
	
	double zita=ex_double[1];
	double work1[mblocksize];
	double Phi[mm];
	double A[mm];
	double Phiinv[mm];
	double Phiinv_A[mm];
	double Phiinv_A_Phiinv[mm];
	
	
	//initialize
	for( i=0; i<mm; i++){
		Phi[i] = 0.0;
		Phiinv[i] = 0.0;
		Phiinv_A[i] = 0.0;
		Phiinv_A_Phiinv[i] = 0.0;
	}
	for( i=0; i<m; i++){
		Phi[i+i*m] = 1.0;
	}
	for( i=0; i<mblocksize; i++){
		work1[i]=0.0;
	}
	for( i=0; i<m; i++){
		ipiv_Phi[i]=0;
	}
	
	count1 =0;
	for( i=0; i<m-1; i++){
		for(j=i+1; j<m; j++){
			Phi[i+j*m]=x[count1];
			Phi[j+i*m]=x[count1];
			count1=count1+1;
		}
	}
	for( i=0; i<mm; i++){
		A[i]=ex_double[i+2];
	}
	
	//Phi^(-1)
	F77_CALL(dcopy)(&mm, Phi, &one, Phiinv, &one);
	F77_CALL(dgetrf)(&m, &m, Phiinv, &m, ipiv_Phi,  &info);
	F77_CALL(dgetri)(&m, Phiinv, &m, ipiv_Phi, work1, &mblocksize, &info);
	//Phi^(-1)A
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &m, &alphaOne, Phiinv, &m, A, &m, &betaZero, Phiinv_A, &m);
	//Phi^(-1)APhi^(-1)
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &m, &alphaOne, Phiinv_A, &m, Phiinv, &m, &betaZero, Phiinv_A_Phiinv, &m);
	
	// first derivative
	count1=0;
	for( i=0; i<m-1; i++){
		for(j=i+1; j<m; j++){
			gr[count1] = 2.0*(1.0-zita)*Phiinv[i+j*m] - 2.0*Phiinv_A_Phiinv[i+j*m];
			count1 = count1 + 1;
		}
	}
	
	
	
}


/* optimize Phi using BFGS (vmmin) */
void vmmin_Phi(int m, double Phi[m*m], double A[m*m], double ex_zita, int ex_maxcount_phi, double ex_tol, int ex_failflag_phi[1])
{
	double Fmin;
	int fncount;
	int maxit=ex_maxcount_phi;
	int trace=0;
	double abstol= ex_tol;
	double reltol= ex_tol;
	int nREPORT=1;
	int grcount;
	int fail;
	int npara = m*(m-1)/2;
	int mask[npara];
	
	int i,j;
	int count1;
	double Phi_nondiag[npara];
	double ex[2+m*m];
	
	char *vmax; //used for freeing the memory
	
	ex[0] = (double) m;
	ex[1] = ex_zita;
	for( i=0; i<m*m; i++){
		ex[i+2]=A[i];
	}
	
	for( i=0; i<npara; i++){
		mask[i] = 1;
	}
	
	count1 =0;
	for( i=0; i<m-1; i++){
		for(j=i+1; j<m; j++){
			Phi_nondiag[count1]= Phi[i+j*m];
			count1=count1+1;
		}
	}
	
	/*vmmin*/
	vmax = vmaxget();
	vmmin(npara, Phi_nondiag, &Fmin, vmmin_Phi_f, vmmin_Phi_gr, maxit, trace, mask,
          abstol, reltol, nREPORT, ex, &fncount, &grcount, &fail);
	vmaxset(vmax);
	if(fail==1) ex_failflag_phi[0]=1;
	count1 =0;
	for( i=0; i<m-1; i++){
		for(j=i+1; j<m; j++){
			Phi[i+j*m] = Phi_nondiag[count1];
			Phi[j+i*m] = Phi_nondiag[count1];
			count1=count1+1;
		}
	}
	
}

/////////////


/*sign function*/
double sign_C(double z)
{
	double ans;
	if(z > 0) ans = 1;
	if(z == 0) ans = 0;
	if(z < 0) ans = -1;
	return(ans);
	
}


/*function of explicit formula of update equation of MC+ for one-dimebsion*/
/*
 double S_func_MC_C(double z, double rho, double gamma)
 {
 
 double ans;
 
 if(fabs(z) <= rho) ans = 0;
 else if(fabs(z) <= rho*gamma) ans = sign_C(z) * (fabs(z) - rho) / (1-1/gamma);
 else ans = z;
 
 return(ans);
 
 
 }
 */

/*function of explicit formula of update equation of MC+ for one-dimebsion*/
double S_func_MC_C(double z, double rho, double gamma)
{
	
	double ans;
	double signz = sign_C(z);
	double delta = signz * (fabs(z) - rho) / (1-1/gamma);
    
	if(gamma>1.0){
		if(fabs(z) > rho*gamma) ans = z;
        else if(fabs(z) <= rho) ans = 0.0;
        else ans = delta;
	}else if(gamma<1.0){
		if(fabs(z) > rho*gamma) ans = z;
        else if(fabs(delta) <= rho*gamma/2) ans = 0.0;
        else ans = signz*rho*gamma;
        //            else ans = 0.0;
	}else{
		ans=0.0;
	}
	
	return(ans);
	
	
}


/*penalty function (MC+)*/
double penaltyfunc_C(double Lambda, double rho, double gamma){
	
	double ans=5.0;
	if(fabs(Lambda) < rho*gamma) ans = rho * (fabs(Lambda) - Lambda*Lambda/(2.0*rho*gamma));
	else ans = rho*rho*gamma/2.0;
	if(rho==0.0) ans=0.0;
	
	return(ans);
	
}



// Update parameters (Lambda, Psi and Phi) using EM with coordinate descent
void LambdaPsiPhiupdate_C(int p, int m, int N, int ex_Npflag, int ex_dimx0, int ex_dimS, double ex_Lambda[p*m], double ex_diagPsi[p], double ex_Phi[m*m], double ex_tol1[1], double ex_tol2[1], double ex_tol3[1], int ex_maxcount1[1], int ex_maxcount2[1], int ex_maxcount_phi[1], double ex_eta[1], double ex_tolPsi[1], double ex_rho[1], double ex_gamma[1], double ex_S[ex_dimS], double ex_diagS[p], double ex_X[ex_dimx0], int ex_flagpenalty[1], int ex_fixindex[p*m], double ans[p*m+p+m*m+3], int corr_factor[1], double ex_zita[1], int ex_failflag[3], int ex_omp[1], int ex_omp_num[1]){
	
	int pm = p*m;
	int mm = m*m;
	int m0;
	int t1=0;
	int blocksize=64;
	int mblocksize=m*blocksize;
	int one=1;
	int info=0;
	int i;
	int j;
	int j0=0;
	int l;
	int l0;
	int zerocolumncheck_flag=0;
	int zerorowcheck_flag=0;
	int flag_phi=0;
	int phi_failflag[1];
	phi_failflag[0]=0;
	
	int m0_vec[m];
	int count2=0;
    
    int n_threads;
    if(ex_omp_num[0]==0){
        n_threads = omp_get_max_threads();
        if(n_threads<1) n_threads=1;
    }else{
        n_threads = ex_omp_num[0];
    }
	
	double sum_Lsabun1=0.0;
	double ALambdanewA_forPsi=0.0;
	double BLambdanew_forPsi=0.0;
	double logF=0.0;
	double penalty0=0.0;
	double sumabsm0=0.0;
	double temp;
	
	
	double *Lambdanew;
	double *Lambdaold;
	double *Lambda0;
	double *Lambda1;
	double *Lambdasabun;
	double *diagPsinew;
	double *diagPsiinv;
	double *diagPsiold;
	double *Psiinv_Lambda;
	double *Psiinv_Lambda_Minv;
	double Phitemp[m*m];
	double Phitemp2[m*m];
	double Phinew[mm];
	double Phiold[mm];
	double Phiinv[mm];
	double A[mm];
	double A0[mm];
	double Atemp[mm];
	double *B;
	double MA[mm];
	double *D;
	double M[mm];
	double Mtemp[mm];
	int ipiv_M[m];
	int ipiv_A[m];
	double Minv[mm];
	double Ainv[mm];
	double *Psiinvhalf_Lambda;
	double diagm[mm];
	double work1[mblocksize];
	double penalized_likelihood[3];
	double eigenMvalues[m];
	double eigenPhivalues[m];
	
	
	//malloc
	Lambdanew = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Lambdanew == NULL ) {
		free( Lambdanew );
		error( "memory allocation is failed\n" );
	}
	Lambdaold = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Lambdaold == NULL ) {
		free( Lambdaold );
		error( "memory allocation is failed\n" );
	}
	Lambda0 = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Lambda0 == NULL ) {
		free( Lambda0 );
		error( "memory allocation is failed\n" );
	}
	Lambda1 = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Lambda1 == NULL ) {
		free( Lambda1 );
		error( "memory allocation is failed\n" );
	}
	Lambdasabun = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Lambdasabun == NULL ) {
		free( Lambdasabun );
		error( "memory allocation is failed\n" );
	}
	diagPsinew = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
	if( diagPsinew == NULL ) {
		free( diagPsinew );
		error( "memory allocation is failed\n" );
	}
	diagPsiinv = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
	if( diagPsiinv == NULL ) {
		free( diagPsiinv );
		error( "memory allocation is failed\n" );
	}
	diagPsiold = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
	if( diagPsiold == NULL ) {
		free( diagPsiold );
		error( "memory allocation is failed\n" );
	}
	Psiinv_Lambda = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Psiinv_Lambda == NULL ) {
		free( Psiinv_Lambda );
		error( "memory allocation is failed\n" );
	}
	Psiinv_Lambda_Minv = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Psiinv_Lambda_Minv == NULL ) {
		free( Psiinv_Lambda_Minv );
		error( "memory allocation is failed\n" );
	}
	B = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( B == NULL ) {
		free( B );
		error( "memory allocation is failed\n" );
	}
	D = (double *)malloc( sizeof( double ) * (N*m) );	// memory allocation
	if( D == NULL ) {
		free( D );
		error( "memory allocation is failed\n" );
	}
	Psiinvhalf_Lambda = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
	if( Psiinvhalf_Lambda == NULL ) {
		free( Psiinvhalf_Lambda );
		error( "memory allocation is failed\n" );
	}
	
	
	
	
	
	
	//initialize
	for( i=0; i<p*m; i++){
		Lambdaold[i] = 10.0;
		Lambdanew[i] = ex_Lambda[i];
		Lambda0[i] = 0.0;
		Lambda1[i] = 0.0;
		Lambdasabun[i] = 0.0;
		Psiinv_Lambda_Minv[i] = 0.0;
		B[i] = 0.0;
		Psiinvhalf_Lambda[i] = 0.0;
		Psiinv_Lambda[i] = 0.0;
	}
	
	
	
	for( i=0; i<mblocksize; i++){
		work1[i] = 0.0;
	}
	
	
	for( i=0; i<m*m; i++){
		A[i] = 0.0;
		A0[i] = 0.0;
		Atemp[i] = 0.0;
		M[i] = 0.0;
		Mtemp[i] = 0.0;
		MA[i] = 0.0;
		Minv[i] = 0.0;
		diagm[i] = 0.0;
		Phitemp[i] = ex_Phi[i];
		Phiold[i] =  ex_Phi[i];
		Phinew[i] =  ex_Phi[i];
		Phitemp2[i] = 0.0;
		Phiinv[i] =  0.0;
	}
	
	for(i=0; i<m; i++){
		diagm[i*m + i] = 1.0;
	}
	
	
	for( i=0; i<m*N; i++){
		D[i] = 0.0;
	}
	
	
	for( i=0; i<p*m; i++){
		Lambdasabun[i] = Lambdanew[i] - Lambdaold[i] ;
	}
	
	for( i=0; i<p; i++){
		diagPsinew[i] = ex_diagPsi[i];
		diagPsiinv[i] = 1.0 / diagPsinew[i];
		diagPsiold[i] = 0.0;
	}
	
	for( i=0; i<p; i++){
		diagPsiinv[i] = 1.0 / ex_diagPsi[i];
	}
	
	for( i=0; i<m; i++){
		ipiv_M[i] = 0;
		ipiv_A[i] = 0;
		eigenMvalues[i] = 0.0;
		eigenPhivalues[i]=0.0;
	}
	
	for( i=0; i<3; i++){
		penalized_likelihood[i] = 0.0;
	}
	
	sum_Lsabun1 = F77_CALL(dnrm2)(&pm, Lambdasabun, &one);
	
	//start
	t1=0;
    //if(ex_omp[0]==1){
    while( sum_Lsabun1 > ex_tol1[0]  && t1 < ex_maxcount1[0]){
        //prof_start(0,"matrix computation");
        t1=t1+1;
        F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambdaold, &one);
        F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda1, &one);
        F77_CALL(dcopy)(&pm, Lambdaold, &one, Psiinvhalf_Lambda, &one);
        F77_CALL(dcopy)(&pm, Lambdaold, &one, Psiinv_Lambda, &one);
        F77_CALL(dcopy)(&p, diagPsinew, &one, diagPsiold, &one);
        F77_CALL(dcopy)(&mm, Phinew, &one, Phiold, &one);
        for( i=0; i<p; i++){
            diagPsiinv[i] = 1.0 / diagPsiold[i];
            for( j=0; j<m; j++){
                Psiinvhalf_Lambda[i+j*p] = Psiinvhalf_Lambda[i+j*p] * sqrt(diagPsiinv[i]);
                Psiinv_Lambda[i+j*p] = Psiinv_Lambda[i+j*p] * diagPsiinv[i];
            }
        }
        
        //M <- t(Lambda) * Psi^(-1) * Lambda + diag(m)
        if(corr_factor[0]==1 && m>1) {
            F77_CALL(dcopy)(&mm, Phiold, &one, Phiinv, &one);
            F77_CALL(dgetrf)(&m, &m, Phiinv, &m, ipiv_M,  &info);
            F77_CALL(dgetri)(&m, Phiinv, &m, ipiv_M, work1, &mblocksize, &info);
            F77_CALL(dcopy)(&mm, Phiinv, &one, M, &one);
        }else{
            F77_CALL(dcopy)(&mm, diagm, &one, M, &one);
        }
        
        F77_CALL(dsyrk)(&UPPER, &TRANST, &m, &p, &alphaOne, Psiinvhalf_Lambda, &p, &betaOne, M, &m);
        
        for( i=0; i<m; i++){
            for( j=0; j<m; j++){
                M[j+i*m] = M[i+j*m];
            }
        }
        
        //M^(-1)
        F77_CALL(dcopy)(&mm, M, &one, Minv, &one);
        F77_CALL(dgetrf)(&m, &m, Minv, &m, ipiv_M,  &info);
        F77_CALL(dgetri)(&m, Minv, &m, ipiv_M, work1, &mblocksize, &info);
        
        //Psiinv_Lambda_Minv
        F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, Psiinv_Lambda, &p,Minv, &m, &betaZero, Psiinv_Lambda_Minv, &p);
        if(ex_Npflag==1){
            //B
            F77_CALL(dgemm)(&TRANST, &TRANSN, &m, &p, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_S, &p, &betaZero, B, &m);
            //A
            F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
            F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &p, &alphaOne, B, &m, Psiinv_Lambda_Minv, &p, &betaOne, A, &m);
        }else{
            //D
            F77_CALL(dgemm)(&TRANST, &TRANST, &m, &N, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_X, &N, &betaZero, D, &m);
            //B
            F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &p, &N, &alphaOne, D, &m, ex_X, &N, &betaZero, B, &m);
            //A
            F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
            F77_CALL(dgemm)(&TRANSN, &TRANST, &m, &m, &N, &alphaOne, D, &m, D, &m, &betaOne, A, &m);
        }
        //prof_stop(0);
        
        
        
        if(ex_flagpenalty[0]==0){
            //No Penalty
            //A^(-1)
            F77_CALL(dcopy)(&mm, A, &one, Ainv, &one);
            F77_CALL(dgetrf)(&m, &m, Ainv, &m, ipiv_A,  &info);
            F77_CALL(dgetri)(&m, Ainv, &m, ipiv_A, work1, &mblocksize, &info);
            //Lambdanew <- B^TA^(-1)
            F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &m, &m, &alphaOne, B, &m, Ainv, &m, &betaZero, Lambdanew, &p);
            //fix lambda which is assumed to be fixed
            for( i=0; i<pm; i++){
                if(ex_fixindex[i]==1) Lambdanew[i] = ex_Lambda[i];
            }
        }
        
        
        if(ex_flagpenalty[0]==1){
            /*Coordinate descent*/
            F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda1, &one);
            //#ifdef _OPENMP
            //		omp_set_num_threads(4);
            //#endif
            //prof_start(1,"coordinate descent");
            
            if(ex_omp[0]==1){
                
//#pragma omp parallel for
#pragma omp parallel for num_threads(n_threads)
                for( j=0; j<p; j++){
                    double sum_Lsabun2=100.0;
                    int t2 = 0;
                    double ALambdanew=0.0;
                    int l_omp;
                    int l0_omp;
                    double weight_cd=1.0;
                    double theta_tilde=0.0;
                    double rho_adj;
                    double gamma_adj;
                    
                    while(sum_Lsabun2 > ex_tol2[0] && t2 < ex_maxcount2[0]){
                        t2=t2+1;
                        for( l_omp=0; l_omp<m; l_omp++){
                            ALambdanew=0.0;
                            for(l0_omp=0; l0_omp<m; l0_omp++){
                                if(l0_omp != l_omp){
                                    ALambdanew = ALambdanew + A[l_omp*m+l0_omp] * Lambdanew[p*l0_omp + j];
                                }
                            }
                            if(ex_fixindex[j+p*l_omp]==0){
                                theta_tilde = (B[j*m+l_omp] - ALambdanew ) / A[l_omp+l_omp*m];
                                weight_cd = diagPsiold[j] / A[l_omp+l_omp*m];
                                rho_adj = ex_rho[0]*weight_cd;
                                gamma_adj = ex_gamma[0]/weight_cd;
                                Lambdanew[j+p*l_omp] = S_func_MC_C(theta_tilde, rho_adj, gamma_adj );
                                /*
                                 if(gamma_adj>1) Lambdanew[j+p*l_omp] = S_func_MC_C(theta_tilde, rho_adj, gamma_adj );
                                 if(gamma_adj<=1){
                                 if(theta_tilde*theta_tilde/2 > penaltyfunc_C(theta_tilde,rho_adj,gamma_adj)) Lambdanew[j+p*l_omp] = theta_tilde;
                                 else Lambdanew[j+p*l_omp]=0.0;
                                 }
                                 */
                                //Lambdanew[j+p*l_omp] = S_func_MC_C(1.0 / diagPsiold[j] * ( B[j*m+l_omp] - ALambdanew ), ex_rho[0] , ex_gamma[0] )  /  (  1.0 / diagPsiold[j] * A[l_omp+l_omp*m] );
                            }
                        }
                        
                        
                        sum_Lsabun2=0.0;
                        for( l_omp=0; l_omp<m; l_omp++){
                            sum_Lsabun2 = sum_Lsabun2+ (Lambdanew[j+p*l_omp] - Lambda1[j+p*l_omp]) * (Lambdanew[j+p*l_omp] - Lambda1[j+p*l_omp]);
                        }
                        sum_Lsabun2=sqrt(sum_Lsabun2);
                        
                        
                        for( l_omp=0; l_omp<m; l_omp++){
                            Lambda1[j+p*l_omp] = Lambdanew[j+p*l_omp];
                        }
                        
                    }
                    
                    if(t2==ex_maxcount2[0]) ex_failflag[1]=1;
                }
            }else{
                
                for( j=0; j<p; j++){
                    double sum_Lsabun2=100.0;
                    int t2 = 0;
                    double ALambdanew=0.0;
                    int l_omp;
                    int l0_omp;
                    double weight_cd=1.0;
                    double theta_tilde=0.0;
                    double rho_adj;
                    double gamma_adj;
                    
                    while(sum_Lsabun2 > ex_tol2[0] && t2 < ex_maxcount2[0]){
                        t2=t2+1;
                        for( l_omp=0; l_omp<m; l_omp++){
                            ALambdanew=0.0;
                            for(l0_omp=0; l0_omp<m; l0_omp++){
                                if(l0_omp != l_omp){
                                    ALambdanew = ALambdanew + A[l_omp*m+l0_omp] * Lambdanew[p*l0_omp + j];
                                }
                            }
                            if(ex_fixindex[j+p*l_omp]==0){
                                theta_tilde = (B[j*m+l_omp] - ALambdanew ) / A[l_omp+l_omp*m];
                                weight_cd = diagPsiold[j] / A[l_omp+l_omp*m];
                                rho_adj = ex_rho[0]*weight_cd;
                                gamma_adj = ex_gamma[0]/weight_cd;
                                Lambdanew[j+p*l_omp] = S_func_MC_C(theta_tilde, rho_adj, gamma_adj );
                                /*
                                 if(gamma_adj>1) Lambdanew[j+p*l_omp] = S_func_MC_C(theta_tilde, rho_adj, gamma_adj );
                                 if(gamma_adj<=1){
                                 if(theta_tilde*theta_tilde/2 > penaltyfunc_C(theta_tilde,rho_adj,gamma_adj)) Lambdanew[j+p*l_omp] = theta_tilde;
                                 else Lambdanew[j+p*l_omp]=0.0;
                                 }
                                 */
                                //Lambdanew[j+p*l_omp] = S_func_MC_C(1.0 / diagPsiold[j] * ( B[j*m+l_omp] - ALambdanew ), ex_rho[0] , ex_gamma[0] )  /  (  1.0 / diagPsiold[j] * A[l_omp+l_omp*m] );
                            }
                        }
                        
                        
                        sum_Lsabun2=0.0;
                        for( l_omp=0; l_omp<m; l_omp++){
                            sum_Lsabun2 = sum_Lsabun2+ (Lambdanew[j+p*l_omp] - Lambda1[j+p*l_omp]) * (Lambdanew[j+p*l_omp] - Lambda1[j+p*l_omp]);
                        }
                        sum_Lsabun2=sqrt(sum_Lsabun2);
                        
                        
                        for( l_omp=0; l_omp<m; l_omp++){
                            Lambda1[j+p*l_omp] = Lambdanew[j+p*l_omp];
                        }
                        
                    }
                    
                    if(t2==ex_maxcount2[0]) ex_failflag[1]=1;
                }
            }
            
            
            //prof_stop(1);
            
        }
        //end coordinate descent
        
        
        
        
        
        
        /*update Psi*/
        for( i=0; i<p; i++){
            ALambdanewA_forPsi=0.0;
            BLambdanew_forPsi=0.0;
            for(l=0; l<m; l++){
                BLambdanew_forPsi=BLambdanew_forPsi+ B[i*m+l] * Lambdanew[p*l + i];
                for(l0=0; l0<m; l0++){
                    ALambdanewA_forPsi = ALambdanewA_forPsi + A[l*m+l0] * Lambdanew[p*l0 + i] *  Lambdanew[p*l + i];
                }
            }
            diagPsinew[i] = ex_diagS[i] - 2 * BLambdanew_forPsi + ALambdanewA_forPsi + ex_diagS[i]*ex_eta[0];
            if(diagPsinew[i] < ex_tolPsi[0]) diagPsinew[i] = ex_tolPsi[0];
        }
        
        
        
        /*update Pji*/
        
        //compute m0
        m0=0;
        for(i=0;i<m;i++){
            m0_vec[i]=0;
        }
        for(j=0;j<m;j++){
            sumabsm0=0.0;
            for(i=0;i<p;i++){
                sumabsm0 = sumabsm0 + fabs(Lambdanew[i+j*p]);
            }
            if(sumabsm0 > ex_tol1[0]){
                m0=m0+1;
                m0_vec[j]=1;
            }
        }
        
        if(m0<m){
            //post-processing: The zero-columns move to right side
            F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda0, &one);
            F77_CALL(dcopy)(&mm, A, &one, A0, &one);
            for(i=0;i<p*m;i++){
                Lambdanew[i]=0.0;
            }
            count2=0;
            for(i=0;i<m;i++){
                if(m0_vec[i]==1){
                    for(j=0;j<p;j++){
                        Lambdanew[j+count2*p]=Lambda0[j+i*p];
                    }
                    count2=count2+1;
                }
            }
        }
        
        
        
        
        flag_phi=0;
        for(i=0;i<m0;i++){
            if(m0_vec[i]!=1){
                flag_phi=1;
            }
        }
        
        /*update Phi if flag_phi=1*/
        
        //initialize Phi
        if(m0>1 && corr_factor[0]==1){
            if(flag_phi==1){
                for( i=0; i<m*m; i++){
                    Phiold[i] = 0.0;
                }
                for( i=0; i<m; i++){
                    Phiold[i+i*m] = 1.0;
                }
            }
        }
        
        
        for(i=0; i<m; i++){
            for(j=0; j<m; j++){
                Phinew[i*m+j]=0.0;
            }
        }
        for(i=0; i<m; i++){
            Phinew[i*m+i]=1.0;
        }
        
        //initialize Atemp
        if(m0>1 && corr_factor[0]==1){
            for( i=0; i<mm; i++){
                Atemp[i]=0.0;
            }
            count2=0;
            for(i=0;i<m;i++){
                if(m0_vec[i]==1){
                    for(j=0;j<m;j++){
                        if(m0_vec[j]==1){
                            Atemp[count2]=A[i*m+j];
                            count2=count2+1;
                        }
                    }
                }
            }
            
            
            
            for( i=0; i<mm; i++){
                Phitemp[i]=0.0;
            }
            for(i=0; i<m0; i++){
                for(j=0; j<m0; j++){
                    Phitemp[i*m0+j]=Phiold[i*m+j];
                }
            }
            
            //update phi using BFGS
            //prof_start(2,"BFGS");
            vmmin_Phi(m0, Phitemp, Atemp, ex_zita[0],ex_maxcount_phi[0], ex_tol3[0], phi_failflag);
            //prof_stop(2);
            if(phi_failflag[0]==1) ex_failflag[2] = phi_failflag[0];
            
            
            
            for(i=0; i<m0; i++){
                for(j=0; j<m0; j++){
                    Phinew[i*m+j]=Phitemp[i*m0+j];
                }
            }
            
        }
        
        //calculate the distance between Lambdanew and Lambdaold
        for( i=0; i<p*m; i++){
            Lambdasabun[i] = Lambdanew[i] - Lambdaold[i];
        }
        sum_Lsabun1 = F77_CALL(dnrm2)(&pm, Lambdasabun, &one);
        //                }
        //end pragmaomp single
    }
    //end EM
    //        }
    //end pragmaomp
    
    //}
    //end if
    
    /*
    if(ex_omp[0]==0){
        while( sum_Lsabun1 > ex_tol1[0]  && t1 < ex_maxcount1[0]){
            //prof_start(0,"matrix computation");
            t1=t1+1;
            F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambdaold, &one);
            F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda1, &one);
            F77_CALL(dcopy)(&pm, Lambdaold, &one, Psiinvhalf_Lambda, &one);
            F77_CALL(dcopy)(&pm, Lambdaold, &one, Psiinv_Lambda, &one);
            F77_CALL(dcopy)(&p, diagPsinew, &one, diagPsiold, &one);
            F77_CALL(dcopy)(&mm, Phinew, &one, Phiold, &one);
            for( i=0; i<p; i++){
                diagPsiinv[i] = 1.0 / diagPsiold[i];
                for( j=0; j<m; j++){
                    Psiinvhalf_Lambda[i+j*p] = Psiinvhalf_Lambda[i+j*p] * sqrt(diagPsiinv[i]);
                    Psiinv_Lambda[i+j*p] = Psiinv_Lambda[i+j*p] * diagPsiinv[i];
                }
            }
            
            //M <- t(Lambda) * Psi^(-1) * Lambda + diag(m)
            if(corr_factor[0]==1 && m>1) {
                F77_CALL(dcopy)(&mm, Phiold, &one, Phiinv, &one);
                F77_CALL(dgetrf)(&m, &m, Phiinv, &m, ipiv_M,  &info);
                F77_CALL(dgetri)(&m, Phiinv, &m, ipiv_M, work1, &mblocksize, &info);
                F77_CALL(dcopy)(&mm, Phiinv, &one, M, &one);
            }else{
                F77_CALL(dcopy)(&mm, diagm, &one, M, &one);
            }
            
            F77_CALL(dsyrk)(&UPPER, &TRANST, &m, &p, &alphaOne, Psiinvhalf_Lambda, &p, &betaOne, M, &m);
            
            for( i=0; i<m; i++){
                for( j=0; j<m; j++){
                    M[j+i*m] = M[i+j*m];
                }
            }
            
            //M^(-1)
            F77_CALL(dcopy)(&mm, M, &one, Minv, &one);
            F77_CALL(dgetrf)(&m, &m, Minv, &m, ipiv_M,  &info);
            F77_CALL(dgetri)(&m, Minv, &m, ipiv_M, work1, &mblocksize, &info);
            
            //Psiinv_Lambda_Minv
            F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, Psiinv_Lambda, &p,Minv, &m, &betaZero, Psiinv_Lambda_Minv, &p);
            if(ex_Npflag==1){
                //B
                F77_CALL(dgemm)(&TRANST, &TRANSN, &m, &p, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_S, &p, &betaZero, B, &m);
                //A
                F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
                F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &p, &alphaOne, B, &m, Psiinv_Lambda_Minv, &p, &betaOne, A, &m);
            }else{
                //D
                F77_CALL(dgemm)(&TRANST, &TRANST, &m, &N, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_X, &N, &betaZero, D, &m);
                //B
                F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &p, &N, &alphaOne, D, &m, ex_X, &N, &betaZero, B, &m);
                //A
                F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
                F77_CALL(dgemm)(&TRANSN, &TRANST, &m, &m, &N, &alphaOne, D, &m, D, &m, &betaOne, A, &m);
            }
            //prof_stop(0);
            
            
            if(ex_flagpenalty[0]==0){
                //No Penalty
                //A^(-1)
                F77_CALL(dcopy)(&mm, A, &one, Ainv, &one);
                F77_CALL(dgetrf)(&m, &m, Ainv, &m, ipiv_A,  &info);
                F77_CALL(dgetri)(&m, Ainv, &m, ipiv_A, work1, &mblocksize, &info);
                //Lambdanew <- B^TA^(-1)
                F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &m, &m, &alphaOne, B, &m, Ainv, &m, &betaZero, Lambdanew, &p);
                //fix lambda which is assumed to be fixed
                for( i=0; i<pm; i++){
                    if(ex_fixindex[i]==1) Lambdanew[i] = ex_Lambda[i];
                }
            }
            
            
            if(ex_flagpenalty[0]==1){
                //Coordinate descent//
                F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda1, &one);
                //#ifdef _OPENMP
                //		omp_set_num_threads(4);
                //#endif
                
                //prof_start(1,"coordinate descent");
                
                for( j=0; j<p; j++){
                    double sum_Lsabun2=100.0;
                    int t2 = 0;
                    double ALambdanew=0.0;
                    int l_omp;
                    int l0_omp;
                    double weight_cd=1.0;
                    double theta_tilde=0.0;
                    double rho_adj;
                    double gamma_adj;
                    
                    while(sum_Lsabun2 > ex_tol2[0] && t2 < ex_maxcount2[0]){
                        t2=t2+1;
                        for( l_omp=0; l_omp<m; l_omp++){
                            ALambdanew=0.0;
                            for(l0_omp=0; l0_omp<m; l0_omp++){
                                if(l0_omp != l_omp){
                                    ALambdanew = ALambdanew + A[l_omp*m+l0_omp] * Lambdanew[p*l0_omp + j];
                                }
                            }
                            if(ex_fixindex[j+p*l_omp]==0){
                                theta_tilde = (B[j*m+l_omp] - ALambdanew ) / A[l_omp+l_omp*m];
                                weight_cd = diagPsiold[j] / A[l_omp+l_omp*m];
                                rho_adj = ex_rho[0]*weight_cd;
                                gamma_adj = ex_gamma[0]/weight_cd;
                                Lambdanew[j+p*l_omp] = S_func_MC_C(theta_tilde, rho_adj, gamma_adj );
                                //Lambdanew[j+p*l_omp] = S_func_MC_C(1.0 / diagPsiold[j] * ( B[j*m+l_omp] - ALambdanew ), ex_rho[0] , ex_gamma[0] )  /  (  1.0 / diagPsiold[j] * A[l_omp+l_omp*m] );
                            }
                        }
                        
                        
                        sum_Lsabun2=0.0;
                        for( l_omp=0; l_omp<m; l_omp++){
                            sum_Lsabun2 = sum_Lsabun2+ (Lambdanew[j+p*l_omp] - Lambda1[j+p*l_omp]) * (Lambdanew[j+p*l_omp] - Lambda1[j+p*l_omp]);
                        }
                        sum_Lsabun2=sqrt(sum_Lsabun2);
                        
                        
                        for( l_omp=0; l_omp<m; l_omp++){
                            Lambda1[j+p*l_omp] = Lambdanew[j+p*l_omp];
                        }
                        
                    }
                    
                    if(t2==ex_maxcount2[0]) ex_failflag[1]=1;
                }
            }
            //prof_stop(1);
            
            
            
            
            
            
            //update Psi//
            for( i=0; i<p; i++){
                ALambdanewA_forPsi=0.0;
                BLambdanew_forPsi=0.0;
                for(l=0; l<m; l++){
                    BLambdanew_forPsi=BLambdanew_forPsi+ B[i*m+l] * Lambdanew[p*l + i];
                    for(l0=0; l0<m; l0++){
                        ALambdanewA_forPsi = ALambdanewA_forPsi + A[l*m+l0] * Lambdanew[p*l0 + i] *  Lambdanew[p*l + i];
                    }
                }
                diagPsinew[i] = ex_diagS[i] - 2 * BLambdanew_forPsi + ALambdanewA_forPsi + ex_diagS[i]*ex_eta[0];
                if(diagPsinew[i] < ex_tolPsi[0]) diagPsinew[i] = ex_tolPsi[0];
            }
            
            
            
            //update Pji//
            
            //compute m0
            m0=0;
            for(i=0;i<m;i++){
                m0_vec[i]=0;
            }
            for(j=0;j<m;j++){
                sumabsm0=0.0;
                for(i=0;i<p;i++){
                    sumabsm0 = sumabsm0 + fabs(Lambdanew[i+j*p]);
                }
                if(sumabsm0 > ex_tol1[0]){
                    m0=m0+1;
                    m0_vec[j]=1;
                }
            }
            
            if(m0<m){
                //post-processing: The zero-columns move to right side
                F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda0, &one);
                F77_CALL(dcopy)(&mm, A, &one, A0, &one);
                for(i=0;i<p*m;i++){
                    Lambdanew[i]=0.0;
                }
                count2=0;
                for(i=0;i<m;i++){
                    if(m0_vec[i]==1){
                        for(j=0;j<p;j++){
                            Lambdanew[j+count2*p]=Lambda0[j+i*p];
                        }
                        count2=count2+1;
                    }
                }
            }
            
            
            
            
            flag_phi=0;
            for(i=0;i<m0;i++){
                if(m0_vec[i]!=1){
                    flag_phi=1;
                }
            }
            
            //update Phi if flag_phi=1//
            
            //initialize Phi
            if(m0>1 && corr_factor[0]==1){
                if(flag_phi==1){
                    for( i=0; i<m*m; i++){
                        Phiold[i] = 0.0;
                    }
                    for( i=0; i<m; i++){
                        Phiold[i+i*m] = 1.0;
                    }
                }
            }
            
            
            for(i=0; i<m; i++){
                for(j=0; j<m; j++){
                    Phinew[i*m+j]=0.0;
                }
            }
            for(i=0; i<m; i++){
                Phinew[i*m+i]=1.0;
            }
            
            //initialize Atemp
            if(m0>1 && corr_factor[0]==1){
                for( i=0; i<mm; i++){
                    Atemp[i]=0.0;
                }
                count2=0;
                for(i=0;i<m;i++){
                    if(m0_vec[i]==1){
                        for(j=0;j<m;j++){
                            if(m0_vec[j]==1){
                                Atemp[count2]=A[i*m+j];
                                count2=count2+1;
                            }
                        }
                    }
                }
                
                
                
                for( i=0; i<mm; i++){
                    Phitemp[i]=0.0;
                }
                for(i=0; i<m0; i++){
                    for(j=0; j<m0; j++){
                        Phitemp[i*m0+j]=Phiold[i*m+j];
                    }
                }
                
                //update phi using BFGS
                //prof_start(2,"BFGS");
                vmmin_Phi(m0, Phitemp, Atemp, ex_zita[0],ex_maxcount_phi[0], ex_tol3[0], phi_failflag);
                //prof_stop(2);
                if(phi_failflag[0]==1) ex_failflag[2] = phi_failflag[0];
                
                
                
                for(i=0; i<m0; i++){
                    for(j=0; j<m0; j++){
                        Phinew[i*m+j]=Phitemp[i*m0+j];
                    }
                }
                
            }
            
            //calculate the distance between Lambdanew and Lambdaold
            for( i=0; i<p*m; i++){
                Lambdasabun[i] = Lambdanew[i] - Lambdaold[i];
            }
            sum_Lsabun1 = F77_CALL(dnrm2)(&pm, Lambdasabun, &one);
        }
        //end EM
    }
    //end if


*/
    
    R_CheckUserInterrupt();
    
    if(t1==ex_maxcount1[0]) ex_failflag[0]=1;
    
    //If we have only one element of zero at each column, we modify it by zero
    for( j=0; j<m; j++){
        for( j0=0; j0<m; j0++){
            if(fabs(Phinew[j0+j*m]) < ex_tol3[0]) Phinew[j0+j*m]=0.0;
        }
    }
    for( j=0; j<m; j++){
        zerocolumncheck_flag=0;
        for(i=0; i<p; i++){
            if(fabs(Lambdanew[i+j*p]) > ex_tol1[0]) zerocolumncheck_flag = zerocolumncheck_flag + 1;
        }
        if(zerocolumncheck_flag==0 || zerocolumncheck_flag==1){
            for( i=0; i<p; i++){
                if(fabs(Lambdanew[i+j*p]) > ex_tol1[0]){
                    temp=0.0;
                    for( j0=0; j0<m; j0++){
                        if(j0!=j) temp = temp + fabs(Phinew[j0+j*m]);
                    }
                    if(temp==0.0){
                        diagPsinew[i]=diagPsinew[i]+Lambdanew[i+j*p]*Lambdanew[i+j*p];
                        Lambdanew[i+j*p]=0.0;
                    }
                }else{
                    Lambdanew[i+j*p]=0.0;
                }
            }
        }
    }
    
    /*The diagonal elements of Psi is diagPsi for row with zero-elements.  */
    for( i=0; i<p; i++){
        zerorowcheck_flag=0;
        for( j=0; j<m; j++){
            if(fabs(Lambdanew[i+j*p]) < ex_tol1[0]) zerorowcheck_flag = zerorowcheck_flag + 1;
        }
        if(zerorowcheck_flag==m){
            diagPsinew[i]=ex_diagS[i];
        }
    }
    
    
    //change columns
    m0=0;
    for(i=0;i<m;i++){
        m0_vec[i]=0;
    }
    for(j=0;j<m;j++){
        sumabsm0=0.0;
        for(i=0;i<p;i++){
            sumabsm0 = sumabsm0 + fabs(Lambdanew[i+j*p]);
        }
        if(sumabsm0 > ex_tol1[0]){
            m0=m0+1;
            m0_vec[j]=1;
        }
    }
    
    //F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda0, &one);
    if(m0<m){
        //post-processing2
        F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda0, &one);
        F77_CALL(dcopy)(&mm, A, &one, A0, &one);
        for(i=0;i<p*m;i++){
            Lambdanew[i]=0.0;
        }
        count2=0;
        for(i=0;i<m;i++){
            if(m0_vec[i]==1){
                for(j=0;j<p;j++){
                    Lambdanew[j+count2*p]=Lambda0[j+i*p];
                }
                count2=count2+1;
            }
        }
    }
    
    
    
    
    
    
    //calculation of log-likelihood
    //M <- t(Lambda) * Psi^(-1) * Lambda + diagm
    F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambdaold, &one);
    F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda1, &one);
    F77_CALL(dcopy)(&pm, Lambdaold, &one, Psiinvhalf_Lambda, &one);
    F77_CALL(dcopy)(&pm, Lambdaold, &one, Psiinv_Lambda, &one);
    for( i=0; i<p; i++){
        diagPsiinv[i] = 1.0 / diagPsinew[i];
        for( j=0; j<m; j++){
            Psiinvhalf_Lambda[i+j*p] = Psiinvhalf_Lambda[i+j*p] * sqrt(diagPsiinv[i]);
            Psiinv_Lambda[i+j*p] = Psiinv_Lambda[i+j*p] * diagPsiinv[i];
        }
    }
    
    if(corr_factor[0]==1 && m>1) {
        F77_CALL(dcopy)(&mm, Phinew, &one, Phiinv, &one);
        F77_CALL(dgetrf)(&m, &m, Phiinv, &m, ipiv_M,  &info);
        F77_CALL(dgetri)(&m, Phiinv, &m, ipiv_M, work1, &mblocksize, &info);
        F77_CALL(dcopy)(&mm, Phiinv, &one, M, &one);
    }else{
        F77_CALL(dcopy)(&mm, diagm, &one, M, &one);
    }
    
    F77_CALL(dsyrk)(&UPPER, &TRANST, &m, &p, &alphaOne, Psiinvhalf_Lambda, &p, &betaOne, M, &m);
    
    for( i=0; i<m; i++){
        for( j=0; j<m; j++){
            M[j+i*m] = M[i+j*m];
        }
    }
    
    //M^(-1)
    F77_CALL(dcopy)(&mm, M, &one, Minv, &one);
    F77_CALL(dgetrf)(&m, &m, Minv, &m, ipiv_M,  &info);
    F77_CALL(dgetri)(&m, Minv, &m, ipiv_M, work1, &mblocksize, &info);
    
    //Psiinv_Lambda_Minv
    F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, Psiinv_Lambda, &p, Minv, &m, &betaZero, Psiinv_Lambda_Minv, &p);
    if(ex_Npflag==1){
        //B
        F77_CALL(dgemm)(&TRANST, &TRANSN, &m, &p, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_S, &p, &betaZero, B, &m);
        //A
        F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
        F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &p, &alphaOne, B, &m, Psiinv_Lambda_Minv, &p, &betaOne, A, &m);
    }else{
        //D
        F77_CALL(dgemm)(&TRANST, &TRANST, &m, &N, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_X, &N, &betaZero, D, &m);
        //B
        F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &p, &N, &alphaOne, D, &m, ex_X, &N, &betaZero, B, &m);
        //A
        F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
        F77_CALL(dgemm)(&TRANSN, &TRANST, &m, &m, &N, &alphaOne, D, &m, D, &m, &betaOne, A, &m);
    }
    
    //M*A
    F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &m, &alphaOne, M, &m, A, &m, &betaZero, MA, &m);
    //det(M)
    F77_CALL(dcopy)(&mm, M, &one, Mtemp, &one);
    F77_CALL(dsyev)(&jobzN, &UPPER, &m, Mtemp, &m, eigenMvalues, work1, &mblocksize, &info);
    //det(Phi)
    if(corr_factor[0]==1){
        F77_CALL(dcopy)(&mm, Phinew, &one, Phitemp2, &one);
        F77_CALL(dsyev)(&jobzN, &UPPER, &m, Phitemp2, &m, eigenPhivalues, work1, &mblocksize, &info);
    }
    logF = 0.0;
    for( i=0; i<p; i++){
        logF = logF + log(diagPsinew[i]) + ex_diagS[i] / diagPsinew[i];
    }
    
    if(corr_factor[0]==1){
        for( i=0; i<m; i++){
            logF = logF + log(eigenMvalues[i]) + log(eigenPhivalues[i])  - ( MA[i+i*m] - 1.0);
        }
    }
    
    if(corr_factor[0]==0){
        for( i=0; i<m; i++){
            logF = logF + log(eigenMvalues[i]) - ( MA[i+i*m] - 1.0);
        }
    }
    
    
    //calculation of penalty
    penalty0=0.0;
    for( i=0; i<p*m; i++){
        penalty0 = penalty0 + penaltyfunc_C(Lambdanew[i],ex_rho[0],ex_gamma[0]);
    }
    for( i=0; i<p; i++){
        penalty0 = penalty0 + ex_eta[0] * ex_diagS[i] * diagPsiinv[i] / 2.0 ;
    }
    
    penalized_likelihood[0] = logF;
    penalized_likelihood[1] = penalty0;
    penalized_likelihood[2] = logF + 2.0*penalty0;
    
    
    for( i=0; i<p*m; i++){
        ans[i] = Lambdanew[i];
    }
    for( i=0; i<p; i++){
        ans[i+p*m] = diagPsinew[i];
    }
    for( i=0; i<m*m; i++){
        ans[i+p*(m+1)] = Phinew[i];
    }
    for( i=0; i<3; i++){
        ans[i+p*(m+1)+m*m] = penalized_likelihood[i];
    }
    
    //free memory
    free(Lambdanew);
    free(Lambdaold);
    free(Lambda0);
    free(Lambda1);
    free(Lambdasabun);
    free(diagPsinew);
    free(diagPsiinv);
    free(diagPsiold);
    free(Psiinv_Lambda);
    free(Psiinv_Lambda_Minv);
    free(B);
    free(D);
    free(Psiinvhalf_Lambda);
    
    
}


void diagPsiupdate_C(int p, int m, int N, int ex_Npflag, int ex_dimx0, int ex_dimS, double ex_Lambda[p*m], double ex_diagPsi[p], double ex_Phi[m*m], double ex_tol1[1],  int ex_maxcount1[1], double ex_tolPsi[1], double ex_eta[1], double ex_S[ex_dimS], double ex_diagS[p], double ex_X[ex_dimx0],double ans[p],double A[m*m],double B[m*p], int corr_factor[1], double ex_zita[1])
//ans:Lambda,Psi,
{
    int pm = p*m;
    int mm = m*m;
    int t1=0;
    int blocksize=64;
    int mblocksize=m*blocksize;
    int one=1;
    int info=0;
    int i=0;
    int j=0;
    int l=0;
    int l0=0;
    
    double sum_diagPsisabun1=0.0;
    double ALambdanewA_forPsi=0.0;
    double BLambdanew_forPsi=0.0;
    double logF=0.0;
    double penalty0=0.0;
    
    double *Lambdanew;
    double *diagPsinew;
    double *diagPsiinv;
    double *diagPsiold;
    double *diagPsisabun;
    double *Psiinv_Lambda;
    double *Psiinv_Lambda_Minv;
    double *D;
    double M[mm];
    double Mtemp[mm];
    double MA[mm];
    int ipiv_M[m];
    double Minv[mm];
    double *Psiinvhalf_Lambda;
    double diagm[mm];
    double work1[mblocksize];
    double eigenMvalues[m];
    double eigenPhivalues[m];
    double Phinew[mm];
    double Phiinv[mm];
    double Phitemp2[mm];
    
    
    
    //malloc
    Lambdanew = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
    if( Lambdanew == NULL ) {
        free( Lambdanew );
        error( "memory allocation is failed\n" );
    }
    diagPsinew = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
    if( diagPsinew == NULL ) {
        free( diagPsinew );
        error( "memory allocation is failed\n" );
    }
    diagPsiinv = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
    if( diagPsiinv == NULL ) {
        free( diagPsiinv );
        error( "memory allocation is failed\n" );
    }
    diagPsiold = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
    if( diagPsiold == NULL ) {
        free( diagPsiold );
        error( "memory allocation is failed\n" );
    }
    diagPsisabun = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
    if( diagPsisabun == NULL ) {
        free( diagPsisabun );
        error( "memory allocation is failed\n" );
    }
    Psiinv_Lambda = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
    if( Psiinv_Lambda == NULL ) {
        free( Psiinv_Lambda );
        error( "memory allocation is failed\n" );
    }
    Psiinv_Lambda_Minv = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
    if( Psiinv_Lambda_Minv == NULL ) {
        free( Psiinv_Lambda_Minv );
        error( "memory allocation is failed\n" );
    }
    Psiinvhalf_Lambda = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
    if( Psiinvhalf_Lambda == NULL ) {
        free( Psiinvhalf_Lambda );
        error( "memory allocation is failed\n" );
    }
    D = (double *)malloc( sizeof( double ) * (N*m) );	// memory allocation
    if( D == NULL ) {
        free( D );
        error( "memory allocation is failed\n" );
    }
    
    
    
    
    //initialize
    
    for( i=0; i<p*m; i++){
        Lambdanew[i] = ex_Lambda[i];
        Psiinv_Lambda_Minv[i] = 0.0;
        B[i] = 0.0;
        Psiinvhalf_Lambda[i] = 0.0;
        Psiinv_Lambda[i] = 0.0;
    }
    
    for( i=0; i<m*m; i++){
        Phinew[i] = ex_Phi[i];
        Phiinv[i] = 0.0;
        Phitemp2[i] = 0.0;
    }
    
    
    
    for( i=0; i<mblocksize; i++){
        work1[i] = 0.0;
    }
    
    
    for( i=0; i<m*m; i++){
        A[i] = 0.0;
        M[i] = 0.0;
        Mtemp[i] = 0.0;
        Minv[i] = 0.0;
        diagm[i] = 0.0;
    }
    
    for(i=0; i<m; i++){
        diagm[i*m + i] = 1.0;
    }
    
    
    for( i=0; i<m*N; i++){
        D[i] = 0.0;
    }
    
    
    for( i=0; i<p; i++){
        diagPsinew[i] = ex_diagPsi[i];
        diagPsiinv[i] = 1.0 / diagPsinew[i];
        diagPsiold[i] = 0.0;
        diagPsisabun[i] = 10.0;
    }
    
    
    for( i=0; i<m; i++){
        ipiv_M[i] = 0;
        eigenPhivalues[i]=0.0;
    }
    
    sum_diagPsisabun1 = F77_CALL(dnrm2)(&p, diagPsisabun, &one);
    
    
    
    //Start
    
    t1=0;
    while(sum_diagPsisabun1 > ex_tol1[0] && t1 < ex_maxcount1[0]){
        t1=t1+1;
        F77_CALL(dcopy)(&p, diagPsinew, &one, diagPsiold, &one);
        F77_CALL(dcopy)(&pm, Lambdanew, &one, Psiinvhalf_Lambda, &one);
        F77_CALL(dcopy)(&pm, Lambdanew, &one, Psiinv_Lambda, &one);
        for( i=0; i<p; i++){
            diagPsiinv[i] = 1.0 / diagPsinew[i];
            for( j=0; j<m; j++){
                Psiinvhalf_Lambda[i+j*p] = Psiinvhalf_Lambda[i+j*p] * sqrt(diagPsiinv[i]);
                Psiinv_Lambda[i+j*p] = Psiinv_Lambda[i+j*p] * diagPsiinv[i];
            }
        }
        
        //M <- t(Lambda) * Psi^(-1) * Lambda + diagm
        if(corr_factor[0]==1 && m>1) {
            F77_CALL(dcopy)(&mm, Phinew, &one, Phiinv, &one);
            F77_CALL(dgetrf)(&m, &m, Phiinv, &m, ipiv_M,  &info);
            F77_CALL(dgetri)(&m, Phiinv, &m, ipiv_M, work1, &mblocksize, &info);
            F77_CALL(dcopy)(&mm, Phiinv, &one, M, &one);
        }else{
            F77_CALL(dcopy)(&mm, diagm, &one, M, &one);
        }
        
        F77_CALL(dsyrk)(&UPPER, &TRANST, &m, &p, &alphaOne, Psiinvhalf_Lambda, &p, &betaOne, M, &m);
        
        for( i=0; i<m; i++){
            for( j=0; j<m; j++){
                M[j+i*m] = M[i+j*m];
            }
        }
        
        //M^(-1)
        F77_CALL(dcopy)(&mm, M, &one, Minv, &one);
        F77_CALL(dgetrf)(&m, &m, Minv, &m, ipiv_M,  &info);
        F77_CALL(dgetri)(&m, Minv, &m, ipiv_M, work1, &mblocksize, &info);
        
        //Psiinv_Lambda_Minv
        F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, Psiinv_Lambda, &p, Minv, &m, &betaZero, Psiinv_Lambda_Minv, &p);
        if(ex_Npflag==1){
            //B
            F77_CALL(dgemm)(&TRANST, &TRANSN, &m, &p, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_S, &p, &betaZero, B, &m);
            //A
            F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
            F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &p, &alphaOne, B, &m, Psiinv_Lambda_Minv, &p, &betaOne, A, &m);
        }else{
            //D
            F77_CALL(dgemm)(&TRANST, &TRANST, &m, &N, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_X, &N, &betaZero, D, &m);
            //B
            F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &p, &N, &alphaOne, D, &m, ex_X, &N, &betaZero, B, &m);
            //A
            F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
            F77_CALL(dgemm)(&TRANSN, &TRANST, &m, &m, &N, &alphaOne, D, &m, D, &m, &betaOne, A, &m);
        }
        
        
        //update Psi
        for( i=0; i<p; i++){
            ALambdanewA_forPsi=0.0;
            BLambdanew_forPsi=0.0;
            for(l=0; l<m; l++){
                BLambdanew_forPsi=BLambdanew_forPsi+ B[i*m+l] * Lambdanew[p*l + i];
                for(l0=0; l0<m; l0++){
                    ALambdanewA_forPsi = ALambdanewA_forPsi + A[l*m+l0] * Lambdanew[p*l0 + i] *  Lambdanew[p*l + i];
                }
            }
            diagPsinew[i] = ex_diagS[i] - 2 * BLambdanew_forPsi + ALambdanewA_forPsi + ex_diagS[i]*ex_eta[0];
            if(diagPsinew[i] < ex_tolPsi[0]) diagPsinew[i] = ex_tolPsi[0];
        }
        
        //calulation of distance between old Psi and new Psi
        for( i=0; i<p; i++){
            diagPsisabun[i] = diagPsinew[i] - diagPsiold[i];
        }
        sum_diagPsisabun1 = F77_CALL(dnrm2)(&p, diagPsisabun, &one);
    }
    
    
    //calulation of log-likelihood
    //M*A
    F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &m, &alphaOne, M, &m, A, &m, &betaZero, MA, &m);
    //det(M)
    F77_CALL(dcopy)(&mm, M, &one, Mtemp, &one);
    F77_CALL(dsyev)(&jobzN, &UPPER, &m, Mtemp, &m, eigenMvalues, work1, &mblocksize, &info);
    //det(Phi)
    if(corr_factor[0]==1){
        F77_CALL(dcopy)(&mm, Phinew, &one, Phitemp2, &one);
        F77_CALL(dsyev)(&jobzN, &UPPER, &m, Phitemp2, &m, eigenPhivalues, work1, &mblocksize, &info);
    }
    
    logF = 0.0;
    for( i=0; i<p; i++){
        logF = logF + log(diagPsinew[i]) + ex_diagS[i] / diagPsinew[i]  ;
    }
    
    if(corr_factor[0]==1){
        for( i=0; i<m; i++){
            logF = logF + log(eigenMvalues[i]) + log(eigenPhivalues[i])  - ( MA[i+i*m] - 1.0) ;
        }
    }
    
    if(corr_factor[0]==0){
        for( i=0; i<m; i++){
            logF = logF + log(eigenMvalues[i]) - ( MA[i+i*m] - 1.0) ;
        }
    }
    
    
    
    penalty0=0.0;
    for( i=0; i<p; i++){
        penalty0 = penalty0 + ex_eta[0] * ex_diagS[i] * diagPsiinv[i] / 2.0 ;
    }
    
    //output answer
    for( i=0; i<p; i++){
        ans[i] = diagPsinew[i];
    }
    
    free(Lambdanew);
    free(diagPsinew);
    free(diagPsiinv);
    free(diagPsiold);
    free(diagPsisabun);
    free(Psiinv_Lambda);
    free(Psiinv_Lambda_Minv);
    free(D);
    free(Psiinvhalf_Lambda);
    
}



/*Calculate GFI (this program works only when p<500)*/
void GFI_C(int p, int m, int ex_dimS, double ex_Lambda[p*m], double ex_diagPsi[p], double ex_Phi[m*m], double ex_S[ex_dimS], double ans[1]){
    
    double *Sigma;
    double *Sigmainv;
    double *SigmainvS;
    int i;
    int pp=p*p;
    int pm=p*m;
    int one=1;
    int *ipiv_Sigma;
    int blocksize=64;
    int pblocksize=p*blocksize;
    int info=0;
    double *LambdaPhi;
    double *work1;
    
    double GFI0=0.0;  //tr(SigmainvS^2)
    double GFI1=0.0;  //tr(SigmainvS)
    
    
    
    //malloc
    Sigma = (double *)malloc( sizeof( double ) * (ex_dimS) );	// memory allocation
    if( Sigma == NULL ) {
        free( Sigma );
        error( "memory allocation is failed\n" );
    }
    
    Sigmainv = (double *)malloc( sizeof( double ) * (ex_dimS) );	// memory allocation
    if( Sigmainv == NULL ) {
        free( Sigmainv );
        error( "memory allocation is failed\n" );
    }
    
    SigmainvS = (double *)malloc( sizeof( double ) * (ex_dimS) );	// memory allocation
    if( SigmainvS == NULL ) {
        free( SigmainvS );
        error( "memory allocation is failed\n" );
    }
    
    ipiv_Sigma = (int *)malloc( sizeof( int ) * (p) );	// memory allocation
    if( ipiv_Sigma == NULL ) {
        free( ipiv_Sigma );
        error( "memory allocation is failed\n" );
    }
    
    LambdaPhi = (double *)malloc( sizeof( double ) * (p*m) );	// memory allocation
    if( LambdaPhi == NULL ) {
        free( LambdaPhi );
        error( "memory allocation is failed\n" );
    }
    
    work1 = (double *)malloc( sizeof( double ) * (pblocksize) );	// memory allocation
    if( work1 == NULL ) {
        free( work1 );
        error( "memory allocation is failed\n" );
    }
    
    
    
    for( i=0; i<ex_dimS; i++){
        Sigma[i] = 0.0;
        Sigmainv[i] = 0.0;
        SigmainvS[i] = 0.0;
    }
    
    
    for( i=0; i<p; i++){
        ipiv_Sigma[i] = 0;
    }
    
    for( i=0; i<pm; i++){
        LambdaPhi[i] = 0.0;
    }
    
    for( i=0; i<pblocksize; i++){
        work1[i] = 0.0;
    }
    
    
    //Sigma=Lambda Phi Lambda^T+Psi
    F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, ex_Lambda, &p, ex_Phi, &m, &betaZero, LambdaPhi, &p);
    F77_CALL(dgemm)(&TRANSN, &TRANST, &p, &p, &m, &alphaOne, LambdaPhi, &p, ex_Lambda, &p, &betaZero, Sigma, &p);
    for( i=0; i<p; i++){
        Sigma[i+i*p] = Sigma[i+i*p] + ex_diagPsi[i];
    }
    
    //Sigmainv
    F77_CALL(dcopy)(&pp, Sigma, &one, Sigmainv, &one);
    F77_CALL(dgetrf)(&p, &p, Sigmainv, &p, ipiv_Sigma,  &info);
    F77_CALL(dgetri)(&p, Sigmainv, &p, ipiv_Sigma, work1, &pblocksize, &info);
    
    //(Sigmainv%*%S)
    F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &p, &p, &alphaOne, Sigmainv, &p, ex_S, &p, &betaZero, SigmainvS, &p);
    
    GFI0=0.0;
    for( i=0; i<ex_dimS; i++){
        GFI0 = GFI0 + SigmainvS[i]*SigmainvS[i];
    }
    GFI1=0.0;
    for( i=0; i<p; i++){
        GFI1 = GFI1 + SigmainvS[i+i*p];
    }
    
    ans[0] = 1.0 - (GFI0 - 2.0 * GFI1 + (double) p) / GFI0;
    //ans = GFI;
    
    
    
    free( Sigma );
    free( Sigmainv );
    free( SigmainvS );
    free( ipiv_Sigma );
    free( LambdaPhi );
    free( work1 );
}


/* Update Lambda, Psi and Phi using the random initial values  */
void LambdaPsiPhiupdate_RANDOM_C(int p, int m, int N, int ex_Npflag, int ex_dimx0, int ex_dimS, int initialiter, double ex_tol1[1], double ex_tol2[1], double ex_tol3[1],  int ex_maxcountinitial[1], int ex_maxcount2[1], int ex_maxcount_phi[1], double ex_eta[1], double ex_tolPsi[1], double ex_rho[1], double ex_gamma[1], double ex_S[ex_dimS], double ex_diagS[p], double ex_X[ex_dimx0], int ex_flagpenalty[1], int ex_fixindex[p*m], int ex_maxcount1[1], double ans[p*m+p+m*m+3], int corr_factor[1], double ex_zita[1], int ex_failflag[3], int ex_omp[1], int ex_omp_num[1])
{
    
    
    int pm = p*m;
    int mm = m*m;
    int i = 0;
    int initiali = 0;
    int one = 1;
    int three = 3;
    double penalized_likelihood_kouho = 1000000000000000.0;
    double temp = 0.0;
    
    
    double *ans0;
    double *Lambda;
    double *Lambda0;
    double *diagPsi;
    double *diagPsi0;
    double Phi[m*m];
    double Phi0[m*m];
    double penalized_likelihood[3];
    double penalized_likelihood0[3];
    
    int flagpenalty2=1;
    int failflag_dummy[3]={0,0,0};
    
    //malloc
    ans0 = (double *)malloc( sizeof( double ) * (p*m+p+m*m+3) );	// memory allocation
    if( ans0 == NULL ) {
        free( ans0 );
        error( "memory allocation is failed\n" );
    }
    
    Lambda = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
    if( Lambda == NULL ) {
        free( Lambda );
        error( "memory allocation is failed\n" );
    }
    
    Lambda0 = (double *)malloc( sizeof( double ) * (pm) );	// memory allocation
    if( Lambda0 == NULL ) {
        free( Lambda0 );
        error( "memory allocation is failed\n" );
    }
    
    diagPsi = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
    if( diagPsi == NULL ) {
        free( diagPsi );
        error( "memory allocation is failed\n" );
    }
    
    diagPsi0 = (double *)malloc( sizeof( double ) * (p) );	// memory allocation
    if( diagPsi0 == NULL ) {
        free( diagPsi0 );
        error( "memory allocation is failed\n" );
    }
    
    
    
    //initialize
    for( i=0; i<p*m; i++){
        Lambda[i] = 0.0;
        Lambda0[i] = 0.0;
    }
    for( i=0; i<m*m; i++){
        Phi[i] = 0.0;
    }
    for( i=0; i<m; i++){
        Phi[i+i*m] = 1.0;
    }
    for( i=0; i<p; i++){
        diagPsi[i] = 1.0;
    }
    for( i=0; i<3; i++){
        penalized_likelihood[i] = 0.0;
        penalized_likelihood0[i] = 0.0;
    }
    
    
    
    // implement EM with coordinate descent for "initialiter" initial values.
    for( initiali=0; initiali<initialiter; initiali++){
        
        for( i=0; i<p*m; i++){
            GetRNGstate();
            Lambda0[i] = unif_rand();
            PutRNGstate();
        }
        for( i=0; i<p; i++){
            GetRNGstate();
            diagPsi0[i] = unif_rand();
            PutRNGstate();
        }
        for( i=0; i<m*m; i++){
            Phi0[i] = 0.0;
        }
        for( i=0; i<m; i++){
            Phi0[i+i*m] = 1.0;
        }
        
        // Update Parameters (if ex_maxcountinitial is small, the number of iteration of EM is small)
        LambdaPsiPhiupdate_C(p,m,N,ex_Npflag, ex_dimx0, ex_dimS, Lambda0, diagPsi0, Phi0, ex_tol1, ex_tol2, ex_tol3, ex_maxcountinitial, ex_maxcount2, ex_maxcount_phi, ex_eta, ex_tolPsi, ex_rho, ex_gamma, ex_S, ex_diagS, ex_X, &flagpenalty2, ex_fixindex, ans0, corr_factor, ex_zita, failflag_dummy, ex_omp, ex_omp_num);
        for(i=0; i<p*m;i++){
            temp = ans0[i];
            Lambda0[i] =temp;
        }
        for(i=0; i<p;i++){
            temp = ans0[i+p*m];
            diagPsi0[i] =temp;
        }
        for(i=0; i<m*m;i++){
            temp = ans0[i+p*m+p];
            Phi0[i]=temp;
        }
        
        /*
         // Update Parameters (if ex_maxcountinitial is small, the number of iteration of EM is small)
         LambdaPsiPhiupdate_C(p,m,N,ex_Npflag, ex_dimx0, ex_dimS, Lambda0, diagPsi0, Phi0, ex_tol1, ex_tol2, ex_tol3, ex_maxcountinitial, ex_maxcount2, ex_maxcount_phi, ex_eta, ex_tolPsi, ex_rho, ex_gamma, ex_S, ex_diagS, ex_X, ex_flagpenalty, ex_fixindex, ans0, corr_factor, ex_zita, failflag_dummy);
         for(i=0; i<p*m;i++){
         temp = ans0[i];
         Lambda0[i] =temp;
         }
         for(i=0; i<p;i++){
         temp = ans0[i+p*m];
         diagPsi0[i] =temp;
         }
         for(i=0; i<m*m;i++){
         temp = ans0[i+p*m+p];
         Phi0[i]=temp;
         }
         */
        
        for(i=0; i<3;i++){
            temp = ans0[i+p*m+p+m*m];
            penalized_likelihood0[i]=temp;
        }
        temp=ans0[2+p*m+p];
        if(penalized_likelihood_kouho > temp){
            penalized_likelihood_kouho = temp;
            F77_CALL(dcopy)(&pm, Lambda0, &one, Lambda, &one);
            F77_CALL(dcopy)(&p, diagPsi0, &one, diagPsi, &one);
            F77_CALL(dcopy)(&mm, Phi0, &one, Phi, &one);
            F77_CALL(dcopy)(&three, penalized_likelihood0, &one, penalized_likelihood, &one);
        }
    }
    
    /* Update Parameters with best initial values */
    LambdaPsiPhiupdate_C(p,m,N,ex_Npflag, ex_dimx0, ex_dimS,Lambda, diagPsi, Phi, ex_tol1, ex_tol2, ex_tol3, ex_maxcount1, ex_maxcount2, ex_maxcount_phi, ex_eta, ex_tolPsi, ex_rho, ex_gamma, ex_S, ex_diagS, ex_X, &flagpenalty2, ex_fixindex,ans,corr_factor, ex_zita, ex_failflag, ex_omp, ex_omp_num);
    
    
    free( ans0 );
    free( Lambda );
    free( Lambda0 );
    free( diagPsi );
    free( diagPsi0 );
}

/* set initial value of Lambda, Psi and rhomax*/
SEXP LamdPsiupdate_forinitial_C(SEXP ex_p, SEXP ex_N, SEXP ex_tol1, SEXP ex_tol2, SEXP ex_initial_iter, SEXP ex_maxcount1, SEXP ex_maxcount2, SEXP ex_max_countphi, SEXP ex_eta, SEXP ex_tolPsi, SEXP ex_gamma, SEXP ex_S, SEXP ex_diagS, SEXP ex_X, SEXP flagpenalty, SEXP ex_Lambdainitialij_for_rhoinitial, SEXP ex_Npflag, SEXP ex_dimx0, SEXP ex_dimS, SEXP ex_omp, SEXP ex_omp_num)
{
    
    SEXP ans0;
    SEXP ans;
    SEXP Lambda;
    SEXP diagPsi;
    SEXP Lambda0;
    SEXP diagPsi0;
    SEXP Phi;
    SEXP Lambda_ans;
    SEXP diagPsi_ans;
    SEXP penalized_likelihood;
    SEXP penalized_likelihood0;
    SEXP A;
    SEXP B;
    SEXP rho;
    
    int p = INTEGER(ex_p)[0];
    int m=1;
    int N = INTEGER(ex_N)[0];
    int initialiter = INTEGER(ex_initial_iter)[0];
    int i = 0;
    int initiali = 0;
    int pm=p*m;
    int three=3;
    int one=1;
    int l,j,l0 = 0;
    int failflag[3]={0,0,0};
    int fixindex[p];
    
    
    int Lambdaidamax=1;
    int corr_factor=0;
    
    
    double penalized_likelihood_kouho = 1000000000000000.0;
    double ALambda = 0.0;
    double initial_check_kouho = 0.0;
    double initial_check = 0.0;
    double x = 0.0;
    double temp = 0.0;
    double rho_zero[1] = {0.0};
    double zita=0.0;
    
    PROTECT(Lambda = allocMatrix(REALSXP,p,m));
    PROTECT(diagPsi = allocVector(REALSXP,p));
    PROTECT(Lambda0 = allocMatrix(REALSXP,p,m));
    PROTECT(diagPsi0 = allocVector(REALSXP,p));
    PROTECT(Phi = allocMatrix(REALSXP,m,m));
    PROTECT(Lambda_ans = allocMatrix(REALSXP,p,m));
    PROTECT(diagPsi_ans = allocVector(REALSXP,p));
    PROTECT(penalized_likelihood = allocVector(REALSXP, 3));
    PROTECT(penalized_likelihood0 = allocVector(REALSXP, 3));
    PROTECT(A = allocMatrix(REALSXP,m,m));
    PROTECT(B = allocMatrix(REALSXP, m,p));
    PROTECT(rho = allocVector(REALSXP, 1));
    PROTECT(ans0 = allocVector(REALSXP, p*m+p+m*m+3));
    PROTECT(ans = allocVector(VECSXP, 5));
    
    for( i=0; i<m*m; i++){
        REAL(A)[i] = 0.0;
        REAL(Phi)[i] = 0.0;
    }
    for( i=0; i<p*m; i++){
        REAL(B)[i] = 0.0;
    }
    
    
    for( i=0; i<p*m; i++){
        REAL(Lambda)[i] = 0.0;
        REAL(Lambda0)[i] = 0.0;
        REAL(Lambda_ans)[i] = 0.0;
    }
    for( i=0; i<p; i++){
        REAL(diagPsi)[i] = 1.0;
        REAL(diagPsi0)[i] = 1.0;
        REAL(diagPsi_ans)[i] = 1.0;
        fixindex[i]=0;
    }
    for( i=0; i<3; i++){
        REAL(penalized_likelihood)[i] = 0.0;
        REAL(penalized_likelihood0)[i] = 0.0;
    }
    
    REAL(rho)[0]=0.0;
    
    
    
    for( initiali=0; initiali<initialiter; initiali++){
        for(i=0; i<p;i++){
            GetRNGstate();
            REAL(Lambda0)[i] = unif_rand();
            PutRNGstate();
            GetRNGstate();
            REAL(diagPsi0)[i] = unif_rand();
            PutRNGstate();
        }
        LambdaPsiPhiupdate_C(p, m, N, INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], REAL(Lambda0), REAL(diagPsi0), REAL(Phi), REAL(ex_tol1), REAL(ex_tol2), REAL(ex_tol2), INTEGER(ex_maxcount1), INTEGER(ex_maxcount2), INTEGER(ex_max_countphi), REAL(ex_eta), REAL(ex_tolPsi), rho_zero, REAL(ex_gamma), REAL(ex_S), REAL(ex_diagS), REAL(ex_X), INTEGER(flagpenalty), fixindex, REAL(ans0), &corr_factor ,&zita, failflag, INTEGER(ex_omp), INTEGER(ex_omp_num));
        for(i=0; i<p*m;i++){
            temp= REAL(ans0)[i];
            REAL(Lambda0)[i] =temp;
        }
        for(i=0; i<p;i++){
            temp=  REAL(ans0)[i+p*m];
            REAL(diagPsi0)[i]=temp;
        }
        for(i=0; i<3;i++){
            temp=REAL(ans0)[i+p*m+p];
            REAL(penalized_likelihood0)[i] = temp;
        }
        temp=REAL(ans0)[2+p*m+p];
        if(penalized_likelihood_kouho > temp){
            penalized_likelihood_kouho = temp;
            F77_CALL(dcopy)(&pm, REAL(Lambda0), &one, REAL(Lambda), &one);
            F77_CALL(dcopy)(&p, REAL(diagPsi0), &one, REAL(diagPsi), &one);
            F77_CALL(dcopy)(&three, REAL(penalized_likelihood0), &one, REAL(penalized_likelihood), &one);
        }
    }
    
    //choose two elements of maximum absolute value of Lambda
    F77_CALL(dcopy)(&pm, REAL(Lambda), &one, REAL(Lambda0), &one);
    F77_CALL(dcopy)(&p, REAL(diagPsi), &one, REAL(diagPsi0), &one);
    F77_CALL(dcopy)(&three, REAL(penalized_likelihood), &one, REAL(penalized_likelihood0), &one);
    
    Lambdaidamax = F77_CALL(idamax)(&pm,REAL(Lambda0),&one);
    Lambdaidamax = Lambdaidamax-1;
    
    for(i=0; i<pm;i++){
        if(i!=Lambdaidamax) REAL(Lambda0)[i]=0.0;
    }
    F77_CALL(dcopy)(&pm, REAL(Lambda0), &one, REAL(Lambda), &one);
    
    //update diagPsi while fixing Lambda
    for( initiali=0; initiali<initialiter; initiali++){
        
        for( j=0; j<p*m; j++){
            REAL(Lambda)[j]=REAL(Lambda0)[j] * REAL(ex_Lambdainitialij_for_rhoinitial)[initiali];
        }
        
        //update psi(diagPsi0:old, diagPsi:new)
        diagPsiupdate_C(p, m, N, INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0],  REAL(Lambda), REAL(diagPsi0), REAL(Phi), REAL(ex_tol1), INTEGER(ex_maxcount1), REAL(ex_tolPsi), REAL(ex_eta), REAL(ex_S), REAL(ex_diagS), REAL(ex_X), REAL(diagPsi), REAL(A), REAL(B), &corr_factor ,&zita);
        
        initial_check_kouho = 0.0;
        for( j=0; j<p; j++){
            for( l=0; l<1; l++){
                ALambda=0.0;
                for(l0=0; l0<m; l0++){
                    if(l0 != l) ALambda = ALambda + REAL(A)[l*m+l0] * REAL(Lambda)[p*l0 + j];
                }
                if(Lambdaidamax!=j){
                    x = fabs( 1.0 / REAL(diagPsi)[j] * ( REAL(B)[j*m+l] - ALambda ) );
                    if(initial_check_kouho < x)  initial_check_kouho = x;
                }
            }
        }
        
        
        if(initial_check < initial_check_kouho){
            initial_check = initial_check_kouho;
            F77_CALL(dcopy)(&pm, REAL(Lambda), &one, REAL(Lambda_ans), &one);
            F77_CALL(dcopy)(&p, REAL(diagPsi), &one, REAL(diagPsi_ans), &one);
            REAL(rho)[0]=initial_check;
        }
    }
    
    SET_VECTOR_ELT(ans, 0,Lambda_ans);
    SET_VECTOR_ELT(ans, 1,diagPsi_ans);
    SET_VECTOR_ELT(ans, 2,rho);
    
    
    UNPROTECT(14);
    
    return(ans);
}



SEXP fancmain(SEXP ex_Lambda_ini, SEXP ex_diagPsi_ini, SEXP ex_N, SEXP ex_tol1, SEXP ex_tol2, SEXP ex_tol3, SEXP ex_initial_iter, SEXP ex_maxcount1, SEXP ex_maxcount2, SEXP ex_maxcount_phi, SEXP ex_maxcountinitial, SEXP ex_eta, SEXP ex_tolPsi, SEXP ex_rhomatrix, SEXP ex_gamma_vec, SEXP ex_X, SEXP ex_diagS, SEXP ex_S, SEXP ex_Npflag, SEXP ex_dimx0, SEXP ex_dimS, SEXP ex_pmax_for_S, SEXP ex_corr_factor, SEXP ex_zita, SEXP ex_all_random_start, SEXP ex_trace, SEXP ex_omp, SEXP ex_omp_num, SEXP pBar){
    SEXP Lambdaold;
    SEXP Lambdanew;
    SEXP diagPsiold;
    SEXP diagPsinew;
    SEXP Phiold;
    SEXP Phinew;
    SEXP Lambda_index;
    SEXP parameterupdate;
    SEXP parameterupdate_new;
    SEXP ans_diagPsi;
    SEXP ans_Phi;
    SEXP ans_AIC;
    SEXP ans_BIC;
    SEXP ans_CAIC;
    SEXP ans_DF;
    SEXP ans_GFI;
    SEXP ans_logF2c;
    SEXP ans_failflag;
    SEXP ans_nonzero;
    SEXP ans;
    
    SEXP utilsPackage, percentComplete, percenttemp;
    int N = INTEGER(ex_N)[0];
    int p = INTEGER(GET_DIM(ex_Lambda_ini))[0];
    int m = INTEGER(GET_DIM(ex_Lambda_ini))[1];
    int i,j,l,i0,j0;
    int pm = p*m;
    int mm = m*m;
    int stepgamma,steprho;
    int STEPrho =  INTEGER(GET_DIM(ex_rhomatrix))[0];
    int STEPgamma =  INTEGER(GET_DIM(ex_rhomatrix))[1];
    int initialiter = INTEGER(ex_initial_iter)[0];
    
    int one=1;
    int m0=1;
    int m0new=1;
    
    
    int flagpenalty=1;
    int flag_m=1;
    int flag_NewOld=1; //1:new
    
    
    
    int corr_factor_zero[1]; // no update of factor correlation
    corr_factor_zero[0]=0;
    
    double sumabsLambda=0.0;
    double sumabsm0=0.0;
    double temp=0.0;
    double tempnew=0.0;
    double DF=0.0;
    double nonzero=0.0;
    double GFI[1];
    double DF_factorloadings_lasso[STEPrho];
    double rho;
    double gamma;
    
    
    
    int unit_sparseloadings=1000;
    int n_sparseloadings=0;
    int n0_sparseloadings=unit_sparseloadings;
    int *sparseloading_i;
    double *sparseloading_x;
    int n_sparsenet_update=0;
    
    //malloc
    //scanf( "%d", &n0_sparseloadings );
    sparseloading_x = (double *)malloc( sizeof( double ) * n0_sparseloadings );	// memory allocation
    sparseloading_i = (int *)malloc( sizeof( int ) * n0_sparseloadings );	// memory allocation
    if( sparseloading_x == NULL ) {
        free( sparseloading_i );
        free( sparseloading_x );
        error( "memory allocation is failed\n" );
    }
    if( sparseloading_i == NULL ) {
        free( sparseloading_i );
        free( sparseloading_x );
        error( "memory allocation is failed\n" );
    }
    
    
    PROTECT(Lambdaold = allocMatrix(REALSXP, p,m));
    PROTECT(Lambdanew = allocMatrix(REALSXP, p,m));
    PROTECT(diagPsiold = allocVector(REALSXP, p));
    PROTECT(diagPsinew = allocVector(REALSXP, p));
    PROTECT(Phiold = allocMatrix(REALSXP, m,m));
    PROTECT(Phinew = allocMatrix(REALSXP, m,m));
    PROTECT(Lambda_index = allocMatrix(INTSXP, p,m));
    PROTECT(parameterupdate = allocVector(REALSXP, p*m+p+m*m+3));
    PROTECT(parameterupdate_new = allocVector(REALSXP, p*m+p+m*m+3));
    PROTECT(ans_AIC = allocMatrix(REALSXP, STEPrho,STEPgamma));
    PROTECT(ans_BIC = allocMatrix(REALSXP, STEPrho,STEPgamma));
    PROTECT(ans_CAIC = allocMatrix(REALSXP, STEPrho,STEPgamma));
    PROTECT(ans_DF = allocMatrix(REALSXP, STEPrho,STEPgamma));
    PROTECT(ans_nonzero = allocMatrix(REALSXP, STEPrho,STEPgamma));
    PROTECT(ans_GFI = allocMatrix(REALSXP, STEPrho,STEPgamma));
    PROTECT(ans_diagPsi = allocMatrix(REALSXP, p,STEPrho*STEPgamma));
    PROTECT(ans_logF2c = allocMatrix(REALSXP, 3,STEPrho*STEPgamma));
    PROTECT(ans_Phi = allocMatrix(REALSXP, m, m*STEPrho*STEPgamma));
    PROTECT(ans_failflag = allocVector(INTSXP, 3)); //em,cd,phi
    PROTECT(ans = allocVector(VECSXP, 12));
    PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
    PROTECT(percentComplete = allocVector(REALSXP,1));
    
    
    //initialize
    //prof_init(10);
    for(i=0;i<p*m;i++){
        REAL(Lambdaold)[i]=REAL(ex_Lambda_ini)[i];
        REAL(Lambdanew)[i]=0.0;
        INTEGER(Lambda_index)[i]=1;
    }
    for(i=0;i<p;i++){
        REAL(diagPsiold)[i]=REAL(ex_diagPsi_ini)[i];
        REAL(diagPsinew)[i]=0.0;
        INTEGER(Lambda_index)[i]=0;
    }
    for(i=0;i<m*m;i++){
        REAL(Phiold)[i]=0.0;
        REAL(Phinew)[i]=0.0;
    }
    for(i=0;i<m;i++){
        REAL(Phiold)[i+i*m]=1.0;
        REAL(Phinew)[i+i*m]=1.0;
    }
    
    for(i=0;i<STEPrho*STEPgamma;i++){
        REAL(ans_AIC)[i]=0.0;
        REAL(ans_BIC)[i]=0.0;
        REAL(ans_CAIC)[i]=0.0;
        REAL(ans_DF)[i]=0.0;
        REAL(ans_nonzero)[i]=0.0;
        REAL(ans_GFI)[i]=0.0;
    }
    for(i=0;i<STEPrho;i++){
        DF_factorloadings_lasso[i] = 0.0;
    }
    for(i=0;i<p*STEPrho*STEPgamma;i++){
        REAL(ans_diagPsi)[i]=0.0;
    }
    for(i=0;i<3*STEPrho*STEPgamma;i++){
        REAL(ans_logF2c)[i]=0.0;
    }
    for(i=0;i<m*m*STEPrho*STEPgamma;i++){
        REAL(ans_Phi)[i]=0.0;
    }
    for(i=0;i< p*m+p+mm+3 ;i++){
        REAL(parameterupdate)[i]=0.0;
        REAL(parameterupdate_new)[i]=0.0;
    }
    for(i=0;i< 3 ;i++){
        INTEGER(ans_failflag)[i]=0;
    }
    
    gamma=0.0;
    rho=0.0;
    
    //start!!
    for(stepgamma=0; stepgamma<STEPgamma; stepgamma++){
        gamma = REAL(ex_gamma_vec)[stepgamma];
        for(steprho=0; steprho<STEPrho; steprho++){
            rho = REAL(ex_rhomatrix)[steprho + stepgamma*STEPrho];
            //initialize Lambda (warm start)
            if(stepgamma==0  && steprho>0){
                for(i=0;i<p*m;i++){
                    REAL(Lambdaold)[i] = REAL(Lambdanew)[i];
                }
                for(i=0;i<p;i++){
                    REAL(diagPsiold)[i] = REAL(ans_diagPsi)[i+(steprho-1)*p];
                }
                for(i=0;i<mm;i++){
                    REAL(Phiold)[i] = REAL(ans_Phi)[i+(steprho-1)*m*m];
                }
            }
            
            if(stepgamma >= 1){
                for(i=0;i<p*m;i++){
                    if(sparseloading_i[n_sparsenet_update]==i+(stepgamma-1)*STEPrho*p*m + steprho*p*m){
                        REAL(Lambdaold)[i] = sparseloading_x[n_sparsenet_update];
                        n_sparsenet_update = n_sparsenet_update + 1;
                    }else{
                        REAL(Lambdaold)[i] = 0.0;
                    }
                }
                for(i=0;i<p;i++){
                    REAL(diagPsiold)[i] = REAL(ans_diagPsi)[i+(stepgamma-1)*STEPrho*p + steprho*p];
                }
                for(i=0;i<m*m;i++){
                    REAL(Phiold)[i] = REAL(ans_Phi)[i+(stepgamma-1)*STEPrho*m*m + steprho*m*m];
                }
            }
            
            
            //if all elements are closed to zero, set initial values (gamma=Inf)
            if(stepgamma == 0){
                sumabsLambda = F77_CALL(dnrm2)(&pm, REAL(Lambdaold), &one);
                if(sumabsLambda < REAL(ex_tol1)[0]){
                    F77_CALL(dcopy)(&pm, REAL(ex_Lambda_ini), &one, REAL(Lambdaold), &one);
                    F77_CALL(dcopy)(&p, REAL(ex_diagPsi_ini), &one, REAL(diagPsiold), &one);
                    for(i=0;i<mm;i++){
                        REAL(Phiold)[i]=0.0;
                    }
                    for(i=0;i<m;i++){
                        REAL(Phiold)[i+i*m]=1.0;
                    }
                }
            }
            
            //update parameter
            //compute m0
            m0=0;
            for(j=0;j<m;j++){
                sumabsm0=0.0;
                for(i=0;i<p;i++){
                    sumabsm0 = sumabsm0 + fabs(REAL(Lambdaold)[i+j*p]);
                }
                if(sumabsm0 > REAL(ex_tol1)[0]) m0=m0+1;
            }
            
            //set Lambda_index
            for(i=0;i<p*m;i++){
                INTEGER(Lambda_index)[i]=1;
            }
            for(i=0;i<p*m0;i++){
                INTEGER(Lambda_index)[i] = 0;
            }
            //Rprintf("%d\n" , m0);
            
            for(i=0;i< p*m+p+mm+3 ;i++){
                REAL(parameterupdate)[i]=0.0;
            }
            if(INTEGER(ex_all_random_start)[0]==1){
                LambdaPsiPhiupdate_RANDOM_C(p,m,N,INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], initialiter, REAL(ex_tol1), REAL(ex_tol2), REAL(ex_tol3), INTEGER(ex_maxcount1), INTEGER(ex_maxcount2) , INTEGER(ex_maxcount_phi), REAL(ex_eta), REAL(ex_tolPsi), &rho, &gamma, REAL(ex_S), REAL(ex_diagS), REAL(ex_X), &flagpenalty, INTEGER(Lambda_index), INTEGER(ex_maxcountinitial),REAL(parameterupdate), INTEGER(ex_corr_factor), REAL(ex_zita), INTEGER(ans_failflag), INTEGER(ex_omp), INTEGER(ex_omp_num));
            }else{
                LambdaPsiPhiupdate_C(p,m,N, INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], REAL(Lambdaold), REAL(diagPsiold), REAL(Phiold), REAL(ex_tol1),REAL(ex_tol2),REAL(ex_tol3),INTEGER(ex_maxcount1),INTEGER(ex_maxcount2), INTEGER(ex_maxcount_phi), REAL(ex_eta), REAL(ex_tolPsi), &rho, &gamma, REAL(ex_S), REAL(ex_diagS), REAL(ex_X), &flagpenalty, INTEGER(Lambda_index),REAL(parameterupdate), INTEGER(ex_corr_factor), REAL(ex_zita), INTEGER(ans_failflag), INTEGER(ex_omp), INTEGER(ex_omp_num));
            }
            //compute m0
            m0=0;
            for(j=0;j<m;j++){
                sumabsm0=0.0;
                for(i=0;i<p;i++){
                    sumabsm0 = sumabsm0 + fabs(REAL(parameterupdate)[i+j*p]);
                }
                if(sumabsm0 > REAL(ex_tol1)[0]) m0=m0+1;
            }
            
            flag_m=0;
            flag_NewOld = 1;
            
            //if(m0==0) flag_m=1;
            if(m0==m) flag_m=1;
            if(stepgamma!=0) flag_m=1;
            
            if(flag_m==1) flag_NewOld = 0;
            
            
            //check if the number of factors should be updated
            if(flag_m == 0){
                
                for(i=0;i<p*m;i++){
                    INTEGER(Lambda_index)[i]=1;
                }
                for(i=0;i<p*m;i++){
                    INTEGER(Lambda_index)[i] = 0;
                }
                
                for(i=0;i<p*m;i++){
                    temp= REAL(parameterupdate)[i];
                    REAL(Lambdaold)[i] =temp;
                }
                for(i=0;i<p;i++){
                    temp=REAL(parameterupdate)[i+p*m];
                    REAL(diagPsiold)[i] =temp;
                }
                for(i=0;i<mm;i++){
                    temp=REAL(parameterupdate)[i+p*m+p];
                    REAL(Phiold)[i] =temp;
                }
                
                for(i=0;i< p*m+p+mm+3 ;i++){
                    REAL(parameterupdate_new)[i]=0.0;
                }
                
                
                //estimation of oblique structure without estimation of orthogonal structure (REJECTION!!)
                //LambdaPsiPhiupdate_RANDOM_C(p,m,N,INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], initialiter,REAL(Lambdaold), REAL(diagPsiold), REAL(Phiold), REAL(ex_tol1), REAL(ex_tol2), REAL(ex_tol3), INTEGER(ex_maxcount1), INTEGER(ex_maxcount2) , INTEGER(ex_maxcount_phi), REAL(ex_eta), REAL(ex_tolPsi), REAL(rho), REAL(gamma), REAL(ex_S), REAL(ex_diagS), REAL(ex_X), INTEGER(flagpenalty), INTEGER(Lambda_index),  INTEGER(m0SEXP), INTEGER(ex_maxcountinitial),REAL(parameterupdate_new), INTEGER(ex_corr_factor), REAL(ex_zita), INTEGER(ans_failflag));
                
                
                //orthogonal -> oblique
                // estimation of orthogonal structure
                LambdaPsiPhiupdate_RANDOM_C(p,m,N,INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], initialiter, REAL(ex_tol1), REAL(ex_tol2), REAL(ex_tol3), INTEGER(ex_maxcount1), INTEGER(ex_maxcount2) , INTEGER(ex_maxcount_phi), REAL(ex_eta), REAL(ex_tolPsi), &rho ,&gamma, REAL(ex_S), REAL(ex_diagS), REAL(ex_X), &flagpenalty, INTEGER(Lambda_index),  INTEGER(ex_maxcountinitial),REAL(parameterupdate_new), corr_factor_zero, REAL(ex_zita), INTEGER(ans_failflag), INTEGER(ex_omp), INTEGER(ex_omp_num));
                
                //compute m0new
                m0new=0;
                for(j=0;j<m;j++){
                    sumabsm0=0.0;
                    for(i=0;i<p;i++){
                        sumabsm0 = sumabsm0 + fabs(REAL(parameterupdate_new)[i+j*p]);
                    }
                    if(sumabsm0 > REAL(ex_tol1)[0]){
                        m0new=m0new+1;
                    }
                }
                
                tempnew=REAL(parameterupdate_new)[p*m+p+m*m+2];
                temp=REAL(parameterupdate)[p*m+p+m*m+2];
                
                // estimation of oblique structure
                if(m0new>m0 && tempnew <temp){
                    if(INTEGER(ex_corr_factor)[0]==1){
                        for(i=0;i<p*m;i++){
                            temp= REAL(parameterupdate_new)[i];
                            REAL(Lambdaold)[i] =temp;
                        }
                        for(i=0;i<p;i++){
                            temp=REAL(parameterupdate_new)[i+p*m];
                            REAL(diagPsiold)[i] =temp;
                        }
                        for(i=0;i<mm;i++){
                            temp=REAL(parameterupdate_new)[i+p*m+p];
                            REAL(Phiold)[i] =temp;
                        }
                        
                        for(i=0;i< p*m+p+mm+3 ;i++){
                            REAL(parameterupdate_new)[i]=0.0;
                        }
                        LambdaPsiPhiupdate_C(p,m,N, INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], REAL(Lambdaold), REAL(diagPsiold), REAL(Phiold), REAL(ex_tol1),REAL(ex_tol2), REAL(ex_tol3), INTEGER(ex_maxcount1),INTEGER(ex_maxcount2),  INTEGER(ex_maxcount_phi), REAL(ex_eta), REAL(ex_tolPsi), &rho ,&gamma , REAL(ex_S), REAL(ex_diagS), REAL(ex_X), &flagpenalty, INTEGER(Lambda_index),REAL(parameterupdate_new), INTEGER(ex_corr_factor), REAL(ex_zita), INTEGER(ans_failflag), INTEGER(ex_omp), INTEGER(ex_omp_num));
                    }
                }
                //orthogonal -> oblique(end)
                
                
                //compute m0new
                m0new=0;
                for(j=0;j<m;j++){
                    sumabsm0=0.0;
                    for(i=0;i<p;i++){
                        sumabsm0 = sumabsm0 + fabs(REAL(parameterupdate_new)[i+j*p]);
                    }
                    if(sumabsm0 > REAL(ex_tol1)[0]){
                        m0new=m0new+1;
                    }
                }
                
                
                //check updating if m0 <- m0new
                tempnew=REAL(parameterupdate_new)[p*m+p+m*m+2];
                temp=REAL(parameterupdate)[p*m+p+m*m+2];
                
                if(m0new == m || m0new <= m0){
                    flag_m = 1;
                    if(tempnew >temp){
                        flag_NewOld = 0;
                    }
                }else{
                    flag_NewOld = 1;
                    if(tempnew >temp){
                        flag_NewOld = 0;
                        flag_m = 1;
                    }
                }
                
                if(flag_m == 0){
                    for(i=0;i<p*m+p+m*m+3;i++){
                        REAL(parameterupdate)[i] = REAL(parameterupdate_new)[i];
                    }
                    //update m0
                    m0=m0new;
                }
            }
            //endwhile
            
            //copy the result
            for(i=0;i<p*m;i++){
                REAL(Lambdanew)[i]=0.0;
            }
            for(i=0;i<p;i++){
                REAL(diagPsinew)[i]=0.0;
            }
            for(i=0;i<mm;i++){
                REAL(Phinew)[i]=0.0;
            }
            
            if(flag_NewOld==1){
                l=0;
                for(j=0;j<m;j++){
                    for(i=0;i<p;i++){
                        i0 = i+l*p+steprho*p*m+STEPrho*stepgamma*p*m;
                        j0 = i+j*p;
                        temp =  REAL(parameterupdate_new)[j0];
                        REAL(Lambdanew)[j0]=temp;
                        
                        //realloc
                        if(temp!=0.0){
                            n_sparseloadings=n_sparseloadings+1;
                            if(n_sparseloadings - n0_sparseloadings==0){
                                unit_sparseloadings=unit_sparseloadings*2;
                                n0_sparseloadings=unit_sparseloadings;
                                //( "%d", &n0_sparseloadings );
                                sparseloading_x = (double *)realloc(sparseloading_x,  sizeof( double ) * n0_sparseloadings );	// memory allocation
                                sparseloading_i = (int *)realloc(sparseloading_i,  sizeof( int ) * n0_sparseloadings );	// memory allocation
                                if( sparseloading_x == NULL ) {
                                    free( sparseloading_i );
                                    free( sparseloading_x );
                                    error( "memory allocation is failed\n" );
                                }
                                if( sparseloading_i == NULL ) {
                                    free( sparseloading_i );
                                    free( sparseloading_x );
                                    error( "memory allocation is failed\n" );
                                }
                            }
                            sparseloading_i[n_sparseloadings-1] = i0;
                            sparseloading_x[n_sparseloadings-1] = temp;
                            
                        }
                        
                        
                        
                        
                    }
                    l=l+1;
                }
                for(i=0;i<p;i++){
                    temp = REAL(parameterupdate_new)[i+p*m];
                    i0=i+steprho*p+STEPrho*stepgamma*p;
                    REAL(ans_diagPsi)[i0]=temp;
                    REAL(diagPsinew)[i]=temp;
                }
                for(i=0;i<mm;i++){
                    temp = REAL(parameterupdate_new)[i+p*m+p];
                    i0=i+steprho*mm+STEPrho*stepgamma*mm;
                    REAL(ans_Phi)[i0]=temp;
                    REAL(Phinew)[i]=temp;
                }
                REAL(ans_logF2c)[0+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate_new)[p*m+p+m*m+0];
                REAL(ans_logF2c)[1+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate_new)[p*m+p+m*m+1];
                REAL(ans_logF2c)[2+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate_new)[p*m+p+m*m+2];
                
                
                //compute m0
                m0=0;
                for(j=0;j<m;j++){
                    sumabsm0=0.0;
                    for(i=0;i<p;i++){
                        sumabsm0 = sumabsm0 + fabs(REAL(parameterupdate_new)[i+j*p]);
                    }
                    if(sumabsm0 > REAL(ex_tol1)[0]) m0=m0+1;
                }
                
                
                if(stepgamma==0){
                    DF=0.0;
                    for(i=0;i<pm;i++){
                        if(fabs(REAL(parameterupdate_new)[i]) > REAL(ex_tol1)[0]) DF=DF+1.0;
                    }
                    DF_factorloadings_lasso[steprho] = DF;
                    DF=DF+ (double) p;
                    if(INTEGER(ex_corr_factor)[0]==1 && m0>1) DF=DF+(double) m0 *((double) m0-1) / 2.0;
                }else{
                    //DF=REAL(ans_DF)[steprho + STEPrho*(stepgamma-1)];
                    DF = DF_factorloadings_lasso[steprho];
                    DF=DF+ (double) p;
                    if(INTEGER(ex_corr_factor)[0]==1 && m0>1) DF=DF+(double) m0 *((double) m0-1) / 2.0;
                }
                
                //calculate nonzero factor loadings
                nonzero=0.0;
                for(i=0;i<pm;i++){
                    if(fabs(REAL(parameterupdate_new)[i]) > REAL(ex_tol1)[0]) nonzero=nonzero+1.0;
                }
                REAL(ans_nonzero)[steprho + STEPrho*stepgamma] = nonzero;
                
                REAL(ans_DF)[steprho + STEPrho*stepgamma]= DF;
                REAL(ans_AIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate_new)[p*m+p+m*m+0] + 2*DF;
                REAL(ans_BIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate_new)[p*m+p+m*m+0] + log(N)*DF;
                REAL(ans_CAIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate_new)[p*m+p+m*m+0] + (log(N)+1)*DF;
                
                GFI[0]=0.0;
                if(p<=INTEGER(ex_pmax_for_S)[0]) GFI_C(p, m, INTEGER(ex_dimS)[0], REAL(Lambdanew), REAL(diagPsinew), REAL(Phinew), REAL(ex_S), GFI);
                REAL(ans_GFI)[steprho + STEPrho*stepgamma]= GFI[0];
            }
            
            
            if(flag_NewOld==0){
                for(j=0;j<m;j++){
                    for(i=0;i<p;i++){
                        i0 = i+j*p+steprho*p*m+STEPrho*stepgamma*p*m;
                        j0 = i+j*p;
                        temp =  REAL(parameterupdate)[j0];
                        REAL(Lambdanew)[j0]=temp;
                        
                        
                        //realloc
                        if(temp!=0.0){
                            n_sparseloadings=n_sparseloadings+1;
                            if(n_sparseloadings - n0_sparseloadings==0){
                                unit_sparseloadings=unit_sparseloadings*2;
                                n0_sparseloadings=unit_sparseloadings;
                                //scanf( "%d", &n0_sparseloadings );
                                sparseloading_x = (double *)realloc(sparseloading_x,  sizeof( double ) * n0_sparseloadings );	// memory allocation
                                sparseloading_i = (int *)realloc(sparseloading_i,  sizeof( int ) * n0_sparseloadings );	//memory allocation
                                if( sparseloading_x == NULL ) {
                                    free( sparseloading_i );
                                    free( sparseloading_x );
                                    error( "memory allocation is failed\n" );
                                }
                                if( sparseloading_i == NULL ) {
                                    free( sparseloading_i );
                                    free( sparseloading_x );
                                    error( "memory allocation is failed\n" );
                                }
                            }
                            sparseloading_i[n_sparseloadings-1] = i0;
                            sparseloading_x[n_sparseloadings-1] = temp;
                            
                        }
                        
                        
                    }
                }
                
                for(i=0;i<p;i++){
                    temp = REAL(parameterupdate)[i+p*m];
                    i0=i+steprho*p+STEPrho*stepgamma*p;
                    REAL(ans_diagPsi)[i0]=temp;
                    REAL(diagPsinew)[i]=temp;
                }
                for(i=0;i<mm;i++){
                    temp = REAL(parameterupdate)[i+p*m+p];
                    i0=i+steprho*mm+STEPrho*stepgamma*mm;
                    REAL(ans_Phi)[i0]=temp;
                    REAL(Phinew)[i]=temp;
                }
                REAL(ans_logF2c)[0+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate)[p*m+p+m*m+0];
                REAL(ans_logF2c)[1+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate)[p*m+p+m*m+1];
                REAL(ans_logF2c)[2+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate)[p*m+p+m*m+2];
                
                /* compute m0*/
                m0=0;
                for(j=0;j<m;j++){
                    sumabsm0=0.0;
                    for(i=0;i<p;i++){
                        sumabsm0 = sumabsm0 + fabs(REAL(parameterupdate)[i+j*p]);
                    }
                    if(sumabsm0 > REAL(ex_tol1)[0]) m0=m0+1;
                }
                
                
                if(stepgamma==0){
                    DF=0.0;
                    for(i=0;i<pm;i++){
                        if(fabs(REAL(parameterupdate)[i]) > REAL(ex_tol1)[0]) DF=DF+1.0;
                    }
                    DF_factorloadings_lasso[steprho] = DF;
                    DF=DF+ (double) p;
                    if(INTEGER(ex_corr_factor)[0]==1 && m0>1) DF=DF+(double) m0 *((double) m0-1) / 2.0;
                }else{
                    //DF=REAL(ans_DF)[steprho + STEPrho*(stepgamma-1)];
                    DF = DF_factorloadings_lasso[steprho];
                    DF=DF+ (double) p;
                    if(INTEGER(ex_corr_factor)[0]==1 && m0>1) DF=DF+(double) m0 *((double) m0-1) / 2.0;
                }
                
                /* calculate nonzero factor loadings */
                nonzero=0.0;
                for(i=0;i<pm;i++){
                    if(fabs(REAL(parameterupdate)[i]) > REAL(ex_tol1)[0]) nonzero=nonzero+1.0;
                }
                REAL(ans_nonzero)[steprho + STEPrho*stepgamma] = nonzero;
                REAL(ans_DF)[steprho + STEPrho*stepgamma]= DF;
                REAL(ans_AIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate)[p*m+p+m*m+0] + 2*DF;
                REAL(ans_BIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate)[p*m+p+m*m+0] + log(N)*DF;
                REAL(ans_CAIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate)[p*m+p+m*m+0] + (log(N)+1)*DF;
                
                
                GFI[0]=0.0;
                if(p<=INTEGER(ex_pmax_for_S)[0]) GFI_C(p, m, INTEGER(ex_dimS)[0], REAL(Lambdanew), REAL(diagPsinew),  REAL(Phinew), REAL(ex_S), GFI);
                REAL(ans_GFI)[steprho + STEPrho*stepgamma]= GFI[0];
                
                
            }
            if(INTEGER(ex_trace)[0]==1){
                //Rprintf("The solution at gamma%drho%d has been computed. (%d/%d)\n" , stepgamma+1, steprho+1,stepgamma*STEPrho+steprho+1, STEPgamma*STEPrho);
                REAL(percentComplete)[0]=(double) (stepgamma*STEPrho+steprho+2);
                PROTECT(percenttemp=eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage));
                UNPROTECT(1);
            }
        }
    }
    
    
    
    SEXP ans_Lambda_i;
    SEXP ans_Lambda_x;
    PROTECT(ans_Lambda_i = allocVector(INTSXP, n_sparseloadings));
    PROTECT(ans_Lambda_x = allocVector(REALSXP, n_sparseloadings));
    for(i=0; i<n_sparseloadings;i++){
        INTEGER(ans_Lambda_i)[i]=sparseloading_i[i];
        REAL(ans_Lambda_x)[i]=sparseloading_x[i];
    }
    
    SET_VECTOR_ELT(ans, 0, ans_Lambda_i);
    SET_VECTOR_ELT(ans, 1, ans_Lambda_x);
    SET_VECTOR_ELT(ans, 2, ans_diagPsi);
    SET_VECTOR_ELT(ans, 3, ans_logF2c);
    SET_VECTOR_ELT(ans, 4, ans_AIC);
    SET_VECTOR_ELT(ans, 5, ans_BIC);
    SET_VECTOR_ELT(ans, 6, ans_CAIC);
    SET_VECTOR_ELT(ans, 7, ans_DF);
    SET_VECTOR_ELT(ans, 8, ans_GFI);
    SET_VECTOR_ELT(ans, 9, ans_Phi);
    SET_VECTOR_ELT(ans, 10, ans_failflag);
    SET_VECTOR_ELT(ans, 11, ans_nonzero);
    
    
    free( sparseloading_i );
    free( sparseloading_x );
    UNPROTECT(24);
    //prof_print();
    return(ans);
}



