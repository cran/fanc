//Ｒ用 Ｃ
//#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Parse.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>


//////////////////////////////////////////////////
//BLAS: y := alpha*A*x + beta*y
int gl_incx = 1;
int gl_incy = 1;
double alphaOne = 1.0;
double alphaminusOne = -1.0;
double betaOne = 1.0;
double betaZero = 0.0;
double Zero = 0.0;
char TRANSN = 'N';
char TRANST = 'T';
char UPPER = 'U';
char jobzN = 'N';



double sign_C(double z)
{
	double ans;
	if(z > 0) ans = 1;
	if(z == 0) ans = 0;
	if(z < 0) ans = -1;
	return(ans);
	
	
}



double S_func_MC_C(double z, double rho, double gamma)
{
	
	double ans;
	
	if(fabs(z) <= rho) ans = 0;
	else if(fabs(z) <= rho*gamma) ans = sign_C(z) * (fabs(z) - rho) / (1-1/gamma);
	else ans = z;
		
	return(ans);
	
	
}


double penaltyfunc_C(double Lambda, double rho, double gamma){
	
	double ans=5.0;
	if(fabs(Lambda) < rho*gamma) ans = rho * (fabs(Lambda) - Lambda*Lambda/(2.0*rho*gamma));
	else ans = rho*rho*gamma/2.0;
	if(rho==0.0) ans=0.0;
	
	return(ans);
	
}



// EM
// 入力：データ(y,X)，step幅（delta_t），ステップ数（STEP），
// 出力：係数の推定値行列(betahat_matrix)，自由度ベクトル（df_zou）






void LamdPsiupdate_C2(int p, int m, int N, int ex_Npflag, int ex_dimx0, int ex_dimS, double ex_Lambda[p*m], double ex_diagPsi[p], double ex_tol1[1], double ex_tol2[1], int ex_maxcount1[1], int ex_maxcount2[1], double ex_eta[1], double ex_tolPsi[1], double ex_rho[1], double ex_gamma[1], double ex_S[ex_dimS], double ex_diagS[p], double ex_X[ex_dimx0], int ex_flagpenalty[1], int ex_fixindex[p*m], double ans[p*m+p+3]){
	
	int pm = p*m;
	int pp = p*p;
	int mm = m*m;
	int count=0;
	int t1=0;
	int t2=0;
	int blocksize=64;
	int pblocksize=p*blocksize;
	int mblocksize=m*blocksize;
	int one=1;
	int zero=0;
	int info=0;
	int i;
	int j;
	int l;
	int l0;
	int goro;
	int zerocolumncheck_flag=0;
	int zerorowcheck_flag=0;

	
	double beta_jstar_previous=0.0;
	double L0norm;
	double sum_Lsabun1=0.0;
	double sum_Lsabun2=0.0;
	double ALambdanew=0.0;
	double ALambdanewA_forPsi=0.0;
	double BLambdanew_forPsi=0.0;
	double logF=0.0;
	double penalty0=0.0;
	
	
	double Lambdanew[pm];
	double Lambdaold[pm];
	double Lambda0[pm];
	double Lambda1[pm];
	double Lambdasabun[pm];
	double Lambdasabun2[m];
	double diagPsinew[p];
	double diagPsiinv[p];
	double diagPsiold[p];
	double Psiinv_Lambda[pm];
	double Psiinv_Lambda_Minv[pm];
	double A[mm];
	double B[pm];
	double MA[mm];
	double D[m*N];
	double M[mm];
	int ipiv_M[m];
	int ipiv_A[m];
	double Minv[mm];
	double Ainv[mm];
	double Psiinvhalf_Lambda[pm];
	double diagm[mm];
	double work1[mblocksize];
	double penalized_likelihood[3];
	double eigenMvalues[m];
	
	//初期化
	
	
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
		M[i] = 0.0;
		MA[i] = 0.0;
		Minv[i] = 0.0;
		diagm[i] = 0.0;
	}
	
	for(i=0; i<m; i++){
		diagm[i*m + i] = 1.0;
		Lambdasabun2[i] = 0.0;
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
	}
	
	for( i=0; i<3; i++){
		penalized_likelihood[i] = 0.0;
	}
	sum_Lsabun1 = F77_CALL(dnrm2)(&pm, Lambdasabun, &one);

	
	//スタート！！

	t1=0;
	while(sum_Lsabun1 > ex_tol1[0] && t1 < ex_maxcount1[0]){
		//for(goro=0; goro<1000; goro++){
		t1=t1+1;
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
		
		//M <- t(Lambda) * Psi^(-1) * Lambda + diagm の計算
		F77_CALL(dcopy)(&mm, diagm, &one, M, &one);
		F77_CALL(dsyrk)(&UPPER, &TRANST, &m, &p, &alphaOne, Psiinvhalf_Lambda, &p, &betaOne, M, &m); 
		
		for( i=0; i<m; i++){	
			for( j=0; j<m; j++){	
				M[j+i*m] = M[i+j*m];
			}
		}
		
		//M^(-1)の計算	
		F77_CALL(dcopy)(&mm, M, &one, Minv, &one);
		F77_CALL(dgetrf)(&m, &m, Minv, &m, ipiv_M,  &info);
		F77_CALL(dgetri)(&m, Minv, &m, ipiv_M, work1, &mblocksize, &info);
		
		//Psiinv_Lambda_Minvの計算
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, Psiinv_Lambda, &p,Minv, &m, &betaZero, Psiinv_Lambda_Minv, &p);
		if(ex_Npflag==1){
			//Bの計算
			F77_CALL(dgemm)(&TRANST, &TRANSN, &m, &p, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_S, &p, &betaZero, B, &m);
			//Aの計算
			F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
			F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &p, &alphaOne, B, &m, Psiinv_Lambda_Minv, &p, &betaOne, A, &m);
		}else{
			//Dの計算
			F77_CALL(dgemm)(&TRANST, &TRANST, &m, &N, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_X, &N, &betaZero, D, &m);
			//Bの計算
			F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &p, &N, &alphaOne, D, &m, ex_X, &N, &betaZero, B, &m);
			//Aの計算
			F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
			F77_CALL(dgemm)(&TRANSN, &TRANST, &m, &m, &N, &alphaOne, D, &m, D, &m, &betaOne, A, &m);
		}
		
		
		
		if(ex_flagpenalty[0]==0){
			//No Penalty（普通にupdate）
			//A^(-1)の計算	
			F77_CALL(dcopy)(&mm, A, &one, Ainv, &one);
			F77_CALL(dgetrf)(&m, &m, Ainv, &m, ipiv_A,  &info);
			F77_CALL(dgetri)(&m, Ainv, &m, ipiv_A, work1, &mblocksize, &info);
			//Lambdanew <- B^TA^(-1)
			F77_CALL(dgemm)(&TRANST, &TRANSN, &p, &m, &m, &alphaOne, B, &m, Ainv, &m, &betaZero, Lambdanew, &p);
			//fixed lambdaはfixする
			for( i=0; i<pm; i++){
				if(ex_fixindex[i]==1) Lambdanew[i] = ex_Lambda[i];
			}
		}
		
		
		
		if(ex_flagpenalty[0]==1){
			
            //Coordinate descent
            F77_CALL(dcopy)(&pm, Lambdanew, &one, Lambda1, &one);
            for( j=0; j<p; j++){	
                sum_Lsabun2 = 100.0;
                t2 = 0;
                while(sum_Lsabun2 > ex_tol2[0] && t2 < ex_maxcount2[0]){
                    t2=t2+1;
					for( l=0; l<m; l++){
						ALambdanew=0.0;
						for(l0=0; l0<m; l0++){
							if(l0 != l) ALambdanew = ALambdanew + A[l*m+l0] * Lambdanew[p*l0 + j];
						}
						//fixedindex1=0のみupdate
						if(ex_fixindex[j+p*l]==0) Lambdanew[j+p*l] = S_func_MC_C(1.0 / diagPsinew[j] * ( B[j*m+l] - ALambdanew ), ex_rho[0] , ex_gamma[0] )  /  (  1.0 / diagPsinew[j] * A[l+l*m] );
					}
                    
                    
                    //距離計算
					for( l=0; l<m; l++){
                        Lambdasabun2[l] = Lambdanew[j+p*l] - Lambda1[j+p*l];
                    }
                    sum_Lsabun2 = F77_CALL(dnrm2)(&m, Lambdasabun2, &one);
                    
                    
                    //更新
					for( l=0; l<m; l++){
                        Lambda1[j+p*l] = Lambdanew[j+p*l];
                    }
                    
                    
				}
				
			}
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
			diagPsinew[i] = ex_diagS[i] - 2 * BLambdanew_forPsi + ALambdanewA_forPsi + ex_eta[0];
			if(diagPsinew[i] < ex_tolPsi[0]) diagPsinew[i] = ex_tolPsi[0];
		}
		
		
		for( i=0; i<p*m; i++){
			Lambdasabun[i] = Lambdanew[i] - Lambdaold[i];
		}
	sum_Lsabun1 = F77_CALL(dnrm2)(&pm, Lambdasabun, &one);
	}
//EM終了

//各列に0でない要素が1つしかなかったらその要素を０に修正（そっちのほうが尤度が大きくなるので）
	for( j=0; j<m; j++){
		zerocolumncheck_flag=0;
		for( i=0; i<p; i++){
			if(fabs(Lambdanew[i+j*p]) > ex_tol1[0]) zerocolumncheck_flag = zerocolumncheck_flag + 1;
		}
		if(zerocolumncheck_flag==0 || zerocolumncheck_flag==1 ){
			for( i=0; i<p; i++){
				diagPsinew[i]=diagPsinew[i]+Lambdanew[i+j*p]*Lambdanew[i+j*p];
				Lambdanew[i+j*p]=0.0;
			}
		}
	}
	
//各行が零ベクトルだったら対応する独自分散をdiagSに置換（そっちのほうが尤度が大きくなるので）
	for( i=0; i<p; i++){
		zerorowcheck_flag=0;
		for( j=0; j<m; j++){
			if(fabs(Lambdanew[i+j*p]) < ex_tol1[0]) zerorowcheck_flag = zerorowcheck_flag + 1;
		}
		if(zerorowcheck_flag==m){
			diagPsinew[i]=ex_diagS[i];
		}
	}
	
	
//log-likelihoodの計算
	
	
	//M <- t(Lambda) * Psi^(-1) * Lambda + diagm の再計算
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
	F77_CALL(dcopy)(&mm, diagm, &one, M, &one);
	F77_CALL(dsyrk)(&UPPER, &TRANST, &m, &p, &alphaOne, Psiinvhalf_Lambda, &p, &betaOne, M, &m); 
	
	for( i=0; i<m; i++){	
		for( j=0; j<m; j++){	
			M[j+i*m] = M[i+j*m];
		}
	}
	
	//M^(-1)のの再計算	
	F77_CALL(dcopy)(&mm, M, &one, Minv, &one);
	F77_CALL(dgetrf)(&m, &m, Minv, &m, ipiv_M,  &info);
	F77_CALL(dgetri)(&m, Minv, &m, ipiv_M, work1, &mblocksize, &info);
	
	//Psiinv_Lambda_Minvの再計算
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, Psiinv_Lambda, &p, Minv, &m, &betaZero, Psiinv_Lambda_Minv, &p);
	if(ex_Npflag==1){
		//Bの再計算
		F77_CALL(dgemm)(&TRANST, &TRANSN, &m, &p, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_S, &p, &betaZero, B, &m);
		//Aの再計算
		F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &p, &alphaOne, B, &m, Psiinv_Lambda_Minv, &p, &betaOne, A, &m);
	}else{
		//Dの再計算
		F77_CALL(dgemm)(&TRANST, &TRANST, &m, &N, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_X, &N, &betaZero, D, &m);
		//Bの再計算
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &p, &N, &alphaOne, D, &m, ex_X, &N, &betaZero, B, &m);
		//Aの再計算
		F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
		F77_CALL(dgemm)(&TRANSN, &TRANST, &m, &m, &N, &alphaOne, D, &m, D, &m, &betaOne, A, &m);
	}
	
	//M*Aの計算
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &m, &alphaOne, M, &m, A, &m, &betaZero, MA, &m);
	//det(M)
	F77_CALL(dsyev)(&jobzN, &UPPER, &m, M, &m, eigenMvalues, work1, &mblocksize, &info);
	logF = 0.0;
	for( i=0; i<p; i++){
		logF = logF + log(diagPsinew[i]) + ex_diagS[i] / diagPsinew[i]  ;
	}
	for( i=0; i<m; i++){
		logF = logF + log(eigenMvalues[i])  - ( MA[i+i*m] - 1.0);
	}
	
	
	//penaltyの計算
	penalty0=0.0;
	for( i=0; i<p*m; i++){
		penalty0 = penalty0 + penaltyfunc_C(Lambdanew[i],ex_rho[0],ex_gamma[0]);
	}
	
	
	//代入r	
	penalized_likelihood[0] = logF;
	penalized_likelihood[1] = penalty0;
	penalized_likelihood[2] = logF + 2.0*penalty0;
	
	
	for( i=0; i<p*m; i++){
		ans[i] = Lambdanew[i];
	}
	for( i=0; i<p; i++){
		ans[i+p*m] = diagPsinew[i];
	}
	for( i=0; i<3; i++){
		ans[i+p*(m+1)] = penalized_likelihood[i];
	}

	


}


void diagPsiupdate_C(int p, int m, int N, int ex_Npflag, int ex_dimx0, int ex_dimS, double ex_Lambda[p*m], double ex_diagPsi[p], double ex_tol1[1],  int ex_maxcount1[1], double ex_tolPsi[1], double ex_eta[1], double ex_S[ex_dimS], double ex_diagS[p], double ex_X[ex_dimx0],double ans[p],double A[m*m],double B[m*p])
//ans:Lambda,Psi,
{
	int pm = p*m;
	int pp = p*p;
	int mm = m*m;
	int count=0;
	int t1=0;
	int t2=0;
	int blocksize=64;
	int pblocksize=p*blocksize;
	int mblocksize=m*blocksize;
	int one=1;
	int zero=0;
	int info=0;
	int i=0;
	int j=0;
	int l=0;
	int l0=0;
	int goro=0;
	
	double beta_jstar_previous=0.0;
	double L0norm=0.0;
	double sum_diagPsisabun1=0.0;
	double ALambdanew=0.0;
	double ALambdanewA_forPsi=0.0;
	double BLambdanew_forPsi=0.0;
	double logF=0.0;
	double penalty0=0.0;
	
	double Lambdanew[pm];
	double diagPsinew[p];
	double diagPsiinv[p];
	double diagPsiold[p];
	double diagPsi0[p];
	double diagPsisabun[p];
	double Psiinv_Lambda[pm];
	double Psiinv_Lambda_Minv[pm];
//	double A[mm];
//	double B[pm];
	double D[N*m];
	double M[mm];
	double MA[mm];
	int ipiv_M[m];
	int ipiv_A[m];
	double Minv[mm];
	double Ainv[mm];
	double Psiinvhalf_Lambda[pm];
	double diagm[mm];
	double work1[mblocksize];
	double penalized_likelihood[3];
	double eigenMvalues[m];
	
	
	
	
	//初期化
	
	
	for( i=0; i<p*m; i++){
		Lambdanew[i] = ex_Lambda[i];
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
		M[i] = 0.0;
		Minv[i] = 0.0;
		Ainv[i] = 0.0;
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
		diagPsi0[i] = 0.0;
		diagPsisabun[i] = 10.0;
	}	
	
	
	for( i=0; i<m; i++){	
		ipiv_M[i] = 0;
		ipiv_A[i] = 0;
	}
	
	sum_diagPsisabun1 = F77_CALL(dnrm2)(&p, diagPsisabun, &one);
	
	
	
	//スタート！！
	
	t1=0;
	while(sum_diagPsisabun1 > ex_tol1[0] && t1 < ex_maxcount1[0]){
		//for(goro=0; goro<1000; goro++){
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
		
		//M <- t(Lambda) * Psi^(-1) * Lambda + diagm の計算
		F77_CALL(dcopy)(&mm, diagm, &one, M, &one);
		F77_CALL(dsyrk)(&UPPER, &TRANST, &m, &p, &alphaOne, Psiinvhalf_Lambda, &p, &betaOne, M, &m); 
		
		for( i=0; i<m; i++){	
			for( j=0; j<m; j++){	
				M[j+i*m] = M[i+j*m];
			}
		}
		
		//M^(-1)の計算	
		F77_CALL(dcopy)(&mm, M, &one, Minv, &one);
		F77_CALL(dgetrf)(&m, &m, Minv, &m, ipiv_M,  &info);
		F77_CALL(dgetri)(&m, Minv, &m, ipiv_M, work1, &mblocksize, &info);
		
		//Psiinv_Lambda_Minvの計算
		F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &m, &m, &alphaOne, Psiinv_Lambda, &p, Minv, &m, &betaZero, Psiinv_Lambda_Minv, &p);
		if(ex_Npflag==1){
			//Bの計算
			F77_CALL(dgemm)(&TRANST, &TRANSN, &m, &p, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_S, &p, &betaZero, B, &m);
			//Aの計算
			F77_CALL(dcopy)(&mm, Minv, &one, A, &one);
			F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &p, &alphaOne, B, &m, Psiinv_Lambda_Minv, &p, &betaOne, A, &m);
		}else{
			//Dの計算
			F77_CALL(dgemm)(&TRANST, &TRANST, &m, &N, &p, &alphaOne, Psiinv_Lambda_Minv, &p, ex_X, &N, &betaZero, D, &m);
			//Bの計算
			F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &p, &N, &alphaOne, D, &m, ex_X, &N, &betaZero, B, &m);
			//Aの計算
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
			diagPsinew[i] = ex_diagS[i] - 2 * BLambdanew_forPsi + ALambdanewA_forPsi + ex_eta[0];
			if(diagPsinew[i] < ex_tolPsi[0]) diagPsinew[i] = ex_tolPsi[0];
		}
		
		//距離計算
		for( i=0; i<p; i++){	
			diagPsisabun[i] = diagPsinew[i] - diagPsiold[i];
		}
		sum_diagPsisabun1 = F77_CALL(dnrm2)(&p, diagPsisabun, &one);
	}
	

	//log-likelihoodの計算
	//M*Aの計算
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &m, &m, &m, &alphaOne, M, &m, A, &m, &betaZero, MA, &m);
	//det(M)
	F77_CALL(dsyev)(&jobzN, &UPPER, &m, M, &m, eigenMvalues, work1, &mblocksize, &info);
	logF = 0.0;
	for( i=0; i<p; i++){
		logF = logF + log(diagPsinew[i]) + ex_diagS[i] / diagPsinew[i]  ;
	}
	for( i=0; i<m; i++){
		logF = logF + log(eigenMvalues[i])  - ( MA[i+i*m] - 1.0);
	}
	
	
	
	//代入r(penalty0=0.0)
	penalized_likelihood[0] = logF;
	penalized_likelihood[1] = penalty0;
	penalized_likelihood[2] = logF + 2.0*penalty0;
	
	
	
	//答えを出力
//	for( i=0; i<p*m; i++){
//		ans[i] = Lambdanew[i];
//	}
	for( i=0; i<p; i++){
		ans[i] = diagPsinew[i];
	}
//	for( i=0; i<3; i++){
//		ans[i+p*(m+1)] = penalized_likelihood[i];
//	}
	 
}




void GFI_C(int p, int m, int ex_dimS, double ex_Lambda[p*m], double ex_diagPsi[p], double ex_S[ex_dimS], double ans[1]){
	
	double Sigma[ex_dimS];
	double Sigmainv[ex_dimS];
	double SigmainvS[ex_dimS];
	int i;
	int pp=p*p;
	int one=1;
	int ipiv_Sigma[p];
	int blocksize=64;
	int pblocksize=p*blocksize;
	int info=0;
	double work1[pblocksize];
	
	double GFI0=0.0;  //tr(SigmainvS^2)
	double GFI1=0.0;  //tr(SigmainvS)

	
	for( i=0; i<ex_dimS; i++){
		Sigma[i] = 0.0;
		Sigmainv[i] = 0.0;
		SigmainvS[i] = 0.0;
	}
	
	
	for( i=0; i<p; i++){	
		ipiv_Sigma[i] = 0;
	}

	for( i=0; i<pblocksize; i++){	
		work1[i] = 0.0;
	}
	
	
	//Sigma=LambdaLambda^T+Psi
	F77_CALL(dgemm)(&TRANSN, &TRANST, &p, &p, &m, &alphaOne, ex_Lambda, &p, ex_Lambda, &p, &betaZero, Sigma, &p);
	for( i=0; i<p; i++){
		Sigma[i+i*p] = Sigma[i+i*p] + ex_diagPsi[i];
	}
	
	//Sigmainv
	F77_CALL(dcopy)(&pp, Sigma, &one, Sigmainv, &one);
	F77_CALL(dgetrf)(&p, &p, Sigmainv, &p, ipiv_Sigma,  &info);
	F77_CALL(dgetri)(&p, Sigmainv, &p, ipiv_Sigma, work1, &pblocksize, &info);
	
	//(Sigmainv%*%S)
	F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &p, &p, &alphaOne, Sigmainv, &p, ex_S, &p, &betaZero, SigmainvS, &p);
	//F77_CALL(dgemm)(&TRANSN, &TRANSN, &p, &p, &p, &alphaOne, SigmainvS, &p, SigmainvS, &p, &betaZero, SigmainvS2, &p);
	
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
}

void LamdPsiupdate_initialver2_C2(int p, int m, int N, int ex_Npflag, int ex_dimx0, int ex_dimS, int initaliter, double ex_Lambda[p*m], double ex_diagPsi[p], double ex_tol1[1], double ex_tol2[1], int ex_maxcountinitial[1], int ex_maxcount2[1], double ex_eta[1], double ex_tolPsi[1], double ex_rho[1], double ex_gamma[1], double ex_S[ex_dimS], double ex_diagS[p], double ex_X[ex_dimx0], int ex_flagpenalty[1], int ex_fixindex[p*m], double ex_Lambdainitial_rand[p*initaliter],int ex_m0[1], int ex_maxcount1[1], double ans[p*m+p+3])
{
	
	
	int pm = p*m;
	int i = 0;
	int initiali = 0;
	int one = 1;
	int three = 3;
	double penalized_likelihood_kouho = 1000000000000000.0;
	double temp = 0.0;

	
	double ans0[p*m+p+3];
	double Lambda[pm];
	double Lambda0[pm];
	double diagPsi[p];
	double diagPsi0[p];
	double penalized_likelihood[3];
	double penalized_likelihood0[3];
	
	int Lambda_index_allfree[pm];
	int flagpenalty2[1];
	
	
	
	
	for( i=0; i<p*m; i++){
		Lambda[i] = ex_Lambda[i];
		Lambda0[i] = ex_Lambda[i];
		Lambda_index_allfree[i] = 0;
	}
	for( i=0; i<p; i++){
		diagPsi[i] = ex_diagPsi[i];
		diagPsi0[i] = ex_diagPsi[i];
	}
	for( i=0; i<3; i++){
		penalized_likelihood[i] = 0.0;
		penalized_likelihood0[i] = 0.0;
	}
	flagpenalty2[0]=1;

	
	
	
	 for( initiali=0; initiali<initaliter; initiali++){
		
		for(i=0; i<p;i++){
			Lambda0[i+p*ex_m0[0]] = ex_Lambdainitial_rand[i+initiali*p];
		}
		LamdPsiupdate_C2(p,m,N, ex_Npflag, ex_dimx0, ex_dimS, Lambda0, diagPsi, ex_tol1, ex_tol2, ex_maxcountinitial, ex_maxcount2, ex_eta, ex_tolPsi, ex_rho, ex_gamma, ex_S, ex_diagS, ex_X, ex_flagpenalty, ex_fixindex, ans0);
		for(i=0; i<p*m;i++){
			temp = ans0[i];
			Lambda0[i]=temp;
		}
		for(i=0; i<p;i++){
			temp = ans0[i+p*m];
			diagPsi0[i] =temp;
		}
		LamdPsiupdate_C2(p,m,N,ex_Npflag, ex_dimx0, ex_dimS, Lambda0, diagPsi0, ex_tol1, ex_tol2, ex_maxcountinitial, ex_maxcount2, ex_eta, ex_tolPsi, ex_rho, ex_gamma, ex_S, ex_diagS, ex_X, flagpenalty2, Lambda_index_allfree, ans0);
		for(i=0; i<p*m;i++){
			temp = ans0[i];
			Lambda0[i] =temp;
		}
		for(i=0; i<p;i++){
			temp = ans0[i+p*m];
			diagPsi0[i] =temp;
		}
		for(i=0; i<3;i++){
			temp = ans0[i+p*m+p];
			penalized_likelihood0[i]=temp;
		}
		temp=ans0[2+p*m+p];
		if(penalized_likelihood_kouho > temp){
			penalized_likelihood_kouho = temp;
			F77_CALL(dcopy)(&pm, Lambda0, &one, Lambda, &one);
			F77_CALL(dcopy)(&p, diagPsi0, &one, diagPsi, &one);
			F77_CALL(dcopy)(&three, penalized_likelihood0, &one, penalized_likelihood, &one);
		}
	}
	
	
	
	LamdPsiupdate_C2(p,m,N,ex_Npflag, ex_dimx0, ex_dimS,Lambda, diagPsi, ex_tol1, ex_tol2, ex_maxcount1, ex_maxcount2, ex_eta, ex_tolPsi, ex_rho, ex_gamma, ex_S, ex_diagS, ex_X, flagpenalty2, Lambda_index_allfree,ans);
	
}


SEXP LamdPsiupdate_initialver1_R(SEXP ex_Lambda, SEXP ex_diagPsi, SEXP ex_N, SEXP ex_tol1, SEXP ex_tol2, SEXP ex_maxcount1, SEXP ex_eta, SEXP ex_tolPsi, SEXP ex_rho, SEXP ex_gamma, SEXP ex_S, SEXP ex_diagS, SEXP ex_X, SEXP flagpenalty, SEXP ex_fixindex, SEXP ex_Lambdainitial_rand, SEXP ex_diagPsiinitial_rand,  SEXP ex_Lambdainitialij_for_rhoinitial, SEXP ex_Npflag, SEXP ex_dimx0, SEXP ex_dimS)
{
	
	SEXP ans0;
	SEXP ans;
	SEXP Lambda;
	SEXP diagPsi;
	SEXP Lambda0;
	SEXP diagPsi0;
	SEXP Lambda_ans;
	SEXP Lambda_initial;
	SEXP diagPsi_ans;
	SEXP diagPsi_initial;
	SEXP penalized_likelihood;
	SEXP penalized_likelihood0;
	SEXP A;
	SEXP B;
	SEXP rho;

	int p = INTEGER(GET_DIM(ex_Lambda))[0];
	int m = INTEGER(GET_DIM(ex_Lambda))[1];
	int N = INTEGER(ex_N)[0];
//	int N = INTEGER(GET_DIM(ex_X))[0];
	int initaliter = INTEGER(GET_DIM(ex_Lambdainitial_rand))[1];
	int i = 0;
	int initiali = 0;
	int pm=p*m;
	int three=3;
	int one=1;
	int l,j,l0 = 0;

	
	int Lambdaidamax=1;;
	int Lambdaidamax2=1;
	double penalized_likelihood_kouho = 1000000000000000.0;
	double ALambda = 0.0;
	double initial_check_kouho = 0.0;
	double initial_check = 0.0;
	double ALambdanew = 0.0;
	double x = 0.0;
	double temp = 0.0;
	
	PROTECT(Lambda = allocMatrix(REALSXP,p,m));
	PROTECT(diagPsi = allocVector(REALSXP,p));
	PROTECT(Lambda0 = allocMatrix(REALSXP,p,m));
	PROTECT(diagPsi0 = allocVector(REALSXP,p));
	PROTECT(Lambda_ans = allocMatrix(REALSXP,p,m));
	PROTECT(Lambda_initial = allocMatrix(REALSXP,p,m));
	PROTECT(diagPsi_ans = allocVector(REALSXP,p));
	PROTECT(diagPsi_initial= allocVector(REALSXP,p));
	PROTECT(penalized_likelihood = allocVector(REALSXP, 3));
	PROTECT(penalized_likelihood0 = allocVector(REALSXP, 3));
	PROTECT(A = allocMatrix(REALSXP,m,m));
	PROTECT(B = allocMatrix(REALSXP, m,p));
	PROTECT(rho = allocVector(REALSXP, 1));
	PROTECT(ans0 = allocVector(REALSXP, p*m+p+3));
	PROTECT(ans = allocVector(VECSXP, 5));
	
	for( i=0; i<m*m; i++){
		REAL(A)[i] = 0.0;
	}
	for( i=0; i<p*m; i++){
		REAL(B)[i] = 0.0;
	}

	
	for( i=0; i<p*m; i++){
		REAL(Lambda)[i] = REAL(ex_Lambda)[i];
		REAL(Lambda0)[i] = REAL(ex_Lambda)[i];
		REAL(Lambda_ans)[i] = REAL(ex_Lambda)[i];
		REAL(Lambda_initial)[i] = REAL(ex_Lambda)[i];
	}
	for( i=0; i<p; i++){
		REAL(diagPsi)[i] = REAL(ex_diagPsi)[i];
		REAL(diagPsi0)[i] = REAL(ex_diagPsi)[i];
		REAL(diagPsi_ans)[i] = REAL(ex_diagPsi)[i];
	}
	for( i=0; i<3; i++){
		REAL(penalized_likelihood)[i] = 0.0;
		REAL(penalized_likelihood0)[i] = 0.0;
	}
	
	REAL(rho)[0]=0.0;
	
	//SET_VECTOR_ELT(ans0, 0,Lambda0);
	//SET_VECTOR_ELT(ans0, 1,diagPsi0);	
	//SET_VECTOR_ELT(ans0, 2,penalized_likelihood0);
	
	
	for( initiali=0; initiali<initaliter; initiali++){
		for(i=0; i<p;i++){
			REAL(Lambda0)[i] = REAL(ex_Lambdainitial_rand)[i+initiali*p];
			REAL(diagPsi0)[i] = REAL(ex_diagPsiinitial_rand)[i+initiali*p];
		}
		LamdPsiupdate_C2(p,m,N, INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], REAL(Lambda0), REAL(diagPsi0), REAL(ex_tol1), REAL(ex_tol2), INTEGER(ex_maxcount1), INTEGER(ex_maxcount1), REAL(ex_eta), REAL(ex_tolPsi), REAL(ex_rho), REAL(ex_gamma), REAL(ex_S), REAL(ex_diagS), REAL(ex_X), INTEGER(flagpenalty), INTEGER(ex_fixindex),REAL(ans0));
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
			F77_CALL(dcopy)(&pm, REAL(Lambda0), &one, REAL(Lambda_initial), &one);
			F77_CALL(dcopy)(&p, REAL(diagPsi0), &one, REAL(diagPsi), &one);
			F77_CALL(dcopy)(&p, REAL(diagPsi0), &one, REAL(diagPsi_initial), &one);
			F77_CALL(dcopy)(&three, REAL(penalized_likelihood0), &one, REAL(penalized_likelihood), &one);
		}
	}
	
	//choose two elements of maximum absolute value of Lambda
	F77_CALL(dcopy)(&pm, REAL(Lambda), &one, REAL(Lambda0), &one);
	F77_CALL(dcopy)(&p, REAL(diagPsi), &one, REAL(diagPsi0), &one);
	F77_CALL(dcopy)(&three, REAL(penalized_likelihood), &one, REAL(penalized_likelihood0), &one);
	
	Lambdaidamax = F77_CALL(idamax)(&pm,REAL(Lambda0),&one);
	Lambdaidamax = Lambdaidamax-1;
	/*
	REAL(Lambda0)[Lambdaidamax]=0.0;
	Lambdaidamax2 = F77_CALL(idamax)(&pm,REAL(Lambda0),&one);
	Lambdaidamax2 = Lambdaidamax2-1;
	REAL(Lambda0)[Lambdaidamax]=REAL(Lambda)[Lambdaidamax];
	*/
	
	for(i=0; i<pm;i++){
		if(i!=Lambdaidamax) REAL(Lambda0)[i]=0.0;
//		if(i!=Lambdaidamax && i!=Lambdaidamax2) REAL(Lambda0)[i]=0.0;
	}
	F77_CALL(dcopy)(&pm, REAL(Lambda0), &one, REAL(Lambda), &one);

	//update diagPsi while fixing Lambda
	for( initiali=0; initiali<initaliter; initiali++){

		for( j=0; j<p*m; j++){	
			REAL(Lambda)[j]=REAL(Lambda0)[j] * REAL(ex_Lambdainitialij_for_rhoinitial)[initiali];
		}
		
		//update psi(diagPsi0:old, diagPsi:new)
		 diagPsiupdate_C(p, m, N, INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0],  REAL(Lambda), REAL(diagPsi0), REAL(ex_tol1), INTEGER(ex_maxcount1), REAL(ex_tolPsi), REAL(ex_eta), REAL(ex_S), REAL(ex_diagS), REAL(ex_X), REAL(diagPsi), REAL(A), REAL(B));

		initial_check_kouho = 0.0;
		for( j=0; j<p; j++){	
			for( l=0; l<1; l++){
				ALambda=0.0;
				for(l0=0; l0<m; l0++){
					if(l0 != l) ALambda = ALambda + REAL(A)[l*m+l0] * REAL(Lambda)[p*l0 + j];
				}
				if(Lambdaidamax!=j){
//				if(Lambdaidamax!=j && Lambdaidamax2!=j){
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
	SET_VECTOR_ELT(ans, 1,Lambda_initial);
	SET_VECTOR_ELT(ans, 2,diagPsi_ans);	
	SET_VECTOR_ELT(ans, 3,diagPsi_initial);	
	SET_VECTOR_ELT(ans, 4,rho);
	
	
	UNPROTECT(15);
	
	return(ans);
}




SEXP fancmain(SEXP ex_Lambda_ini, SEXP ex_diagPsi_ini, SEXP ex_N, SEXP ex_tol1, SEXP ex_tol2, SEXP ex_maxcount1, SEXP ex_maxcount2, SEXP ex_maxcountinitial, SEXP ex_eta, SEXP ex_tolPsi, SEXP ex_rhomatrix, SEXP ex_gamma_vec, SEXP ex_X, SEXP ex_diagS, SEXP ex_S, SEXP ex_Lambdainitial_rand , SEXP ex_Npflag, SEXP ex_dimx0, SEXP ex_dimS){
	SEXP Lambda_old;
	SEXP Lambda_new;
	SEXP Lambda_temp;
	SEXP diagPsi_old;
	SEXP diagPsi_new;
	SEXP diagPsi_temp;
	SEXP Lambda_index_allfree;
	SEXP Lambda_index00;

	SEXP parameterupdate;
	SEXP parameterupdate_new;
	SEXP penalized_likelihood;
	SEXP penalized_likelihood_new;
	SEXP penalized_likelihood_temp;
	
	SEXP flagpenalty;

	SEXP rho;
	SEXP gamma;
	SEXP m0SEXP;

	SEXP sumabsm0new_flagvec;

	SEXP ans_Lambda;
	SEXP ans_diagPsi;
	SEXP ans_AIC;
	SEXP ans_BIC;
	SEXP ans_CAIC;
	SEXP ans_DF;
	SEXP ans_GFI;
	SEXP ans_logF2c;
	SEXP ans;
	
	//int N = INTEGER(GET_DIM(ex_X))[0];
	int N = INTEGER(ex_N)[0];
	int p = INTEGER(GET_DIM(ex_Lambda_ini))[0];
	int m = INTEGER(GET_DIM(ex_Lambda_ini))[1];
	int i,j,l,i0,j0;
	int pm = p*m;
	int stepgamma,steprho;
	int STEPrho =  INTEGER(GET_DIM(ex_rhomatrix))[0];
	int STEPgamma =  INTEGER(GET_DIM(ex_rhomatrix))[1];
	int initaliter = INTEGER(GET_DIM(ex_Lambdainitial_rand))[1];
	
	int one=1;
	int m0=1;
	int m0new=1;
	

	int flag_m=1;
	int flag_m_FIRST=1;
	int flag_NewOld=1; //1:new

	double sumabsLambda=0.0;
	double sumabsm0=0.0;
	double logF_Lambda0=0.0;
	double temp=0.0;
	double tempnew=0.0;
	double DF=0.0;
	double GFI[1];
	


	
	PROTECT(Lambda_old = allocMatrix(REALSXP, p,m));
	PROTECT(Lambda_new = allocMatrix(REALSXP, p,m));
	PROTECT(Lambda_temp = allocMatrix(REALSXP, p,m));
	PROTECT(diagPsi_old = allocVector(REALSXP, p));
	PROTECT(diagPsi_new = allocVector(REALSXP, p));
	PROTECT(diagPsi_temp = allocVector(REALSXP, p));
	PROTECT(Lambda_index_allfree = allocMatrix(INTSXP, p,m));
	PROTECT(Lambda_index00 = allocMatrix(INTSXP, p,m));
	PROTECT(penalized_likelihood = allocVector(REALSXP, 3));
	PROTECT(penalized_likelihood_new = allocVector(REALSXP, 3));
	PROTECT(penalized_likelihood_temp = allocVector(REALSXP, 3));
	PROTECT(parameterupdate = allocVector(REALSXP, p*m+p+3));
	PROTECT(parameterupdate_new = allocVector(REALSXP, p*m+p+3));
	
	PROTECT(flagpenalty = allocVector(INTSXP, 1));
	PROTECT(rho = allocVector(REALSXP, 1));
	PROTECT(gamma = allocVector(REALSXP, 1));
	PROTECT(m0SEXP = allocVector(INTSXP, 1));
	PROTECT(sumabsm0new_flagvec = allocVector(INTSXP, m));

	
	PROTECT(ans_AIC = allocMatrix(REALSXP, STEPrho,STEPgamma));
	PROTECT(ans_BIC = allocMatrix(REALSXP, STEPrho,STEPgamma));
	PROTECT(ans_CAIC = allocMatrix(REALSXP, STEPrho,STEPgamma));
	PROTECT(ans_DF = allocMatrix(REALSXP, STEPrho,STEPgamma));
	PROTECT(ans_GFI = allocMatrix(REALSXP, STEPrho,STEPgamma));
	PROTECT(ans_Lambda = allocMatrix(REALSXP, p,m*STEPrho*STEPgamma));
	PROTECT(ans_diagPsi = allocMatrix(REALSXP, p,STEPrho*STEPgamma));
	PROTECT(ans_logF2c = allocMatrix(REALSXP, 3,STEPrho*STEPgamma));
	PROTECT(ans = allocVector(VECSXP, 8));

	
//初期化
	for(i=0;i<p*m;i++){
		REAL(Lambda_old)[i]=REAL(ex_Lambda_ini)[i];
		REAL(Lambda_new)[i]=0.0;
		REAL(Lambda_temp)[i]=0.0;
		INTEGER(Lambda_index_allfree)[i]=0;
		INTEGER(Lambda_index00)[i]=0;
	}
	for(i=0;i<p;i++){
		REAL(diagPsi_old)[i]=REAL(ex_diagPsi_ini)[i];
		REAL(diagPsi_new)[i]=0.0;
		REAL(diagPsi_temp)[i]=0.0;
	}
	for(i=0;i<STEPrho*STEPgamma;i++){
		REAL(ans_AIC)[i]=0.0;
		REAL(ans_BIC)[i]=0.0;
		REAL(ans_CAIC)[i]=0.0;
		REAL(ans_DF)[i]=0.0;
		REAL(ans_GFI)[i]=0.0;
	}
	for(i=0;i<m;i++){
		INTEGER(sumabsm0new_flagvec)[i]=0;
	}
	for(i=0;i<p*m*STEPrho*STEPgamma;i++){
		REAL(ans_Lambda)[i]=0.0;
	}
	for(i=0;i<p*STEPrho*STEPgamma;i++){
		REAL(ans_diagPsi)[i]=0.0;
	}
	for(i=0;i<3*STEPrho*STEPgamma;i++){
		REAL(ans_logF2c)[i]=0.0;
	}
	for(i=0;i<3;i++){
		REAL(penalized_likelihood)[i]=0.0;
		REAL(penalized_likelihood_new)[i]=0.0;
		REAL(penalized_likelihood_temp)[i]=0.0;
	}
	for(i=0;i< p*m+p+3 ;i++){
		REAL(parameterupdate)[i]=0.0;
		REAL(parameterupdate_new)[i]=0.0;
	}	
	
	REAL(gamma)[0]=0.0;
	REAL(rho)[0]=0.0;
	INTEGER(flagpenalty)[0] = 1;
	INTEGER(m0SEXP)[0] = m0;
	
	//start!!
	for(stepgamma=0; stepgamma<STEPgamma; stepgamma++){
		REAL(gamma)[0] = REAL(ex_gamma_vec)[stepgamma];
		for(steprho=0; steprho<STEPrho; steprho++){
			REAL(rho)[0] = REAL(ex_rhomatrix)[steprho + stepgamma*STEPrho];
				
				//sparsenetの更新に従ってみる
			if(stepgamma==0  && steprho>0){
				for(i=0;i<p*m;i++){
					REAL(Lambda_old)[i] = REAL(ans_Lambda)[i+(steprho-1)*p*m];
				}
				for(i=0;i<p;i++){
					REAL(diagPsi_old)[i] = REAL(ans_diagPsi)[i+(steprho-1)*p];
				}
			}
				
			if(stepgamma >= 1){
				for(i=0;i<p*m;i++){
					REAL(Lambda_old)[i] = REAL(ans_Lambda)[i+(stepgamma-1)*STEPrho*p*m + steprho*p*m];
				}
				for(i=0;i<p;i++){
					REAL(diagPsi_old)[i] = REAL(ans_diagPsi)[i+(stepgamma-1)*STEPrho*p + steprho*p];
				}
			}

			//絶対値が小さかったらinitialに変更
			sumabsLambda = F77_CALL(dnrm2)(&pm, REAL(Lambda_old), &one);
			if(sumabsLambda < REAL(ex_tol1)[0]){
				F77_CALL(dcopy)(&pm, REAL(ex_Lambda_ini), &one, REAL(Lambda_old), &one);
				F77_CALL(dcopy)(&p, REAL(ex_diagPsi_ini), &one, REAL(diagPsi_old), &one);
			}

			

			//update parameter
			INTEGER(flagpenalty)[0] = 1;
			LamdPsiupdate_C2(p,m,N, INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], REAL(Lambda_old), REAL(diagPsi_old),REAL(ex_tol1),REAL(ex_tol2),INTEGER(ex_maxcount1),INTEGER(ex_maxcount2),REAL(ex_eta), REAL(ex_tolPsi), REAL(rho),REAL(gamma), REAL(ex_S), REAL(ex_diagS), REAL(ex_X), INTEGER(flagpenalty), INTEGER(Lambda_index_allfree),REAL(parameterupdate));
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
			
			if(m0==0) flag_m=1;						
			if(m0==m) flag_m=1;	
			
			if(flag_m==1) flag_NewOld = 0;

			//check if the number of factors should be updated	
			while(flag_m == 0){
				flag_m=1;//必ず消すこと！
				for(i=0;i<p*m;i++){
					INTEGER(Lambda_index00)[i]=INTEGER(Lambda_index_allfree)[i];
				}	
				for(i=0;i<p*m0;i++){
					INTEGER(Lambda_index00)[i] = 1;
				}
					
				INTEGER(flagpenalty)[0] = 0;
				INTEGER(m0SEXP)[0] = m0;
				
				
				for(i=0;i<p*m;i++){
					temp= REAL(parameterupdate)[i];
					REAL(Lambda_old)[i] =temp;
				}
				for(i=0;i<p;i++){
					temp=REAL(parameterupdate)[i+p*m];
					REAL(diagPsi_old)[i] =temp;
				}
				LamdPsiupdate_initialver2_C2(p,m,N,INTEGER(ex_Npflag)[0], INTEGER(ex_dimx0)[0], INTEGER(ex_dimS)[0], initaliter,REAL(Lambda_old), REAL(diagPsi_old), REAL(ex_tol1), REAL(ex_tol2), INTEGER(ex_maxcountinitial), INTEGER(ex_maxcount2) ,REAL(ex_eta), REAL(ex_tolPsi), REAL(rho), REAL(gamma), REAL(ex_S), REAL(ex_diagS), REAL(ex_X), INTEGER(flagpenalty), INTEGER(Lambda_index00), REAL(ex_Lambdainitial_rand), INTEGER(m0SEXP), INTEGER(ex_maxcount1),REAL(parameterupdate_new));
				//compute m0new
				m0new=0;
				for(j=0;j<m;j++){
					sumabsm0=0.0;
					for(i=0;i<p;i++){
						sumabsm0 = sumabsm0 + fabs(REAL(parameterupdate_new)[i+j*p]);
					}
					INTEGER(sumabsm0new_flagvec)[j] = 0;
					if(sumabsm0 > REAL(ex_tol1)[0]){
						m0new=m0new+1;
						INTEGER(sumabsm0new_flagvec)[j] = 1;
					}
				}
				//check updating if m0 <- m0new
				tempnew=REAL(parameterupdate_new)[p*m+p+2];
				temp=REAL(parameterupdate)[p*m+p+2];
				if(m0new > m0){
					if(tempnew >temp){
						flag_NewOld = 0;
						flag_m = 1;
					}
				}
				if(m0new == 0){
					flag_NewOld = 0;
					flag_m = 1;
				}
					
				if(m0new == m || m0new <= m0){
					flag_m = 1;
					if(tempnew >temp){
						flag_NewOld = 0;
					}
				}
					
				if(flag_m == 0){
					
					for(j=0;j<m;j++){
						for(i=0;i<p;i++){
							temp = REAL(parameterupdate_new)[i+j*p];
							REAL(Lambda_old)[i+j*p]=temp;
						}
					}
					for(i=0;i<p;i++){
						 temp= REAL(parameterupdate_new)[i+p*m];
						REAL(diagPsi_old)[i]=temp;
					}
					
					//update m0
					m0=m0new;
				}
			
				
			}
			//endwhile
			//結果を代入
			for(i=0;i<p*m;i++){
				REAL(Lambda_new)[i]=0.0;
			}			
			for(i=0;i<p;i++){
				REAL(diagPsi_new)[i]=0.0;
			}			
			
			if(flag_NewOld==1){
				l=0;
				for(j=0;j<m;j++){
					if(INTEGER(sumabsm0new_flagvec)[j] == 1){
						for(i=0;i<p;i++){
							i0 = i+l*p+steprho*p*m+STEPrho*stepgamma*p*m;
							j0 = i+j*p;
							temp =  REAL(parameterupdate_new)[j0];
							REAL(ans_Lambda)[i0] = temp;
							REAL(Lambda_new)[j0]=temp;
						}
						l=l+1;
					}
				}
				for(i=0;i<p;i++){
					temp = REAL(parameterupdate_new)[i+p*m];
					i0=i+steprho*p+STEPrho*stepgamma*p;
					REAL(ans_diagPsi)[i0]=temp;
					REAL(diagPsi_new)[i]=temp;
				}
				REAL(ans_logF2c)[0+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate_new)[p*m+p+0];
				REAL(ans_logF2c)[1+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate_new)[p*m+p+1];
				REAL(ans_logF2c)[2+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate_new)[p*m+p+2];
				
				if(stepgamma==0){
					DF=0.0;
					for(i=0;i<pm;i++){
						if(fabs(REAL(parameterupdate_new)[i]) > REAL(ex_tol1)[0]) DF=DF+1.0;
					}
					DF=DF+p;
				}else{
					DF=REAL(ans_DF)[steprho + STEPrho*(stepgamma-1)];
				}
				
				REAL(ans_DF)[steprho + STEPrho*stepgamma]= DF;
				REAL(ans_AIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate_new)[p*m+p+0] + 2*DF;
				REAL(ans_BIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate_new)[p*m+p+0] + log(N)*DF;
				REAL(ans_CAIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate_new)[p*m+p+0] + (log(N)+1)*DF;
				
				GFI[0]=0.0;
				if(INTEGER(ex_Npflag)[0]==1) GFI_C(p, m, INTEGER(ex_dimS)[0], REAL(Lambda_new), REAL(diagPsi_new), REAL(ex_S), GFI);
				REAL(ans_GFI)[steprho + STEPrho*stepgamma]= GFI[0];		
			}
			
			
			if(flag_NewOld==0){
				for(j=0;j<m;j++){
					for(i=0;i<p;i++){
						i0 = i+j*p+steprho*p*m+STEPrho*stepgamma*p*m;
						j0 = i+j*p;
						temp =  REAL(parameterupdate)[j0];
						REAL(ans_Lambda)[i0] = temp;
						REAL(Lambda_new)[j0]=temp;
					}
				}
				
				for(i=0;i<p;i++){
					temp = REAL(parameterupdate)[i+p*m];
					i0=i+steprho*p+STEPrho*stepgamma*p;
					REAL(ans_diagPsi)[i0]=temp;
					REAL(diagPsi_new)[i]=temp;
				}
				REAL(ans_logF2c)[0+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate)[p*m+p+0];
				REAL(ans_logF2c)[1+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate)[p*m+p+1];
				REAL(ans_logF2c)[2+steprho*3+STEPrho*stepgamma*3]=REAL(parameterupdate)[p*m+p+2];
				
				if(stepgamma==0){
					DF=0.0;
					for(i=0;i<pm;i++){
						if(fabs(REAL(parameterupdate)[i]) > REAL(ex_tol1)[0]) DF=DF+1.0;
					}
					DF=DF+p;
				}else{
					DF=REAL(ans_DF)[steprho + STEPrho*(stepgamma-1)];
				}
				
				REAL(ans_DF)[steprho + STEPrho*stepgamma]= DF;
				REAL(ans_AIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate)[p*m+p+0] + 2*DF;
				REAL(ans_BIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate)[p*m+p+0] + log(N)*DF;
				REAL(ans_CAIC)[steprho + STEPrho*stepgamma]= N * REAL(parameterupdate)[p*m+p+0] + (log(N)+1)*DF;
				

				GFI[0]=0.0;
				if(INTEGER(ex_Npflag)[0]==1) GFI_C(p, m, INTEGER(ex_dimS)[0], REAL(Lambda_new), REAL(diagPsi_new), REAL(ex_S), GFI);
				REAL(ans_GFI)[steprho + STEPrho*stepgamma]= GFI[0];		
			}
		}
	}
	
	SET_VECTOR_ELT(ans, 0, ans_Lambda);
	SET_VECTOR_ELT(ans, 1, ans_diagPsi);
	SET_VECTOR_ELT(ans, 2, ans_logF2c);
	SET_VECTOR_ELT(ans, 3, ans_AIC);
	SET_VECTOR_ELT(ans, 4, ans_BIC);
	SET_VECTOR_ELT(ans, 5, ans_CAIC);
	SET_VECTOR_ELT(ans, 6, ans_DF);
	SET_VECTOR_ELT(ans, 7, ans_GFI);
	
	
	UNPROTECT(27);

	return(ans);
}



