#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>


void sphgrad(double *grad, double *X, double *r, double *D, double *ldm, double *udm, double *ldiss, double *udiss, int *N, int *P)
	{
		int i, j, n, p, s;
		double tmpdij;
		n = *N;
		p = *P;
		//
		for(i = 0; i <n; i++){// start i
			grad[n*p + i] = 0;
			for(j = 1; j < n; j++){//start for j
				if(i != j){
					tmpdij = D[i+j*n];
					if(tmpdij==0){
						tmpdij = 0.001;
					}
					grad[n*p + i] += -4*( udiss[i+j*n] -udm[i+j*n])*(r[i]);
					if(ldm[i +j*n]>0 ){
						grad[n*p + i] += +4*( ldiss[i+j*n] -ldm[i+j*n])*(r[i]);
						for(s=0; s<p;s++){//start for s
							grad[i+n*s] += -2*( ldiss[i+j*n] - ldm[i+j*n])*(X[i+n*s]-X[j+n*s])/tmpdij;
						}//end for s
					}
					for(s=0; s<p;s++){//start for s
						grad[i+n*s] += -2*( udiss[i+j*n] -udm[i+j*n])*(X[i+n*s]-X[j+n*s])/tmpdij;
					}//end for s
				}//end for j
			}
		}//end for i
	}
	
	
void mmSph(double *X, double *r, double *str, double *ldist, double *udist, double *ldiss, double *udiss, double *eps, int *relax, int *dilation, int *N, int *P, int *maxit, int *report)
	{//start void	Rprintf("initial stress= %f \n", str[0]);
		/*start: def var*/
		int i, j, iter, k, n, s, p, rp, MAXIT, info, lwork;
		double tmp=0,tmpri=0, tmprj=0, tmpcd=0, tmpstr =0, EPS=0, eta=0, dil =0;
		double *CD, *Y, *q, *alp1,*alp2,*bet1,*bet2, *gam1, *gam2, *gam3, *gam4, *del1, *del2, *del3 , *A1, *A2, *B, *b;
		//def var for SVD
		double *singular_value, *work, *U_svd, *Vt_svd;
		char jobu = 'A', jobvt = 'A';
		//def var for pinv
		double *bw1, *bw2;
		
		n = *N;
		p = *P;
		MAXIT = *maxit;
		lwork = 5*n;
		rp = *report;
		EPS = *eps;
		CD=(double *)calloc(n*n,sizeof(double));
		//Old variables
		Y =(double *)calloc(n*p,sizeof(double));
		q = (double *)calloc(n,sizeof(double));
		//coefficients for MM
		alp1 =(double *)calloc(n*n,sizeof(double));
		alp2 =(double *)calloc(n*n,sizeof(double));
		bet1 =(double *)calloc(n*n,sizeof(double));
		bet2 =(double *)calloc(n*n,sizeof(double));
		gam1 =(double *)calloc(n*n,sizeof(double));
		gam2 =(double *)calloc(n*n,sizeof(double));
		gam3 =(double *)calloc(n*n,sizeof(double));
		gam4 =(double *)calloc(n*n,sizeof(double));
		del1 =(double *)calloc(n*n,sizeof(double));
		del2 =(double *)calloc(n*n,sizeof(double));
		del3 =(double *)calloc(n*n,sizeof(double));
		A1 =(double *)calloc(n*n,sizeof(double));
		A2 =(double *)calloc(n,sizeof(double));
		B =(double *)calloc(n*n,sizeof(double));
		b =(double *)calloc(n,sizeof(double));
		//var for svd
		singular_value = (double *)calloc(n,sizeof(double));
		U_svd = (double *)calloc(n*n,sizeof(double));
		Vt_svd = (double *)calloc(n*n,sizeof(double));
		work = (double *)calloc(lwork,sizeof(double));
		//var for pinv
		bw1=(double *)calloc(n,sizeof(double));
		bw2=(double *)calloc(n,sizeof(double));
		/*end: def var*/
		
		/*start: compute centre distances*/
		for(i=0; i<(n-1); i++){
			for(j =i+1; j <n; j++){
				eta += pow(ldiss[i+j*n],2)+ pow(udiss[i+j*n],2);
				tmpcd =0;
				for(s=0; s<p; s++){
					tmpcd += pow(X[i+s*n]-X[j+s*n],2);
				}
				CD[i+j*n] = sqrt(tmpcd);
				CD[j+i*n] = CD[i+j*n];
			}
		}
		/*end: compute centre dist.*/
		/*start: compute interval dist. and stress*/
		tmpstr = 0;
		for(i = 0; i <(n-1); i++){
			for(j = i+1; j < n; j++){
				ldist[i + j*n] = CD[i + j*n] - r[i] - r[j];
				if(ldist[i + j*n]<0){
					ldist[i + j*n] = 0;
				}
				tmpstr += pow(ldiss[i+j*n]-ldist[i+j*n],2);
				ldist[j+i*n] = ldist[i + j*n];
				udist[i+j*n] = CD[i+j*n] + r[i] + r[j];
				tmpstr += pow(udiss[i+j*n]-udist[i+j*n],2);
			}
		}
		/*end: compute interval dist. and stress*/
		str[0] = tmpstr; 
		Rprintf("initial value %f \n", str[0]);
		/*Main loop===============================*/
		for(iter=1;iter<MAXIT+1;iter++){/*start for r*/
		/*r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r*/
			for(i = 0; i <(n-1); i++){//i
				for(j = i+1; j < n; j++){//j
					tmpri = r[i];
					tmprj = r[j];
					tmpcd = CD[i+j*n];
					if(tmpri==0){
						tmpri = 0.000001;
					}
					if(tmprj==0){
						tmprj = 0.000001;
					}
					if(tmpcd==0){
						tmpcd = 0.000001;
					}
					alp1[i+j*n] = 1+((r[i]+r[j])/tmpcd);
					alp1[j+i*n] = alp1[i+j*n];
					alp2[i+j*n] = (r[i]+r[j]+CD[i+j*n])/tmpri;
					alp2[j+i*n] = (r[i]+r[j]+CD[i+j*n])/tmprj;
					gam1[i+j*n] = 2;
					gam1[j+i*n] =2;
					gam2[i+j*n] = 2*(1+(r[j]/tmpri));
					gam2[j+i*n] = 2*(1+(r[i]/tmprj));
					bet1[i+j*n] = 0;
					bet2[i+j*n] = udiss[i+j*n];
					bet2[j+i*n] = udiss[i+j*n];
					gam3[i+j*n] = 0;
					del1[i+j*n] =0;
					if(CD[i+j*n]>0){
						bet1[i+j*n] = udiss[i+j*n]/CD[i+j*n];
						if(CD[i+j*n]>=r[i]+r[j]){
							gam3[i+j*n] = (r[i] + r[j] +CD[i+j*n])/CD[i+j*n];
							del1[i+j*n] = ldiss[i+j*n]/CD[i+j*n];
						}else{
							gam3[i+j*n] = 2;
						}
					}
					bet1[j+i*n] = bet1[i+j*n];
					gam3[j+i*n] = gam3[i+j*n];
					del1[j+i*n] = del1[i+j*n];
					gam4[i+j*n] =0;
					del2[i+j*n]=0;
					del2[j+i*n]=0;
					del3[j+i*n]=0;
					del3[j+i*n] = 0;
					if(CD[i+j*n]>=r[i]+r[j]){
						gam4[i+j*n] = r[i] + r[j] +CD[i+j*n];
						del2[i+j*n] = ldiss[i+j*n]/tmpri;
						del2[j+i*n] = ldiss[i+j*n]/tmprj;
						del3[i+j*n] = ldiss[i+j*n]*(r[i] + r[j]);
						del3[j+i*n] = del3[i+j*n];
					}else{
						gam4[i+j*n] = 2*(r[i]+r[j]);
					}
					gam4[j+i*n] = gam4[i+j*n];
				}//i
			}//j
			//Comp. A1, A2, B, b
			for(i = 0; i <n; i++){
				b[i] =0;
				A1[i+i*n] = 0;
				A2[i] = 0;
				B[i+i*n] =0;
				alp1[i+i*n] = 0;
				alp2[i+i*n]=0;
				bet1[i+i*n] = 0;
				bet2[i+i*n] = 0;
				gam1[i+i*n] =0;
				gam2[i+i*n] =0;
				gam3[i+i*n] =0;
				gam4[i+i*n] = 0;
				del1[i+i*n] =0;
				del2[i+i*n] =0;
				for(k = 0;k<n;k++){
					b[i] += bet2[i+k*n] + gam4[i+k*n];
					A1[i+i*n] += alp1[i+k*n] +gam1[i+k*n]; 
					A2[i] += alp2[i+k*n] + gam2[i+k*n]+del2[i+k*n]; 
					B[i+i*n] += bet1[i+k*n] + gam3[i+k*n] + del1[i+k*n]; 
				}
				for(j = i+1; j < n; j++){
					A1[i+j*n] = -(alp1[i+j*n]+gam1[i+j*n]);
					A1[j+i*n] = -(alp1[j+i*n]+gam1[j+i*n]);
					B[i+j*n] = -(bet1[i+j*n] + gam3[i+j*n] + del1[i+j*n]);
					B[j+i*n] = -(bet1[j+i*n] + gam3[j+i*n] + del1[j+i*n]);
				}
			}
			/*Start: Renew X*/
			//XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
			dgesvd_(&jobu, &jobvt, &n, &n, A1, &n, singular_value, U_svd, &n, Vt_svd, &n, work, &lwork, &info);
			if(info != 0){
				Rprintf("Error (SVD) code %d \n",info);
				break;
			}
			for(s =0; s<p;s++){//s
				//bw1 := B%*%y_s
				for(i=0;i<n;i++){
					bw1[i] = 0;
					for(j=0;j<n;j++){
						bw1[i] += B[i+j*n]*X[j+s*n];
					}
				}
				//bw2 := U^T%*%bw1 := U^T%*%(B%*%y_s)
				for(i=0;i<n;i++){
					bw2[i] = 0;
					for(j=0;j<n;j++){
						bw2[i] += U_svd[j+i*n]*bw1[j];
					}
				}
				//bw1 := S^(-1)%*%bw2
				for(i=0;i<n;i++){
					bw1[i] = 0;
					if(singular_value[i]>0.00000001){
						bw1[i] = bw2[i]/singular_value[i];
					}
				}
				//Rprintf("minsv= %f \n",singular_value[n-2]);
				//X[,s] := Vt^T%*%bw1
				for(i=0;i<n;i++){
					Y[i+s*n] = X[i+s*n];
					X[i+s*n] = 0;
					for(j=0;j<n;j++){
						X[i+s*n] += Vt_svd[j+i*n]*bw1[j];
					}
				}
			}//s
			//XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
			/*End: Renew X*/
			/*Start: Renew r*/
			for(i=0;i<n;i++){
				q[i] = r[i];
				r[i] = b[i]/A2[i];
			}
			/*End: Renew r*/
			/*Do accelaration by relaxed update*/
			if(*relax==1){
				//Relaxed update
				for(i=0;i<n;i++){
					tmp = r[i];
					if(2*tmp-q[i]>0){
						r[i] = 2*tmp-q[i];
					}
					for(s=0;s<p;s++){
						tmp =X[i+s*n];
						X[i+s*n] = 2*tmp - Y[i+s*n];
					}
				}
			}
			if(*dilation==1){
				/*Start comp. stress*/
				//compute centre distances
				for(i=0; i<(n-1); i++){
					for(j =i+1; j <n; j++){
						tmpcd =0;
						for(s=0; s<p; s++){
							tmpcd += pow(X[i+s*n]-X[j+s*n],2);
						}
						CD[i+j*n] = sqrt(tmpcd);
						CD[j+i*n] = CD[i+j*n];
					}
				}
				tmp = 0;
				//tmpstr = 0;
				/*compute interval dist. and dil. const.*/
				eta = 0;
				for(i = 0; i <(n-1); i++){
					for(j = i+1; j < n; j++){
						ldist[i + j*n] = CD[i + j*n] - r[i] - r[j];
						if(ldist[i + j*n]<0){
							ldist[i + j*n] = 0;
						}
						ldist[j+i*n] = ldist[i + j*n];
						udist[i+j*n] = CD[i+j*n] + r[i] + r[j];
						udist[j+i*n] = udist[i+j*n]; 
						eta += pow(ldist[i+j*n],2)+pow(udist[i+j*n],2);
						//tmpstr += pow(udiss[i+j*n]-udist[i+j*n],2) + pow(ldiss[i+j*n]-ldist[i+j*n],2);
						tmp += ldist[i + j*n]*ldiss[i+j*n] +  udist[i + j*n]*udiss[i+j*n];
					}
				}
				dil = tmp/eta;
				/*end: compute interval dist. and dil. const.*/
				/*Update*/
				for(i=0;i<n;i++){
					tmp = r[i];
					r[i] = dil*tmp;
					for(s=0;s<p;s++){
						tmp =X[i+s*n];
						X[i+s*n] = dil*tmp;
					}
				}
				//compute interval dist. and stress
				tmpstr = 0;
				for(i = 0; i <(n-1); i++){
					for(j = i+1; j < n; j++){
						tmpcd = dil*CD[i+j*n];
						CD[i+j*n] = tmpcd;
						CD[j+i*n] = CD[i+j*n];
						tmp = dil*ldist[i+j*n];
						ldist[i + j*n] = tmp;
						ldist[j+i*n] = ldist[i + j*n];
						tmpstr += pow(ldiss[i+j*n]-ldist[i+j*n],2);
						tmp = dil*udist[i+j*n];
						udist[i+j*n] = tmp;
						udist[j+i*n] = udist[i+j*n]; 
						tmpstr += pow(udiss[i+j*n]-udist[i+j*n],2);
					}
				}
				/*end: compute interval dist. and stress*/
			}else{
				/*Start comp. stress*/
				//compute centre distances
				for(i=0; i<(n-1); i++){
					for(j =i+1; j <n; j++){
						tmpcd =0;
						for(s=0; s<p; s++){
							tmpcd += pow(X[i+s*n]-X[j+s*n],2);
						}
						CD[i+j*n] = sqrt(tmpcd);
						CD[j+i*n] = CD[i+j*n];
					}
				}
				/*compute interval dist. and stress*/
				tmpstr = 0;
				for(i = 0; i <(n-1); i++){
					for(j = i+1; j < n; j++){
						ldist[i + j*n] = CD[i + j*n] - r[i] - r[j];
						if(ldist[i + j*n]<0){
							ldist[i + j*n] = 0;
						}
						tmpstr += pow(ldiss[i+j*n]-ldist[i+j*n],2);
						ldist[j+i*n] = ldist[i + j*n];
						udist[i+j*n] = CD[i+j*n] + r[i] + r[j];
						udist[j+i*n] = udist[i+j*n]; 
						tmpstr += pow(udiss[i+j*n]-udist[i+j*n],2);
					}
				}
				/*end: compute interval dist. and stress*/
			}
			str[iter] = tmpstr; 
			if(iter % rp==0){
				Rprintf("iter %d stress = %f \n",iter, str[iter]);
			}
			/*if(str[iter-1]+1<str[iter]){
				Rprintf("iter (final) %d stress = %f \n",iter, str[iter]);
				Rprintf("Error\n");
				break;
			}*/
			if(fabs(str[iter-1]-str[iter])<EPS){
				Rprintf("iter (final) %d stress = %f \n",iter, str[iter]);
				Rprintf("converged\n");
				break;
			}
			if(iter==MAXIT){
				Rprintf("final (iter %d) stress = %f \n",iter, str[iter]);
				Rprintf("stopped after %d iterations \n",iter);
				break;
			}
		/*r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r*/
		}/*end for r*/
		/*Main loop===============================*/
		free(Y);
		free(q);
		free(CD);
		free(alp1);
		free(alp2);
		free(bet1);
		free(bet2);
		free(gam1);
		free(gam2);
		free(gam3);
		free(gam4);
		free(del1);
		free(del2);
		free(del3);
		free(A1);
		free(A2);
		free(B);
		free(b);
		//var for svd
		free(singular_value);
		free(U_svd);
		free(Vt_svd);
		free(work);
		//var for pinv
		free(bw1);
		free(bw2);
	}//end void
	
	void bidist(double *X, double *R, double *ldist, double *udist, int *N, int *P)
	{//start void	Rprintf("initial stress= %f \n", str[0]);
		/*start: def var*/
		int i, j, n, s, p;
		double tmp=0, tmpl=0, tmpu=0;
		n = *N;
		p = *P;
		/*start: compute interval dist.*/
		for(i = 0; i <(n-1); i++){
			for(j = i+1; j < n; j++){
				tmp =0;
				tmpl = 0;
				tmpu =0;
				for(s =0; s<p; s++){//s
					tmp = fabs(X[i+s*n]-X[j+s*n])-(R[i+s*n]+R[j+s*n]);
					if(tmp >0){
						tmpl += pow(tmp,2);
					}
					tmpu +=  pow(fabs(X[i+s*n]-X[j+s*n])+(R[i+s*n]+R[j+s*n]),2);
				}//s
				ldist[i+j*n] = sqrt(tmpl);
				udist[i+j*n] =sqrt(tmpu);
				ldist[j+i*n] = ldist[i+j*n];
				udist[j+i*n] =udist[i+j*n];
			}
		}
		/*end: compute interval dist.*/
	}//end void


void boxobj(double *str, double *X, double *R, double *ldiss, double *udiss, int *N, int *P)
	{//start void	Rprintf("initial stress= %f \n", str[0]);
		/*start: def var*/
		int i, j, n, s, p;
		double tmp=0, tmpl=0, tmpu=0;
		double *ldist, *udist;
		
		n = *N;
		p = *P;
		
		ldist = (double *)calloc(n*n,sizeof(double));
		udist = (double *)calloc(n*n,sizeof(double));
		
		str[0] = 0;
		/*start: compute interval dist.*/
		for(i = 0; i <(n-1); i++){
			for(j = i+1; j < n; j++){
				tmp =0;
				tmpl = 0;
				tmpu =0;
				for(s =0; s<p; s++){//s
					tmp = fabs(X[i+s*n]-X[j+s*n])-(R[i+s*n]+R[j+s*n]);
					if(tmp >0){
						tmpl += pow(tmp,2);
					}
					tmpu +=  pow(fabs(X[i+s*n]-X[j+s*n])+(R[i+s*n]+R[j+s*n]),2);
				}//s
				ldist[i+j*n] = sqrt(tmpl);
				udist[i+j*n] =sqrt(tmpu);
				ldist[j+i*n] = ldist[i+j*n];
				udist[j+i*n] =udist[i+j*n];
				str[0] += pow(ldist[i+j*n]-ldiss[i+j*n],2) + pow(udist[i+j*n]-udiss[i+j*n],2);
			}
		}
		/*end: compute interval dist.*/
		free(ldist);
		free(udist);
	}//end void
	

void boxgrad(double *grad, double *X, double *R, double *ldiss, double *udiss, int *N, int *P)
	{//start void	Rprintf("initial stress= %f \n", str[0]);
		/*start: def var*/
		int i, j, n, s, p;
		double tmp=0, tmpl=0, tmpu=0, aijs=0, bijs=0;
		double *ldist, *udist;
		
		n = *N;
		p = *P;
		
		ldist = (double *)calloc(n*n,sizeof(double));
		udist = (double *)calloc(n*n,sizeof(double));
		
		/*start: compute interval dist.*/
		for(i = 0; i <(n-1); i++){
			for(j = i+1; j < n; j++){
				tmp =0;
				tmpl = 0;
				tmpu =0;
				for(s =0; s<p; s++){//s
					tmp = fabs(X[i+s*n]-X[j+s*n])-(pow(R[i+s*n],2)+pow(R[j+s*n],2));
					if(tmp >0){
						tmpl += pow(tmp,2);
					}
					tmpu +=  pow(fabs(X[i+s*n]-X[j+s*n])+(pow(R[i+s*n],2)+pow(R[j+s*n],2)),2);
				}//s
				ldist[i+j*n] = sqrt(tmpl);
				udist[i+j*n] =sqrt(tmpu);
				ldist[j+i*n] = ldist[i+j*n];
				udist[j+i*n] =udist[i+j*n];
			}
		}
		/*end: compute interval dist.*/
		/*start: compute grad.*/
		for(s =0; s<p; s++){//s
			for(i = 0; i <n; i++){//i
				grad[i+s*n] = 0;
				grad[n*p + i+s*n] =0;
				for(j = 1; j < n; j++){//j
					if(i!=j){
						tmpl = ldist[i+j*n];
						tmpu =udist[i+j*n];
						if(ldist[i+j*n]==0){
							tmpl=0.001;
						}
						if(udist[i+j*n]==0){
							tmpu=0.001;
						}
						aijs = fabs(X[i+s*n]-X[j+s*n]) - (pow(R[i+s*n],2) + pow(R[j+s*n],2));
						bijs = fabs(X[i+s*n]-X[j+s*n]) + (pow(R[i+s*n],2) + pow(R[j+s*n],2));
						if(aijs>0){
							grad[n*p + i+s*n] += 2*(ldiss[i+j*n]-ldist[i+j*n])*aijs*R[i+s*n]/tmpl;
							if(X[i+s*n]-X[j+s*n]>0){
								grad[i+s*n] += -(ldiss[i+j*n]-ldist[i+j*n])*aijs/tmpl;
							}else if(X[i+s*n]-X[j+s*n] < 0){
								grad[i+s*n] += (ldiss[i+j*n]-ldist[i+j*n])*aijs/tmpl;
							}
						}
						grad[n*p + i+s*n] += -2*(udiss[i+j*n]-udist[i+j*n])*bijs*R[i+s*n]/tmpu;
						if(X[i+s*n]-X[j+s*n]>0){
							grad[i+s*n] += -(udiss[i+j*n]-udist[i+j*n])*bijs/tmpu;
						}else if(X[i+s*n]-X[j+s*n] < 0){
							grad[i+s*n] += (udiss[i+j*n]-udist[i+j*n])*bijs/tmpu;
						}
					}
				}//j
				tmp = grad[i+s*n];
				grad[i+s*n] =2*tmp;
				tmp = grad[n*p+i+s*n];
				grad[n*p+i+s*n] =2*tmp;
			}//i
		}//s
		/*end: compute grad.*/
		free(ldist);
		free(udist);
	}//end void

void mmBox(double *X, double *R, double *str, double *ldist, double *udist, double *ldiss, double *udiss, double *eps, int *relax, int *dilation, int *N, int *P, int *maxit, int *report)
	{//start void	Rprintf("initial stress= %f \n", str[0]);
		/*start: def var*/
		int i, j, iter, k, n, s, p, rp, MAXIT, info, lwork;
		double tmpris=0, tmprjs=0, tmpads=0, tmp=0, tmpstr =0, tmpl=0, tmpu=0, EPS=0, eta =0, dil =0;
		double *Y, *Q, *alp1,*alp2, *alp3, *alp4, *alp5,*bet1,*bet2, *bet3,*bet4,*bet5, *A1, *A2, *B, *b;
		//def var for SVD
		double *singular_value, *work, *U_svd, *Vt_svd;
		char jobu = 'A', jobvt = 'A';
		//def var for pinv
		double *bw1, *bw2;
		
		n = *N;
		p = *P;
		MAXIT = *maxit;
		lwork = 5*n;
		rp = *report;
		EPS = *eps;
		//Old variables
		Y =(double *)calloc(n*p,sizeof(double));
		Q =(double *)calloc(n*p,sizeof(double));
		//coefficients for MM
		alp1 =(double *)calloc(n*n,sizeof(double));
		alp2 =(double *)calloc(n*n,sizeof(double));
		alp3 =(double *)calloc(n*n,sizeof(double));
		alp4 =(double *)calloc(n*n,sizeof(double));
		alp5 =(double *)calloc(n*n,sizeof(double));
		bet1 =(double *)calloc(n*n,sizeof(double));
		bet2 =(double *)calloc(n*n,sizeof(double));
		bet3 =(double *)calloc(n*n,sizeof(double));
		bet4 =(double *)calloc(n*n,sizeof(double));
		bet5 =(double *)calloc(n*n,sizeof(double));
		A1 =(double *)calloc(n*n,sizeof(double));
		A2 =(double *)calloc(n,sizeof(double));
		B =(double *)calloc(n*n,sizeof(double));
		b =(double *)calloc(n,sizeof(double));
		//var for svd
		singular_value = (double *)calloc(n,sizeof(double));
		U_svd = (double *)calloc(n*n,sizeof(double));
		Vt_svd = (double *)calloc(n*n,sizeof(double));
		work = (double *)calloc(lwork,sizeof(double));
		//var for pinv
		bw1=(double *)calloc(n,sizeof(double));
		bw2=(double *)calloc(n,sizeof(double));
		/*end: def var*/
		
		/*start: compute interval dist. and stress*/
		tmpstr = 0;
		for(i = 0; i <(n-1); i++){
			for(j = i+1; j < n; j++){
				eta += pow(ldiss[i+j*n],2)+ pow(udiss[i+j*n],2);
				tmp =0;
				tmpl = 0;
				tmpu =0;
				for(s =0; s<p; s++){//s
					tmp = fabs(X[i+s*n]-X[j+s*n])-(R[i+s*n]+R[j+s*n]);
					if(tmp >0){
						tmpl += pow(tmp,2);
					}
					tmpu +=  pow(fabs(X[i+s*n]-X[j+s*n])+(R[i+s*n]+R[j+s*n]),2);
				}//s
				ldist[i+j*n] = sqrt(tmpl);
				udist[i+j*n] =sqrt(tmpu);
				ldist[j+i*n] = ldist[i+j*n];
				udist[j+i*n] =udist[i+j*n];
				tmpstr += pow(ldist[i+j*n]-ldiss[i+j*n],2) + pow(udist[i+j*n]-udiss[i+j*n],2);
			}
		}
		/*end: compute interval dist. and stress*/
		str[0] = tmpstr; 
		Rprintf("initial value %f \n", str[0]);
		/*Main loop===============================*/
		for(iter=1;iter<MAXIT+1;iter++){/*start for r*/
		/*r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r*/
			for(s =0; s<p;s++){//s
				for(i = 0; i <(n-1); i++){//i
					for(j = i+1; j < n; j++){//j
						tmpris = R[i+s*n];
						tmprjs = R[j+s*n];
						tmpads = fabs(X[i+s*n]-X[j+s*n]);
						tmp = tmpads;
						if(tmpris==0){
							tmpris = 0.000001;
						}
						if(tmprjs==0){
							tmprjs = 0.000001;
						}
						if(tmpads==0){
							tmpads = 0.000001;
						}
						alp1[i+j*n] = 1+(R[i+s*n]+R[j+s*n])/tmpads;
						alp1[j+i*n] = alp1[i+j*n];
						alp2[i+j*n] = (tmp + R[i+s*n]+R[j+s*n])/R[i+s*n];
						alp2[j+i*n] = (tmp + R[i+s*n]+R[j+s*n])/R[j+s*n];
						alp3[i+j*n] = 2;
						alp3[j+i*n] = 2;
						alp4[i+j*n] = 2*(1+R[j+s*n]/tmpris);
						alp4[j+i*n] = 2*(1+R[i+s*n]/tmprjs);
						bet1[i+j*n] = 0;
						bet2[i+j*n] = 0;
						if(udist[i+j*n]>0){
							bet2[i+j*n] = udiss[i+j*n]*(tmp+R[i+s*n]+R[j+s*n])/(udist[i+j*n]);
							if(tmp>0){
								bet1[i+j*n] = udiss[i+j*n]*(tmp+R[i+s*n]+R[j+s*n])/(tmp*udist[i+j*n]);
							}
						}
						bet1[j+i*n] = bet1[i+j*n];
						bet2[j+i*n] = bet2[i+j*n];
						bet3[i+j*n] = 0;
						alp5[i+j*n] = 0;
						alp5[j+i*n] = 0;
						bet5[i+j*n] = 0;
						if(tmp>R[i+s*n]+R[j+s*n]){
							bet4[i+j*n] = tmp + R[i+s*n]+R[j+s*n];
							if(tmp>0){
								bet3[i+j*n] = (tmp + R[i+s*n]+R[j+s*n])/tmp;
							}
							if(ldist[i+j*n]>0){
								alp5[i+j*n] = ldiss[i+j*n]*fmax(0,tmp-(R[i+s*n]+R[j+s*n]))/(tmpris*ldist[i+j*n]);
								alp5[j+i*n] = ldiss[i+j*n]*fmax(0,tmp-(R[i+s*n]+R[j+s*n]))/(tmprjs*ldist[i+j*n]);
								if(tmp>0){
									bet5[i+j*n] = ldiss[i+j*n]*fmax(0,tmp-(R[i+s*n]+R[j+s*n]))/(tmp*ldist[i+j*n]);
								}
							}
						}else{
							bet4[i+j*n] = 2*(R[i+s*n]+R[j+s*n]);
							if(tmp>0){
								bet3[i+j*n] = 2;
							}
						}
						bet3[j+i*n] = bet3[i+j*n];
						bet4[j+i*n] = bet4[i+j*n];
						bet5[j+i*n] = bet5[i+j*n];
					}//i
				}//j
				//Comp. A1, A2, B, b
				for(i = 0; i <n; i++){
					b[i] =0;
					A1[i+i*n] = 0;
					A2[i] = 0;
					B[i+i*n] =0;
					alp1[i+i*n] = 0;
					alp2[i+i*n] = 0;
					alp3[i+i*n] = 0;
					alp4[i+i*n] = 0;
					alp5[i+i*n] = 0;
					bet1[i+i*n] = 0;
					bet2[i+i*n] = 0;
					bet3[i+i*n] = 0;
					bet4[i+i*n] = 0;
					bet5[i+i*n] = 0;
					for(k = 0;k<n;k++){
						b[i] += bet2[i+k*n] + bet4[i+k*n];
						A1[i+i*n] += alp1[i+k*n] +alp3[i+k*n]; 
						A2[i] += alp2[i+k*n] + alp4[i+k*n]+alp5[i+k*n]; 
						B[i+i*n] += bet1[i+k*n] + bet3[i+k*n] + bet5[i+k*n];
					}
					for(j = i+1; j < n; j++){
						A1[i+j*n] = -(alp1[i+j*n]+alp3[i+j*n]);
						A1[j+i*n] = -(alp1[j+i*n]+alp3[j+i*n]);
						B[i+j*n] = -(bet1[i+j*n] + bet3[i+j*n] + bet5[i+j*n]);
						B[j+i*n] = -(bet1[j+i*n] + bet3[j+i*n] + bet5[j+i*n]);
					}
				}
				/*Start: Renew X[,s]*/
				//XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
				dgesvd_(&jobu, &jobvt, &n, &n, A1, &n, singular_value, U_svd, &n, Vt_svd, &n, work, &lwork, &info);
				if(info != 0){
					//Rprintf("Error (SVD) code %d \n",info);
					break;
				}
				//bw1 := B%*%y_s
				for(i=0;i<n;i++){
					bw1[i] = 0;
					for(j=0;j<n;j++){
						bw1[i] += B[i+j*n]*X[j+s*n];
					}
				}
				//bw2 := U^T%*%bw1 := U^T%*%(B%*%y_s)
				for(i=0;i<n;i++){
					bw2[i] = 0;
					for(j=0;j<n;j++){
						bw2[i] += U_svd[j+i*n]*bw1[j];
					}
				}
				//bw1 := S^(-1)%*%bw2
				for(i=0;i<n;i++){
					bw1[i] = 0;
					if(singular_value[i]>0.00001){
						bw1[i] = bw2[i]/singular_value[i];
					}
				}
				//Rprintf("minsv= %f \n",singular_value[n-2]);
				//X[,s] := Vt^T%*%bw1
				for(i=0;i<n;i++){
					Y[i+s*n] = X[i+s*n];
					X[i+s*n] = 0;
					for(j=0;j<n;j++){
						X[i+s*n] += Vt_svd[j+i*n]*bw1[j];
					}
				}
				/*End: Renew X[,s]*/
				/*Start: Renew R[,s]*/
				for(i=0;i<n;i++){
					Q[i+s*n] = R[i+s*n];
					R[i+s*n] = b[i]/A2[i];
				}
				/*End: Renew R[,s]*/
			}//s
			if(info != 0){
				Rprintf("Error (SVD) code %d \n",info);
				break;
			}
			/*Do accelaration by relaxed update*/
			if(*relax==1){
				//Relaxed update
				for(i=0;i<n;i++){
					for(s=0;s<p;s++){
						tmp =X[i+s*n];
						X[i+s*n] = 2*tmp - Y[i+s*n];
						tmp = R[i+s*n];
						if(2*tmp-Q[i+s*n]>0){
							R[i+s*n] = 2*tmp-Q[i+s*n];
						}
					}
				}
			}
			if(*dilation==1){
				//Dil. update
				/*start: compute interval dist. and stress*/
				dil = 0;
				eta = 0;
				for(i = 0; i <(n-1); i++){
					for(j = i+1; j < n; j++){
						tmp =0;
						tmpl = 0;
						tmpu =0;
						for(s =0; s<p; s++){//s
							tmp = fabs(X[i+s*n]-X[j+s*n])-(R[i+s*n]+R[j+s*n]);
							if(tmp >0){
								tmpl += pow(tmp,2);
							}
							tmpu +=  pow(fabs(X[i+s*n]-X[j+s*n])+(R[i+s*n]+R[j+s*n]),2);
						}//s
						ldist[i+j*n] = sqrt(tmpl);
						udist[i+j*n] =sqrt(tmpu);
						ldist[j+i*n] = ldist[i+j*n];
						udist[j+i*n] =udist[i+j*n];
						eta += pow(ldist[i+j*n],2)+pow(udist[i+j*n],2);
						dil += udiss[i+j*n]*udist[i+j*n] + ldiss[i+j*n]*ldist[i+j*n];
					}
				}
				tmp = dil/eta;
				dil = tmp;
				for(i=0;i<n;i++){
					for(s=0;s<p;s++){
						tmp =X[i+s*n];
						X[i+s*n] = dil*tmp;
						tmp = R[i+s*n];
						R[i+s*n] = dil*tmp;
					}
				}
				tmpstr = 0;
				for(i = 0; i <(n-1); i++){
					for(j = i+1; j < n; j++){
						tmpl = dil*ldist[i+j*n];
						ldist[i+j*n] = tmpl;
						tmpu = dil*udist[i+j*n];
						udist[i+j*n] = tmpu;
						ldist[j+i*n] = ldist[i+j*n];
						udist[j+i*n] =udist[i+j*n];
						tmpstr += pow(ldist[i+j*n]-ldiss[i+j*n],2) + pow(udist[i+j*n]-udiss[i+j*n],2);
					}
				}
				/*end: compute interval dist. and stress*/
			}else{
				/*Start comp. stress*/
				/*start: compute interval dist. and stress*/
				tmpstr = 0;
				for(i = 0; i <(n-1); i++){
					for(j = i+1; j < n; j++){
						tmp =0;
						tmpl = 0;
						tmpu =0;
						for(s =0; s<p; s++){//s
							tmp = fabs(X[i+s*n]-X[j+s*n])-(R[i+s*n]+R[j+s*n]);
							if(tmp >0){
								tmpl += pow(tmp,2);
							}
							tmpu +=  pow(fabs(X[i+s*n]-X[j+s*n])+(R[i+s*n]+R[j+s*n]),2);
						}//s
						ldist[i+j*n] = sqrt(tmpl);
						udist[i+j*n] =sqrt(tmpu);
						ldist[j+i*n] = ldist[i+j*n];
						udist[j+i*n] =udist[i+j*n];
						tmpstr += pow(ldist[i+j*n]-ldiss[i+j*n],2) + pow(udist[i+j*n]-udiss[i+j*n],2);
					}
				}
				/*end: compute interval dist. and stress*/
			}/*end if dilation==1 else*/
			//XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
			str[iter] = tmpstr; 
			if(iter % rp==0){
				Rprintf("iter %d stress = %f \n",iter, str[iter]);
			}
			if(fabs(str[iter-1]-str[iter])<EPS){
				Rprintf("iter (final) %d stress = %f \n",iter, str[iter]);
				Rprintf("converged\n");
				break;
			}
			if(iter==MAXIT){
				Rprintf("final (iter %d) stress = %f \n",iter, str[iter]);
				Rprintf("stopped after %d iterations \n",iter);
				break;
			}
		/*r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r-r*/
		}/*end for r*/
		/*Main loop===============================*/
		free(Y);
		free(Q);
		free(alp1);
		free(alp2);
		free(alp3);
		free(alp4);
		free(alp5);
		free(bet1);
		free(bet2);
		free(bet3);
		free(bet4);
		free(bet5);
		free(A1);
		free(A2);
		free(B);
		free(b);
		//var for svd
		free(singular_value);
		free(U_svd);
		free(Vt_svd);
		free(work);
		//var for pinv
		free(bw1);
		free(bw2);
	}//end void
