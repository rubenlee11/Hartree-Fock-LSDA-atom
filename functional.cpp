#include<math.h>
#include<stdio.h>
#include"common.h"
/****************************************************
at present, only PW92 LSDA is supported.
****************************************************/  

static double Cal_xc_density(double ns,double zeta);
static void Cal_xc_potential(double *v_xc,double rs,double zeta);
static double Cal_f(double zeta);
static double Cal_fprim(double zeta);
static double Cal_G(int flag,int idx,double rs);
static double Cal_e_x_0(double rs);

double Cal_xc_energy(){
	/****************************************************
	calculate xc energy
	****************************************************/ 
	int idx_cheb,idx_leb;
	double **F_xc,n_up,n_down,ns,rs,zeta;	
	F_xc = (double **)malloc(sizeof(double*)*N_cheb);
	for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
		F_xc[idx_cheb] = (double*)malloc(sizeof(double)*N_leb);
	}
	
	for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
		for(idx_leb = 0;idx_leb<N_leb;idx_leb++){
			n_up = rho[0][idx_cheb][idx_leb];
			n_down = rho[1][idx_cheb][idx_leb];
			if(n_up>1e-15||n_down>1e-15){
				ns = n_up+n_down;
				zeta = (n_up-n_down)/ns;				
				F_xc[idx_cheb][idx_leb] = Cal_xc_density(ns,zeta);
			}
			else{
				F_xc[idx_cheb][idx_leb] = 0.0;
			}
		}
	}	
	return Quadrature(F_xc);
}

void Cal_xc_matrix(double *result,int mu,int nu){
	/****************************************************
	calculate xc potential in basis using Beck grid
	mu: orbital phi_mu......
	
	<mu|v_xc^up|nu>,<mu|v_xc^down|nu>
	****************************************************/  
	int idx_cheb,idx_leb;
	double ***F_xc,n_up,n_down,ns,rs,zeta;
	double *v_xc;
	
	F_xc = (double ***)malloc(sizeof(double**)*2);	
	F_xc[0] = (double **)malloc(sizeof(double*)*N_cheb);	
	F_xc[1] = (double **)malloc(sizeof(double*)*N_cheb);	
	for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
		F_xc[0][idx_cheb] = (double *)malloc(sizeof(double)*N_leb);	
		F_xc[1][idx_cheb] = (double *)malloc(sizeof(double)*N_leb);	
	}
	
	v_xc = (double *)malloc(sizeof(double)*2);	
	
	for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
		for(idx_leb = 0;idx_leb<N_leb;idx_leb++){
			n_up = rho[0][idx_cheb][idx_leb];
			n_down = rho[1][idx_cheb][idx_leb];
			if(n_up>1e-15||n_down>1e-15){
				ns = n_up+n_down;
				zeta = (n_up-n_down)/ns;
				rs = pow(3.0/(4.0*PI*ns),1.0/3.0);
				Cal_xc_potential(v_xc,rs,zeta);								
				F_xc[0][idx_cheb][idx_leb] = v_xc[0]*phi[mu][idx_cheb][idx_leb]*phi[nu][idx_cheb][idx_leb];
				F_xc[1][idx_cheb][idx_leb] = v_xc[1]*phi[mu][idx_cheb][idx_leb]*phi[nu][idx_cheb][idx_leb];				
			}
			else{
				F_xc[0][idx_cheb][idx_leb] = F_xc[1][idx_cheb][idx_leb] = 0.0;
			}
		}		
	}
	result[0] = Quadrature(F_xc[0]);
	result[1] = Quadrature(F_xc[1]);
	free(v_xc);
	
	for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
		free(F_xc[0][idx_cheb]);	
		free(F_xc[1][idx_cheb]);	
	}		
	free(F_xc[0]);	
	free(F_xc[1]);
	free(F_xc);
}

static double Cal_xc_density(double ns,double zeta){
	/****************************************************
	calculate xc energy density given total density rs and polarization zeta
	e_x_0: e(ns,0)......
	****************************************************/  	
	double e_x,e_x_0,e_x_1;
	double e_c,e_c_0,e_c_1,alpha;
	double rs,f;
	rs = pow(3.0/(4.0*PI*ns),1.0/3.0);
	f = Cal_f(zeta);
	// exchange energy per electron
	e_x_0 = Cal_e_x_0(rs);
	e_x_1 = pow(2.0,1.0/3.0)*e_x_0;
	e_x = e_x_0+(e_x_1-e_x_0)*f;
	
	// correlation energy per electron
	e_c_0 = Cal_G(0,0,rs);
	e_c_1 = Cal_G(0,1,rs);
	alpha = -Cal_G(0,2,rs);

	double z4 = pow(zeta,4);
	
	e_c = e_c_0+alpha*f/1.709921*(1.0-z4)+(e_c_1-e_c_0)*f*z4;
	
	return (e_x+e_c)*ns;
}	

static void Cal_xc_potential(double *v_xc,double rs,double zeta){
	/****************************************************
	calculate xc potential given total density rs and polarization zeta
	****************************************************/  		
	double e_x,e_x_0,e_x_1;
	double e_c,e_c_0,e_c_1,alpha;
	double f,fprim;
	double *v_x,*v_c;
	double PexPrs,PexPz,PecPrs,PecPz;
	double z4 = pow(zeta,4);;
	
	
	v_x = (double *)malloc(sizeof(double)*2);
	v_c = (double *)malloc(sizeof(double)*2);
		
	f = Cal_f(zeta);
	fprim = Cal_fprim(zeta);
	// exchange potential
	e_x_0 = Cal_e_x_0(rs);
	e_x_1 = pow(2.0,1.0/3.0)*e_x_0;
	e_x = e_x_0+(e_x_1-e_x_0)*f;

	PexPrs = -1.0*e_x/rs;
	PexPz = (e_x_1-e_x_0)*fprim;
	
	v_x[0] = e_x - rs/3.0*PexPrs-(zeta-1.0)*PexPz;
	v_x[1] = e_x - rs/3.0*PexPrs-(zeta+1.0)*PexPz;
	// correlation potential	

	e_c_0 = Cal_G(0,0,rs);
	e_c_1 = Cal_G(0,1,rs);
	alpha = -Cal_G(0,2,rs);
	e_c = e_c_0+alpha*f/1.709921*(1.0-z4)+(e_c_1-e_c_0)*f*z4;	
	
	PecPrs = Cal_G(1,0,rs)*(1-f*z4) + Cal_G(1,1,rs)*f*z4 - Cal_G(1,2,rs)*f/1.709921*(1-z4);
	PecPz = 4.0*pow(zeta,3)*f*(e_c_1-e_c_0-alpha/1.709921)+fprim*(z4*(e_c_1-e_c_0)+(1-z4)*alpha/1.709921);
	
	v_c[0] = e_c - rs/3.0*PecPrs-(zeta-1.0)*PecPz;
	v_c[1] = e_c - rs/3.0*PecPrs-(zeta+1.0)*PecPz;
	
	v_xc[0] = v_x[0]+v_c[0];
	v_xc[1] = v_x[1]+v_c[1];
}

static double Cal_f(double zeta){
	return (pow(1.0+zeta,4.0/3.0)+pow(1-zeta,4.0/3.0)-2.0)/(pow(2.0,4.0/3.0)-2.0);
}

static double Cal_fprim(double zeta){
	return 4.0/3.0*(pow(1+zeta,1.0/3.0)-pow(1-zeta,1.0/3.0))/(pow(2.0,4.0/3.0)-2.0);
}

static double Cal_G(int flag,int idx,double rs){
	/****************************************************
	Calculate the G function in PW92 functional.
	flag = 0: G
		   1: dG/drs
		   
	idx = 0: e_c(r,0)
		1: e_c(r,1)
		2: -alpha_c(r)
	****************************************************/  
	double p[3],A[3],alpha_1[3],beta_1[3],beta_2[3],beta_3[3],beta_4[3],c_1[3],c_2[3],c_3[3],d_0[3],d_1[3];
	p[0] = 1.00;p[1] = 1.00;p[2] = 1.0;
	A[0] = 0.031091; A[1] = 0.015545;A[2] = 0.016887;
	alpha_1[0] = 0.21370;alpha_1[1] = 0.20548;alpha_1[2] = 0.11125;
	beta_1[0] = 7.5957;beta_1[1] = 14.1189;beta_1[2] = 10.357;
	beta_2[0] = 3.5876;beta_2[1] = 6.1977;beta_2[2] = 3.6231;
	beta_3[0] = 1.6382;beta_3[1] = 3.3662;beta_3[2] = 0.88026;
	beta_4[0] = 0.49294;beta_4[1] = 0.62517;beta_4[2] = 0.49671;
	c_1[0] = 0.046644;c_1[1] = 0.025599;c_1[2] = 0.035475;
	c_2[0] = 0.00664;c_2[1] = 0.00319;c_2[2] = 0.00188;
	c_3[0] = 0.01043;c_3[1] = 0.00384;c_3[2] = 0.00521;
	d_0[0] = 0.4335;d_0[1] = 0.3287;d_0[2] = 0.2240;
	d_1[0] = 1.4408;d_1[1] = 1.7697;d_1[2] = 0.3969;
	
	if(flag == 0){
		double Q1;
		
		Q1 = 2.0*A[idx]*(beta_1[idx]*pow(rs,0.5)+beta_2[idx]*rs+beta_3[idx]*pow(rs,1.5)+beta_4[idx]*pow(rs,p[idx]+1.0));	
		
		return -2.0*A[idx]*(1.0+alpha_1[idx]*rs)*log(1.0+1.0/Q1);	
	}
	else if(flag == 1){
		double Q0,Q1,Q1p;
		
		Q0 = -2.0*A[idx]*(1.0+alpha_1[idx]*rs);
		Q1 = 2.0*A[idx]*(beta_1[idx]*pow(rs,0.5)+beta_2[idx]*rs+beta_3[idx]*pow(rs,1.5)+beta_4[idx]*pow(rs,p[idx]+1.0));
		Q1p = A[idx]*(beta_1[idx]*pow(rs,-0.5)+2.0*beta_2[idx]+3.0*beta_3[idx]*pow(rs,0.5)+2.0*(p[idx]+1.0)*beta_4[idx]*pow(rs,p[idx]));
		return -2.0*A[idx]*alpha_1[idx]*log(1.0+1.0/Q1)-Q0*Q1p/(Q1*Q1+Q1);
	}	
	else{
		return -1.0;
	}
}

static double Cal_e_x_0(double rs){
	return -3.0/(4.0*PI*rs)*pow(9.0*PI/4.0,1.0/3.0);
}