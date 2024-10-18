#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

/****************************************************
reference: 
    lebedev rule:
    https://people.sc.fsu.edu/~jburkardt/cpp_src/sphere_lebedev_rule/sphere_lebedev_rule.html
    chebyshev rule:
    https://people.math.sc.edu/Burkardt/cpp_src/chebyshev2_rule/chebyshev2_rule.html
****************************************************/  

static void Cal_phi();
static int gen_oh ( int code, double a, double b, double v, double *x, double *y, double *z, double *w );
static void ld0026 ( double *x, double *y, double *z, double *w );

int N_leb = 26;
int N_cheb = 25;
double *r_cheb,*x_cheb,*w_cheb;  
double *x_leb,*y_leb,*z_leb,*w_leb;        
double rm = 1.0;

void Init_Quadrature(){
    /****************************************************
    Init chebyshev grid for radial integration, lebedev grid for sphere integration
    ****************************************************/ 

	int idx_cheb,idx_leb;

    // allocate memory for grid
    r_cheb = (double *)malloc(sizeof(double)*N_cheb);
    x_cheb = (double *)malloc(sizeof(double)*N_cheb);
    w_cheb = (double *)malloc(sizeof(double)*N_cheb);

    x_leb = (double *)malloc(sizeof(double)*N_leb);
    y_leb = (double *)malloc(sizeof(double)*N_leb);
    z_leb = (double *)malloc(sizeof(double)*N_leb);
    w_leb = (double *)malloc(sizeof(double)*N_leb);

    rho = (double ***)malloc(sizeof(double**)*2);
    rho[0] = (double**)malloc(sizeof(double*)*N_cheb);
    rho[1] = (double**)malloc(sizeof(double*)*N_cheb);        
    for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
        rho[0][idx_cheb] = (double*)malloc(sizeof(double)*N_leb);
        rho[1][idx_cheb] = (double*)malloc(sizeof(double)*N_leb);
    }

    vxc = (double ***)malloc(sizeof(double**)*2);
    vxc[0] = (double**)malloc(sizeof(double*)*N_cheb);
    vxc[1] = (double**)malloc(sizeof(double*)*N_cheb);        
    for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
        vxc[0][idx_cheb] = (double*)malloc(sizeof(double)*N_leb);
        vxc[1][idx_cheb] = (double*)malloc(sizeof(double)*N_leb);
    }    

    phi = (double ***)malloc(sizeof(double**)*Aorbnum);
    for(int k = 0;k<Aorbnum;k++){
        phi[k] = (double **)malloc(sizeof(double*)*N_cheb);
        for(idx_cheb = 0;idx_cheb<N_cheb;idx_cheb++){
            phi[k][idx_cheb] = (double*)malloc(sizeof(double)*N_leb);     
        }  
    }
  
    // chebyshev rule of order 25 
    w_cheb[0]=0.00175555717117274;
    w_cheb[1]=0.00692020208309326;
    w_cheb[2]=0.0151937843462012;
    w_cheb[3]=0.0260954734326124;
    w_cheb[4]=0.0389917027986525;
    w_cheb[5]=0.0531329904694531;
    w_cheb[6]=0.0676974962070771;
    w_cheb[7]=0.0818387838778785;
    w_cheb[8]=0.0947350132439175;
    w_cheb[9]=0.105636702330329;
    w_cheb[10]=0.113910284593437;
    w_cheb[11]=0.119074929505357;
    w_cheb[12]=0.12083048667653;
    w_cheb[13]=0.119074929505358;
    w_cheb[14]=0.113910284593436;
    w_cheb[15]=0.105636702330329;
    w_cheb[16]=0.0947350132439181;
    w_cheb[17]=0.0818387838778785;
    w_cheb[18]=0.0676974962070765;
    w_cheb[19]=0.053132990469454;
    w_cheb[20]=0.0389917027986522;
    w_cheb[21]=0.0260954734326127;
    w_cheb[22]=0.015193784346201;
    w_cheb[23]=0.00692020208309326;
    w_cheb[24]=0.00175555717117274;

    x_cheb[0]=-0.992708874098054;
    x_cheb[1]=-0.970941817426052;
    x_cheb[2]=-0.935016242685414;
    x_cheb[3]=-0.885456025653209;
    x_cheb[4]=-0.822983865893656;
    x_cheb[5]=-0.7485107481711;
    x_cheb[6]=-0.663122658240795;
    x_cheb[7]=-0.568064746731155;
    x_cheb[8]=-0.464723172043768;
    x_cheb[9]=-0.354604887042535;
    x_cheb[10]=-0.239315664287557;
    x_cheb[11]=-0.120536680255322;
    x_cheb[12]=-0.000;
    x_cheb[13]=0.120536680255323;
    x_cheb[14]=0.239315664287557;
    x_cheb[15]=0.354604887042536;
    x_cheb[16]=0.464723172043768;
    x_cheb[17]=0.568064746731156;
    x_cheb[18]=0.663122658240795;
    x_cheb[19]=0.748510748171101;
    x_cheb[20]=0.822983865893656;
    x_cheb[21]=0.885456025653209;
    x_cheb[22]=0.935016242685414;
    x_cheb[23]=0.970941817426051;
    x_cheb[24]=0.992708874098054;    

    for(idx_cheb = 0;idx_cheb <N_cheb;idx_cheb++){
        r_cheb[idx_cheb] = rm*(x_cheb[idx_cheb]+1)/(1-x_cheb[idx_cheb]);
    }
    // lebedev rule of order 26
    ld0026(x_leb,y_leb,z_leb,w_leb);	
	Cal_phi();
}

double Quadrature(double **F){
    /****************************************************
    calculate the integral of the given function F over Beck grid
    ****************************************************/
	
	double value,x,r;
	double Int_sphere = 0.0;
    int idx_cheb,idx_leb;
	
	value = 0.0;
    for(idx_cheb = 0;idx_cheb < N_cheb;idx_cheb++){
		Int_sphere = 0.0;
		for(idx_leb = 0; idx_leb < N_leb;idx_leb++){
			Int_sphere += w_leb[idx_leb]*F[idx_cheb][idx_leb];		
		}
		Int_sphere = 4.0*PI*Int_sphere;
		
		x = x_cheb[idx_cheb];
		r = r_cheb[idx_cheb];
		value += w_cheb[idx_cheb]*Int_sphere*pow(r,2)*2.0*rm/(pow(1-x,2)*sqrt(1-x*x));	
	}	
	
	return value;
}

static int gen_oh ( int code, double a, double b, double v, double *x, double *y, double *z, double *w )
	/******************************************************************************/
	/*
	  Purpose:

		GEN_OH generates points under OH symmetry.

	  Discussion:

		Given a point on a sphere, specified by A and B, this routine generates
		all the equivalent points under OH symmetry, making grid points with
		weight V.

		The variable NUM is increased by the number of different points
		generated.

		Depending on CODE, there are from 6 to 48 different but equivalent
		points that are generated:

		  CODE=1:   (0,0,1) etc                                (  6 points)
		  CODE=2:   (0,A,A) etc, A=1/sqrt(2)                   ( 12 points)
		  CODE=3:   (A,A,A) etc, A=1/sqrt(3)                   (  8 points)
		  CODE=4:   (A,A,B) etc, B=sqrt(1-2 A^2)               ( 24 points)
		  CODE=5:   (A,B,0) etc, B=sqrt(1-A^2), A input        ( 24 points)
		  CODE=6:   (A,B,C) etc, C=sqrt(1-A^2-B^2), A, B input ( 48 points)

	  Modified:

		11 September 2010

	  Author:

		Dmitri Laikov

	  Reference:

		Vyacheslav Lebedev, Dmitri Laikov,
		A quadrature formula for the sphere of the 131st
		algebraic order of accuracy,
		Russian Academy of Sciences Doklady Mathematics,
		Volume 59, Number 3, 1999, pages 477-481.

	  Parameters:

		Input, int CODE, selects the symmetry group.

		Input, double A, B, information that may be needed to
		generate the coordinates of the points (for code = 5 or 6 only).

		Input, double V, the weight to be assigned the points.

		Output, double X[NUM], Y[NUM], Z[NUM], W[NUM], the coordinates
		and weights of the symmetric points generated on this call.

		Output, int GEN_OH, the number of points generated by this call.
	*/
	{
	  double c;
	  int n;

	  if ( code == 1 )
	  {
		a = 1.0;
		x[0] =   a; y[0] = 0.0; z[0] = 0.0; w[0] = v;
		x[1] =  -a; y[1] = 0.0; z[1] = 0.0; w[1] = v;
		x[2] = 0.0; y[2] =   a; z[2] = 0.0; w[2] = v;
		x[3] = 0.0; y[3] =  -a; z[3] = 0.0; w[3] = v;
		x[4] = 0.0; y[4] = 0.0; z[4] =   a; w[4] = v;
		x[5] = 0.0; y[5] = 0.0; z[5] =  -a; w[5] = v;
		n = 6;
	  }
	  else if ( code == 2 )
	  {
		a = 1.0 / sqrt ( 2.0 );
		x[ 0] = 0.0; y[ 0] =   a; z[ 0] =   a; w[ 0] = v;
		x[ 1] = 0.0; y[ 1] =   a; z[ 1] =  -a; w[ 1] = v;
		x[ 2] = 0.0; y[ 2] =  -a; z[ 2] =   a; w[ 2] = v;
		x[ 3] = 0.0; y[ 3] =  -a; z[ 3] =  -a; w[ 3] = v;
		x[ 4] =   a; y[ 4] = 0.0; z[ 4] =   a; w[ 4] = v;
		x[ 5] =   a; y[ 5] = 0.0; z[ 5] =  -a; w[ 5] = v;
		x[ 6] =  -a; y[ 6] = 0.0; z[ 6] =   a; w[ 6] = v;
		x[ 7] =  -a; y[ 7] = 0.0; z[ 7] =  -a; w[ 7] = v;
		x[ 8] =   a; y[ 8] =   a; z[ 8] = 0.0; w[ 8] = v;
		x[ 9] =   a; y[ 9] =  -a; z[ 9] = 0.0; w[ 9] = v;
		x[10] =  -a; y[10] =   a; z[10] = 0.0; w[10] = v;
		x[11] =  -a; y[11] =  -a; z[11] = 0.0; w[11] = v;
		n = 12;
	  }
	  else if ( code == 3 )
	  {
		a = 1.0 / sqrt ( 3.0 );
		x[0] =   a; y[0] =   a; z[0] =   a; w[0] = v;
		x[1] =   a; y[1] =   a; z[1] =  -a; w[1] = v;
		x[2] =   a; y[2] =  -a; z[2] =   a; w[2] = v;
		x[3] =   a; y[3] =  -a; z[3] =  -a; w[3] = v;
		x[4] =  -a; y[4] =   a; z[4] =   a; w[4] = v;
		x[5] =  -a; y[5] =   a; z[5] =  -a; w[5] = v;
		x[6] =  -a; y[6] =  -a; z[6] =   a; w[6] = v;
		x[7] =  -a; y[7] =  -a; z[7] =  -a; w[7] = v;
		n = 8;
	  }
	  else if ( code == 4 )
	  {
		b = sqrt ( 1.0 - 2.0 * a * a );
		x[ 0] =   a; y[ 0] =   a; z[ 0] =   b; w[ 0] = v;
		x[ 1] =   a; y[ 1] =   a; z[ 1] =  -b; w[ 1] = v;
		x[ 2] =   a; y[ 2] =  -a; z[ 2] =   b; w[ 2] = v;
		x[ 3] =   a; y[ 3] =  -a; z[ 3] =  -b; w[ 3] = v;
		x[ 4] =  -a; y[ 4] =   a; z[ 4] =   b; w[ 4] = v;
		x[ 5] =  -a; y[ 5] =   a; z[ 5] =  -b; w[ 5] = v;
		x[ 6] =  -a; y[ 6] =  -a; z[ 6] =   b; w[ 6] = v;
		x[ 7] =  -a; y[ 7] =  -a; z[ 7] =  -b; w[ 7] = v;
		x[ 8] =   a; y[ 8] =   b; z[ 8] =   a; w[ 8] = v;
		x[ 9] =   a; y[ 9] =  -b; z[ 9] =   a; w[ 9] = v;
		x[10] =   a; y[10] =   b; z[10] =  -a; w[10] = v;
		x[11] =   a; y[11] =  -b; z[11] =  -a; w[11] = v;
		x[12] =  -a; y[12] =   b; z[12] =   a; w[12] = v;
		x[13] =  -a; y[13] =  -b; z[13] =   a; w[13] = v;
		x[14] =  -a; y[14] =   b; z[14] =  -a; w[14] = v;
		x[15] =  -a; y[15] =  -b; z[15] =  -a; w[15] = v;
		x[16] =   b; y[16] =   a; z[16] =   a; w[16] = v;
		x[17] =  -b; y[17] =   a; z[17] =   a; w[17] = v;
		x[18] =   b; y[18] =   a; z[18] =  -a; w[18] = v;
		x[19] =  -b; y[19] =   a; z[19] =  -a; w[19] = v;
		x[20] =   b; y[20] =  -a; z[20] =   a; w[20] = v;
		x[21] =  -b; y[21] =  -a; z[21] =   a; w[21] = v;
		x[22] =   b; y[22] =  -a; z[22] =  -a; w[22] = v;
		x[23] =  -b; y[23] =  -a; z[23] =  -a; w[23] = v;
		n = 24;
	  }
	  else if ( code == 5 )
	  {
		b = sqrt ( 1.0 - a * a );
		x[ 0] =   a; y[ 0] =   b; z[ 0] = 0.0; w[ 0] = v;
		x[ 1] =   a; y[ 1] =  -b; z[ 1] = 0.0; w[ 1] = v;
		x[ 2] =  -a; y[ 2] =   b; z[ 2] = 0.0; w[ 2] = v;
		x[ 3] =  -a; y[ 3] =  -b; z[ 3] = 0.0; w[ 3] = v;
		x[ 4] =   b; y[ 4] =   a; z[ 4] = 0.0; w[ 4] = v;
		x[ 5] =   b; y[ 5] =  -a; z[ 5] = 0.0; w[ 5] = v;
		x[ 6] =  -b; y[ 6] =   a; z[ 6] = 0.0; w[ 6] = v;
		x[ 7] =  -b; y[ 7] =  -a; z[ 7] = 0.0; w[ 7] = v;
		x[ 8] =   a; y[ 8] = 0.0; z[ 8] =   b; w[ 8] = v;
		x[ 9] =   a; y[ 9] = 0.0; z[ 9] =  -b; w[ 9] = v;
		x[10] =  -a; y[10] = 0.0; z[10] =   b; w[10] = v;
		x[11] =  -a; y[11] = 0.0; z[11] =  -b; w[11] = v;
		x[12] =   b; y[12] = 0.0; z[12] =   a; w[12] = v;
		x[13] =   b; y[13] = 0.0; z[13] =  -a; w[13] = v;
		x[14] =  -b; y[14] = 0.0; z[14] =   a; w[14] = v;
		x[15] =  -b; y[15] = 0.0; z[15] =  -a; w[15] = v;
		x[16] = 0.0; y[16] =   a; z[16] =   b; w[16] = v;
		x[17] = 0.0; y[17] =   a; z[17] =  -b; w[17] = v;
		x[18] = 0.0; y[18] =  -a; z[18] =   b; w[18] = v;
		x[19] = 0.0; y[19] =  -a; z[19] =  -b; w[19] = v;
		x[20] = 0.0; y[20] =   b; z[20] =   a; w[20] = v;
		x[21] = 0.0; y[21] =   b; z[21] =  -a; w[21] = v;
		x[22] = 0.0; y[22] =  -b; z[22] =   a; w[22] = v;
		x[23] = 0.0; y[23] =  -b; z[23] =  -a; w[23] = v;
		n = 24;
	  }
	  else if ( code == 6 )
	  {
		c = sqrt ( 1.0 - a * a - b * b );
		x[ 0] =   a; y[ 0] =   b; z[ 0] =   c; w[ 0] = v;
		x[ 1] =   a; y[ 1] =   b; z[ 1] =  -c; w[ 1] = v;
		x[ 2] =   a; y[ 2] =  -b; z[ 2] =   c; w[ 2] = v;
		x[ 3] =   a; y[ 3] =  -b; z[ 3] =  -c; w[ 3] = v;
		x[ 4] =  -a; y[ 4] =   b; z[ 4] =   c; w[ 4] = v;
		x[ 5] =  -a; y[ 5] =   b; z[ 5] =  -c; w[ 5] = v;
		x[ 6] =  -a; y[ 6] =  -b; z[ 6] =   c; w[ 6] = v;
		x[ 7] =  -a; y[ 7] =  -b; z[ 7] =  -c; w[ 7] = v;
		x[ 8] =   a; y[ 8] =   c; z[ 8] =   b; w[ 8] = v;
		x[ 9] =   a; y[ 9] =   c; z[ 9] =  -b; w[ 9] = v;
		x[10] =   a; y[10] =  -c; z[10] =   b; w[10] = v;
		x[11] =   a; y[11] =  -c; z[11] =  -b; w[11] = v;
		x[12] =  -a; y[12] =   c; z[12] =   b; w[12] = v;
		x[13] =  -a; y[13] =   c; z[13] =  -b; w[13] = v;
		x[14] =  -a; y[14] =  -c; z[14] =   b; w[14] = v;
		x[15] =  -a; y[15] =  -c; z[15] =  -b; w[15] = v;
		x[16] =   b; y[16] =   a; z[16] =   c; w[16] = v;
		x[17] =   b; y[17] =   a; z[17] =  -c; w[17] = v;
		x[18] =   b; y[18] =  -a; z[18] =   c; w[18] = v;
		x[19] =   b; y[19] =  -a; z[19] =  -c; w[19] = v;
		x[20] =  -b; y[20] =   a; z[20] =   c; w[20] = v;
		x[21] =  -b; y[21] =   a; z[21] =  -c; w[21] = v;
		x[22] =  -b; y[22] =  -a; z[22] =   c; w[22] = v;
		x[23] =  -b; y[23] =  -a; z[23] =  -c; w[23] = v;
		x[24] =   b; y[24] =   c; z[24] =   a; w[24] = v;
		x[25] =   b; y[25] =   c; z[25] =  -a; w[25] = v;
		x[26] =   b; y[26] =  -c; z[26] =   a; w[26] = v;
		x[27] =   b; y[27] =  -c; z[27] =  -a; w[27] = v;
		x[28] =  -b; y[28] =   c; z[28] =   a; w[28] = v;
		x[29] =  -b; y[29] =   c; z[29] =  -a; w[29] = v;
		x[30] =  -b; y[30] =  -c; z[30] =   a; w[30] = v;
		x[31] =  -b; y[31] =  -c; z[31] =  -a; w[31] = v;
		x[32] =   c; y[32] =   a; z[32] =   b; w[32] = v;
		x[33] =   c; y[33] =   a; z[33] =  -b; w[33] = v;
		x[34] =   c; y[34] =  -a; z[34] =   b; w[34] = v;
		x[35] =   c; y[35] =  -a; z[35] =  -b; w[35] = v;
		x[36] =  -c; y[36] =   a; z[36] =   b; w[36] = v;
		x[37] =  -c; y[37] =   a; z[37] =  -b; w[37] = v;
		x[38] =  -c; y[38] =  -a; z[38] =   b; w[38] = v;
		x[39] =  -c; y[39] =  -a; z[39] =  -b; w[39] = v;
		x[40] =   c; y[40] =   b; z[40] =   a; w[40] = v;
		x[41] =   c; y[41] =   b; z[41] =  -a; w[41] = v;
		x[42] =   c; y[42] =  -b; z[42] =   a; w[42] = v;
		x[43] =   c; y[43] =  -b; z[43] =  -a; w[43] = v;
		x[44] =  -c; y[44] =   b; z[44] =   a; w[44] = v;
		x[45] =  -c; y[45] =   b; z[45] =  -a; w[45] = v;
		x[46] =  -c; y[46] =  -b; z[46] =   a; w[46] = v;
		x[47] =  -c; y[47] =  -b; z[47] =  -a; w[47] = v;
		n = 48;
	  }
	  else
	  {
		fprintf ( stderr, "\n" );
		fprintf ( stderr, "GEN_OH - Fatal error!\n" );
		fprintf ( stderr, "  Illegal value of code.\n" );
		exit ( 1 );
	  }
	  return n;
}
/******************************************************************************/

static void ld0026 ( double *x, double *y, double *z, double *w )

	/******************************************************************************/
	/*
	Purpose:

		LD0026 computes the 26 point Lebedev angular grid.

	Modified:

		12 September 2010

	Author:

		Dmitri Laikov

	Reference:

		Vyacheslav Lebedev, Dmitri Laikov,
		A quadrature formula for the sphere of the 131st
		algebraic order of accuracy,
		Russian Academy of Sciences Doklady Mathematics,
		Volume 59, Number 3, 1999, pages 477-481.

	Parameters:

		Output, double X[N], Y[N], Z[N], W[N], the coordinates
		and weights of the points.
	*/
{
  double a = 0.0;
  double b = 0.0;
  int n;
  double v;

  n = 0;
  v = 0.4761904761904762e-1;
  n = n + gen_oh ( 1, a, b, v, x + n, y + n, z + n, w + n );
  v = 0.3809523809523810e-1;
  n = n + gen_oh ( 2, a, b, v, x + n, y + n, z + n, w + n );
  v = 0.3214285714285714e-1;
  n = n + gen_oh ( 3, a, b, v, x + n, y + n, z + n, w + n );
  n = n - 1;

  return;
}

void Cal_rho(){
    int idx_cheb,idx_leb,mu,nu;
    double rho_0,rho_1;
    for(idx_cheb = 0; idx_cheb < N_cheb; idx_cheb++){        
        for(idx_leb = 0; idx_leb < N_leb;idx_leb++){ 
            rho_0 = rho_1 = 0.0;           
            for(mu = 0; mu < Aorbnum;mu++){
                for(nu = 0; nu < Aorbnum;nu++){
                    rho_0 += P[0][mu][nu]*phi[mu][idx_cheb][idx_leb]*phi[nu][idx_cheb][idx_leb];
                    rho_1 += P[1][mu][nu]*phi[mu][idx_cheb][idx_leb]*phi[nu][idx_cheb][idx_leb];
                }
            }      
            rho[0][idx_cheb][idx_leb] = rho_0;
            rho[1][idx_cheb][idx_leb] = rho_1;  
        }
    }
}

static void Cal_phi(){
    int k,mu,idx_cheb,idx_leb;
    int L,mx,my,mz;
    double value,r,x,y,z;
    for(k = 0; k < Aorbnum;k++){        
        L = basis[k].L;
        mx = basis[k].mx;
        my = basis[k].my;
        mz = basis[k].mz;
        for(idx_cheb=0;idx_cheb<N_cheb;idx_cheb++){
            for(idx_leb = 0;idx_leb<N_leb;idx_leb++){
                value = 0;
                r = r_cheb[idx_cheb];
                x = x_leb[idx_leb];
                y = y_leb[idx_leb];
                z = z_leb[idx_leb];
                for(mu = 0;mu<L;mu++){                
                    value += basis[k].coeff[mu]*exp(-1.0*basis[k].alpha[mu]*r*r);
                }
                value = value*pow(r*x,mx)*pow(r*y,my)*pow(r*z,mz);
                phi[k][idx_cheb][idx_leb] = value;            
            }
        }
    }
}
