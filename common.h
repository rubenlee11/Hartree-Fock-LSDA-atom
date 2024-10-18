#define PI              3.1415926535897932384626

#ifndef COMMON_H
#define COMMON_H

/* typedef struct {
    double alpha;
    int x,y,z;
} Prim_Gaussian; */

typedef struct {
    int L;    
    double *coeff;
    double *alpha;
    int mx,my,mz;
} Basis;

#endif


extern Basis *basis;

/* 
Most names of variables follow those of in Szabo 1996, Dover.
 */

/*****************************************************************************
                             system related variables
*****************************************************************************/

/****************************************************
int Z: 
    atomic number
****************************************************/
extern int Z;

/****************************************************
int Enum: 
    number of electrons
****************************************************/
extern int Enum;

/****************************************************
int N_up: 
    number of electrons of spin up
****************************************************/
extern int N_up;

/****************************************************
int N_down: 
    number of electrons of spin down
****************************************************/
extern int N_down;

/****************************************************
int Orbnum: 
    number of physical orbitals, for example, 
    1s,2s,2px,2py,2pz
****************************************************/
extern int Orbnum;

/****************************************************
int *STOnum: 
    number of STO for each physical orbital
    size:
        STOnum[Orbnum]
****************************************************/
extern int *STOnum;

/****************************************************
int Aorbnum:
    number of total orbitals,equals to sum(STOnum)
****************************************************/
extern int Aorbnum;

/****************************************************
int *Ct_L:
    contraction length for each STO
    size:
        Ct_L[AOrbnum]
****************************************************/
extern int *Ct_L;

/****************************************************
double **Ct_D:
    contraction coefficient for each primitive
    gaussian orbital in STO
    size:
        Ct_D[Aorbnum][*Ct_L]
****************************************************/
extern double **Ct_D;

/****************************************************
double **Ct_A:
    orbital exponents for each primitive
    gaussian orbital in STO
    size:
        Ct_A[Aorbnum][*Ct_L]
****************************************************/
extern double **Ct_A;

/****************************************************
int **Ct_M:
    orbital angular momentum for each primitive
    gaussian orbital in STO
    size:
        Ct_M[Orbnum][3]
****************************************************/
extern int **Ct_M;

/****************************************************
int N_cheb:
    number of chebyshev radial grid points
****************************************************/
extern int N_cheb;

/****************************************************
int N_leb:
    number of lebedev sphere grid points
****************************************************/
extern int N_leb;

/****************************************************
double ***rho:
    electron density values at Beck grid
    size:
        rho[2][N_cheb][N_leb]
****************************************************/
extern double ***rho;

/****************************************************
double ***phi:
    basis function values at Beck grid
    size:
        phi[Aorbnum][N_cheb][N_leb]
****************************************************/
extern double ***phi;

extern double rm;

/****************************************************
double ***vxc:
    xc potential values at Beck grid
    size:
        vxc[2][N_cheb][N_leb]
****************************************************/
extern double ***vxc;

extern double *r_cheb,*x_cheb,*w_cheb;
extern double *x_leb,*y_leb,*z_leb,*w_leb;

/****************************************************
int Max_SCF: 
    maximum interation number
****************************************************/
extern int Max_SCF;

/****************************************************
double E_Conv:
    convergence criteria
****************************************************/
extern double E_Conv;

/****************************************************
int Xc_flag:
    type of exchange correlation. 
    0: Hartree Fock
    1: LSDA
****************************************************/
extern int Xc_flag;

/*****************************************************************************
                             SCF related variables
*****************************************************************************/

/****************************************************
double *Etot:
    total energy at each SCF step
    size:
        Etot[Max_SCF]
****************************************************/
extern double *Etot;

/****************************************************
double Fermi_Energy:
    Fermi_Energy at each SCF step
****************************************************/
extern double Fermi_Energy;

/****************************************************
double **E_tot:
    total energy matrix at each SCF step, see Exercise 3.40
    size:
        E_tot[Aorbnum][Aorbnum]
****************************************************/
extern double **E_tot;


/****************************************************
double ****Int_two:
    two body integral
    size:
        Int_two[Aorbnum][Aorbnum][Aorbnum][Aorbnum]
****************************************************/
extern double ****Int_two;

/****************************************************
double ***H:
    Hamiltonian matrix. the first index points to spin.
    size: 
        H[2][Aorbnum][Aorbnum]
****************************************************/
extern double ***H;

/****************************************************
double **H_kin:
    kinetic energy matrix.
    size: 
        H_kin[Aorbnum][Aorbnum]
****************************************************/
extern double **H_kin;

/****************************************************
double **H_ec:
    eletron-core interaction matrix.
    size: 
        H_ec[Aorbnum][Aorbnum]
****************************************************/
extern double **H_ec;

/****************************************************
double **H_core:
    one body matrix.
    size: 
        H_core[Aorbnum][Aorbnum]
****************************************************/
extern double **H_core;

/****************************************************
double **H_ee:
    eletron-electron interaction matrix.
    size: 
        H_ee[Aorbnum][Aorbnum]
****************************************************/
extern double **H_ee;

/****************************************************
double ***H_xc:
    Hamiltonian matrix. the first index points to spin.
    size: 
        H_xc[2][Aorbnum][Aorbnum]
****************************************************/
extern double ***H_xc;

/****************************************************
double **S:
    overlap matrix.
    size: 
        S[Aorbnum][Aorbnum]
****************************************************/
extern double **S;

/****************************************************
double *S_data:
    overlap matrix, continuous memory.
    size: 
        S[Aorbnum*Aorbnum]
****************************************************/
extern double *S_data;

/****************************************************
double ***C:
    expansion coefficients of orbitals 
    over STOs. the first index points to spin.
    size: 
        C[2][Aorbnum][Aorbnum]
****************************************************/
extern double ***C;

/****************************************************
double ***P_old:
    expansion coefficients of orbitals 
    over STOs. the first index points to spin.
    size: 
        P_old[2][Aorbnum][Aorbnum]
****************************************************/
extern double ***P_old;

/****************************************************
double **C_data:
    expansion coefficients of orbitals, continuous memory.
    size: 
        C[2][Aorbnum*Aorbnum]
****************************************************/
extern double **C_data;

/****************************************************
double ***P:
    spin polarized density matrix.
    size:
        P[2][Aorbnum][Aorbnum]
****************************************************/
extern double ***P;

/****************************************************
double **P_tot:
    total density matrix.
    size:
        P[Aorbnum][Aorbnum]
****************************************************/
extern double **P_tot;

/****************************************************
double **Eorb:
    orbital energy
    size:
        Eorb[2][Aorbnum]
****************************************************/
extern double **Eorb;

/****************************************************
int SCF_Step:
    global variables.
****************************************************/
extern int SCF_Step;

/****************************************************
double Mixing:
    global variables.
****************************************************/
extern double Mixing;

/*****************************************************************************
                             Main sub routines of DFT(HF) calculation.
*****************************************************************************/

extern void Init_Basis();
extern void Init_Arrays();
extern void Init_Grid();
extern void Integral_C();
extern void Integral_H();
extern void Mix();
extern void SCF();
extern void Total_Energy();
extern void Wave_function();

extern void Cal_rho();
extern void Init_Quadrature();
extern double Quadrature(double **F);

extern double Cal_xc_energy();
extern void Cal_xc_matrix(double *result,int mu,int nu);

/*****************************************************************************
                             useful sub routines
*****************************************************************************/

/****************************************************
void Addition:
    Add N*N square matrix B to A.
****************************************************/
extern void Addition(int N, double **A, double **B);

/****************************************************
void Trace:
    Compute the trace of N*N matrix A times B.
****************************************************/
extern double Trace(int N, double **A, double **B);

/****************************************************
void Set_Zero:
    Set matrix elements of A to 0.
****************************************************/
extern void Set_Zero(int N, double **A);

extern void Print_Matrix(int N1, int N2, double **A);
