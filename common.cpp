#define PI              3.1415926535897932384626
#include<stdio.h>
#include"common.h"
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
int Z;

/****************************************************
int Enum: 
    number of electrons
****************************************************/
int Enum;

/****************************************************
int N_up: 
    number of electrons of spin up
****************************************************/
int N_up;

/****************************************************
int N_down: 
    number of electrons of spin down
****************************************************/
int N_down;

/****************************************************
int Orbnum: 
    number of physical orbitals, for example, 
    1s,2s,2px,2py,2pz
****************************************************/
int Orbnum;

/****************************************************
int *STOnum: 
    number of STO for each physical orbital
    size:
        STOnum[Orbnum]
****************************************************/
int *STOnum;

/****************************************************
int Aorbnum:
    number of total orbitals,equals to sum(STOnum)
****************************************************/
int Aorbnum;

/****************************************************
Basis *basis:
    #Aorbnum Basis structs.
****************************************************/
Basis *basis;

/****************************************************
int *Ct_L:
    contraction length for each STO
    size:
        Ct_L[AOrbnum]
****************************************************/
int *Ct_L;

/****************************************************
double **Ct_D:
    contraction coefficient for each primitive
    gaussian orbital in STO
    size:
        Ct_D[Aorbnum][*Ct_L]
****************************************************/
double **Ct_D;

/****************************************************
double **Ct_A:
    orbital exponents for each primitive
    gaussian orbital in STO
    size:
        Ct_A[Aorbnum][*Ct_L]
****************************************************/
double **Ct_A;

/****************************************************
int **Ct_M:
    orbital angular momentum for each primitive
    gaussian orbital in STO
    size:
        Ct_M[Orbnum][3]
****************************************************/
int **Ct_M;

/****************************************************
double ***rho:
    electron density values at Beck grid
    size:
        rho[2][N_cheb][N_leb]
****************************************************/
double ***rho;

/****************************************************
double ***phi:
    basis function values at Beck grid
    size:
        phi[Aorbnum][N_cheb][N_leb]
****************************************************/
double ***phi;

/****************************************************
double ***vxc:
    xc potential values at Beck grid
    size:
        vxc[2][N_cheb][N_leb]
****************************************************/
double ***vxc;

/****************************************************
int Max_SCF: 
    maximum interation number
****************************************************/
int Max_SCF;

/****************************************************
double E_Conv:
    convergence criteria
****************************************************/
double E_Conv;

/****************************************************
int Xc_flag:
    type of exchange correlation. 
    0: Hartree Fock
    1: LSDA
****************************************************/
int Xc_flag;

/*****************************************************************************
                             SCF related variables
*****************************************************************************/

/****************************************************
double *Etot:
    total energy at each SCF step
    size:
        Etot[Max_SCF]
****************************************************/
double *Etot;

/****************************************************
double Fermi_Energy:
    Fermi_Energy at each SCF step
****************************************************/
double Fermi_Energy;

/****************************************************
double **E_tot:
    total energy matrix at each SCF step, see Exercise 3.40
    size:
        E_tot[Aorbnum][Aorbnum]
****************************************************/
double **E_tot;


/****************************************************
double ****Int_two:
    two body integral
    size:
        Int_two[Aorbnum][Aorbnum][Aorbnum][Aorbnum]
****************************************************/
double ****Int_two;

/****************************************************
double ***H:
    Hamiltonian matrix. the first index points to spin.
    size: 
        H[2][Aorbnum][Aorbnum]
****************************************************/
double ***H;

/****************************************************
double **H_kin:
    kinetic energy matrix.
    size: 
        H_kin[Aorbnum][Aorbnum]
****************************************************/
double **H_kin;

/****************************************************
double **H_ec:
    eletron-core interaction matrix.
    size: 
        H_ec[Aorbnum][Aorbnum]
****************************************************/
double **H_ec;

/****************************************************
double **H_core:
    one body matrix.
    size: 
        H_core[Aorbnum][Aorbnum]
****************************************************/
double **H_core;

/****************************************************
double **H_ee:
    eletron-electron interaction matrix.
    size: 
        H_ee[Aorbnum][Aorbnum]
****************************************************/
double **H_ee;

/****************************************************
double ***H_xc:
    Hamiltonian matrix. the first index points to spin.
    size: 
        H_xc[2][Aorbnum][Aorbnum]
****************************************************/
double ***H_xc;

/****************************************************
double **S:
    overlap matrix.
    size: 
        S[Aorbnum][Aorbnum]
****************************************************/
double **S;

/****************************************************
double *S_data:
    overlap matrix, continuous memory.
    size: 
        S[Aorbnum*Aorbnum]
****************************************************/
double *S_data;

/****************************************************
double ***C:
    expansion coefficients of orbitals 
    over STOs. the first index points to spin.
    size: 
        C[2][Aorbnum][Aorbnum]
****************************************************/
double ***C;

/****************************************************
double ***P_old:
    expansion coefficients of orbitals 
    over STOs. the first index points to spin.
    size: 
        P_old[2][Aorbnum][Aorbnum]
****************************************************/
double ***P_old;

/****************************************************
double **C_data:
    expansion coefficients of orbitals, continuous memory.
    size: 
        C[2][Aorbnum*Aorbnum]
****************************************************/
double **C_data;

/****************************************************
double ***P:
    spin polarized density matrix.
    size:
        P[2][Aorbnum][Aorbnum]
****************************************************/
double ***P;

/****************************************************
double **P_tot:
    total density matrix.
    size:
        P[Aorbnum][Aorbnum]
****************************************************/
double **P_tot;

/****************************************************
double **Eorb:
    orbital energy
    size:
        Eorb[2][Aorbnum]
****************************************************/
double **Eorb;

/****************************************************
int SCF_Step:
    global variables.
****************************************************/
int SCF_Step;

/****************************************************
double Mixing:
    global variables.
****************************************************/
double Mixing;


/*****************************************************************************
                             useful sub routines
*****************************************************************************/

void Addition(int N, double **A, double **B){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            A[i][j]+=B[i][j];
        }
    }
}

void Set_Zero(int N, double **A){
    for(int i = 0;i<N;i++){
        for (int j = 0;j<N;j++){
            A[i][j] = 0.0;
        }
    }
}

double Trace(int N, double **A, double **B){
    double t = 0.0;
    for(int i = 0;i<N;i++){
        for(int j = 0;j<N;j++){
            t+=A[i][j]*B[j][i];
        }
    }
    return t;
}

void Print_Matrix(int N1, int N2, double **A){
    for (int i = 0;i<N1;i++){
        for(int j = 0;j<N2;j++){
            printf("%.10f  ",A[i][j]);
        }
        printf("\n");
    }
}