#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>
#include "common.h"

/****************************************************
static int Find_Value:
    for carbon atom, given ERI index, find the 
    sorted ERI index, while ERI.
    the sorting rule is: s->px->py->pz 
****************************************************/  
static int Find_Value(int *idx);

/****************************************************
static void Sort_Integer:
    sort: A<B
****************************************************/ 
static void Sort_Integer(int *A,int *B);

/****************************************************
static int Find_M:
    for carbon atom, find the angular momentum of
    the input primitive index.
****************************************************/  
static int Find_M(int A);

/****************************************************
static int Theta:
    theta function
****************************************************/  
static int Theta(int x);

/****************************************************
static int Trans_idx:
    flag: the target angular momentum 
    x: input primitive index

    for carbon atom, find the target primitive index
    while fixing ERI.
****************************************************/  
static int Trans_idx(int flag,int x);

/****************************************************
static int start_C:
    for carbon atom, given the basis index, 
    return the primitive index starting value. 
****************************************************/  
static int start_C(int k);

void Init_Basis(){
    /****************************************************
    allocation of calculating parameters
    ****************************************************/ 
    int i,mu;
    double zeta;       
    if(Z == 6){        
        Enum = 6;
        N_up = 4;
        N_down = Enum - N_up;
        Orbnum = 5;    
        Aorbnum = 9;  
        /****************************************************
        6-31G basis for C
        ****************************************************/      
        STOnum = (int*)malloc(sizeof(int)*Orbnum); 
        STOnum[0] = 1;STOnum[1] = 2;STOnum[2] = 2;STOnum[3] = 2;STOnum[4] = 2;   

        Ct_L = (int*)malloc(sizeof(int)*Aorbnum);

        Ct_L[0] = 6;
        Ct_L[1] = 3;
        Ct_L[2] = 1;
        Ct_L[3] = 3;
        Ct_L[4] = 1;
        Ct_L[5] = 3;
        Ct_L[6] = 1;
        Ct_L[7] = 3;
        Ct_L[8] = 1;

        Ct_D = (double**)malloc(sizeof(double*)*Aorbnum);
        for(int i=0; i<Aorbnum;i++){
            Ct_D[i] = (double*)malloc(sizeof(double)*Ct_L[i]);
        }

        Ct_A = (double**)malloc(sizeof(double*)*Aorbnum);
        for(int i=0; i<Aorbnum;i++){
            Ct_A[i] = (double*)malloc(sizeof(double)*Ct_L[i]);
        }

        Ct_A[0][0] = .3047524880e4;Ct_D[0][0] = .1834737130e-2;
        Ct_A[0][1] = .4573695180e3;Ct_D[0][1] = .1403732280e-1;
        Ct_A[0][2] = .1039486850e3;Ct_D[0][2] = .06884262220;
        Ct_A[0][3] = .2921015530e2;Ct_D[0][3] = .2321844430;
        Ct_A[0][4] = .9286662960e1;Ct_D[0][4] = .4679413480;
        Ct_A[0][5] = .3163926960e1;Ct_D[0][5] = .3623119850;

        Ct_A[1][0] = Ct_A[3][0] = Ct_A[5][0] = Ct_A[7][0] = .7868272350e+01;Ct_D[1][0] = -.1193324200;Ct_D[3][0] = Ct_D[5][0] = Ct_D[7][0] = .6899906660e-1;
        Ct_A[1][1] = Ct_A[3][1] = Ct_A[5][1] = Ct_A[7][1] = .1881288540e+01;Ct_D[1][1] = -.1608541520;Ct_D[3][1] = Ct_D[5][1] = Ct_D[7][1] = .3164239610;
        Ct_A[1][2] = Ct_A[3][2] = Ct_A[5][2] = Ct_A[7][2] = .5442492580;Ct_D[1][2] = .114345644e1;Ct_D[3][2] = Ct_D[5][2] = Ct_D[7][2] = .7443082910;

        Ct_A[2][0] = Ct_A[4][0] = Ct_A[6][0] = Ct_A[8][0] = .1687144782;Ct_D[2][0] = 1.000000000;Ct_D[4][0] = Ct_D[6][0] = Ct_D[8][0] = 1.000000000;  

        //normalization s orbital primitive coefficients
        for(i = 0; i < 3;i++){
            for(mu = 0;mu<Ct_L[i];mu++){
                zeta = 2*Ct_A[i][mu];
                Ct_D[i][mu] = Ct_D[i][mu]/sqrt(pow(sqrt(PI/(zeta)),3));
            }        
        }         
        //normalization p orbital primitive coefficients
        for(i = 3; i < Aorbnum;i++){
            for(mu = 0;mu<Ct_L[i];mu++){
                zeta = 2*Ct_A[i][mu];
                Ct_D[i][mu] = Ct_D[i][mu]/sqrt(pow(sqrt(PI/(zeta)),3)/(2.0*zeta));
            }        
        } 
        // Basis structs
        basis = (Basis*)malloc(sizeof(Basis)*Aorbnum);
        for(i = 0; i < Aorbnum;i++){
            basis[i].L = Ct_L[i];
            basis[i].coeff = (double*)malloc(sizeof(double)*basis[i].L);
            basis[i].alpha = (double*)malloc(sizeof(double)*basis[i].L);

            for(mu = 0; mu < basis[i].L;mu++){
                basis[i].coeff[mu] = Ct_D[i][mu];
                basis[i].alpha[mu] = Ct_A[i][mu];
            }

            if(i<3){
                basis[i].mx = basis[i].my = basis[i].mz = 0;
            }
            else if(i>=3&&i<5){
                basis[i].mx = 1;
                basis[i].my = basis[i].mz = 0;
            }
            else if(i>=5&&i<7){
                basis[i].my = 1;
                basis[i].mx = basis[i].mz = 0;
            }
            else{
                basis[i].mz = 1;
                basis[i].mx = basis[i].my = 0;
            }
        }

    }
    else if(Z == 1){
        // 3-21G basis
        Enum = 1;
        N_up = 1;
        N_down = Enum - N_up;
        Orbnum = 2;    
        Aorbnum = 2;

        STOnum = (int*)malloc(sizeof(int)*Orbnum); 
        STOnum[0] = 2;STOnum[1] = 1;      

        Ct_L = (int*)malloc(sizeof(int)*Aorbnum);

        Ct_L[0] = 2;
        Ct_L[1] = 1;  

        Ct_D = (double**)malloc(sizeof(double*)*Aorbnum);
        for(int i=0; i<Aorbnum;i++){
            Ct_D[i] = (double*)malloc(sizeof(double)*Ct_L[i]);
        }

        Ct_A = (double**)malloc(sizeof(double*)*Aorbnum);
        for(int i=0; i<Aorbnum;i++){
            Ct_A[i] = (double*)malloc(sizeof(double)*Ct_L[i]);
        }        

        Ct_A[0][0] = 0.5447178000e1; Ct_D[0][0] = 0.1562849787;
        Ct_A[0][1] = 0.8245472400; Ct_D[0][1] = 0.9046908767;

        Ct_A[1][0] = 0.1831915800; Ct_D[1][0] =1.0000000;   
		
		//normalization
		for(i = 0; i <Aorbnum;i++){
			for(mu = 0;mu<Ct_L[i];mu++){
				zeta = 2*Ct_A[i][mu];
				Ct_D[i][mu] = Ct_D[i][mu]/sqrt(pow(sqrt(PI/(zeta)),3));
			}        
		}				
        // Basis structs
        basis = (Basis*)malloc(sizeof(Basis)*Aorbnum);
        for(i = 0; i < Aorbnum;i++){
            basis[i].L = Ct_L[i];
            basis[i].coeff = (double*)malloc(sizeof(double)*basis[i].L);
            basis[i].alpha = (double*)malloc(sizeof(double)*basis[i].L);

            for(mu = 0; mu < basis[i].L;mu++){
                basis[i].coeff[mu] = Ct_D[i][mu];
                basis[i].alpha[mu] = Ct_A[i][mu];
            }

            basis[i].mx = basis[i].my = basis[i].mz = 0;
        }		
    }
}

void Init_Arrays(){
    /****************************************************
    allocation of memory for arrays
    ****************************************************/  
    H = (double***)malloc(sizeof(double**)*(2));
    for (int spin=0;spin<=1;spin++){
        H[spin] = (double**)malloc(sizeof(double*)*Aorbnum);
        for (int j=0;j<Aorbnum;j++){
            H[spin][j] = (double*)malloc(sizeof(double)*Aorbnum);
        }
    }   

    S_data = (double*)malloc(sizeof(double)*Aorbnum*Aorbnum);
    S = (double**)malloc(sizeof(double*)*Aorbnum);
    for (int j=0;j<Aorbnum;j++){
        S[j] = (double*)malloc(sizeof(double)*Aorbnum);
    }
    Set_Zero(Aorbnum,S);

    H_kin = (double**)malloc(sizeof(double*)*Aorbnum);
    for (int j=0;j<Aorbnum;j++){
        H_kin[j] = (double*)malloc(sizeof(double)*Aorbnum);
    }
    Set_Zero(Aorbnum,H_kin);

    H_ec = (double**)malloc(sizeof(double*)*Aorbnum);
    for (int j=0;j<Aorbnum;j++){
        H_ec[j] = (double*)malloc(sizeof(double)*Aorbnum);
    }    
    Set_Zero(Aorbnum,H_ec);

    H_core = (double**)malloc(sizeof(double*)*Aorbnum);
    for (int j=0;j<Aorbnum;j++){
        H_core[j] = (double*)malloc(sizeof(double)*Aorbnum);
    }   
    Set_Zero(Aorbnum,H_core);  

    Int_two = (double****)malloc(sizeof(double***)*Aorbnum);
    for (int i = 0; i < Aorbnum; i++){
        Int_two[i] = (double***)malloc(sizeof(double**)*Aorbnum);
        for (int j = 0; j < Aorbnum; j++){
            Int_two[i][j] = (double**)malloc(sizeof(double*)*Aorbnum);
            for (int k = 0; k < Aorbnum; k++){
                Int_two[i][j][k] = (double*)malloc(sizeof(double)*Aorbnum);
            }
        }
    }
    /*     for(int i = 0; i < Aorbnum;i++){
        for(int j = 0; j < Aorbnum; j++){
            Set_Zero(Aorbnum,Int_two[i][j]);
        }
    } */

    H_ee = (double**)malloc(sizeof(double*)*Aorbnum);
    for (int j=0;j<Aorbnum;j++){
        H_ee[j] = (double*)malloc(sizeof(double)*Aorbnum);
    }    

    H_xc = (double***)malloc(sizeof(double**)*(2));
    for (int spin=0;spin<=1;spin++){
        H_xc[spin] = (double**)malloc(sizeof(double*)*Aorbnum);
        for (int j=0;j<Aorbnum;j++){
            H_xc[spin][j] = (double*)malloc(sizeof(double)*Aorbnum);
        }
    }

    C_data = (double**)malloc(sizeof(double*)*2);
    C_data[0] = (double*)malloc(sizeof(double)*Aorbnum*Aorbnum);
    C_data[1] = (double*)malloc(sizeof(double)*Aorbnum*Aorbnum);
    C = (double***)malloc(sizeof(double**)*2);
    for (int spin=0;spin<=1;spin++){
        C[spin] = (double**)malloc(sizeof(double*)*Aorbnum);
        for (int j=0;j<Aorbnum;j++){
            C[spin][j] = C_data[spin]+j*Aorbnum;
        }
    }
    for(int i=0;i<Aorbnum;i++){
        for(int j=0;j<Aorbnum;j++){
            C_data[0][Aorbnum*i+j] = 0.0;
            C_data[1][Aorbnum*i+j] = 0.0;
        }
    } 

    P = (double***)malloc(sizeof(double**)*(2));
    for (int spin=0;spin<2;spin++){
        P[spin] = (double**)malloc(sizeof(double*)*Aorbnum);
        for (int j=0;j<Aorbnum;j++){
            P[spin][j] = (double*)malloc(sizeof(double)*Aorbnum);
        }
    } 

    P_tot = (double**)malloc(sizeof(double*)*Aorbnum);
    for (int j=0;j<Aorbnum;j++){
        P_tot[j] = (double*)malloc(sizeof(double)*Aorbnum);
    }       

    Eorb = (double**)malloc(sizeof(double*)*2);
    for (int j=0;j<Aorbnum;j++){
        Eorb[j] = (double*)malloc(sizeof(double)*Aorbnum);
    }    

    Etot = (double*)malloc(sizeof(double)*(Max_SCF+1));
}

void Init_Grid(){}

void Integral_C(){
    /****************************************************
    Only works for 6-31G basis with 1s,2s,2p orbital!
    ****************************************************/
    // the 6-31G basis of carbon consists of 6 1s exponents, 4 2s-2p shared exponents.
    int N_1s = 6;int N_2sp = 4;int N_prim;    
    N_prim = N_1s + N_2sp;
    double zeta,xi,rho,ita;
    int i,j;
    int k1,k2,k3,k4;
    int mu,nu,lambda,sigma;
    int idx_1 = 0;int idx_2 = 0;int idx_3 = 0;int idx_4 = 0;    
    // exponents of s shell, sp shell primitive Gaussians. 
    double Alpha[N_prim];      
    for (int i = 0; i < N_1s;i++){
        Alpha[i] = Ct_A[0][i];
    }
    for (int i = 0; i < 3;i++){
        Alpha[i+6] = Ct_A[1][i];
    }
    Alpha[9] = Ct_A[2][0];    
    /****************************************************
    overlap primitive matrix, s-s, px-px
    kinetic primitive matrix, s-s, px-px
    electron-core primitive matrix, s-s, px-px
    ****************************************************/    
    double Olp_Prim_ss[N_prim][N_prim] = {0.0};
    double Olp_Prim_pp[N_2sp][N_2sp] = {0.0};
    double T_Prim_ss[N_prim][N_prim] = {0.0};
    double T_Prim_pp[N_2sp][N_2sp] = {0.0};
    double Ec_Prim_ss[N_prim][N_prim] = {0.0};
    double Ec_Prim_pp[N_2sp][N_2sp] = {0.0};   

    /****************************************************
    s-s primitive matrix
    ****************************************************/  
    for(i = 0; i <N_prim;i++){
        for(j = 0; j <N_prim;j++){
            zeta = Alpha[i]+Alpha[j];
            xi = Alpha[i]*Alpha[j]/zeta;
            Olp_Prim_ss[i][j] = pow(sqrt(PI/zeta),3);
            T_Prim_ss[i][j] = 3*xi*Olp_Prim_ss[i][j];
            Ec_Prim_ss[i][j] = -Z*2.*PI/zeta;
        }
    }     

    /****************************************************
    px-px primitive matrix
    ****************************************************/  
    for(i = 0; i < N_2sp;i++){
        for(j = 0; j < N_2sp;j++){
            zeta = Alpha[i+6]+Alpha[j+6];
            xi = Alpha[i+6]*Alpha[j+6]/zeta;            
            Olp_Prim_pp[i][j] = 0.5*Olp_Prim_ss[i+6][j+6]/zeta;
            T_Prim_pp[i][j] = 5*xi*Olp_Prim_pp[i][j];
            Ec_Prim_pp[i][j] = -Z*4./3.*sqrt(zeta/PI)*Olp_Prim_pp[i][j];
        }
    }      
    /****************************************************
    primitive matrix to basis matrix
    ****************************************************/  
    // s-s 
     for (k1 = 0; k1 < 3;k1++){
        for(k2 = k1; k2 < 3;k2++){            
            for(mu = 0;mu<Ct_L[k1];mu++){
                for(nu = 0;nu<Ct_L[k2];nu++){
                    S[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*Olp_Prim_ss[start_C(k1)+mu][start_C(k2)+nu];
                    H_kin[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*T_Prim_ss[start_C(k1)+mu][start_C(k2)+nu];
                    H_ec[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ec_Prim_ss[start_C(k1)+mu][start_C(k2)+nu];
                }
            }   
        }        
    } 
    // px-px 
     for (k1 = 3; k1 < 5;k1++){
        for(k2 = k1; k2 < 5;k2++){            
            for(mu = 0;mu<Ct_L[k1];mu++){
                for(nu = 0;nu<Ct_L[k2];nu++){
                    S[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*Olp_Prim_pp[start_C(k1)+mu][start_C(k2)+nu];
                    H_kin[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*T_Prim_pp[start_C(k1)+mu][start_C(k2)+nu];
                    H_ec[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ec_Prim_pp[start_C(k1)+mu][start_C(k2)+nu];
                }
            }   
        }        
    }     
    // py-py,pz-pz
    S[5][5] = S[7][7] = S[3][3];S[6][6] = S[8][8] = S[4][4];S[7][8] = S[5][6] = S[3][4];
    H_kin[5][5] = H_kin[7][7] = H_kin[3][3];H_kin[6][6] = H_kin[8][8] = H_kin[4][4];H_kin[7][8] = H_kin[5][6] = H_kin[3][4];
    H_ec[5][5] = H_ec[7][7] = H_ec[3][3];H_ec[6][6] = H_ec[8][8] = H_ec[4][4];H_ec[7][8] = H_ec[5][6] = H_ec[3][4];
    // symmetrize
    for (k1 = 0; k1 < Aorbnum;k1++){
        for (k2 = 0; k2 < k1; k2++){
            S[k1][k2] = S[k2][k1];
            H_kin[k1][k2] = H_kin[k2][k1];
            H_ec[k1][k2] = H_ec[k2][k1];
        }
    }

    Addition(Aorbnum,H_core,H_kin);
    Addition(Aorbnum,H_core,H_ec);

    /****************************************************
    electron repulsion integrals
    (ss,ss),(sp,sp),(ss,pp),(pp,pp)    

    ERI_pppp_iijj means (pi pi, pj pj), and no matter i!=j takes whatever value, the integrals remain the same.
    ERT_pppp_ijij means (pi pj, pi pj)
    ****************************************************/    
    double ERI_ssss[N_prim][N_prim][N_prim][N_prim] = {0.0};
    double ERI_spsp[N_prim][N_2sp][N_prim][N_2sp] = {0.0};
    double ERI_sspp[N_prim][N_prim][N_2sp][N_2sp] = {0.0};
    double ERI_pppp_iijj[N_2sp][N_2sp][N_2sp][N_2sp] = {0.0};  
    double ERI_pppp_ijij[N_2sp][N_2sp][N_2sp][N_2sp] = {0.0};
    double ERI_pppp_iiii[N_2sp][N_2sp][N_2sp][N_2sp] = {0.0};

    for (int i = 0; i < N_prim;i++){
        for (int j = 0; j < N_prim;j++){
            for (int k = 0; k < N_prim;k++){
                for (int l = 0; l < N_prim;l++){
                    zeta = Alpha[i] + Alpha[j];
                    ita = Alpha[k] + Alpha[l];
                    ERI_ssss[i][j][k][l] = 2.0*pow(sqrt(PI),5)/(zeta*ita*sqrt(zeta+ita));                                                    
                }
            }
        }
    }  

    for (int i = 0; i < N_prim;i++){
        for (int j = 0; j < N_2sp;j++){
            for (int k = 0; k < N_prim;k++){
                for (int l = 0; l < N_2sp;l++){
                    zeta = Alpha[i] + Alpha[j+6];
                    ita = Alpha[k] + Alpha[l+6];
                    ERI_spsp[i][j][k][l] = ERI_ssss[i][j+6][k][l+6]/(6.0*(zeta+ita));       
                }
            }
        }
    }  

    for (int i = 0; i < N_prim;i++){
        for (int j = 0; j < N_prim;j++){
            for (int k = 0; k < N_2sp;k++){
                for (int l = 0; l < N_2sp;l++){
                    zeta = Alpha[i] + Alpha[j];
                    ita = Alpha[k+6] + Alpha[l+6];
                    rho = zeta*ita/(zeta+ita);
                    ERI_sspp[i][j][k][l] = ERI_ssss[i][j][k+6][l+6]*(1.0-rho/(3.0*ita))/(2.0*ita);                   
                }
            }
        }
    }   

    for (int i = 0; i < N_2sp;i++){
        for (int j = 0; j < N_2sp;j++){
            for (int k = 0; k < N_2sp;k++){
                for (int l = 0; l < N_2sp;l++){
                    zeta = Alpha[i+6] + Alpha[j+6];
                    ita = Alpha[k+6] + Alpha[l+6];
                    rho = zeta*ita/(zeta+ita);
                    ERI_pppp_iijj[i][j][k][l] = (2.0/3.0+rho*rho/(5.0*zeta*ita))/(4.0*zeta*ita)*ERI_ssss[i+6][j+6][k+6][l+6];                           
                    ERI_pppp_ijij[i][j][k][l] = ERI_ssss[i+6][j+6][k+6][l+6]/(20.0*pow(zeta+ita,2));    
                    ERI_pppp_iiii[i][j][k][l] = 2.0*ERI_pppp_ijij[i][j][k][l]+ERI_pppp_iijj[i][j][k][l];
                }
            }
        }
    }       

    /****************************************************
    primitive matrix to basis matrix
    ****************************************************/  
    // (s s,s s)
    for (k1 = 0; k1 < 3;k1++){
        for (k2 = 0; k2 < 3;k2++){
            for (k3 = 0; k3 < 3;k3++){
                for (k4 = 0; k4 < 3;k4++){   
                    Int_two[k1][k2][k3][k4] = 0.0;                 
                    for(mu = 0;mu<Ct_L[k1];mu++){
                        for(nu = 0;nu<Ct_L[k2];nu++){
                            for(lambda = 0;lambda<Ct_L[k3];lambda++){
                                for(sigma = 0;sigma<Ct_L[k4];sigma++){
                                    Int_two[k1][k2][k3][k4] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ct_D[k3][lambda]*Ct_D[k4][sigma]*ERI_ssss[start_C(k1)+mu][start_C(k2)+nu][start_C(k3)+lambda][start_C(k4)+sigma];
                                }
                            }
                        }
                    }                  
                }
            }
        }
    }

    // (s px,s px)
    for (k1 = 0; k1 < 3;k1++){               
        for (k2 = 3; k2 < 5;k2++){        
            for (k3 = 0; k3 < 3;k3++){
                for (k4 = 3; k4 < 5;k4++){                    
                    for(mu = 0;mu<Ct_L[k1];mu++){
                        for(nu = 0;nu<Ct_L[k2];nu++){
                            for(lambda = 0;lambda<Ct_L[k3];lambda++){
                                for(sigma = 0;sigma<Ct_L[k4];sigma++){
                                    Int_two[k1][k2][k3][k4] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ct_D[k3][lambda]*Ct_D[k4][sigma]*ERI_spsp[start_C(k1)+mu][start_C(k2)+nu][start_C(k3)+lambda][start_C(k4)+sigma];
                                }
                            }
                        }
                    }                   
                }           
            }
        }
    }    

    // (s s,px px)
    for (k1 = 0; k1 < 3;k1++){               
        for (k2 = 0; k2 < 3;k2++){        
            for (k3 = 3; k3 < 5;k3++){
                for (k4 = 3; k4 < 5;k4++){                    
                    for(mu = 0;mu<Ct_L[k1];mu++){
                        for(nu = 0;nu<Ct_L[k2];nu++){
                            for(lambda = 0;lambda<Ct_L[k3];lambda++){
                                for(sigma = 0;sigma<Ct_L[k4];sigma++){
                                    Int_two[k1][k2][k3][k4] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ct_D[k3][lambda]*Ct_D[k4][sigma]*ERI_sspp[start_C(k1)+mu][start_C(k2)+nu][start_C(k3)+lambda][start_C(k4)+sigma];
                                }
                            }
                        }
                    }                   
                }           
            }
        }
    }  

    // (px px,py py)
    for (k1 = 3; k1 < 5;k1++){               
        for (k2 = 3; k2 < 5;k2++){         
            for (k3 = 5; k3 < 7;k3++){
                for (k4 = 5; k4 < 7;k4++){                    
                    for(mu = 0;mu<Ct_L[k1];mu++){
                        for(nu = 0;nu<Ct_L[k2];nu++){
                            for(lambda = 0;lambda<Ct_L[k3];lambda++){
                                for(sigma = 0;sigma<Ct_L[k4];sigma++){
                                    Int_two[k1][k2][k3][k4] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ct_D[k3][lambda]*Ct_D[k4][sigma]*ERI_pppp_iijj[start_C(k1)+mu][start_C(k2)+nu][start_C(k3)+lambda][start_C(k4)+sigma];
                                }
                            }
                        }
                    }                  
                }           
            }
        }
    } 

    // (px py,px py) 
    for (k1 = 3; k1 < 5;k1++){               
        for (k2 = 5; k2 < 7;k2++){         
            for (k3 = 3; k3 < 5;k3++){
                for (k4 = 5; k4 < 7;k4++){                    
                    for(mu = 0;mu<Ct_L[k1];mu++){
                        for(nu = 0;nu<Ct_L[k2];nu++){
                            for(lambda = 0;lambda<Ct_L[k3];lambda++){
                                for(sigma = 0;sigma<Ct_L[k4];sigma++){
                                    Int_two[k1][k2][k3][k4] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ct_D[k3][lambda]*Ct_D[k4][sigma]*ERI_pppp_ijij[start_C(k1)+mu][start_C(k2)+nu][start_C(k3)+lambda][start_C(k4)+sigma];
                                }
                            }
                        }
                    }                  
                }           
            }
        }
    }     

    
    // (px px,px px)
    for (k1 = 3; k1 < 5;k1++){               
        for (k2 = 3; k2 < 5;k2++){       
            for (k3 = 3; k3 < 5;k3++){
                for (k4 = 3; k4 < 5;k4++){                    
                    for(mu = 0;mu<Ct_L[k1];mu++){
                        for(nu = 0;nu<Ct_L[k2];nu++){
                            for(lambda = 0;lambda<Ct_L[k3];lambda++){
                                for(sigma = 0;sigma<Ct_L[k4];sigma++){
                                    Int_two[k1][k2][k3][k4] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ct_D[k3][lambda]*Ct_D[k4][sigma]*ERI_pppp_iiii[start_C(k1)+mu][start_C(k2)+nu][start_C(k3)+lambda][start_C(k4)+sigma];                                    
                                }
                            }
                        }
                    }                    
                }              
            }
        }
    }           

    // fill the rest part of Int_two
    int *idx;
    idx = (int*)malloc(sizeof(int)*4);
    for (k1 = 0; k1 < Aorbnum; k1++){
        for (k2 = 0; k2 < Aorbnum; k2++){
            for (k3 = 0; k3 < Aorbnum;k3++){
                for (k4 = 0; k4 < Aorbnum; k4++){
                    idx[0] = k1;idx[1] = k2;idx[2] = k3;idx[3] = k4;
                    int flag = Find_Value(idx);
                    if (flag){
                        Int_two[k1][k2][k3][k4] = Int_two[idx[0]][idx[1]][idx[2]][idx[3]];
                    }
                    else{
                        Int_two[k1][k2][k3][k4] = 0.0;
                    }
                }
            }
        }
    }

}

void Integral_H(){
    /****************************************************
    Only works for 3-21G basis with 1s orbital!
    ****************************************************/
    int N_prim = 3;    
    double zeta,xi,rho,ita;
    int i,j;
    int k1,k2,k3,k4;
    int mu,nu,lambda,sigma;
    int idx_1 = 0;int idx_2 = 0;int idx_3 = 0;int idx_4 = 0;  

    double Alpha[N_prim];

    Alpha[0] = Ct_A[0][0];Alpha[1] = Ct_A[0][1];Alpha[2] = Ct_A[1][0];

    /****************************************************
    overlap primitive matrix, s-s
    kinetic primitive matrix, s-s
    electron-core primitive matrix, s-s
    ****************************************************/    
    double Olp_Prim_ss[N_prim][N_prim] = {0.0};
    double T_Prim_ss[N_prim][N_prim] = {0.0};
    double Ec_Prim_ss[N_prim][N_prim] = {0.0};

    for(i = 0; i <N_prim;i++){
        for(j = 0; j <N_prim;j++){
            zeta = Alpha[i]+Alpha[j];
            xi = Alpha[i]*Alpha[j]/zeta;
            Olp_Prim_ss[i][j] = pow(sqrt(PI/zeta),3);
            T_Prim_ss[i][j] = 3*xi*Olp_Prim_ss[i][j];
            Ec_Prim_ss[i][j] = -Z*2.0*sqrt(zeta/PI)*Olp_Prim_ss[i][j];
        }
    }      

    /****************************************************
    primitive matrix to basis matrix
    ****************************************************/  
    // s-s 
     for (k1 = 0; k1 < 2;k1++){
        idx_2 = idx_1;
        for(k2 = k1; k2 < 2;k2++){            
            for(mu = 0;mu<Ct_L[k1];mu++){
                for(nu = 0;nu<Ct_L[k2];nu++){
                    S[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*Olp_Prim_ss[idx_1+mu][idx_2+nu];
                    H_kin[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*T_Prim_ss[idx_1+mu][idx_2+nu];
                    H_ec[k1][k2] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ec_Prim_ss[idx_1+mu][idx_2+nu];
                }
            }   
            idx_2 += Ct_L[k2];
        }        
        idx_1 += Ct_L[k1];
    }    

    // symmetrize
    for (k1 = 0; k1 < Aorbnum;k1++){
        for (k2 = 0; k2 < k1; k2++){
            S[k1][k2] = S[k2][k1];
            H_kin[k1][k2] = H_kin[k2][k1];
            H_ec[k1][k2] = H_ec[k2][k1];
        }
    }

    Addition(Aorbnum,H_core,H_kin);
    Addition(Aorbnum,H_core,H_ec);        

    /****************************************************
    electron repulsion integrals
    (ss,ss)
    ****************************************************/    
    double ERI_ssss[N_prim][N_prim][N_prim][N_prim] = {0.0};   

    for (int i = 0; i < N_prim;i++){
        for (int j = 0; j < N_prim;j++){
            for (int k = 0; k < N_prim;k++){
                for (int l = 0; l < N_prim;l++){
                    zeta = Alpha[i] + Alpha[j];
                    ita = Alpha[k] + Alpha[l];
                    ERI_ssss[i][j][k][l] = 2.0*pow(sqrt(PI),5)/(zeta*ita*sqrt(zeta+ita));                                                    
                }
            }
        }
    }    

    /****************************************************
    primitive matrix to basis matrix
    ****************************************************/  
    // (s s,s s)
    idx_1 = 0;
    for (k1 = 0; k1 < 2;k1++){
        idx_2 = 0;
        for (k2 = 0; k2 < 2;k2++){
            idx_3 = 0;
            for (k3 = 0; k3 < 2;k3++){
                idx_4 = 0;
                for (k4 = 0; k4 < 2;k4++){   
                    Int_two[k1][k2][k3][k4] = 0.0;                 
                    for(mu = 0;mu<Ct_L[k1];mu++){
                        for(nu = 0;nu<Ct_L[k2];nu++){
                            for(lambda = 0;lambda<Ct_L[k3];lambda++){
                                for(sigma = 0;sigma<Ct_L[k4];sigma++){
                                    Int_two[k1][k2][k3][k4] += Ct_D[k1][mu]*Ct_D[k2][nu]*Ct_D[k3][lambda]*Ct_D[k4][sigma]*ERI_ssss[idx_1+mu][idx_2+nu][idx_3+lambda][idx_4+sigma];
                                }
                            }
                        }
                    } 
                    idx_4+=Ct_L[k4];                   
                }
                idx_3+=Ct_L[k3]; 
            }
            idx_2+=Ct_L[k2]; 
        }
        idx_1+=Ct_L[k1]; 
    }       

    // fill the rest part of Int_two
    int *idx;
    idx = (int*)malloc(sizeof(int)*4);
    for (k1 = 0; k1 < Aorbnum; k1++){
        for (k2 = 0; k2 < Aorbnum; k2++){
            for (k3 = 0; k3 < Aorbnum;k3++){
                for (k4 = 0; k4 < Aorbnum; k4++){
                    idx[0] = k1;idx[1] = k2;idx[2] = k3;idx[3] = k4;
                    int flag = Find_Value(idx);
                    if (flag){
                        Int_two[k1][k2][k3][k4] = Int_two[idx[0]][idx[1]][idx[2]][idx[3]];
                    }
                }
            }
        }
    }    

}

static int Find_Value(int *idx){
    // find ERI index
    Sort_Integer(idx,idx+1);
    Sort_Integer(idx+2,idx+3);
    int i;int sign[4];
    // sort
    if(idx[0]>idx[2]){
        int a,b;
        a = idx[0];b = idx[1];
        idx[0] = idx[2];
        idx[1] = idx[3];
        idx[2] = a;
        idx[3] = b;
    }
    // check parity
    for (i = 0; i < 4; i++){
        sign[i] = 1;       
    }    
    for (i = 0; i < 4; i++){
        sign[Find_M(idx[i])] =  sign[Find_M(idx[i])]*(-1);       
    }

    
    for (i = 0; i < 4; i++){
        if(sign[i] == -1){
            return 0;
        }
    }
    // find corresponding ERI index  
    int flag = Find_M(idx[0]);
    if(flag == 0){
        idx[0] = Trans_idx(1,idx[0]);
        idx[1] = Trans_idx(1,idx[1]);
        idx[2] = Trans_idx(1,idx[2]);
        idx[3] = Trans_idx(1,idx[3]);
    }
    else if (flag == 1){
        if (Find_M(idx[1]) == 3){
            idx[1] = Trans_idx(2,idx[1]);
            idx[3] = Trans_idx(2,idx[3]);            
        }
        else if (Find_M(idx[2]) ==3){
            idx[2] = Trans_idx(2,idx[2]);
            idx[3] = Trans_idx(2,idx[3]);  
        }
    }
    else if (flag ==2){
        if(Find_M(idx[1]) == 3){
            // (py pz, py pz)->(px py, px py)
            idx[0] = Trans_idx(1,idx[0]);
            idx[1] = Trans_idx(2,idx[1]);
            idx[2] = Trans_idx(1,idx[2]);
            idx[3] = Trans_idx(2,idx[3]);            
        }
        else if(Find_M(idx[2]) == 2){
            // (py py, py py)->(px px, px px)
            idx[0] = Trans_idx(1,idx[0]);
            idx[1] = Trans_idx(1,idx[1]);
            idx[2] = Trans_idx(1,idx[2]);
            idx[3] = Trans_idx(1,idx[3]);             
        }
        else{
            // (py py, pz pz)->(px px, py py)
            idx[0] = Trans_idx(1,idx[0]);
            idx[1] = Trans_idx(1,idx[1]);
            idx[2] = Trans_idx(2,idx[2]);
            idx[3] = Trans_idx(2,idx[3]);               
        }
    }
    else if (flag == 3){
        // (pz pz, pz pz)->(px px, py py)
        idx[0] = Trans_idx(1,idx[0]);
        idx[1] = Trans_idx(1,idx[1]);
        idx[2] = Trans_idx(1,idx[2]);
        idx[3] = Trans_idx(1,idx[3]);         
    }
    else{}
    return 1;
}

static int Theta(int x){
    if(x>=0){
        return 1;
    }
    else{
        return 0;
    }
}

static int Trans_idx(int flag,int x){
    if (flag==1){
        return (x-((x-3)/2)*2*Theta(x-3));
    }
    else if (flag==2){
        return (x-((x-5)/2)*2*Theta(x-5));
    }    
    else if (flag==3){
        return x;
    }
    else{
        return -1;
    }
}

static int Find_M(int A){
    // find angular momentum
    if(A>=0&&A<3){
        return 0;
    }
    else if(A>=3&&A<5){
        return 1;
    }
    else if(A>=5&&A<7){
        return 2;
    }
    else if(A>=7&&A<9){
        return 3;
    }
    else{
        return -1;
    }
}

static void Sort_Integer(int *A,int *B){
    int C;
    if(*A > *B){
        C = *A;
        *A = *B;
        *B = C;
    }
}

static int start_C(int k){
    if(k == 0){
        return 0;
    }
    else if(k == 1){
        return 6;
    }
    else if(k == 2){
        return 9;
    }
    else if(k == 3){
        return 0;
    }
    else if(k == 4){
        return 3;
    }
    else if(k == 5){
        return 0;
    } 
    else if(k == 6){
        return 3;
    }      
    else{
        return -1;
    }     
}

static int start_H(int k){
    if(k == 0){
        return 0;
    }
    else if(k == 1){
        return 2;
    }
    else{
        return -1;
    }
}
