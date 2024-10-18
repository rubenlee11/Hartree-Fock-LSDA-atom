#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>
#include "common.h"

void Mix(){
    /****************************************************
    Mix the density matrix with the previous one 
    ****************************************************/        
	double P1,P2;
	double E_xc0,E_xc1,E_ee;
	// form new density matrix
	for(int mu=0; mu<Aorbnum;mu++){
		for(int nu=0; nu<Aorbnum;nu++){
			P1 = P2 = 0.0;
			for(int i=0; i<N_up;i++){
				P1 += C[0][mu][i]*C[0][nu][i];
			}
			for(int i=0; i<N_down;i++){
				P2 += C[1][mu][i]*C[1][nu][i];
			}
			P[0][mu][nu] = Mixing*P1 + (1.0-Mixing)*P[0][mu][nu];
			P[1][mu][nu] = Mixing*P2 + (1.0-Mixing)*P[1][mu][nu];
			P_tot[mu][nu] = P[0][mu][nu] + P[1][mu][nu];                        
		}
	}    
	// form new eletron electron matrix
	for(int mu=0; mu<Aorbnum;mu++){
		for(int nu=0; nu<Aorbnum;nu++){
			E_ee = 0.0;
			for(int lambda=0; lambda<Aorbnum;lambda++){
				for(int sigma=0; sigma<Aorbnum;sigma++){
					E_ee += P_tot[lambda][sigma]*Int_two[mu][nu][lambda][sigma];      
				}
			}
			H_ee[mu][nu] = E_ee;                  
		}
	}   	
	// form new xc interaction matrix
	if(Xc_flag == 0){			
		for(int mu=0; mu<Aorbnum;mu++){
			for(int nu=0; nu<Aorbnum;nu++){
				E_xc0 = E_xc1 = 0.0;
				for(int lambda=0; lambda<Aorbnum;lambda++){
					for(int sigma=0; sigma<Aorbnum;sigma++){    
						E_xc0 += -1.*P[0][lambda][sigma]*Int_two[mu][lambda][sigma][nu];
						E_xc1 += -1.*P[1][lambda][sigma]*Int_two[mu][lambda][sigma][nu];
					}
				}
				H_xc[0][mu][nu] = E_xc0;
				H_xc[1][mu][nu] = E_xc1;                     
			}
		}   
	}
	else if(Xc_flag == 1 && SCF_Step>0){        
		Cal_rho();
        printf("Total charge: %f\n",Quadrature(rho[0])+Quadrature(rho[1]));
		double *result;
		result = (double *)malloc(sizeof(double)*2);
		for(int mu=0; mu<Aorbnum;mu++){
			for(int nu=0; nu<Aorbnum;nu++){                
				Cal_xc_matrix(result,mu,nu);
				H_xc[0][mu][nu] = result[0];
				H_xc[1][mu][nu] = result[1];                                        
			}
		}
   		free(result);
	}
	// form new Hamiltonian matrix
	Set_Zero(Aorbnum,H[0]);
	Set_Zero(Aorbnum,H[1]);   
	Addition(Aorbnum,H[0],H_core);Addition(Aorbnum,H[0],H_ee);Addition(Aorbnum,H[0],H_xc[0]);
	Addition(Aorbnum,H[1],H_core);Addition(Aorbnum,H[1],H_ee);Addition(Aorbnum,H[1],H_xc[1]);
}

void Total_Energy(){
    double Energy = 0.0; 
    if(Xc_flag == 0){
        for(int mu = 0; mu < Aorbnum; mu++){
            for(int nu = 0; nu < Aorbnum; nu++){
                Energy+=0.5*(P_tot[mu][nu]*H_core[mu][nu]+H[0][mu][nu]*P[0][mu][nu]+H[1][mu][nu]*P[1][mu][nu]);
            }
        }
        Etot[SCF_Step] = Energy;
    } 
    else if(Xc_flag == 1){		
		double Energy_one = 0.0;
        // one body energy
        for(int mu = 0; mu < Aorbnum; mu++){
            for(int nu = 0; nu < Aorbnum; nu++){
				Energy_one = 0.0;
				for(int lambda = 0; lambda< Aorbnum;lambda++){
					for(int sigma = 0; sigma< Aorbnum; sigma++){
						Energy_one+=P_tot[lambda][sigma]*Int_two[mu][nu][lambda][sigma];
					}
				}
                Energy+=P_tot[mu][nu]*H_core[mu][nu]+0.5*Energy_one*P_tot[mu][nu];
            }
        }
		// xc energy
		if(SCF_Step>0){
			Cal_rho();
			Energy += Cal_xc_energy();
			Etot[SCF_Step] = Energy;
		}
		else{
			Etot[SCF_Step] = Energy;
		}
    }
}

void SCF(){    
    int info_0,info_1;
    double D_Etot;

    // SCF 0
    Mix();
    for(int mu = 0; mu < Aorbnum; mu++){
        for(int nu = 0; nu < Aorbnum; nu++){
            C[0][mu][nu] = H[0][mu][nu];
            C[1][mu][nu] = H[1][mu][nu];
        }
    }
    for(int i=0; i<Aorbnum;i++){
        for(int j=0; j<Aorbnum; j++){
            S_data[i*Aorbnum+j] = S[i][j];
        }
    }    
    info_0 = LAPACKE_dsygvd(LAPACK_ROW_MAJOR,1, 'V', 'U', Aorbnum, C_data[0], Aorbnum, S_data,Aorbnum, Eorb[0]);
    for(int i=0; i<Aorbnum;i++){
        for(int j=0; j<Aorbnum; j++){
            S_data[i*Aorbnum+j] = S[i][j];
        }
    }        
    info_1 = LAPACKE_dsygvd(LAPACK_ROW_MAJOR,1, 'V', 'U', Aorbnum, C_data[1], Aorbnum, S_data,Aorbnum, Eorb[1]);     
    // calculate the total energy   
    SCF_Step = 0;
    Total_Energy();            
    printf("E = %.7f\n",Etot[0]); 

    for(SCF_Step=1;SCF_Step<=Max_SCF;SCF_Step++){  
        printf("**********************SCF = %-2d*************************\n",SCF_Step);                      
        Mix();
        /*         
         dsygvd write the eigenvectors over the input matrix A 
         the Cholesky factorization over B
         In order to keep the value of H and S unchanged 
         I copy H to C and using C as input, and copy S to S_dummy
         */        
        for(int mu = 0; mu < Aorbnum; mu++){
            for(int nu = 0; nu < Aorbnum; nu++){
                C_data[0][mu*Aorbnum+nu] = H[0][mu][nu];
                C_data[1][mu*Aorbnum+nu] = H[1][mu][nu];
            }
        }
        for(int i=0; i<Aorbnum;i++){
            for(int j=0; j<Aorbnum; j++){
                S_data[i*Aorbnum+j] = S[i][j];
            }
        }                    
                         
        // solve the generalized eigenvalue problem
        info_0 = LAPACKE_dsygvd(LAPACK_ROW_MAJOR,1, 'V', 'U', Aorbnum, C_data[0], Aorbnum, S_data,Aorbnum, Eorb[0]);
        for(int i=0; i<Aorbnum;i++){
            for(int j=0; j<Aorbnum; j++){
                S_data[i*Aorbnum+j] = S[i][j];
            }
        }        
        info_1 = LAPACKE_dsygvd(LAPACK_ROW_MAJOR,1, 'V', 'U', Aorbnum, C_data[1], Aorbnum, S_data,Aorbnum, Eorb[1]);            

        // find the occupation states
        N_up = N_down = 0;        
        for(int i=0;i<Enum;i++){
            if (Eorb[0][N_up]<=Eorb[1][N_down]){            
                Fermi_Energy = Eorb[0][N_up];
                N_up += 1;
            }
            else{            
                Fermi_Energy = Eorb[1][N_down];  
                N_down += 1;      
            }        
        } 
        printf("energy of HOMO: %.8f\n",Fermi_Energy);
        printf("up = %-2d, down = %-2d\n",N_up,N_down);
        // calculate the total energy   
        Total_Energy();   
        if(SCF_Step > 1){
            D_Etot = Etot[SCF_Step] - Etot[SCF_Step-1];
            // convergence test
            if(abs(D_Etot)<=E_Conv){
                printf("\nConvergence reached at SCF = %d\n",SCF_Step);
                printf("E = %.10f, dE = %.10f\n",Etot[SCF_Step], D_Etot);   
                printf("orbital energy\n");
                Print_Matrix(2,Aorbnum,Eorb);   
                printf("spin up eigen vector\n");
                Print_Matrix(Aorbnum,Aorbnum,C[0]); 
                printf("spin down eigen vector\n");
                Print_Matrix(Aorbnum,Aorbnum,C[1]);                                               
                break;
            }   
            else{
                printf("E = %.10f, dE = %.10f\n",Etot[SCF_Step], D_Etot);
            }              
        }       
        else{
            D_Etot = 0.0;
            printf("E = %.10f, dE = %.10f\n",Etot[SCF_Step], D_Etot);
        }                     
    }
}

void Wave_function(){}
