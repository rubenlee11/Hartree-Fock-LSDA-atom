/****************************************************
Features:
        single atom SCF solver with spin polarization 
    Functional: 
        Unrestricted Hartree-Fock:
            Szabo, Attila, and Neil S. Ostlund. Modern quantum chemistry: introduction to advanced electronic structure theory. Courier Corporation, 1996. 
        LSDA:
            Perdew, J. P., & Wang, Y. (1992). Accurate and simple analytic representation of the electron-gas correlation energy. Physical review B, 45(23), 13244.

    Basis: 
        https://www.basissetexchange.org

    Gaussian integral:
        Obara, Shigeru, and A. Saika. "Efficient recursive computation of molecular integrals over Cartesian Gaussian functions." The Journal of chemical physics 84.7 (1986): 3963-3974.    

    Integration Grid: 
        scheme:
            A.Becke, "A multicenter numerical integration scheme for polyatomic molecules", J.Chem.Phys. 88, 2547 (1988)
        lebedev rule:
            https://people.sc.fsu.edu/~jburkardt/cpp_src/sphere_lebedev_rule/sphere_lebedev_rule.html
        chebyshev rule:
            https://people.math.sc.edu/Burkardt/cpp_src/chebyshev2_rule/chebyshev2_rule.html    

Installation:
    see the makefile

Test system: 
    H, C
****************************************************/

#include "common.h"
#include<stdio.h>
#include<time.h>

static double start_time;
static void Record_time();

int main(){
    // setup in atomic units.
    start_time = clock();
    Z = 6;
    Xc_flag = 0;
    E_Conv = 0.000001;
    Max_SCF = 40;
    Mixing = 0.9;

    // initialize variables
    printf("Init parameters...\n");
    Init_Basis();  
    Init_Arrays();  

    Record_time();
    // integrals calculation
    printf("calculating integrals...\n");  
    if(Z == 6){
        Integral_C();
    }
    else if(Z == 1){
        Integral_H();
    }
    else{
        printf("At present we only support H and C......\n");
        return 0;
    }
	
	// Init real space grid for LSDA calculation
	if(Xc_flag == 1){
        printf("Init integration grid...\n");
		Init_Quadrature();
	}	
	Record_time();

    // SCF calculation
    if(Xc_flag == 0){
        printf("SCF: Unrestricted Hartree Fock\n");
    }
    else if(Xc_flag == 1){
        printf("SCF: Unrestricted DFT using PW92 functional\n");
    }
    SCF();
    Record_time();

    // postprocess
    //Wave_function();
};

static void Record_time(){
    printf("The program has beening running for %f seconds\n",(clock()-start_time)/CLOCKS_PER_SEC);
}
