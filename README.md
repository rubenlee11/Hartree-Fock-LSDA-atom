# Hartree-Fock-LSDA-atom

## Features

Homework for Theory of Eletronic Structure and Ab Initio Calculation
Only supports H and C atom.

For H, I used 3-21G basis
For C, I used 6-31G basis

performs both Unrestricted Hartree Fock and DFT LSDA calculation.For LSDA, I used PW92 functional.

For Gaussian integral, I performed an over simplified OS algorithm.

For numerical integration, I used Chebyshev and Lebelev grid.

## Usage
edit the main file and makefile, compile.

## reference

[1] Szabo, Attila, and Neil S. Ostlund. Modern quantum chemistry: introduction to advanced electronic structure theory. Courier Corporation, 1996. 

[2] Perdew, J. P., & Wang, Y. (1992). Accurate and simple analytic representation of the electron-gas correlation energy. Physical review B, 45(23), 13244.

[3] Obara, Shigeru, and A. Saika. "Efficient recursive computation of molecular integrals over Cartesian Gaussian functions." The Journal of chemical physics 84.7 (1986): 3963-3974. 

[4] A.Becke, "A multicenter numerical integration scheme for polyatomic molecules", J.Chem.Phys. 88, 2547 (1988)

[5] https://people.sc.fsu.edu/~jburkardt/cpp_src/sphere_lebedev_rule/sphere_lebedev_rule.html

[6] https://people.math.sc.edu/Burkardt/cpp_src/chebyshev2_rule/chebyshev2_rule.html  

[7] https://www.basissetexchange.org
