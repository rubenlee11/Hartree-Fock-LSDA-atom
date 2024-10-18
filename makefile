# edit lapack,openblas library and include directories before compiling
CC = gcc -Wall -I/usr/local/opt/lapack/include -I/usr/local/opt/openblas/include
LIB = -L/usr/local/opt/lapack/lib -L/usr/local/opt/openblas/lib -L/usr/local/opt/libc/lib -llapacke -llapack -lblas

scf : main.o init.o scf.o common.o quadrature.o functional.o
	$(CC) $(LIB) -o scf main.o init.o scf.o common.o quadrature.o functional.o

main.o : main.cpp common.h
	$(CC) -c main.cpp

init.o : init.cpp common.h
	$(CC) $(LIB) -c init.cpp

scf.o : scf.cpp common.h
	$(CC) $(LIB) -c scf.cpp

common.o : common.cpp
	$(CC) $(LIB) -c common.cpp

quadrature.o: quadrature.cpp common.h
	$(CC) $(LIB) -c quadrature.cpp 

functional.o: functional.cpp
	$(CC) $(LIB) -c functional.cpp

clean: 
	rm -f *.o