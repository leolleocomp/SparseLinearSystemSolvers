#MakeFile, this is fun, bro!
CC=g++
CFLAGS=-c

all:novos_testes

novos_testes: ConjugateGradient.o FCC.o SparseMatrix.o VectorOpr.o SystemOpr.o
	$(CC)  main.cpp ConjugateGradient.o FCC.o SparseMatrix.o VectorOpr.o SystemOpr.o -o novos_testes

ConjugateGradient.o: ConjugateGradient.cpp
	$(CC) $(CFLAGS) ConjugateGradient.cpp  

FCC.o: FCC.cpp
	$(CC) $(CFLAGS) FCC.cpp

SparseMatrix.o: SparseMatrix.cpp
	$(CC) $(CFLAGS) SparseMatrix.cpp

VectorOpr.o: VectorOpr.cpp
	$(CC) $(CFLAGS) VectorOpr.cpp

SystemOpr.o: SystemOpr.cpp
	$(CC) $(CFLAGS) SystemOpr.cpp
