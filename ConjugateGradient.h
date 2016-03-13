#ifndef ConjugateGradient_H
#define ConjugateGradient_H

#include "SparseMatrix.h"
#include "VectorOpr.h"
#include "SystemOpr.h"
#include "FCC.h"
#include <vector>
#include <ctime>
#define TOL 1E-6
using namespace std;

vector<double> GC_precon(SparseMatrix &M, SparseMatrix& A, vector<double> &b, int &iter, clock_t &clk);
vector<double> modify(vector<double> &b, vector<double> &diag);

#endif
