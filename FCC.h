#ifndef FCC_H
#define FCC_H

#include "SparseMatrix.h"
#define IT_MAX 25

bool comp_aux(pair<double, int> a, pair<double, int>b);
SparseMatrix FCC(const SparseMatrix &A, int eta);

#endif
