#ifndef SYSTEM_OPR_H
#define SYSTEM_OPR_H

/**
 * - A linear sytem operator
 *	
 *   intended to perform final operations towards solving
 *   a given linear system (ex:. diagonal, triangular sytems)
 */

#include "SparseMatrix.h"
#include "VectorOpr.h"

enum SystemType { DIAGONAL, TRIANGULAR };

class SystemOpr
{
public:

	void solve(const SparseMatrix &A, vector<double> &x, const vector<double> &b);
	void operator()(const SparseMatrix &A, vector<double> &x, const vector<double> &b);
	void setType(int type);

private:

	int type;

	void solveTriang(const SparseMatrix &A, vector<double> &x, const vector<double> &b);
	void    solveDiag(const SparseMatrix &A, vector<double> &x, const vector<double> &b);
};

#endif
