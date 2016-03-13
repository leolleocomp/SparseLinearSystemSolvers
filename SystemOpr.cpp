#include "SystemOpr.h"
#include <cstdio>

void SystemOpr::solve(const SparseMatrix &A, vector<double> &x, const vector<double> &b)
{
	if (type == DIAGONAL)
		solveDiag(A, x, b);
	else if (type == TRIANGULAR)
		solveTriang(A, x, b);
	else
		printf("not a existent type of system operator...\n");
}

void SystemOpr::operator()(const SparseMatrix &A, vector<double> &x, const vector<double> &b)
{
	solve(A, x, b);
}

void SystemOpr::setType(int t)
{
	if (t == DIAGONAL || t == TRIANGULAR)
		type = t;
	else
		printf("bad type setting, the desired system type doesn't exist\n");
}	

/**
 * - SystemOpr::solveDiag()
 *
 *   assumes A as a diagonal coeficient matrix with a not null diagonal,
 *   as well as pre-initialized params.
 *
 */

void SystemOpr::solveDiag(const SparseMatrix &A, vector<double> &x, const vector<double> &b)
{
	for (int i = 0; i < A.l_size(); ++i)
		x[i] = b[i] / A[A.get(i, i)].v;
}

void SystemOpr::solveTriang(const SparseMatrix &A, vector<double> &x, const vector<double> &b)
{
	vector<double> y(b);

	for (int i = 0; i < A.nnz();++i)
	{
		if (A[i].r > A[i].c)
			y[A[i].r] -= A[i].v * y[A[i].c];
		else if (A[i].r == A[i].c)
			y[A[i].r] /= A[i].v;
	}

	x.assign(y.begin(), y.end());

	for (int i = A.nnz() - 1; i >= 0; --i)
	{
		if (A[i].r < A[i].c)
			x[A[i].r] -= A[i].v * x[A[i].c];
		else if (A[i].r == A[i].c)
			x[A[i].r] /= A[i].v;
	}
}

