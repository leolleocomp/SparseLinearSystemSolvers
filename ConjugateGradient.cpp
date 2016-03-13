#include "ConjugateGradient.h"
#include <cmath>

vector<double> modify(vector<double> &b, vector<double> &diag)
{
	vector<double> ans;

	for (int i = 0; i < b.size(); ++i)
		ans.push_back(b[i] / sqrt(diag[i]));

	return ans;
}

vector<double> GC_precon(SparseMatrix &M, SparseMatrix &A, vector<double> &b, int &iter, clock_t &clk)
{
	vector<double> x, r, y, z, v, diag;
	int k, max = 2 * A.l_size();
	double aux, tmp, s, m;
	SystemOpr solve;

	iter = 0;

	solve.setType(TRIANGULAR);
	
	clk = clock();
	x.assign(A.l_size(), 0);
	diag = A.getDiagonal();
	r = b;

	solve(M, v, modify(b, diag));
	v = modify(v, diag);

	solve(M, y, modify(r, diag));
	y = modify(y, diag);

	aux = y * r;
	for (k = 0; k < max; ++k)
	{
		z = A * v;
		s = aux / (v * z);
		x = x + s * v;
		r = r - s * z;
		tmp = (r * r);

		solve(M, y, modify(r, diag));
		y = modify(y, diag);

		if (islessequal(tmp, TOL)) break;

		tmp = y * r;
		m = tmp / aux;
		aux = tmp;
		v = y + m * v;
		iter++;
	}

	clk = clock() - clk;
	return x;
}
