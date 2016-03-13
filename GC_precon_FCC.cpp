#include <cstdio>
#include <iostream>
#include <cmath>
#include <sstream>
#include "SparseMatrix.h"
#include "VectorOpr.h"
#include "SystemOpr.h"
#include "FCC.h"
#define TOL 1E-6

using std::cout;

void modify(vector<double> &b, const vector<double> &diag)
{
	for (int i = 0; i < b.size(); ++i)
		b[i] = b[i] / sqrt(diag[i]);
}

void _modify(vector<double> &b, const vector<double> &diag)
{
	for (int i = 0; i < b.size(); ++i)
		b[i] = b[i] * sqrt(diag[i]);
}

int main()
{
	int N, M, num_elem;

	SparseMatrix A, precon;
	SystemOpr solve;

	vector<double> b, x, r, v, z, y, diag;

	solve.setType(TRIANGULAR);

	scanf("%d %d %d\n", &N, &M, &num_elem);

	A.init(N, M);
	b.assign(N, 1);
	
	for (int i = 0; i < num_elem; ++i)
	{
		double v;
		int l, c;

		scanf("%d %d %lf \n", &l, &c, &v);

		l--; c--;

		A.put(l, c, v);

		if (l != c)
			A.put(c, l, v);
	}

	A.sortRow();

	for (int i = 0; i < A.nnz(); ++i)
		if (A[i].r == A[i].c)
			diag.push_back(A[i].v);

	cout << "eta;iterações" << endl;

	for (int eta = -N; eta <= N; eta += 10)
	{
		precon = FCC(A, eta);
		b.assign(N, 1);
		//b = A * b;
		x.assign(N, 0);

		int k, max = 2 * N;
		double aux, tmp, s, m;

		r = b;

		modify(b, diag);
		solve(precon, v, b);
		modify(v, diag);

		modify(r, diag);
		solve(precon, y, r);
		modify(y, diag);
		_modify(r, diag);

		aux = y * r;

		for (k = 0; k < max; ++k)
		{
			z = A * v;
	                s = aux / (v * z);
	                x = x + s * v;
	                r = r - s * z;
	                tmp = (r * r);

			modify(r, diag);
			solve(precon, y, r);
			modify(y, diag);
			_modify(r, diag);

		//	for (int i = 0; i < y.size(); ++i)
		//		fprintf(stderr, "%lf\n", y[i]);
		//	printf("%d;%E\n", k + 1, tmp);

	                if (islessequal(tmp, TOL)) break;

			tmp = y * r;
	                m = tmp / aux;
	                aux = tmp;
	                v = y + m * v;
		}

		cout << "ans: " << eta << ";" << k + 1 << endl;
	}

	//for (int i = 0; i < x.size(); ++i)
	//	fprintf(stderr, "%lf\n", x[i]);
}
