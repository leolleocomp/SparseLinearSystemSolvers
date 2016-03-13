#include "FCC.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cfloat>

using namespace std;

typedef vector< vector< pair<int , double> > > SparseTmp;

bool comp_aux(pair<double, int> a, pair<double, int>b)
{
	return fabs(a.first) < fabs(b.first);
}

SparseMatrix FCC(const SparseMatrix &A, int eta)
{
	SparseMatrix G(A.l_size(), A.c_size());
	vector<double> diag;
	
	for (int i = 0; i < A.nnz(); ++i)
		if (A[i].c == A[i].r)
			diag.push_back(A[i].v);

	for (int i = 0; i < A.nnz(); ++i)
		G.put(A[i].r, A[i].c, A[i].v / sqrt(diag[A[i].c] * diag[A[i].r]));

	G.sortRow();

	SparseTmp L;
	vector< pair<double, int> > col_elem;

	int l, mj;
	bool fail;
	double sigma;

	l = 0;
	sigma = 0.0;

	do
	{
		fail = false;

		diag.clear();
		L.assign(G.l_size(), vector< pair<int, double> > ());

		for (int j = 0; j < G.c_size(); ++j)
		{
			mj = 0;
			int gj = 0;

			col_elem.clear();

			for (int i = j + 1; i < G.l_size(); ++i)
			{
				double sum = 0.0;
				int k1, k2; 

				k1 = k2 = 0;

				while ((k1 < L[i].size()) && (k2 < L[j].size()) && (L[i][k1].first < j - 1) && (L[j][k2].first < j - 1))
				{
					if (L[i][k1].first == L[j][k2].first)
						sum += L[i][k1].second * diag[L[i][k1].first] * L[j][k2].second;

					if (L[i][k1].first < L[j][k2].first)
						k1++;
					else if (L[j][k2].first < L[i][k1].first)
						k2++;	
					else
					{
						k1++;
						k2++;
					}
				}

				int pos = G.get(i, j);	


				if (pos >= 0)
				{
					col_elem.push_back(make_pair(G[pos].v - sum, i));
					mj++;
					if (fabs(G[pos].v - sum) > DBL_EPSILON) gj++;
				}
				else 
				{
					col_elem.push_back(make_pair(-sum, i));
					if (fabs(sum) > DBL_EPSILON) gj++; 	
				}
			}

			sort(col_elem.rbegin(), col_elem.rend(), comp_aux);

			double sum = 0.0;

			for (int k = 0; k < L[j].size() && L[j][k].first < j - 1; ++k)
				sum += diag[L[j][k].first] * L[j][k].second * L[j][k].second;

			diag.push_back(1 - sum + sigma);

			if (diag[j] < DBL_EPSILON)
			{
				l++;
				sigma = 5E-4 * pow(2, l);
				fail = true;
				fprintf(stderr, "%d %lf<<\n", l, sigma); 
				break;
			}

			int limite_superior = min(G.l_size() - (j + 1), mj + eta);

			for (int i = 0; i < limite_superior; ++i)
				if (fabs(col_elem[i].first) > DBL_EPSILON)
				{
					L[col_elem[i].second].push_back(make_pair(j, col_elem[i].first / diag[j]));				
				}

		}

	} while (fail && l < IT_MAX);

	if (!fail && l < IT_MAX)
	{
		G.init(A.l_size(), A.c_size());

		for (int i = 0; i < L.size(); ++i)
		{
			G.put(i, i, sqrt(diag[i]));

			for (int j = 0; j < L[i].size(); ++j)
			{
				G.put(i, L[i][j].first, L[i][j].second * sqrt(diag[L[i][j].first]));
				G.put(L[i][j].first, i, L[i][j].second * sqrt(diag[L[i][j].first]));
			}
		}

		G.sortRow();
	} else
		fprintf(stderr, "FCC ERROR!\n");

	return G;
}
