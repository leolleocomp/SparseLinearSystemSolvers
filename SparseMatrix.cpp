#include "SparseMatrix.h"
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <cfloat>
#include <cmath>
using namespace std;

SparseMatrix::SparseMatrix()
{
	init(0, 0);
}

SparseMatrix::SparseMatrix(int N, int M)
{
	init(N, M);
}

void SparseMatrix::assign(const SparseMatrix &toCpy)
{
	matrix.assign(toCpy.matrix.begin(), toCpy.matrix.end());	
	nnl_count.assign(toCpy.nnl_count.begin(), toCpy.nnl_count.begin());	
	init(toCpy.N, toCpy.M);
}

void SparseMatrix::init(int N, int M)
{
	if (N < 0 || M < 0)
	{
		printf("bad SparseMatrix initialization!\n");
		exit(EXIT_FAILURE);
	}

	this->N = N;
	this->M = M;

	matrix.clear();
	nnl_count.assign(N, 0);
}

void SparseMatrix::put(int i, int j, double val)
{
	if (i < 0 || i >= N || j < 0 || j >= M)
	{
		printf("bad SparseMatrix element insertion!\n");
		exit(EXIT_FAILURE);
	}

	if (fabs(val) <= DBL_EPSILON) return;

	element toAdd;

	toAdd.r = i;
	toAdd.c = j;
	toAdd.v = val;

	matrix.push_back(toAdd);
	nnl_count[i]++;
}

int SparseMatrix::get(int i, int j) const
{
	if (i < 0 || i >= N || j < 0 || j >= M)
	{
		printf("bad SparseMatrix get operation!\n");
		exit(EXIT_SUCCESS);
	}

	element toSearch;

	toSearch.r = i;
	toSearch.c = j;
	toSearch.v = 0.0;

	int pos = lower_bound(matrix.begin(), matrix.end(), toSearch) - matrix.begin();

	if (pos != matrix.size() && matrix[pos].r == i && matrix[pos].c == j) 
		return pos;

	return -1;
}

int SparseMatrix::l_size() const
{
	return N;
}

int SparseMatrix::c_size() const
{
	return M;
}

element SparseMatrix::operator[](int i) const
{
	if (i < 0 || i >= matrix.size())
	{
		printf("bad SparseMatrix operator[] usage!\n");
		exit(EXIT_FAILURE);
	}

	return matrix[i];
}

int SparseMatrix::nnz() const
{
	return matrix.size();
}

int SparseMatrix::nnl(int l) const
{
	if (l < 0 || l >= N)
	{
		printf("bad SparseMatrix nnl method call!\n");
		exit(EXIT_FAILURE);
	}

	return nnl_count[l];
}

void SparseMatrix::sortRow()
{
	sort(matrix.begin(), matrix.end());
}

vector<double> SparseMatrix::getDiagonal()
{
	vector<double> diag;

	for (int i = 0; i < nnz(); ++i)
		if (matrix[i].r == matrix[i].c)
			diag.push_back(matrix[i].v);

	return diag;
}

void SparseMatrix::printSparse() const
{
	for (int i = 0; i < nnz(); ++i)
		printf("%d %d %E\n", matrix[i].r, matrix[i].c, matrix[i].v);
}

void SparseMatrix::printDense() const
{
	//toComplete!
}

vector<double> operator*(const SparseMatrix &A, const vector<double> &v)
{
	if (v.size() != A.M)
	{
		printf("bad SparseMatrix vector multiplication!\n");
		exit(EXIT_FAILURE);
	}

	vector<double> out;

	out.assign(A.N, 0.0);

	for (int i = 0; i < A.nnz(); ++i)
	{
		out[A[i].r] += A[i].v * v[A[i].c]; 
	}

	return out;
}

vector<double> operator*(const vector<double> &v, const SparseMatrix &A)
{
	if (v.size() != A.N)
	{
		printf("bad vector x SparseMatrix multiplication!\n");
		exit(EXIT_FAILURE);
	}

	vector<double> out;

	out.assign(A.M, 0.0);

	for (int i = 0; i < A.nnz(); ++i)
	{
		out[A[i].c] += A[i].v * v[A[i].c];
	}

	return out;
}

SparseMatrix operator*(const SparseMatrix &A, const SparseMatrix &B)
{
	if (A.M != B.N)
	{
		printf("bad SparseMatrix x SparseMatrix multiplication!\n");
		exit(EXIT_FAILURE);
	}

	// toImplement!
}

SparseMatrix operator*(double k, const SparseMatrix &B)
{
	SparseMatrix C(B.N, B.M);
	
	for (int i = 0; i < B.nnz(); ++i)
		C.put(B[i].r, B[i].c, B[i].v * k);

	return C;
}

SparseMatrix operator+(const SparseMatrix &A, const SparseMatrix &B)
{
	if (A.N != B.N || A.M != B.M)
	{
		printf("bad SparseMatrix sum, incompatible matrices!\n");
		exit(EXIT_FAILURE);
	}

	SparseMatrix C;
	int it1, it2;

	C.init(A.N, B.M);

	it1 = it2 = 0;

	while (it1 < A.nnz() || it2 < B.nnz())
	{
		if (it1 < A.nnz() && it2 < B.nnz())
		{
			if (A[it1] < B[it2])
			{
				C.put(A[it1].r, A[it1].c, A[it1].v);
				it1++;
			}
			else if (B[it2] < A[it1])
			{
				C.put(B[it2].r, B[it2].c, B[it2].v); 
				it2++;
			}
			else
			{
				C.put(A[it1].r, A[it1].c, A[it1].v + B[it2].v);
				it1++;
				it2++;
			}
		} 
		else if (it1 < A.nnz())
		{
			C.put(A[it1].r, A[it1].c, A[it1].v);
			it1++;
		}
		else
		{
			C.put(B[it2].r, B[it2].c, B[it2].v);
			it2++;
		}
	}

	return C;
}

SparseMatrix operator-(const SparseMatrix &A, const SparseMatrix &B)
{
	if (A.N != B.N || A.M != B.M)
	{
		printf("bad SparseMatrix sum, incompatible matrices!\n");
		exit(EXIT_FAILURE);
	}

	SparseMatrix C;
	int it1, it2;

	C.init(A.N, B.M);

	it1 = it2 = 0;

	while (it1 < A.nnz() || it2 < B.nnz())
	{
		if (it1 < A.nnz() && it2 < B.nnz())
		{
			if (A[it1] < B[it2])
			{
				C.put(A[it1].r, A[it1].c, A[it1].v);
				it1++;
			}
			else if (B[it2] < A[it1])
			{
				C.put(B[it2].r, B[it2].c, B[it2].v); 
				it2++;
			}
			else
			{
				C.put(A[it1].r, A[it1].c, A[it1].v - B[it2].v);
				it1++;
				it2++;
			}
		} 
		else if (it1 < A.nnz())
		{
			C.put(A[it1].r, A[it1].c, A[it1].v);
			it1++;
		}
		else
		{
			C.put(B[it2].r, B[it2].c, B[it2].v);
			it2++;
		}
	}

	return C;
}

SparseMatrix operator~(const SparseMatrix &A)
{
	vector<int> nnc;

	nnc.assign(A.M, 0);

	for (int i = 0; i < A.nnz(); ++i)
		nnc[A[i].c]++;

	vector<int> acm_freq;
	acm_freq.assign(A.M, 0);
	acm_freq[0] = nnc[0];

	for (int i = 1; i < acm_freq.size(); ++i)
		acm_freq[i] = nnc[i] + acm_freq[i-1];

	int col[A.nnz()], row[A.nnz()], val[A.nnz()];

	nnc.assign(A.M, 0);

	for (int i = 0; i < A.nnz(); ++i)
	{
		int index;

		if (A[i].c > 0)
			index = acm_freq[A[i].c - 1] + nnc[A[i].c];
		else
			index = A[i].c + nnc[A[i].c];
			
		col[index]  = A[i].c;
		row[index]  = A[i].r;
		val[index]  = A[i].v;

		nnc[A[i].c]++;
	}

	SparseMatrix C(A.M, A.N);

	for (int i = 0; i < A.nnz(); ++i)
	{
		C.put(col[i], row[i], val[i]);	
	}

	return C;
}
