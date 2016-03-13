#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H
#include <vector>
using namespace std;

/**
 * - Sparse Matrix data structure implementation
 *
 *   @author leolleo
 */

/**
 * - struct element: an element of the sparse matrix
 *
 *   @elem1: 	row index
 *   @elem2:	column index
 *   @elem3:	element value 	
 */

typedef struct element
{
	int r, c;
	double v;

	bool operator<(const struct element &b) const {

		if (r < b.r)
			return true;	
		else if (r == b.r) {
			return c < b.c;
		} else
			return false;
	}

} element;

/**
 * - SparseMatrix: A sparse matrix data structure implementation
 *   
 *   @matrix:		vector containing the nonzero matrix elements.
 *   @nnl_count: 	vector containing the number of nonzero elements
 *   	         	on each line of the matrix.
 *
 *   @N: 	number of rows.
 *   @M:	number of columns.
 */

class SparseMatrix
{
	// operators overload
	friend vector<double> operator*(const SparseMatrix &A, const vector<double> &v);
	friend vector<double> operator*(const vector<double> &v, const SparseMatrix &A);
	friend SparseMatrix operator*(const SparseMatrix &A, const SparseMatrix &B);
	friend SparseMatrix operator*(double k, const SparseMatrix &B);
	friend SparseMatrix operator+(const SparseMatrix &A, const SparseMatrix &B);
	friend SparseMatrix operator-(const SparseMatrix &A, const SparseMatrix &B);
	friend SparseMatrix operator~(const SparseMatrix &A);

public:
	// constructors
	SparseMatrix();
	SparseMatrix(int N, int M);
	void assign(const SparseMatrix &toCpy);

	// initialization methods
	void init(int N, int M);

	// access methods
	void    put(int i, int j, double val);
	int     get(int i, int j) const;
	int     l_size() const;
	int     c_size() const;
	element operator[](int i) const;

	// generic methods
	int  nnz()  	const;
	int  nnl(int l) const;
	void sortRow();
	vector<double> getDiagonal();

	// print methods
	void printSparse() const;
	void printDense()  const;

private:
	vector<element>	matrix;
	vector<int> nnl_count;
	int N, M;
};

#endif
