#ifndef VECTOR_OPR_H
#define VECTOR_OPR_H

#include <vector>
using namespace std;

/**
 *	-- Vector Operation definitions
 *	
 *	@author: leolleo
 */

// vector sum
vector<double> operator+(const vector<double> a, const vector<double> b);

// vector sub
vector<double> operator-(const vector<double> a, const vector<double> b);

// dot product
double operator*(const vector<double> a, const vector<double> b);

// scaling
vector<double> operator*(double k, const vector<double> v);

#endif
