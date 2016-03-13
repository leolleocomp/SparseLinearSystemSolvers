#include "VectorOpr.h"
#include <cstdio>
#include <cstdlib>

vector<double> operator+(const vector<double> a, const vector<double> b)
{
	if (a.size() != b.size())
	{
		printf("bad vector summation!\n");
		exit(EXIT_FAILURE);
	}

	vector<double> c;

	for (int i = 0; i < a.size(); ++i)
		c.push_back(a[i] + b[i]);

	return c;
}


vector<double> operator-(const vector<double> a, const vector<double> b)
{
	if (a.size() != b.size())
	{
		printf("bad vector subtraction!\n");
		exit(EXIT_FAILURE);
	}

	vector<double> c;

	for (int i = 0; i < a.size(); ++i)
		c.push_back(a[i] - b[i]);

	return c;
}


double operator*(const vector<double> a, const vector<double> b)
{
	if (a.size() != b.size())
	{
		printf("bad vector dot product!\n");
		exit(EXIT_FAILURE);
	}

	double sum = 0.0;
	
	for (int i = 0; i < a.size(); ++i)
		sum += a[i] * b[i];

	return sum;
}

vector<double> operator*(double k, const vector<double> v)
{
	vector<double> ans;

	for (int i = 0; i < v.size(); ++i)
		ans.push_back(k * v[i]);

	return ans;
}
