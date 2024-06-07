/*
 * utilities.h
 *
 *  Created on: Feb 2, 2015
 *      Author: jpb
 */

// ANSI color codes
#define RESET "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <string>

using namespace std;

namespace util
{
	void print(MatDoub mat, string symbol = "")
	{
		if (symbol.compare(""))
			cout << YELLOW << symbol << "	Matrix " << mat.nrows() << "x" << mat.ncols() << ":" << RESET << endl;

		for (int m = 0; m < mat.nrows(); m++)
		{
			for (int n = 0; n < mat.ncols(); n++)
			{
				cout << setw(15) << mat[m][n] << "\t";
			}
			cout << endl;
		}
		cout << endl;
		cout << "\n";
	}

	void printDiag(MatDoub mat, string symbol = "")
	{
		if (symbol.compare(""))
			cout << YELLOW << symbol << "	MatrixDiag " << mat.nrows() << "x" << mat.ncols() << ":" << RESET << endl;
		double nmax = mat.nrows() < mat.nrows() ? mat.nrows() : mat.nrows();
		for (int n = 0; n < nmax; n++)
		{
			cout << setw(15) << mat[n][n] << "\t";
		}
		cout << endl;
	}

	MatDoub diag(VecDoub &V)
	{
		double m = V.size();
		MatDoub M(m, m);

		for (int i = 0; i < m; i++)
			for (int j = 0; j < m; j++)
				M[i][j] = 0;

		for (int i = 0; i < m; i++)
			M[i][i] = V[i];

		return M;
	}

	void print(VecDoub vec, string symbol = "")
	{
		if (symbol.compare(""))
			cout << YELLOW << symbol << "	Vector " << vec.size() << "D:" << RESET << endl;

		for (int m = 0; m < vec.size(); m++)
		{
			cout << setw(15) << vec[m] << endl;
		}
		cout << endl;
	}

	MatDoub Transpose(const MatDoub &Mat)
	{
		MatDoub res(Mat.ncols(), Mat.nrows());
		for (int n = 0; n < res.nrows(); n++)
		{
			for (int m = 0; m < res.ncols(); m++)
			{
				res[n][m] = Mat[m][n];
			}
		}
		return res;
	}

	MatDoub T(const MatDoub &Mat)
	{
		return Transpose(Mat);
	}

}

MatDoub operator*(const MatDoub &A1, const MatDoub &A2)
{
	if (A1.ncols() != A2.nrows())
	{
		cerr << "in prod: the number of rows in A1 is not equal to the number of cols in A2" << endl;
	}

	MatDoub res(A1.nrows(), A2.ncols());
	for (int n = 0; n < A1.nrows(); n++)
	{
		for (int m = 0; m < A2.ncols(); m++)
		{
			double temp = 0;
			for (int i = 0; i < A1.ncols(); i++)
			{
				temp += A1[n][i] * A2[i][m];
			}
			res[n][m] = temp;
		}
	}
	return res;
}

VecDoub operator*(const MatDoub &A, const VecDoub &b)
{
	if (A.ncols() != b.size())
	{
		cerr << "in prod: the number of rows in A is not equal to the size of vector b" << endl;
	}

	VecDoub res(A.nrows());
	for (int n = 0; n < A.nrows(); n++)
	{
		double temp = 0;
		for (int m = 0; m < A.ncols(); m++)
		{
			temp += A[n][m] * b[m];
		}
		res[n] = temp;
	}
	return res;
}

VecDoub operator-(const VecDoub &a, const VecDoub &b)
{
	if (a.size() != b.size())
	{
		cerr << "in prod: the size of vector a is not equal to the size of vector b" << endl;
	}

	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = a[n] - b[n];
	}
	return res;
}

VecDoub operator-(const VecDoub &a, const double &b)
{
	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = a[n] - b;
	}
	return res;
}

// function to add a scalar to a vector
VecDoub operator+(const VecDoub &a, const double &b)
{
	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = a[n] + b;
	}
	return res;
}

// function to add a vector to a vector
VecDoub operator+(const VecDoub &a, const VecDoub &b)
{
	if (a.size() != b.size())
	{
		cerr << "in prod: the size of vector a is not equal to the size of vector b" << endl;
	}

	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = a[n] + b[n];
	}
	return res;
}

// function to divide a vector by a vector
VecDoub operator/(const VecDoub &a, const VecDoub &b)
{
	if (a.size() != b.size())
	{
		cerr << "in prod: the size of vector a is not equal to the size of vector b" << endl;
	}

	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = a[n] / b[n];
	}
	return res;
}

// function to divide a vector by a scalar
VecDoub operator/(const VecDoub &a, const Doub &b)
{
	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = a[n] / b;
	}
	return res;
}

// function to multiply a vector by a scalar
VecDoub operator*(const VecDoub &a, const Doub &b)
{
	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = a[n] * b;
	}
	return res;
}

// function to raise a scalar to a vector
VecDoub operator^(const Doub &b, const VecDoub &a)
{
	VecDoub res(a.size());
	for (int n = 0; n < a.size(); n++)
	{
		res[n] = pow(b, a[n]);
	}
	return res;
}

#endif /* UTILITIES_H_ */
