#pragma once
#include "svd.h"
#include <assert.h>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

namespace util {
void print(MatDoub mat, string symbol = "") {
  if (symbol.compare(""))
    cout << symbol << "	Matrix " << mat.nrows() << "x" << mat.ncols() << ":"
         << endl;

  for (int m = 0; m < mat.nrows(); m++) {
    for (int n = 0; n < mat.ncols(); n++) {
      cout << setw(15) << mat[m][n] << "\t";
    }
    cout << endl;
  }
  cout << endl;
}

void printDiag(MatDoub mat, string symbol = "") {
  if (symbol.compare(""))
    cout << symbol << "	MatrixDiag " << mat.nrows() << "x" << mat.ncols() << ":"
         << endl;
  double nmax = mat.nrows() < mat.nrows() ? mat.nrows() : mat.nrows();
  for (int n = 0; n < nmax; n++) {
    cout << setw(15) << mat[n][n] << "\t";
  }
  cout << endl;
}

MatDoub diag(VecDoub &V) {
  double m = V.size();
  MatDoub M(m, m);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      M[i][j] = 0;

  for (int i = 0; i < m; i++)
    M[i][i] = V[i];

  return M;
}

void print(VecDoub vec, string symbol = "") {
  if (symbol.compare(""))
    cout << symbol << "	Vector " << vec.size() << "D:" << endl;

  for (int m = 0; m < vec.size(); m++) {
    cout << setw(15) << vec[m];
  }
  cout << endl;
}

MatDoub Transpose(const MatDoub &Mat) {
  MatDoub res(Mat.ncols(), Mat.nrows());
  for (int n = 0; n < res.nrows(); n++) {
    for (int m = 0; m < res.ncols(); m++) {
      res[n][m] = Mat[m][n];
    }
  }
  return res;
}

MatDoub T(const MatDoub &Mat) { return Transpose(Mat); }

double norm(const VecDoub &A) {
  double res = 0;
  for (int i = 0; i < A.size(); i++) {
    res += A[i] * A[i];
  }
  return sqrt(res);
}

VecDoub SVDerror(const SVD &svd) {
  const int n = svd.v.ncols();
  VecDoub res(n);
  for (int j = 0; j < n; j++) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
      sum += std::pow(svd.v[j][i] / svd.w[i], 2);
    }
    res[j] = sqrt(sum);
  }
  return res;
}

// Parses an M x N matrix from a file with the following format:
//   M
//   N
// Row1
//  ...
// RowM
MatDoub parseFile(std::string filename) {
  std::ifstream ex1a(filename);
  int A_rows, A_cols;
  ex1a >> A_rows;
  ex1a >> A_cols;
  MatDoub A(A_rows, A_cols);

  for (int r = 0; r < A.nrows(); r++) {
    for (int c = 0; c < A.ncols(); c++) {
      double a;
      ex1a >> a;
      A[r][c] = a;
    }
  }
  return A;
}

// Converts a matrix consisting of a single column to a vector
// Can be used with parseFile() to parse vectors from files like so:
// VecDoub b = toVec(parseFile("Ex1b.dat"));
VecDoub toVec(MatDoub A) {
  assert(A.ncols() == 1);
  VecDoub vec(A.nrows());
  for (int i = 0; i < A.nrows(); i++) {
    vec[i] = A[i][0];
  }
  return vec;
}

VecDoub vec2(Doub a, Doub b)
{
  VecDoub vec(2);
  vec[0] = a;
  vec[1] = b;
  return vec;
}

VecDoub vec3(double a, double b, double c) {
  VecDoub res(3);
  res[0] = a;
  res[1] = b;
  res[2] = c;
  return res;
}

VecDoub vec4(double a, double b, double c, double d) {
  VecDoub res(4);
  res[0] = a;
  res[1] = b;
  res[2] = c;
  res[3] = d;
  return res;
}

VecDoub randVec(int size) {
  VecDoub v(size);
  for (int i = 0; i < size; i++) {
    v[i] = float(rand()) / float(RAND_MAX);
  }
  return v;
}

Doub dot(const VecDoub &A, const VecDoub &B) {
  Doub sum = 0;
  for (int i = 0; i < A.size(); i++) {
    sum += A[i] * B[i];
  }
  return sum;
}
} // namespace util

MatDoub operator*(const MatDoub &A1, const MatDoub &A2) {
  if (A1.ncols() != A2.nrows()) {
    cerr << "in prod: the number of rows in A1 is not equal to the number of "
            "cols in A2"
         << endl;
  }

  MatDoub res(A1.nrows(), A2.ncols());
  for (int n = 0; n < A1.nrows(); n++) {
    for (int m = 0; m < A2.ncols(); m++) {
      double temp = 0;
      for (int i = 0; i < A1.ncols(); i++) {
        temp += A1[n][i] * A2[i][m];
      }
      res[n][m] = temp;
    }
  }
  return res;
}

VecDoub operator*(const MatDoub &A, const VecDoub &b) {
  if (A.ncols() != b.size()) {
    cerr << "in prod: the number of rows in A is not equal to the size of "
            "vector b"
         << endl;
  }

  VecDoub res(A.nrows());
  for (int n = 0; n < A.nrows(); n++) {
    double temp = 0;
    for (int m = 0; m < A.ncols(); m++) {
      temp += A[n][m] * b[m];
    }
    res[n] = temp;
  }
  return res;
}
VecDoub operator-(const VecDoub &A, const VecDoub &b) {
  if (A.size() != b.size()) {
    cerr << "Not the same size" << endl;
  }

  VecDoub res(A.size());
  for (int i = 0; i < A.size(); i++) {
    res[i] = A[i] - b[i];
  }
  return res;
}

VecDoub operator+(const VecDoub &A, const VecDoub &b) {
  if (A.size() != b.size()) {
    cerr << "Not the same size" << endl;
  }

  VecDoub res(A.size());
  for (int i = 0; i < A.size(); i++) {
    res[i] = A[i] + b[i];
  }
  return res;
}

VecDoub operator*(const VecDoub &A, double c) {
  VecDoub res(A.size());
  for (int i = 0; i < A.size(); i++) {
    res[i] = c * A[i];
  }
  return res;
}

VecDoub operator*(double c, const VecDoub &A) { return A * c; }

VecDoub operator/(const VecDoub &A, Doub c) { return A * (1. / c); }

VecDoub vecZeros(int size) {
  VecDoub res(size);
  for (int row = 0; row < size; row++) {
    res[row] = 0.;
  }
  return res;
}

VecDoub residuals(const MatDoub &A, const VecDoub &x, const VecDoub &b) {
  VecDoub residuals = A * x - b;
  return residuals;
}

Doub residualError(const VecDoub &residuals, const VecDoub &b) {
  Doub residual_error = util::norm(residuals) / util::norm(b);
  return residual_error;
}

Doub isRandom(const Doub &residual_error, const MatDoub &A) {
  Doub comparison = sqrt( Doub((A.nrows() - A.ncols())) / Doub(A.nrows()));
  Doub bullshit_ratio = residual_error / comparison;
  return bullshit_ratio;
}