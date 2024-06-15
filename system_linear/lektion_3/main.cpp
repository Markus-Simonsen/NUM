#include <iostream>
#include <fstream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../utility_code/utilities.h"

// include library for SVD
#include "../source_code/Numerical-Recipes-master/svd.h"

using namespace std;

bool read_data_2_columns(const string &filename, VecDoub_IO &x, VecDoub_IO &y)
{
	ifstream file(filename);
	if (!file.is_open())
	{
		cerr << "Could not open file " << filename << endl;
		return false;
	}
	double x_val, y_val;
	int i = 0;
	while (file >> x_val >> y_val)
	{
		x[i] = x_val;
		y[i] = y_val;
		i++;
	}
	file.close();
	return true;
}

bool polynomial_Amatrix(MatDoub &A, VecDoub &x, int n)
{
	int m = x.size();
	A.resize(m, n);
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = pow(x[i], j);
		}
	}
	return true;
}

int main()
{
	VecDoub xFilip(82);
	VecDoub yFilip(82);
	read_data_2_columns("FilipData.dat", xFilip, yFilip);

	VecDoub xPont(40);
	VecDoub yPont(40);
	read_data_2_columns("PontData.dat", xPont, yPont);

	// create design matrix
	MatDoub APont(40, 3);
	for (int i = 0; i < 40; i++)
	{
		APont[i][0] = 1;
		APont[i][1] = xPont[i];
		APont[i][2] = xPont[i] * xPont[i];
	}

	// create vector b
	VecDoub BPont = yPont;

	//  solve by LU decomposition
	auto C = util::Transpose(APont) * APont;
	auto c = util::Transpose(APont) * BPont;

	LUdcmp alu(C);
	VecDoub x_LU(3);
	alu.solve(c, x_LU);

	util::print(x_LU, "a");

	// find variance of the parameters
	alu.inverse(C);
	util::printDiag(C, "Variance of estimated parameters: ");

	// solve by Cholesky decomposition
	Cholesky chol(C);
	VecDoub aChol(3);
	chol.solve(c, aChol);

	util::print(aChol, "aChol");

	/* ---------------------------------- svd --------------------------------- */
	// solve by SVD decomposition
	SVD svd_pont(APont);
	VecDoub xSvd_pont(3);
	svd_pont.solve(BPont, xSvd_pont);
	util::print(xSvd_pont, "xSVD");

	// residuals
	VecDoub residuals = APont * xSvd_pont - BPont;

	util::print(residuals, "residuals");

	/* ---------------------------------- filip --------------------------------- */
	util::title("Filip");

	// create design matrix
	MatDoub AFilip(82, 11);
	for (int i = 0; i < 82; i++)
	{
		for (size_t j = 0; j < 11; j++)
		{
			AFilip[i][j] = pow(xFilip[i], j);
		}
	}

	// create vector b
	VecDoub BFilip = yFilip;

	//  solve by LU decomposition
	auto CFilip = util::Transpose(AFilip) * AFilip;
	auto cFilip = util::Transpose(AFilip) * BFilip;

	LUdcmp aluFilip(CFilip);
	VecDoub aFilip(11);
	aluFilip.solve(cFilip, aFilip);

	// find variance of the parameters
	MatDoub C2Filip;
	aluFilip.inverse(C2Filip);
	util::printDiag(C2Filip, "Variance of estimated parameters: ");
	cout << endl;

	// solve by Cholesky decomposition
	Cholesky cholFilip(CFilip);
	VecDoub aCholFilip(11);
	cholFilip.solve(cFilip, aCholFilip);

	/* ---------------------------------- svd --------------------------------- */
	SVD svd_filip(AFilip);
	VecDoub xSvd_Filip(3);
	svd_filip.solve(BFilip, xSvd_Filip);

	// residuals
	VecDoub residualsFilip = AFilip * xSvd_Filip - BFilip;

	util::print(residualsFilip, "residuals");

	util::print(xSvd_Filip, "xSVD");
	return 0;
}
