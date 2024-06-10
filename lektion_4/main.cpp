#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../source_code/Numerical-Recipes-master/svd.h"
#include "../utility_code/utilities.h"

#include <fstream>

// ANSI color codes
#define RESET "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"

#define THRESHOLD 1.0e-10

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

int main()
{
	/* -------------------------------- read data ------------------------------- */
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

	/* --------------------------------- pontius -------------------------------- */
	util::title("Pontius");
	/* ---------------------------- LU decomposition ---------------------------- */
	//  solve by LU decomposition
	auto C = util::Transpose(APont) * APont;
	auto c = util::Transpose(APont) * BPont;

	LUdcmp alu(C);
	VecDoub x_LU(3);
	alu.solve(c, x_LU);

	util::print(x_LU, "x_LU");

	// find variance of the parameters
	alu.inverse(C);
	util::printDiag(C, "Variance of estimated parameters: ");

	/* -------------------------------- cholesky -------------------------------- */
	// solve by cholosky decomposition
	Cholesky chol(C);
	VecDoub aChol(3);
	chol.solve(c, aChol);

	util::print(aChol, "x_Chol");

	/* ---------------------------------- svd --------------------------------- */
	// solve by SVD
	SVD svd_pont(APont);
	VecDoub xSvd_pont(3);
	svd_pont.solve(BPont, xSvd_pont, THRESHOLD);
	util::print(xSvd_pont, "x_SVD");

	// errors on parameters
	VecDoub errors_pont = svd_pont.errors();
	util::print(errors_pont, "errors");

	// residuals
	VecDoub residuals = APont * xSvd_pont - BPont;

	util::print(residuals, "residuals");

	/* ---------------------------------- filip --------------------------------- */
	util::title("Filip");
	int filip_params = 11;

	/* ---------------------------- LU decomposition ---------------------------- */

	// create design matrix
	MatDoub AFilip(82, filip_params);
	for (int i = 0; i < 82; i++)
	{
		AFilip[i][0] = 1;
		AFilip[i][1] = xFilip[i];
		AFilip[i][2] = pow(xFilip[i], 2);
		AFilip[i][3] = pow(xFilip[i], 3);
		AFilip[i][4] = pow(xFilip[i], 4);
		AFilip[i][5] = pow(xFilip[i], 5);
		AFilip[i][6] = pow(xFilip[i], 6);
		AFilip[i][7] = pow(xFilip[i], 7);
		AFilip[i][8] = pow(xFilip[i], 8);
		AFilip[i][9] = pow(xFilip[i], 9);
		AFilip[i][10] = pow(xFilip[i], 10);
	}

	// create vector b
	VecDoub BFilip = yFilip;

	auto CFilip = util::Transpose(AFilip) * AFilip;
	auto cFilip = util::Transpose(AFilip) * BFilip;

	LUdcmp aluFilip(CFilip);
	VecDoub aFilip(filip_params);
	aluFilip.solve(cFilip, aFilip);

	util::print(aFilip, "x_LU");

	// find variance of the parameters
	MatDoub C2Filip;
	aluFilip.inverse(C2Filip);
	util::printDiag(C2Filip, "Variance of estimated parameters: ");

	/* -------------------------------- cholesky -------------------------------- */
	cerr << "Cholesky fails because the matrix is not positive semi definite" << endl;

	// Cholesky cholFilip(CFilip);
	// VecDoub xCholFilip(11);
	// cholFilip.solve(cFilip, xCholFilip);

	/* ---------------------------------- svd --------------------------------- */
	SVD svd_filip(AFilip);
	VecDoub xSvd_Filip(filip_params);
	svd_filip.solve(BFilip, xSvd_Filip, THRESHOLD);

	// residuals
	VecDoub residualsFilip = AFilip * xSvd_Filip - BFilip;

	util::print(residualsFilip, "residuals");

	util::print(xSvd_Filip, "xSVD");

	VecDoub errors = svd_filip.errors();
	util::print(errors, "errors");
	return 0;
}
