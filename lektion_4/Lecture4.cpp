#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../source_code/Numerical-Recipes-master/svd.h"
#include "utilities.h"

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

int main()
{
	/* -------------------------------- read data ------------------------------- */
	VecDoub xFilip(82);
	VecDoub yFilip(82);
	ifstream Filip("FilipData.dat");
	for (int i = 0; i < 82; i++)
	{
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}

	VecDoub xPont(40);
	VecDoub yPont(40);
	ifstream Pont("PontiusData.dat");
	for (int i = 0; i < 40; i++)
	{
		Pont >> yPont[i];
		Pont >> xPont[i];
	}

	// create design matrix
	MatDoub APont(40, 3);
	for (int i = 0; i < 40; i++)
	{
		APont[i][0] = 1;
		APont[i][1] = xPont[i];
		APont[i][2] = xPont[i] * xPont[i];
	}

	// create vector b
	VecDoub BPont(40);
	for (int i = 0; i < 40; i++)
	{
		BPont[i] = yPont[i];
	}

	/* --------------------------------- pontius -------------------------------- */

	cout << GREEN << "Pontius" << RESET << endl;
	/* ---------------------------- LU decomposition ---------------------------- */
	//  solve by LU decomposition
	auto C = util::Transpose(APont) * APont;
	// util::print(C, "C");

	auto c = util::Transpose(APont) * BPont;

	LUdcmp alu(C);
	VecDoub x_LU(3);
	alu.solve(c, x_LU);

	util::print(x_LU, "x_LU");
	cout << endl;

	// find variance of the parameters
	alu.inverse(C);
	util::printDiag(C, "Variance of estimated parameters: ");
	cout << endl;

	/* -------------------------------- cholesky -------------------------------- */
	// solve by cholosky decomposition
	Cholesky chol(C);
	VecDoub aChol(3);
	chol.solve(c, aChol);

	util::print(aChol, "x_Chol");
	cout << endl;

	/* ---------------------------------- svd --------------------------------- */
	// solve by SVD
	SVD svd_pont(APont);
	VecDoub xSvd_pont(3);
	svd_pont.solve(BPont, xSvd_pont, THRESHOLD);
	util::print(xSvd_pont, "x_SVD");
	cout << endl;

	// errors on parameters
	VecDoub errors_pont = svd_pont.errors();
	util::print(errors_pont, "errors");
	std::cout << std::endl;

	// residuals
	VecDoub residuals = APont * xSvd_pont - BPont;

	util::print(residuals, "residuals");
	cout << endl;

	/* ---------------------------------- filip --------------------------------- */
	std::cout << GREEN << "Filip" << RESET << std::endl;
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
	VecDoub BFilip(82);
	for (int i = 0; i < 82; i++)
	{
		BFilip[i] = yFilip[i];
	}
	auto CFilip = util::Transpose(AFilip) * AFilip;
	// util::print(CFilip, "CFilip");

	auto cFilip = util::Transpose(AFilip) * BFilip;

	LUdcmp aluFilip(CFilip);
	VecDoub aFilip(filip_params);
	aluFilip.solve(cFilip, aFilip);

	util::print(aFilip, "x_LU");
	cout << endl;

	// find variance of the parameters
	MatDoub C2Filip;
	aluFilip.inverse(C2Filip);
	util::printDiag(C2Filip, "Variance of estimated parameters: ");
	cout << endl;

	/* -------------------------------- cholesky -------------------------------- */
	cout << RED << "Cholesky fails because the matrix is not positive semi definite" << RESET << endl;
	cout << endl;
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
	cout << endl;

	util::print(xSvd_Filip, "xSVD");

	VecDoub errors = svd_filip.errors();
	util::print(errors, "errors");
	return 0;
}
