#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "utilities.h"

// include library for svd
#include "../source_code/Numerical-Recipes-master/svd.h"

#include <fstream>

using namespace std;

int main()
{
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

	// util::print(xFilip, "xFilip");
	// util::print(yFilip, "yFilip");

	// util::print(xPont, "xPont");
	// util::print(yPont, "yPont");

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

	// util::print(A, "A");
	//  solve by LU decomposition
	auto C = util::Transpose(APont) * APont;
	// util::print(C, "C");

	auto c = util::Transpose(APont) * BPont;

	LUdcmp alu(C);
	VecDoub x_LU(3);
	alu.solve(c, x_LU);

	util::print(x_LU, "a");

	// find variance of the parameters
	alu.inverse(C);
	util::printDiag(C, "Variance of estimated parameters: ");

	// solve by cholosky decomposition
	Cholesky chol(C);
	VecDoub aChol(3);
	chol.solve(c, aChol);

	util::print(aChol, "aChol");

	/* ---------------------------------- svd --------------------------------- */
	// util::print(ASVD, "ASVD");
	// solve by LU decomposition
	SVD svd_pont(APont);
	VecDoub xSvd_pont(3);
	svd_pont.solve(BPont, xSvd_pont);
	util::print(xSvd_pont, "xSVD");

	// residuals
	VecDoub residuals = APont * xSvd_pont - BPont;

	util::print(residuals, "residuals");

	std::cout << "Filip" << std::endl;

	/* ---------------------------------- filip --------------------------------- */
	// create design matrix
	MatDoub AFilip(82, 11);
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

	// util::print(AFilip, "AFilip");
	//  solve by LU decomposition
	auto CFilip = util::Transpose(AFilip) * AFilip;
	// util::print(CFilip, "CFilip");

	auto cFilip = util::Transpose(AFilip) * BFilip;

	LUdcmp aluFilip(CFilip);
	VecDoub aFilip(11);
	aluFilip.solve(cFilip, aFilip);

	// util::print(aFilip, "aFilip");

	// find variance of the parameters
	MatDoub C2Filip;
	aluFilip.inverse(C2Filip);
	util::printDiag(C2Filip, "Variance of estimated parameters: ");
	cout << endl;

	// solve by cholosky decomposition
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
