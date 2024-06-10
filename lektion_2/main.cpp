#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../utility_code/utilities.h"

#include <fstream>

using namespace std;

// read in the data
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
	// read in the data
	VecDoub x_data_Filip(82);
	VecDoub y_data_Filip(82);
	read_data_2_columns("FilipData.dat", x_data_Filip, y_data_Filip);
	VecDoub x_data_Pont(40);
	VecDoub y_data_Pont(40);
	read_data_2_columns("PontiusData.dat", x_data_Pont, y_data_Pont);
	/* --------------------------------- pontius -------------------------------- */
	// create design matrix
	MatDoub A_filip(40, 3);
	for (int i = 0; i < 40; i++)
	{
		A_filip[i][0] = 1;
		A_filip[i][1] = x_data_Pont[i];
		A_filip[i][2] = x_data_Pont[i] * x_data_Pont[i];
	}

	// create vector b_filip
	VecDoub b_filip = y_data_Pont;

	util::print(A_filip, "A_filip");
	// solve by LU decomposition
	auto C = util::Transpose(A_filip) * A_filip;
	auto c = util::Transpose(A_filip) * b_filip;

	LUdcmp alu(C);
	VecDoub a(3);
	alu.solve(c, a);

	util::print(a, "a");

	// find variance of the parameters
	alu.inverse(C);
	util::printDiag(C, "Variance of estimated parameters: ");

	// solve by cholesky decomposition
	Cholesky chol(C);
	VecDoub aChol(3);
	chol.solve(c, aChol);

	util::print(aChol, "aChol");

	/* ---------------------------------- filip --------------------------------- */
	// create design matrix
	MatDoub AFilip(82, 11);
	for (size_t i = 0; i < 82; i++)
	{
		for (size_t j = 0; j < 11; j++)
		{
			AFilip[i][j] = pow(x_data_Filip[i], j);
		}
	}

	// create vector b_filip
	VecDoub bFilip = y_data_Filip;

	//  solve by LU decomposition
	auto CFilip = util::Transpose(AFilip) * AFilip;
	auto cFilip = util::Transpose(AFilip) * bFilip;

	LUdcmp aluFilip(CFilip);
	VecDoub aFilip(11);
	aluFilip.solve(cFilip, aFilip);

	// find variance of the parameters
	MatDoub C2Filip;
	aluFilip.inverse(C2Filip);
	util::printDiag(C2Filip, "Variance of estimated parameters: ");

	// solve by cholosky decomposition
	Cholesky cholFilip(CFilip);
	VecDoub aCholFilip(11);
	cholFilip.solve(cFilip, aCholFilip);

	return 0;
}
