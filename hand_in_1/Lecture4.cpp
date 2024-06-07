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

// find errors on parameters
VecDoub errors(SVD svd)
{
	VecDoub errors(svd.n);
	for (int j = 0; j < svd.n; j++)
	{
		double error = 0;
		for (int i = 0; i < svd.n; i++)
		{
			error += (svd.v[j][i] / svd.w[i]) * (svd.v[j][i] / svd.w[i]);
		}
		errors[j] = sqrt(error);
	}
	return errors;
}

VecDoub sigma_i(double delta, VecDoub residuals)
{
	VecDoub sigma_i(residuals.size());
	for (int i = 0; i < residuals.size(); i++)
	{
		sigma_i[i] = abs(residuals[i]) < delta ? delta : abs(residuals[i]);
	}

	return sigma_i;
}

int main()
{

	/* -------------------------------- hand in 1 ------------------------------- */

	/* ------------------------------- exercise 1 ------------------------------- */
	std::cout
		<< GREEN << "Exercise 1" << RESET << std::endl;

	/* ------------------------------- exercise 1:i ----------------------------- */
	std::cout << GREEN << "i" << RESET << std::endl;
	// read in data
	ifstream i_A_data("Ex1A.dat");
	// read out M and N
	int M, N;
	i_A_data >> M;
	i_A_data >> N;
	std::cout << "M: " << M << " N: " << N << std::endl;

	VecDoub i_matrix_1(40);
	VecDoub i_matrix_2(40);
	VecDoub i_matrix_3(40);
	VecDoub i_matrix_4(40);
	VecDoub i_matrix_5(40);
	VecDoub i_matrix_6(40);
	for (int i = 0; i < 40; i++)
	{
		i_A_data >> i_matrix_1[i];
		i_A_data >> i_matrix_2[i];
		i_A_data >> i_matrix_3[i];
		i_A_data >> i_matrix_4[i];
		i_A_data >> i_matrix_5[i];
		i_A_data >> i_matrix_6[i];
	}

	MatDoub A(40, 6);
	for (int i = 0; i < 40; i++)
	{
		A[i][0] = i_matrix_1[i];
		A[i][1] = i_matrix_2[i];
		A[i][2] = i_matrix_3[i];
		A[i][3] = i_matrix_4[i];
		A[i][4] = i_matrix_5[i];
		A[i][5] = i_matrix_6[i];
	}
	util::print(A, "A");

	// read data
	ifstream i_b_data("Ex1B.dat");
	// read out M and N
	i_b_data >> M;
	i_b_data >> N;
	std::cout << "M: " << M << " N: " << N << std::endl;
	VecDoub b(M);
	for (int i = 0; i < M; i++)
	{
		i_b_data >> b[i];
	}
	util::print(b, "b");

	// make SVD
	SVD svd_exi(A);
	VecDoub xSvd_exi(6);
	svd_exi.solve(b, xSvd_exi, THRESHOLD);
	cout << endl;
	util::print(svd_exi.w, "svd_w_exi");
	cout << endl;

	/* ------------------------------- exercise 1:ii ----------------------------- */
	std::cout << GREEN << "Exercise ii" << RESET << std::endl;
	std::cout << "\n";
	util::print(xSvd_exi, "svd_X_exi");

	cout << endl;

	/* ------------------------------- exercise 1:iii ---------------------------- */
	std::cout << GREEN << "Exercise iii" << RESET << std::endl;
	std::cout << "\n";
	VecDoub errors_exi = svd_exi.errors();
	util::print(errors_exi, "errors");
	cout << endl;

	/* ------------------------------- exercise 1:iv ---------------------------- */
	std::cout << GREEN << "Exercise iv" << RESET << std::endl;
	std::cout << "\n";

	VecDoub residuals_exi = A * xSvd_exi - b;
	util::print(residuals_exi, "residuals");

	/* ------------------------------- exercise 1:v ---------------------------- */
	std::cout << GREEN << "Exercise v" << RESET << std::endl;
	std::cout << "\n";

	VecDoub sigma = sigma_i(1, residuals_exi);
	util::print(sigma, "sigma");

	// create new A
	MatDoub A_sigma = A;
	for (int i = 0; i < sigma.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			A_sigma[i][j] /= sigma[i];
		}
	}
	util::print(A_sigma, "A_sigma");

	// create new b
	VecDoub b_sigma = b;
	for (int i = 0; i < sigma.size(); i++)
	{
		b_sigma[i] /= sigma[i];
	}
	util::print(b_sigma, "b_sigma");

	std::cout << "b_sigma[6]:" << b_sigma[6] << endl;
	std::cout << "A_sigma[0][0]:" << A_sigma[0][0] << endl;

	/* ------------------------------- exercise 1:vi ---------------------------- */
	std::cout << GREEN << "Exercise vi" << RESET << "\n"
			  << std::endl;

	// make SVD on new values
	SVD svd_sigma_exi(A_sigma);
	VecDoub xSvd_sigma_exi(6);
	svd_sigma_exi.solve(b_sigma, xSvd_sigma_exi, THRESHOLD);
	cout << endl;
	util::print(xSvd_sigma_exi, "x_svd_exi");
	cout << endl;

	return 0;
}
