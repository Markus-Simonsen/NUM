#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/qrdcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../source_code/Numerical-Recipes-master/svd.h"
#include "../source_code/Numerical-Recipes-master/roots.h"
#include "../source_code/Numerical-Recipes-master/roots_multidim.h"
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

double eq_1(double x_0, double x_1)
{
	double result = x_0 + 2 * sin(x_1 - x_0) - exp(-sin(x_1 + x_0));
	return result;
}

double eq_2(double x_0, double x_1)
{
	return x_0 * cos(x_1) + sin(x_0) - 1;
}

VecDoub vecfunc(VecDoub_I x)
{
	VecDoub result(2);
	// result[0] = eq_1(x[0], x[1]);
	// result[1] = eq_2(x[0], x[1]);
	result[0] = x[0] + 2 * sin(x[1] - x[0]) - exp(-sin(x[1] + x[0]));
	result[1] = x[0] * cos(x[1]) + sin(x[0]) - 1;
	return result;
}

void print_header()
{
	std::cout << YELLOW << std::setw(5) << "i" << std::setw(15) << "x_i" << std::setw(15) << "d_k" << std::setw(15) << "C" << std::setw(15) << "e_est" << std::setw(15) << "e_true" << RESET << std::endl;
}

void print_table(double i, double x_i, double x_im1, double x_im2, double k, std::string root_finder)
{
	double d_k = fabs(x_i - x_im1);
	double d_km1 = fabs(x_im1 - x_im2);
	double C = fabs(x_i - x_im1) / pow(fabs(x_im1 - x_im2), k);
	double e_est;
	if (root_finder == "bisection")
	{
		e_est = fabs(d_k);
	}
	else if (root_finder == "ridders")
	{
		e_est = C * pow(fabs(d_k), k);
	}
	else if (root_finder == "secant")
	{
		e_est = C * pow(fabs(d_k), k);
	}
	else if (root_finder == "newton")
	{
		e_est = C * pow(fabs(d_k), k);
	}
	double x_true = 0.73908513321516064165531208767387;
	double e_true = fabs(x_i - x_true);

	std::cout << std::setw(5) << i << std::setw(15) << x_i << std::setw(15) << d_k << std::setw(15) << C << std::setw(15) << e_est << std::setw(15) << e_true << std::endl;
}

int main()
{
	util::title("equations with x_0 and x_1 = 0");

	double x_0 = 1;
	double x_1 = 1;

	std::cout << "eq_1: " << eq_1(x_0, x_1) << std::endl;
	std::cout << "eq_2: " << eq_2(x_0, x_1) << std::endl;

	VecDoub_IO x(2);
	x[0] = 1;
	x[1] = 2;

	bool check = true;

	util::title("newt");
	newt(x, check, vecfunc);

	cout << "solution" << endl;

	std::cout << "x_0: " << x[0] << std::endl;
	std::cout << "x_1: " << x[1] << std::endl;

	return 0;
}