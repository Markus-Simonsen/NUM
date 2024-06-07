#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/qrdcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../source_code/Numerical-Recipes-master/svd.h"
#include "../source_code/Numerical-Recipes-master/roots.h"
#include "../source_code/Numerical-Recipes-master/roots_multidim.h"
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

double f_1(double x)
{
	return cos(pow(x, 2)) * exp(-x);
}

double f_2(double x)
{
	return sqrt(x) * cos(pow(x, 2)) * exp(-x);
}

double f_3(double x)
{
	return 1 / sqrt(x) * cos(pow(x, 2)) * exp(-x);
}

double f_4(double x)
{
	return 1000 * exp(-1 / x) * exp(-1 / (1 - x));
}

double extended_midpoint(double f(double x), double a, double b)
{

	double h_old = 0;
	double integral;
	double integral_m1 = 0;
	double integral_m2 = 0;
	double alp_k = 0;
	double error = 0;
	int N = 1;
	int i = 1;
	int f_eval;
	int max_iter = 10;
	std::cout << setw(5) << "i" << setw(15) << "A(hi)" << setw(15) << "A(hi-1)-A(hi)" << setw(15) << "rich-alp^k" << setw(15) << "rich-fejl" << setw(15) << "f-beregninger" << endl;
	while (i < max_iter) // while number of steps is less than 200
	{
		f_eval = 0;
		double h = (b - a) / N;
		integral = 0;
		for (double j = a + 0.5 * h; j < b; j += h)
		{
			integral += f(j) * h;
			f_eval++;
		}
		if (i > 2) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (integral_m2 - integral_m1) / (integral_m1 - integral);
			error = (integral - integral_m1) / (pow(2, alp_k) - 1);
		}
		std::cout << setw(5) << i << setw(15) << integral << setw(15) << integral_m1 - integral << setw(15) << alp_k << setw(15) << error << setw(15) << f_eval << endl;
		integral_m2 = integral_m1;
		integral_m1 = integral;
		h_old = h;
		N *= 2;
		i++;
	}
	return integral;
}

double trapz(double f(double x), double a, double b)
{

	double h_old = 0;
	double integral = 0;
	double integral_m1 = 0;
	double integral_m2 = 0;
	double alp_k = 0;
	double error = 0;
	int N = 1;
	int i = 1;
	int max_iter = 10;
	int f_eval;
	std::cout << setw(5) << "i" << setw(15) << "A(hi)" << setw(15) << "A(hi-1)-A(hi)" << setw(15) << "rich-alp^k" << setw(15) << "rich-fejl" << setw(15) << "f-beregninger" << endl;
	while (i < max_iter)
	{
		f_eval = 0;
		double h = (b - a) / N;
		integral = 0;
		for (double j = a; j < b; j += h)
		{
			integral += (f(j) + f(j + h)) * h / 2;
			f_eval++;
		}

		if (i > 2)
		{
			alp_k = (integral_m2 - integral_m1) / (integral_m1 - integral);
			error = (integral - integral_m1) / (pow(2, alp_k) - 1);
		}

		std::cout << setw(5) << i << setw(15) << integral << setw(15) << integral_m1 - integral << setw(15) << alp_k << setw(15) << error << setw(15) << f_eval << endl;

		integral_m2 = integral_m1;
		integral_m1 = integral;
		h_old = h;
		N *= 2;
		i++;
	}
	return integral;
}

double simpsons(double f(double x), double a, double b)
{

	double h_old = 0;
	double integral = 0;
	double integral_m1 = 0;
	double integral_m2 = 0;
	double alp_k = 0;
	double error = 0;
	int N = 1;
	int i = 1;
	int max_iter = 10;
	int f_eval;

	std::cout << setw(5) << "i" << setw(15) << "A(hi)" << setw(15) << "A(hi-1)-A(hi)" << setw(15) << "rich-alp^k" << setw(15) << "rich-fejl" << setw(15) << "f-beregninger" << endl;
	while (i < max_iter)
	{
		f_eval = 0;
		double h = (b - a) / N; // step-size
		integral = 0;
		for (double j = a; j < b; j += 2 * h)
		{
			integral += (f(j) + 4 * f(j + h) + f(j + h * 2)) * h / 3;
			f_eval += 3;
		}

		if (i > 2)
		{
			alp_k = log2((integral_m2 - integral_m1) / (integral_m1 - integral));
			error = (integral - integral_m1) / (pow(2, alp_k) - 1);
		}

		std::cout << setw(5) << i << setw(15) << integral << setw(15) << integral_m1 - integral << setw(15) << alp_k << setw(15) << error << setw(15) << f_eval << endl;

		integral_m2 = integral_m1;
		integral_m1 = integral;
		h_old = h;
		N *= 2;
		i++;
	}
	return integral;
}

int main()
{
	std::cout << GREEN << "Integrate from a = 0, b = 1: f(x) = cos(x^2)*exp(-x)" << RESET << std::endl;

	double a = 0;
	double b = 1;

	std::cout
		<< GREEN << "Extended midpoint" << RESET << std::endl;
	double integral = extended_midpoint(f_1, a, b);
	std::cout << "integral = " << integral << std::endl;

	std::cout
		<< GREEN << "Trapezoidal rule" << RESET << std::endl;
	integral = trapz(f_1, a, b);
	std::cout << "integral = " << integral << std::endl;

	std::cout
		<< GREEN << "Simpson's rule" << RESET << std::endl;
	integral = simpsons(f_1, a, b);
	std::cout << "integral = " << integral << std::endl;

	std::cout << GREEN << "\nIntegrate from a = 0, b = 1: f(x) = sqrt(x)*cos(x^2)*exp(-x)" << RESET << std::endl;

	// use simpsons
	std::cout
		<< GREEN << "Simpson's rule" << RESET << std::endl;
	integral = simpsons(f_2, a, b);
	std::cout << "integral = " << integral << std::endl;

	/* ------------------------------- function 3 ------------------------------- */
	std::cout << GREEN << "\nIntegrate from a = 0, b = 1: f(x) = 1/sqrt(x)*cos(x^2)*exp(-x)" << RESET << std::endl;

	// midpoint
	std::cout
		<< GREEN << "Extended midpoint" << RESET << std::endl;
	integral = extended_midpoint(f_3, a, b);
	std::cout << "integral = " << integral << std::endl;

	/* ------------------------------- function 4 ------------------------------- */
	std::cout << GREEN << "\nIntegrate from a = 0, b = 1: f(x) = 1000*exp(-1/x)*exp(-1/(1-x))" << RESET << std::endl;

	// trapezoidal
	std::cout
		<< GREEN << "Trapezoidal rule" << RESET << std::endl;
	integral = trapz(f_4, a, b);
	std::cout << "integral = " << integral << std::endl;

	return 0;
}