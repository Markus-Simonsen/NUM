#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/qrdcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../source_code/Numerical-Recipes-master/svd.h"
#include "../source_code/Numerical-Recipes-master/roots.h"
#include "../source_code/Numerical-Recipes-master/roots_multidim.h"
#include "../source_code/Numerical-Recipes-master/quadrature.h"
#include "../source_code/Numerical-Recipes-master/derule.h"

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

void title(string title)
{
	cout << endl
		 << GREEN << title << RESET << endl;
}

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

double f_1(double x, double delta)
{
	return cos(pow(x, 2)) * exp(-x);
}

double f_2(double x, double delta)
{
	if (x < 0.1)
		return sqrt(delta) * cos(pow(delta, 2)) * exp(-delta);
	else
		return sqrt(x) * cos(pow(x, 2)) * exp(-x);
}

double f_3(double x, double delta)
{
	if (x < 0.1)
		return 1 / sqrt(delta) * cos(pow(delta, 2)) * exp(-delta);
	else
		return 1 / sqrt(x) * cos(pow(x, 2)) * exp(-x);
}

double f_4(double x, double delta)
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
		// f_eval = 0;
		double h = (b - a) / N;
		integral = 0;
		for (double j = a + 0.5 * h; j < b; j += h)
		{
			integral += f(j) * h;
			// f_eval++;
		}
		if (i > 2) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (integral_m2 - integral_m1) / (integral_m1 - integral);
			error = (integral - integral_m1) / (pow(2, alp_k) - 1);
		}
		f_eval = 1 + 2 * pow(2, i);
		std::cout
			<< setw(5) << i << setw(15) << integral << setw(15) << integral_m1 - integral << setw(15) << alp_k << setw(15) << error << setw(15) << f_eval << endl;
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
		// f_eval = 0;
		double h = (b - a) / N;
		integral = 0;
		for (double j = a; j < b; j += h)
		{
			integral += (f(j) + f(j + h)) * h / 2;
			// f_eval++;
		}

		if (i > 2)
		{
			alp_k = (integral_m2 - integral_m1) / (integral_m1 - integral);
			error = (integral - integral_m1) / (pow(2, alp_k) - 1);
		}

		f_eval = 1 + 2 * pow(2, i);
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

void print_de_rule(DErule<double(Doub, Doub)> de)
{
	int its = 6;
	// print header
	double A, A_m1 = 0, A_m2 = 0, alp_k = 0, error = 0, f_eval = 1;

	std::cout << setw(5) << "i" << setw(15) << "A(hi)" << setw(15) << "A(hi-1)-A(hi)" << setw(15) << "rich-alp^k" << setw(15) << "rich-fejl" << setw(15) << "f-beregninger" << endl;
	for (size_t i = 0; i < its; i++)
	{
		A = de.next();
		if (i > 1)
		{
			alp_k = (A_m2 - A_m1) / (A_m1 - A);
			error = (A - A_m1) / (alp_k - 1);
		}
		f_eval = 1 + 2 * pow(2, i);
		std::cout << setw(5) << i << setw(15) << A << setw(15) << A_m1 - A << setw(15) << alp_k << setw(15) << error << setw(15) << f_eval << endl;
		A_m2 = A_m1;
		A_m1 = A;
	}
}

int main()
{
	double a = 0, b = 1;

	/* ------------------------------- function 1 ------------------------------- */

	title("Integrate from a = 0, b = 1: f(x) = cos(x^2)*exp(-x)");

	// title("Extended midpoint");
	double integral = extended_midpoint(f_1, a, b);
	// std::cout << "integral = " << integral << std::endl;

	// title("Trapezoidal rule");
	// integral = trapz(f_1, a, b);
	// std::cout << "integral = " << integral << std::endl;

	// title("Simpson's rule");
	// integral = simpsons(f_1, a, b);
	// std::cout << "integral = " << integral << std::endl;

	// use DErule
	title("DErule");
	DErule<double(Doub, Doub)> de(f_1, a, b);
	print_de_rule(de);

	/* ------------------------------- function 2 ------------------------------- */

	title("\nIntegrate from a = 0, b = 1: f(x) = sqrt(x)*cos(x^2)*exp(-x)");
	// // use simpsons
	// title("Simpson's rule");
	// integral = simpsons(f_2, a, b);
	// std::cout << "integral = " << integral << std::endl;

	// use DErule
	title("DErule");
	DErule<double(Doub, Doub)> de2(f_2, a, b);
	print_de_rule(de2);

	// /* ------------------------------- function 3 ------------------------------- */
	title("\nIntegrate from a = 0, b = 1: f(x) = 1/sqrt(x)*cos(x^2)*exp(-x)");
	// // midpoint
	// title("Extended midpoint");
	// integral = extended_midpoint(f_3, a, b);
	// std::cout << "integral = " << integral << std::endl;

	// use DErule
	title("DErule");
	DErule<double(Doub, Doub)> de3(f_3, a, b);
	print_de_rule(de3);

	// /* ------------------------------- function 4 ------------------------------- */
	title("\nIntegrate from a = 0, b = 1: f(x) = 1000*exp(-1/x)*exp(-1/(1-x))");
	// // trapezoidal
	// title("Trapezoidal rule");
	// integral = trapz(f_4, a, b);
	// std::cout << "integral = " << integral << std::endl;

	// use DErule
	title("DErule");
	DErule<double(Doub, Doub)> de4(f_4, a, b);
	print_de_rule(de4);

	return 0;
}