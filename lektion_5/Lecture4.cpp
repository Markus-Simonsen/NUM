#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../source_code/Numerical-Recipes-master/svd.h"
#include "../source_code/Numerical-Recipes-master/roots.h"
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

double funk(double input)
{
	return input - cos(input);
}

double funk_prime(double input)
{
	return 1 + sin(input);
}

int sign(double x)
{
	if (x > 0)
	{
		return 1;
	}
	else if (x < 0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
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

double bisection(double x_l, double x_h, double tolerance, double (*f)(double))
{

	int alpha = 1;
	// Check for bracketing
	if (f(x_l) * f(x_h) > 0)
	{
		std::cout << "Error: No bracketing" << std::endl;
		return 0;
	}

	// convert to minus and positive values
	double x_m = 0;
	double x_p = 0;

	if (f(x_l) < 0)
	{
		x_m = x_l;
		x_p = x_h;
	}
	else
	{
		x_m = x_h;
		x_p = x_l;
	}

	double x_i = 0;
	double x_im1 = 0;
	double x_im2 = 0;
	int i = 0;

	print_header();
	while (fabs(x_m - x_p) > tolerance)
	{
		x_i = (x_m + x_p) / 2;
		if (f(x_i) < 0)
		{
			x_m = x_i;
		}
		else
		{
			x_p = x_i;
		}
		print_table(i, x_i, x_im1, x_im2, alpha, "bisection");

		i++;
		x_im2 = x_im1;
		x_im1 = x_i;
	}
	return (x_m + x_p) / 2;
}

double ridders(double x_l, double x_h, double tolerance, double (*f)(double))
{
	// Check for bracketing
	if (f(x_l) * f(x_h) > 0)
	{
		std::cout << "Error: No bracketing" << std::endl;
		return 0;
	}

	int alpha = 3;
	double x_m = 0;
	double x_p = 0;
	double z_i = 0; // middle point

	if (f(x_l) < 0)
	{
		x_m = x_l;
		x_p = x_h;
	}
	else
	{
		x_m = x_h;
		x_p = x_l;
	}

	double x_i = 0;
	double x_im1 = 0;
	double x_im2 = 0;
	int i = 0;
	while (fabs(x_m - x_p) > tolerance)
	{
		z_i = (x_m + x_p) / 2;
		x_i = z_i + (z_i - x_p) * f(z_i) * sign(f(x_p) - f(x_m)) / sqrt(pow(f(z_i), 2) - f(x_m) * f(x_p));
		if (f(x_i) < 0)
		{
			x_m = x_i;
		}
		else
		{
			x_p = x_i;
		}
		print_table(i, x_i, x_im1, x_im2, alpha, "ridders");

		i++;
		x_im2 = x_im1;
		x_im1 = x_i;
	}
	return (x_m + x_p) / 2;
}

double secant(double x_l, double x_h, double tolerance, double (*f)(double))
{

	int alpha = 1.62;
	double x_i = x_h;
	double x_im1 = x_l;
	double x_im2 = 0;
	int i = 0;

	print_header();
	while (fabs(x_i - x_im1) > tolerance)
	{
		double f_i = f(x_i);
		double f_im1 = f(x_im1);
		double x_ip1 = x_i - f_i * (x_i - x_im1) / (f_i - f_im1);

		x_im2 = x_im1;
		x_im1 = x_i;
		x_i = x_ip1;

		print_table(i, x_i, x_im1, x_im2, alpha, "secant");

		i++;
	}

	return x_i;
}

double newton(double x_0, double tolerance, double (*f)(double), double (*f_prime)(double))
{
	int alpha = 2;
	double x_i = x_0;
	double x_im1 = 0;
	double x_im2 = 0;
	int i = 0;

	print_header();
	while (fabs(x_i - x_im1) > tolerance)
	{
		double f_i = f(x_i);
		double f_prime_i = f_prime(x_i);
		double x_ip1 = x_i - f_i / f_prime_i;

		x_im2 = x_im1;
		x_im1 = x_i;
		x_i = x_ip1;

		print_table(i, x_i, x_im1, x_im2, alpha, "newton");

		i++;
	}

	return x_i;
}

int main()
{
	/* ------------------------------- exercise 1 ------------------------------- */
	std::cout
		<< GREEN << "Exercise 1" << RESET << std::endl;

	/* -------------------------------- bisection ------------------------------- */
	std::cout
		<< GREEN << "Bisection" << RESET << std::endl;

	std::cout << "x_l = 0, x_h = pi/2, Tolerance = 10^-8"
			  << "\n"
			  << endl;

	double x_l = 0;
	double x_h = M_PI / 2;
	double tolerance = 1.0e-8;

	double root = bisection(x_l, x_h, tolerance, funk);
	std::cout << "Root = " << root << std::endl;

	/* -------------------------------- newton ------------------------------- */
	std::cout
		<< GREEN << "Newton" << RESET << std::endl;
	double x_0 = 0.2;
	std::cout << "x_0 = 0.2, Tolerance = 10^-8"
			  << "\n"
			  << endl;

	root = newton(x_0, tolerance, funk, funk_prime);
	std::cout << "Root = " << root << std::endl;

	/* -------------------------------- secant ------------------------------- */
	std::cout
		<< GREEN << "Secant" << RESET << std::endl;

	std::cout << "x_i = 0.2, x_im1 = 0 , Tolerance = 10^-8"
			  << "\n"
			  << endl;

	root = secant(x_l, x_h, tolerance, funk);
	std::cout << "Root = " << root << std::endl;

	/* -------------------------------- ridders ------------------------------- */
	std::cout
		<< GREEN << "Ridders" << RESET << std::endl;

	std::cout << "x_l = 0, x_h = pi/2, Tolerance = 10^-8"
			  << "\n"
			  << endl;

	print_header();
	root = ridders(x_l, x_h, tolerance, funk);
	std::cout << "Root = " << root << std::endl;

	return 0;
}