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
	return (cos(pow(x, 3)) * exp(-x)) / sqrt(x);
}

double f_1(double x, double delta)
{
	if (x < 0.1)
		return (cos(pow(delta, 3)) * exp(-delta)) / sqrt(delta);
	else
		return (cos(pow(x, 3)) * exp(-x)) / sqrt(x);
}

double extended_midpoint(double f(double x), double a, double b, double threshold = 1e-3)
{

	double alp_k = 0;
	double error = 100;
	double h_old = 0;
	double integral;
	double integral_m1 = 0;
	double integral_m2 = 0;
	int N = 1;
	int f_eval;
	int i = 1;
	int max_iter = 40;
	std::cout << setw(5) << "i" << setw(15) << "A(hi)" << setw(15) << "A(hi-1)-A(hi)" << setw(15) << "alp^k" << setw(15) << "rich-comp" << setw(15) << "f-comp" << endl;
	while (abs(error) > threshold && i < max_iter)
	{
		double h = (b - a) / N;
		integral = 0;
		for (double j = a + 0.5 * h; j < b; j += h)
		{
			integral += f(j) * h;
		}
		if (i > 2) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (integral_m2 - integral_m1) / (integral_m1 - integral);
			error = (integral - integral_m1) / (alp_k - 1);
		}
		f_eval = pow(2, i - 1);
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

void print_de_rule(DErule<double(Doub, Doub)> de, double threshold = 1e-3)
{
	int its = 60;
	double A, A_m1 = 0, A_m2 = 0, alp_k = 0, error = 0, f_eval = 1;

	std::cout << setw(5) << "i" << setw(15) << "A(hi)" << setw(15) << "A(hi-1)-A(hi)" << setw(15) << "alp^k" << setw(15) << "error" << setw(15) << "f-comp" << endl;
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
		if (abs(A_m1 - A) < threshold)
			break;
		A_m2 = A_m1;
		A_m1 = A;
	}
}

int main()
{
	double a = 0, b = 4;

	/* ------------------------------- function 1 ------------------------------- */

	title("Integrate from a = 0, b = 1: f(x) = cos(x^2)*exp(-x)");

	double integral = extended_midpoint(f_1, a, b, 1e-3);

	// use DErule
	title("DErule");
	DErule<double(Doub, Doub)> de(f_1, a, b);
	print_de_rule(de);

	return 0;
}