#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "../source_code/Numerical-Recipes-master/qrdcmp.h"
#include "../source_code/Numerical-Recipes-master/cholesky.h"
#include "../source_code/Numerical-Recipes-master/svd.h"			// Singular Value Decomposition
#include "../source_code/Numerical-Recipes-master/roots.h"			// Root finding ie newt
#include "../source_code/Numerical-Recipes-master/roots_multidim.h" // Root finding ie newt also newt_mine
#include "../source_code/Numerical-Recipes-master/quadrature.h"
#include "../source_code/Numerical-Recipes-master/derule.h"
#include "../source_code/Numerical-Recipes-master/rk4.h" // Runge-Kutta 4th order method
#include <functional>

#include "../utility_code/utilities.h"

#include <fstream>

#define THRESHOLD 1.0e-10

using namespace std;

VecDoub diff_eqs(Doub x, VecDoub_I y)
{
	VecDoub dydx(3);
	dydx[0] = exp(-x) * cos(y[1]) + pow(y[2], 2) - y[0];
	dydx[1] = cos(pow(y[2], 2)) - y[1];
	dydx[2] = cos(x) * exp(-pow(y[0], 2)) - y[2];
	return dydx;
}

struct newt_struct
{
	VecDoub y_init;
	Doub x;
	Doub h;
	VecDoub trapz_eqs(VecDoub_I y)
	{
		return y - y_init - 0.5 * h * (diff_eqs(x, y_init) + diff_eqs(x + h, y));
	}
};

void print_header_trapz()
{
	cout << setw(5) << "n" << setw(15) << "v1" << setw(15) << "v2" << setw(15) << "v3" << setw(15) << "error1" << setw(15) << "error2" << setw(15) << "error3" << endl;
}

void trapezoidal(VecDoub_IO &y, double x_start, const double x_des, const double h)
{
	newt_struct my_newt_struct;
	my_newt_struct.y_init = y;
	my_newt_struct.x = x_start;
	my_newt_struct.h = h;
	VecDoub y_star(3);
	VecDoub dxdy(3);
	bool check = true;
	while (my_newt_struct.x < x_des)
	{
		auto bound_newt_eqs = bind(&newt_struct::trapz_eqs, my_newt_struct, placeholders::_1);
		// euler step
		dxdy = diff_eqs(my_newt_struct.x, my_newt_struct.y_init);
		y_star = my_newt_struct.y_init + dxdy * h;

		// trapezoidal step
		newt(y_star, check, bound_newt_eqs);
		my_newt_struct.x += h;
		my_newt_struct.y_init = y_star;
		y = y_star;
	}
}

// this is made for 3 coupled ODEs
void trapezoidal_print(VecDoub_I y_start, const double x_start, const double x_des, const double tol = 1e-10)
{
	VecDoub y(3, INFINITY);
	double h;
	double n = 50;
	int max_iter = 5;
	int i = 1;

	// for calculating alp_k and error
	VecDoub y_m1(3, INFINITY);
	VecDoub y_m2(3, INFINITY);
	VecDoub alp_k(3, INFINITY);
	VecDoub error(3, INFINITY);

	// print header
	print_header_trapz();
	while (i <= max_iter && (error[0] > tol || error[1] > tol || error[2] || tol))
	{
		y = y_start;
		h = (x_des - x_start) / n;
		trapezoidal(y, x_start, x_des, h);

		if (i > 3) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (y_m2 - y_m1) / (y_m1 - y);
			error = (y - y_m1) / (alp_k - 1.0);
		}
		cout << setw(5) << n << setw(15) << y[0] << setw(15)
			 << setw(15) << y[1] << setw(15) << y[2] << setw(15) << error[0] << setw(15) << error[1] << setw(15) << error[2] << endl;
		y_m2 = y_m1;
		y_m1 = y;
		n *= 2;
		i++;

		// if error is nan set to infinity
		nan_to_inf(error);
	}
}

int main()
{
	/* ------------------------------- exercise i ------------------------------- */

	util::title("exercise i");
	// start conditions
	double v1_0 = 1, v2_0 = 2, v3_0 = 3;
	VecDoub y(3);
	y[0] = v1_0;
	y[1] = v2_0;
	y[2] = v3_0;

	// get values of primes at start conditions
	double x = 0;
	VecDoub dxdy(3);
	dxdy = diff_eqs(x, y);
	cout << "dxdy: " << dxdy[0] << " " << dxdy[1] << " " << dxdy[2] << endl;

	/* ------------------------------- exercise ii ------------------------------ */

	util::title("exercise ii");

	double x_start = 0;
	double x_des = 5;
	double tol = 1e-40;
	trapezoidal_print(y, x_start, x_des, tol);

	/* int N = 50;
	double h = (x_des - x_start) / N;
	trapezoidal(y, x_start, x_des, h); */
	cout << "y: " << y[0] << " " << y[1] << " " << y[2] << endl;

	/* ------------------------------ exercise iii ------------------------------ */

	// util::title("exercise iii");

	return 0;
}