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

void v_diff_eqs(Doub x, VecDoub_I &y, VecDoub_O &dydx)
{
	dydx[0] = y[0] * y[1];
	dydx[1] = -pow(y[0], 3);
}

VecDoub diff_eqs(Doub x, VecDoub_I &y)
{
	VecDoub dydx(2);
	dydx[0] = y[0] * y[1];
	dydx[1] = -pow(y[0], 3);
	return dydx;
}

struct newt_struct
{
	VecDoub y_init;
	Doub x;
	Doub h;
	VecDoub trapz_eqs(const VecDoub_IO y)
	{
		VecDoub res(2);
		res[0] = y[0] - y_init[0] - 0.5 * h * (diff_eqs(x, y_init)[0] + diff_eqs(x + h, y)[0]);
		res[1] = y[1] - y_init[1] - 0.5 * h * (diff_eqs(x, y_init)[1] + diff_eqs(x + h, y)[1]);
		return res;
	}
};

template <class T>
void euler(VecDoub_IO &y, const double x_start, const double x_des, T &vecfunc, const double tol = 1e-10)
{
	VecDoub y_start = y;
	double h;
	double x;
	double n = 5;
	VecDoub y_m1(2, INFINITY);
	VecDoub y_m2(2, INFINITY);
	VecDoub dxdy(2);

	VecDoub alp_k(2, INFINITY);
	VecDoub error(2, INFINITY);

	int max_iter = 20;
	int i = 1;

	// print header
	cout << setw(5) << "n" << setw(15) << "u" << setw(15) << "u - u_prev" << setw(15) << "alp_k" << setw(15) << "error" << setw(15) << "v" << setw(15) << "v - v_prev" << setw(15) << "alp_k" << setw(15) << "error" << setw(15) << "u^2+v^2" << endl;
	while (i < max_iter && (error[0] > tol || error[1] > tol))
	{
		x = x_start;
		y = y_start;
		h = x_des / n;
		while (x < x_des)
		{
			// compute slopes
			vecfunc(x, y, dxdy);
			// update y
			y = y + dxdy * h;
			x += h;
		}
		if (i > 3) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (y_m2 - y_m1) / (y_m1 - y);
			error = (y - y_m1) / ((2 ^ alp_k) - 1.0);
		}
		cout << setw(5) << n << setw(15) << y[0] << setw(15) << y[0] - y_m1[0] << setw(15) << alp_k[0] << setw(15) << error[0]
			 << setw(15) << y[1] << setw(15) << y[1] - y_m1[1] << setw(15) << alp_k[1] << setw(15) << error[1] << setw(15) << y[0] * y[0] + y[1] * y[1] << endl;
		y_m2 = y_m1;
		y_m1 = y;
		n *= 2;
		i++;

		// if error is inf or nan set to 99
		if (isnan(error[0]) || isnan(error[1]))
		{
			error[0] = INFINITY;
			error[1] = INFINITY;
		}
	}
}

template <class T>
void print_rk4(VecDoub_IO &y, const Doub x_start, const Doub x_des, T &vecfunc, const Doub tol = 1e-10)
{
	VecDoub y_start = y;
	Doub h;
	Doub x;
	Doub n = 5;
	VecDoub dxdy(2);
	VecDoub y_out;
	VecDoub_O dydx(2);
	VecDoub y_m1(2, INFINITY);
	VecDoub y_m2(2, INFINITY);
	VecDoub alp_k(2, INFINITY);
	VecDoub error(2, INFINITY);

	int max_iter = 10;
	int i = 1;

	// print header
	cout << setw(5) << "n" << setw(15) << "u" << setw(15) << "u - u_prev" << setw(15) << "alp_k" << setw(15) << "error" << setw(15) << "v" << setw(15) << "v - v_prev" << setw(15) << "alp_k" << setw(15) << "error" << endl;

	while (i < max_iter)
	{
		x = x_start;
		y = y_start;
		h = x_des / n;
		while (x < x_des)
		{
			rk4(y, dydx, x, h, y, vecfunc);
			x += h;
		}
		if (i > 2) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (y_m2 - y_m1) / (y_m1 - y);
			error = (y - y_m1) / ((2 ^ alp_k) - 1.0);
		}
		cout << setw(5) << n << setw(15) << y[0] << setw(15) << y[0] - y_m1[0] << setw(15) << alp_k[0] << setw(15) << error[0] << setw(15) << y[1] << setw(15) << y[1] - y_m1[1] << setw(15) << alp_k[1] << setw(15) << error[1] << endl;
		y_m2 = y_m1;
		y_m1 = y;

		n *= 2;
		i++;
	}
}

template <class T>
void leap_frog(VecDoub_IO &y, const double x_start, const double x_des, T &vecfunc, const double tol = 1e-10)
{
	VecDoub y_start = y;
	double h;
	double x;
	double n = 5;
	VecDoub y_m1(2, INFINITY);
	VecDoub y_m2(2, INFINITY);
	VecDoub dxdy(2);

	VecDoub alp_k(2, INFINITY);
	VecDoub error(2, INFINITY);

	int max_iter = 20;
	int i = 1;

	// print header
	cout << setw(5) << "n" << setw(15) << "u" << setw(15) << "u - u_prev" << setw(15) << "alp_k" << setw(15) << "error" << setw(15)
		 << "v" << setw(15) << "v - v_prev" << setw(15) << "alp_k" << setw(15) << "error" << setw(15) << "u^2+v^2" << endl;
	while (i < max_iter && (error[0] > tol || error[1] > tol))
	{
		x = x_start;
		y = y_start;
		h = x_des / n;
		while (x < x_des)
		{
			// compute slopes
			vecfunc(x, y, dxdy);
			// update y
			y[0] += h * dxdy[0];
			y[1] += h * dxdy[1];
			x += h;
		}
		if (i > 3) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (y_m2 - y_m1) / (y_m1 - y);
			error = (y - y_m1) / ((2 ^ alp_k) - 1.0);
		}
		cout << setw(5) << n << setw(15) << y[0] << setw(15) << y[0] - y_m1[0] << setw(15) << alp_k[0] << setw(15) << error[0]
			 << setw(15) << y[1] << setw(15) << y[0] - y_m1[1] << setw(15) << alp_k[1] << setw(15) << error[1] << setw(15) << y[0] * y[0] + y[1] * y[1] << endl;
		y_m2 = y_m1;
		y_m1 = y;
		n *= 2;
		i++;

		// if error is nan set to 99
		if (isnan(error[0]) || isnan(error[1]))
		{
			error[0] = INFINITY;
			error[1] = INFINITY;
		}
	}
}

template <class T>
void trapezoidal(VecDoub_IO &y, const double x_start, const double x_des, T &vecfunc, const double tol = 1e-10)
{
	VecDoub y_start = y;
	double h;
	double x;
	double n = 5;
	VecDoub y_m1(2, INFINITY);
	VecDoub y_m2(2, INFINITY);
	VecDoub dxdy(2);

	VecDoub alp_k(2, INFINITY);
	VecDoub error(2, INFINITY);

	int max_iter = 20;
	int i = 1;

	// print header
	cout << setw(5) << "n" << setw(15) << "u" << setw(15) << "u - u_prev" << setw(15) << "alp_k" << setw(15) << "error" << setw(15)
		 << "v" << setw(15) << "v - v_prev" << setw(15) << "alp_k" << setw(15) << "error" << setw(15) << "u^2+v^2" << endl;
	while (i < max_iter || (error[0] > tol || error[1] > tol))
	{
		newt_struct my_newt_struct; // y_init, current x, stepsize h

		my_newt_struct.x = x_start;
		my_newt_struct.h = x_des / n;

		// perform initial euler step
		vecfunc(my_newt_struct.x, y, dxdy);
		y = y + dxdy[0] + h;

		my_newt_struct.y_init = y;
		bool check = true;
		y = y_start;
		while (my_newt_struct.x < x_des)
		{

			auto bound_newt_eqs = bind(&newt_struct::trapz_eqs, my_newt_struct, placeholders::_1);

			newt(y, check, bound_newt_eqs); // middle parameter is just for fun
			my_newt_struct.x += my_newt_struct.h;
			my_newt_struct.y_init = y;
		}
		if (i > 3) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (y_m2 - y_m1) / (y_m1 - y);
			error = (y - y_m1) / ((2 ^ alp_k) - 1.0);
		}
		cout << setw(5) << n << setw(15) << y[0] << setw(15) << y[0] - y_m1[0] << setw(15) << alp_k[0] << setw(15) << error[0]
			 << setw(15) << y[1] << setw(15) << y[0] - y_m1[1] << setw(15) << alp_k[1] << setw(15) << error[1] << setw(15) << y[0] * y[0] + y[1] * y[1] << endl;
		y_m2 = y_m1;
		y_m1 = y;
		n *= 2;
		i++;

		// if error is nan set to 99
		if (isnan(error[0]) || isnan(error[1]))
		{
			error[0] = INFINITY;
			error[1] = INFINITY;
		}
	}
}

int main()
{
	double u_0 = 1, v_0 = 1;

	VecDoub_IO y(2);
	y[0] = u_0;
	y[1] = v_0;
	/* ------------------------------- function 1 ------------------------------- */

	title("Differential equation");

	// title("Euler method");
	// euler(y, 0, 10, v_diff_eqs, 1e-3);

	// title("runge kutta 4th order");
	// y[0] = u_0;
	// y[1] = v_0;
	// print_rk4(y, 0, 10, v_diff_eqs, 1e-3);

	// title("Leap-frog method");
	// y[0] = u_0;
	// y[1] = v_0;
	// leap_frog(y, 0, 10, v_diff_eqs, 1e-3);

	title("Trapezoidal method");
	y[0] = u_0;
	y[1] = v_0;
	trapezoidal(y, 0, 10, v_diff_eqs, 1e-90);

	return 0;
}