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
#include "../source_code/Numerical-Recipes-master/tridag.h"

#include "../utility_code/utilities.h"

#include <fstream>

#define THRESHOLD 1.0e-10

using namespace std;

struct start_cond
{
	Doub a;
	Doub b;
	Doub alpha;
	Doub beta;
	Doub N;
	Doub h = (b - a) / N;
};

Doub F(const Doub y_d, const Doub y, const Doub x)
{
	return 2 * x + sin(y_d) - cos(y);
}

Doub F_dy(const Doub y_d, const Doub y, const Doub x)
{
	return sin(y);
}

Doub F_ddy(const Doub y_d, const Doub y, const Doub x)
{
	return cos(y_d);
}

void calcJacobian(VecDoub_IO &a, VecDoub_IO &b, VecDoub_IO &c, const Doub h, const Doub alpha, const Doub beta, const VecDoub_I y, const VecDoub_I x)
{
	a[0] = 0;
	b[0] = 2 + h * h * F_dy((y[2] - alpha) / (2 * h), y[1], x[1]);
	c[0] = -1 + h / 2 * F_ddy((y[2] - alpha) / (2 * h), y[1], x[1]);

	for (size_t i = 1; i < b.size() - 1; i++)
	{
		a[i] = -1 - h / 2 * F_ddy((y[i + 1] - y[i - 1]) / (2 * h), y[i], x[i]);
		b[i] = 2 + h * h * F_dy((y[i + 1] - y[i - 1]) / (2 * h), y[i], x[i]);
		c[i] = -1 + h / 2 * F_ddy((y[i + 1] - y[i - 1]) / (2 * h), y[i], x[i]);
	}

	a[a.size() - 1] = -1 - h / 2 * F_ddy((beta - y[y.size() - 1]) / (2 * h), y[a.size()], x[x.size()]);
	b[b.size() - 1] = 2 + h * h * F_dy((beta - y[y.size() - 1]) / (2 * h), y[y.size()], x[x.size()]);
	c[c.size() - 1] = 0;
}

void calcPhi(VecDoub_IO &phi, const Doub alpha, const Doub beta, const Doub h, VecDoub_I y, VecDoub_I x)
{
	phi[0] = -alpha + 2 * y[0] - y[1] + h * h * F((y[1] - alpha) / (2 * h), y[0], x[0]);
	for (size_t i = 1; i < phi.size() - 1; i++)
	{
		phi[i] = -y[i - 1] + 2 * y[i] - y[i + 1] + h * h * F((y[i + 1] - y[i - 1]) / (2 * h), y[i], x[i]);
	}
	phi[phi.size() - 1] = -y[y.size() - 2] + 2 * y[y.size() - 1] - beta + h * h * F((beta - y[y.size() - 2]) / (2 * h), y[y.size() - 1], x[x.size() - 1]);
}

// Finite Difference Method
void FDM_1it(VecDoub_IO &y, const Doub a, const Doub b, const Doub alpha, const Doub beta, const int N, const Doub tol = 1.0e-4)
{
	int N_points = N - 1;
	VecDoub a_vec(N_points);
	VecDoub b_vec(N_points);
	VecDoub c_vec(N_points);
	VecDoub phi_vec(N_points);
	VecDoub y_delta(N_points);

	Doub h = (b - a) / N;
	// cout << "h: " << h << endl;

	// linearly interpolate
	VecDoub x(N_points);
	for (size_t i = 0; i < x.size(); i++)
	{
		x[i] = a + (i + 1) * h;
	}
	// util::print(x, "x");

	for (size_t i = 0; i < 3; i++)
	{
		// calculate Jacobian aka. tridiagonal matrix or coefficient matrix
		calcJacobian(a_vec, b_vec, c_vec, h, alpha, beta, y, x);

		// calculate phi aka. right hand side
		calcPhi(phi_vec, alpha, beta, h, y, x);
		// util::print(phi_vec, "phi");

		// solve the system of equations
		tridag(a_vec, b_vec, c_vec, -1 * phi_vec, y_delta);
		// perform newton step
		y = y + y_delta;

		// util::print(y, "y");
	}
}

// double vector size and linearly interpolate between previous points
VecDoub doubleVectorLinearize(VecDoub_I input_vec, const Doub alpha, const Doub beta)
{
	VecDoub output_vec(2 * (input_vec.size() + 1) - 1);
	output_vec[0] = (input_vec[0] + alpha) / 2;
	for (size_t i = 0; i < input_vec.size(); i++)
	{
		output_vec[2 * i + 1] = input_vec[i];
		output_vec[2 * i + 2] = (input_vec[i] + input_vec[i + 1]) / 2;
	}
	output_vec[output_vec.size() - 1] = (input_vec[input_vec.size() - 1] + beta) / 2;
	return output_vec;
}

void FDM(const Doub a, const Doub b, const Doub alpha, const Doub beta, const int N_start, const Doub tol = 1.0e-4)
{
	int N = N_start; // y(1) = y[N/2]
	int N_points = N - 1;
	Doub error = 99.0;
	Doub y1_m1 = 0.0;
	Doub y1_m2 = 0.0;
	Doub alp_k = 0.0;

	VecDoub y(N_points, 0.0);
	// linearly interpolate
	for (size_t i = 0; i < y.size(); i++)
	{
		y[i] = alpha + (i + 1) * (beta - alpha) / N_points;
	}

	VecDoub y_new(2 * (y.size() + 1) - 1);

	cout << setw(5) << "N" << setw(15) << "y(1)" << setw(15) << "y(1)-y_prev(1)" << setw(15) << "alp_k" << setw(15) << "error" << setw(22) << "order estimate" << endl;
	for (size_t i = 0; i < 10; i++)
	{
		// util::print(y, "y_before");
		FDM_1it(y, a, b, alpha, beta, N, tol);
		// util::print(y, "y_after");
		Doub y1 = y[N_points / 2];
		if (i > 2) // we need at least 3 iterations to calculate alp_k
		{
			alp_k = (y1_m2 - y1_m1) / (y1_m1 - y1);
			error = (y1 - y1_m1) / (pow(2, alp_k) - 1.0);
		}

		cout << setw(5) << N << setw(15) << y1 << setw(15) << y1 - y1_m1 << setw(15) << alp_k << setw(15) << error << setw(15) << alp_k << endl;
		y1_m2 = y1_m1;
		y1_m1 = y1;

		N *= 2;
		N_points = N - 1;
		// linearly interpolate
		y_new = doubleVectorLinearize(y, alpha, beta);
		y.resize(y_new.size());
		y = y_new;
		y_new.resize(2 * (y.size() + 1) - 1);
	}
}

int main()
{
	/* -------------------- two point boundary value problem -------------------- */

	util::title("finite difference method for two point boundary value problem");

	// a is beginning of interval, b is end of interval in terms of x
	double a = 0.0;
	double b = 2.0;

	// starting guess for y(a) and y(b)
	double alpha = 0.0;
	double beta = 1.0;

	// number of points
	int N = 2;

	FDM(a, b, alpha, beta, N);

	return 0;
}