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

Doub F(const Doub x, const Doub t)
{
	return x * (1 - x) * cos(t) * exp(-t / 10);
}

Doub a(const Doub t)
{
	return 0;
}

Doub b(const Doub t)
{
	return 1;
}

Doub g(const Doub x)
{
	return pow(x, 4);
}

Doub dudt(const Doub dududxdx, const Doub alpha, const Doub x, const Doub t)
{
	return alpha * dududxdx + F(x, t);
}

void calc_rhs(VecDoub_IO &rhs, VecDoub_I u_current, const Doub r, const Doub t, const Doub dt)
{
	rhs[0] = 1 / 2 * r * (a(t) + a(t + dt)) + (1 - r) * u_current[0] + 1 / 2 * r * u_current[1] + 1 / 2 * dt * (F(u_current[0], t) + F(u_current[0], t + dt));

	for (size_t i = 1; i < rhs.size() - 1; i++)
	{
		rhs[i] = 1 / 2 * r * (u_current[i - 1] + u_current[i + 1]) + (1 - r) * u_current[i] + 1 / 2 * dt * (F(u_current[i], t) + F(u_current[i], t + dt));
	}

	rhs[rhs.size() - 1] = 1 / 2 * r * (b(t) + b(t + dt)) + (1 - r) * u_current[u_current.size() - 1] + 1 / 2 * r * u_current[u_current.size() - 2] + 1 / 2 * dt * (F(u_current[u_current.size() - 1], t) + F(u_current[u_current.size() - 1], t + dt));
}

void crank_nicolson(VecDoub_IO &u, const Doub alpha, const Doub N, const Doub t_des, const Doub k, const Doub x_start, const Doub x_end)
{
	const Doub dx = (x_end - x_start) / N;
	const Doub dt = k * dx;
	Doub t = 0;

	Doub r = alpha * dt / (dx * dx);

	VecDoub RHS(u.size());

	// calculating coefficient matrix
	VecDoub a(u.size(), -r / 2);
	a[0] = 0;
	VecDoub b(u.size(), 1 + r);
	VecDoub c(u.size(), -r / 2);
	c[c.size() - 1] = 0;

	// calculating right hand side
	while (t < t_des)
	{
		calc_rhs(RHS, u, r, t, dt);

		tridag(a, b, c, RHS, u);

		t += dt;
		cout << "t: " << t << endl;
		cout << "u(0.5): " << u[u.size() / 2] << endl;
	}
}

int main()
{
	/* -------------------- two point boundary value problem -------------------- */

	util::title("Crank Nicholson Method");

	const Doub alpha = 1;

	const Doub N = 10;
	const Doub x_start = 0;
	const Doub x_end = 1;
	const Doub k = 0.5;
	const Doub t_des = 1;

	// initial condition
	VecDoub u(N - 1);
	for (size_t i = 0; i < u.size(); i++)
	{
		u[i] = g(x_start + (i + 1) * (x_end - x_start) / N);
	}

	crank_nicolson(u, alpha, N, t_des, k, x_start, x_end);
	return 0;
}