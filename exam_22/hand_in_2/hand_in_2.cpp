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

double n_var;

template <class T>
void newto(VecDoub_IO &x, Bool &check, T &vecfunc)
{
    int show_its = 6;
    const Int MAXITS = 2000;
    const Doub TOLF = 1.0e-8, TOLMIN = 1.0e-12, STPMX = 100.0;
    const Doub TOLX = numeric_limits<Doub>::epsilon();
    Int i, j, its, n = x.size();
    Doub den, f, fold, stpmax, sum, temp, test;
    VecDoub g(n), p(n), xold(n);
    MatDoub fjac(n, n);
    NRfmin<T> fmin(vecfunc);
    NRfdjac<T> fdjac(vecfunc);
    VecDoub &fvec = fmin.fvec;
    f = fmin(x);
    test = 0.0;
    // print header
    //  x consists of {L_0, L, p, x, theta, phi, a, H}
    std::cout << YELLOW << "i" << std::setw(15) << "L_0" << std::setw(15) << "L" << std::setw(15) << "p" << std::setw(15) << "x" << std::setw(15) << "theta" << std::setw(15) << "phi" << std::setw(15) << "a" << std::setw(15) << "H" << RESET << std::endl;
    // alternatively if you dont know the variable names
    //  std::cout << YELLOW << "i" << std::setw(15) << "x[0]" << std::setw(15) << "x[1]" << std::setw(15) << "x[2]" << std::setw(15) << "x[3]" << std::setw(15) << "x[4]" << std::setw(15) << "x[5]" << std::setw(15) << "x[6]" << std::setw(15) << "x[7]" << RESET << std::endl;

    for (i = 0; i < n; i++)
        if (abs(fvec[i]) > test)
            test = abs(fvec[i]);
    if (test < 0.01 * TOLF)
    {
        check = false;
        return;
    }
    sum = 0.0;
    for (i = 0; i < n; i++)
        sum += SQR(x[i]);
    stpmax = STPMX * MAX(sqrt(sum), Doub(n));
    for (its = 0; its < MAXITS; its++)
    {
        fjac = fdjac(x, fvec);
        for (i = 0; i < n; i++)
        {
            sum = 0.0;
            for (j = 0; j < n; j++)
                sum += fjac[j][i] * fvec[j];
            g[i] = sum;
        }
        if (its < show_its)
        {
            int markus; // std::cout << its << setw(15) << x[0] << setw(15) << x[1] << setw(15) << x[2] << setw(15) << x[3] << setw(15) << x[4] << setw(15) << x[5] << setw(15) << x[6] << setw(15) << x[7] << std::endl;
        }
        for (i = 0; i < n; i++)
            xold[i] = x[i];
        fold = f;
        for (i = 0; i < n; i++)
            p[i] = -fvec[i];
        LUdcmp alu(fjac);
        alu.solve(p, p);
        lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
        test = 0.0;
        for (i = 0; i < n; i++)
            if (abs(fvec[i]) > test)
                test = abs(fvec[i]);
        if (test < TOLF)
        {
            check = false;
            return;
        }
        if (check)
        {
            test = 0.0;
            den = MAX(f, 0.5 * n);
            for (i = 0; i < n; i++)
            {
                temp = abs(g[i]) * MAX(abs(x[i]), 1.0) / den;
                if (temp > test)
                    test = temp;
            }
            check = (test < TOLMIN);
            return;
        }
        test = 0.0;
        for (i = 0; i < n; i++)
        {
            temp = (abs(x[i] - xold[i])) / MAX(abs(x[i]), 1.0);
            if (temp > test)
                test = temp;
        }
        if (test < TOLX)
            return;
    }
    throw("MAXITS exceeded in newt");
}

VecDoub vecfunc(VecDoub_I x)
{
    // d is held constant
    const double d = 30;

    //  x consists of {L_0, L, p, x, theta, phi, a, H}
    double L_0 = x[0], L = x[1], p = x[2], x_ = x[3], theta = x[4], phi = x[5], a = x[6], H = x[7];

    // material constants
    const double v = 120, k = 2.5, w = 4.0, alpha = 2 * 10e-7;

    VecDoub result(8);
    result[0] = a * (cosh(x_ / a) - 1) - p;
    result[1] = 2 * a * sinh(x_ / a) - L;
    result[2] = 2 * x_ + 2 * k * cos(theta) - d;
    result[3] = p + k * sin(theta) - n_var;
    result[4] = sinh(x_ / a) - tan(phi);
    result[5] = (1 + v / (w * L_0)) * tan(phi) - tan(theta);
    result[6] = L_0 * (1 + alpha * H) - L;
    result[7] = (w * L_0) / (2 * sin(phi)) - H;

    return result;
}

int main()
{
    // matrix for storing L_0 and H for each n
    MatDoub L_0_H(6, 2);

    // initial guess;
    // x consists of {L_0, L, p, x, theta, phi, a, H}
    VecDoub_IO x_guess(8);
    x_guess[0] = 24;
    x_guess[1] = 24;
    x_guess[2] = 1;
    x_guess[3] = 12;
    x_guess[4] = 0.1;
    x_guess[5] = 0.1;
    x_guess[6] = 100;
    x_guess[7] = 100;

    std::vector<double> n_vec = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0};

    bool check = true;
    /* --------------------------------- newtons method -------------------------------- */
    for (size_t i = 0; i < n_vec.size(); i++)
    {
        n_var = n_vec[i];
        cout << GREEN << "n = " << n_vec[i] << RESET << endl;
        VecDoub_IO x = x_guess;
        newto(x, check, vecfunc);

        util::print(x, "solution_x");
        L_0_H[i][0] = x[0];
        L_0_H[i][1] = x[7];
    }

    // loop through the matrix and print the values with corresponding n
    cout << YELLOW << setw(5) << "n" << setw(15) << "L_0" << setw(15) << "H" << RESET << endl;
    for (size_t i = 0; i < L_0_H.nrows(); i++)
    {
        cout << setw(5) << n_vec[i] << setw(15) << L_0_H[i][0] << setw(15) << L_0_H[i][1] << endl;
    }

    return 0;
}