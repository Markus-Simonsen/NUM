#include "nr3.h"
#include "utilities.h"
#include "helper_functions.h"
#include "cholesky.h"
#include <cmath>

Doub F4(Doub yd, Doub y, Doub x)
{
    Doub x1 = pow(y, 3) + 2 * pow(x, 3) * yd;
    Doub x2 = 2 * pow(x, 2) - pow(y, 2) * pow(yd, 2);
    Doub x3 = 1 + 64 * pow(x, 6) + 16 * pow(y, 6);
    return 48 * x1 * x2 / x3;
}

Doub F4y(Doub yd, Doub y, Doub x)
{
    Doub x1  = 16 * pow(y,10) * pow(yd, 2);
    Doub x2  = 96 * pow(y, 8) * pow(x, 2);
    Doub x3  = 128 * pow(y, 7) * pow(x, 3) * pow(yd, 3);
    Doub x4  = 384 * pow(y, 5) * pow(x, 5) * yd;
    Doub x5  = 5 * pow(y, 4) * pow(yd, 2);
    Doub x6  = 320 * pow(y, 4) * pow(x, 6) * pow(yd, 2);
    Doub x7  = 384 * pow(y, 2) * pow(x, 8);
    Doub x8  = 6 * pow(y, 2) * pow(x, 2);
    Doub x9  = 4 * y * pow(x, 3) * pow(yd, 3);
    Doub x10 = 256 * y * pow(x, 9) * pow(yd, 3);
    Doub y1  = 64 * pow(x, 6);
    Doub y2  = 16 * pow(y, 6);
    return 48 * (x1 - x2 + x3 - x4 - x5 - x6 + x7 + x8 - x9 - x10) / pow(1 + y1 + y2, 2);
}

Doub F4yd(Doub yd, Doub y, Doub x)
{
    Doub x1 = 2 * pow(y, 5) * yd;
    Doub x2 = 6 * pow(y, 2) * pow(x, 3) * pow(yd, 2);
    Doub x3 = 4 * pow(x, 5);
    Doub y1  = 64 * pow(x, 6);
    Doub y2  = 16 * pow(y, 6);
    return 48 * (-x1 - x2 + x3) / (1 + y1 + y2);
}

void Exercise1()
{
    MatDoub A = util::parseFile("Ex1A.dat");
    VecDoub b = util::toVec(util::parseFile("Ex1b.dat"));

    Cholesky cho(A);
    util::printDiag(cho.el, "Diagonal elements of L");
    std::cout << std::endl;

    VecDoub x(6);
    cho.solve(b, x);
    util::print(x, "The solution x to Ax=b");
    std::cout << std::endl;

    util::print(residuals(A, x, b), "Residuals");
    std::cout << "Residual error: " << residualError(residuals(A, x, b), b) << std::endl;
}

void Exercise2()
{
    Doub a1 = 1.1, a2 = 2.1, a3 = 0.8;
    Doub b1 = 0.4, b2 = 1.3, b3 = 0.5;
    Doub derror = 1e-8;

    auto f0 = [=](Doub Ua, Doub Ub)
    {
        VecDoub rA(3);
        rA[0] = a1 * pow(cos(1 + Ua), 3);
        rA[1] = a2 * pow(Ua, 2);
        rA[2] = a3 * Ua * sin(Ua);

        VecDoub rB(3);
        rB[0] = b1 * (Ub + exp(- pow(Ub, 2)));
        rB[1] = b2 * pow(Ub, 3);
        rB[2] = b3 * cos(Ub);

        VecDoub rA_(3);
        rA_[0] = a1 * pow(cos(1 + Ua + derror), 3);
        rA_[1] = a2 * pow(Ua + derror, 2);
        rA_[2] = a3 * (Ua + derror) * sin(Ua + derror);

        VecDoub rB_(3);
        rB_[0] = b1 * ((Ub + derror) + exp(- pow(Ub + derror, 2)));
        rB_[1] = b2 * pow(Ub + derror, 3);
        rB_[2] = b3 * cos(Ub + derror);

        VecDoub rAd = (rA_ - rA) / derror;
        VecDoub rBd = (rB_ - rB) / derror;

        return util::dot(rAd, rA - rB);
    };

    auto f1 = [=](Doub Ua, Doub Ub)
    {
        VecDoub rA(3);
        rA[0] = a1 * pow(cos(1 + Ua), 3);
        rA[1] = a2 * pow(Ua, 2);
        rA[2] = a3 * Ua * sin(Ua);

        VecDoub rB(3);
        rB[0] = b1 * (Ub + exp(- pow(Ub, 2)));
        rB[1] = b2 * pow(Ub, 3);
        rB[2] = b3 * cos(Ub);

        VecDoub rA_(3);
        rA_[0] = a1 * pow(cos(1 + Ua + derror), 3);
        rA_[1] = a2 * pow(Ua + derror, 2);
        rA_[2] = a3 * (Ua + derror) * sin(Ua + derror);

        VecDoub rB_(3);
        rB_[0] = b1 * ((Ub + derror) + exp(- pow(Ub + derror, 2)));
        rB_[1] = b2 * pow(Ub + derror, 3);
        rB_[2] = b3 * cos(Ub + derror);

        VecDoub rAd = (rA_ - rA) / derror;
        VecDoub rBd = (rB_ - rB) / derror;

        return util::dot(rBd, rB - rA);
    };

    Doub Ua = 1;
    Doub Ub = -1;
    std::cout << fixed;
    std::setprecision(6);
    std::cout << "f0: " << f0(Ua, Ub) << std::endl;
    std::cout << "f1: " << f1(Ua, Ub) << std::endl;
    std::cout.unsetf(std::ios_base::floatfield);

    auto f = [=](VecDoub x)
    {
        VecDoub out(2);

        Doub Ua = x[0];
        Doub Ub = x[1];

        out[0] = f0(Ua, Ub);
        out[1] = f1(Ua, Ub);

        return out;
    };

    bool check;
    VecDoub vars(2);
    vars[0] = Ua;
    vars[1] = Ub;
    std::cout << "NEWTONS METHOD" << std::endl;
    NewtonTable(vars, check, f);
    std::cout << std::endl;
}

void Exercise3()
{
    Doub a1 = 1, b1 = 1;
    Doub a2 = 0.25, b2 = 0.25;

    VecDoub Y(2), g(2);
    Y[0] = 0, Y[1] = 0;
    g[0] = 2, g[1] = 1;

    VecDoub init(5);
    init[0] = Y[0];
    init[1] = 0;
    init[2] = Y[1];
    init[3] = 0;
    init[4] = 1;

    Doub accuracy = 1e-4;

    auto F = [=](Doub t, VecDoub x)
    {
        Doub u1 = x[0], v1 = x[1], u2 = x[2], v2 = x[3], x_ = x[4], g1 = g[0], g2 = g[1];
        VecDoub res(5);
        
        res[0] = v1;
        res[1] = a1 * (b1 * (g1 - u1) - v1) + x_;
        res[2] = v2;
        res[3] = a2 * (b2 * (g2 - u2) - v2) + (3 * x_ * (1 - pow(x_, 2)));
        res[4] = -x_;

        return res;
    };

    MidpointPlot(F, 0, 20, init);
    // Use the Octave below to plot the trajectory:
    // data = csvread("Trajectory.csv"); x = data(:,1); y = data(:,2); plot(x,y)

    Doub start = 0;
    Doub end   = 5;
    int N = 5;

    std::cout << std::endl;
    std::cout << "MIDPOINT" << std::endl;
    MidpointODE(F, start, end, init, N, accuracy);
    std::cout << std::endl;

}

void Exercise4()
{
    Doub a = -0.9, b = 0.8;
    Doub alpha = -0.85, beta = -0.9;
    Doub accuracy = 1e-6;
    std::cout << "FINITE DIFFERENCE METHOD" << std::endl;
    FiniteDiff(a, b, alpha, beta, F4, F4y, F4yd, accuracy);
    std::cout << std::endl;
}

void Exercise5()
{
    Doub a = 0, b = 4, acc = 1e-3;

    auto F = [&](Doub x)
    {
        return cos(x*x*x) * exp(-x) / sqrt(x);
    };

    std::cout << "EXTENDED MIDPOINT METHOD" << std::endl;
    Midpoint(F, a, b, acc);
    std::cout << std::endl;

    std::cout << "DE-RULE WITH TRAPEZOIDAL" << std::endl;
    DE(F, a, b, acc);
    std::cout << std::endl;
}

Int main()
{

    Exercise1();

    Exercise2();

    Exercise3();
    // As explained in the Exercise3() to see the plot,
    // you would have to have octave and input the following command:
    // data = csvread("Trajectory.csv"); x = data(:,1); y = data(:,2); plot(x,y)

    Exercise4();

    Exercise5();

}