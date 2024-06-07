#pragma once
#include "nr3.h"
#include "utilities.h"
#include "tridag.h"
#include "quadrature.h"
#include "roots_multidim.h"
#include <fstream>
#include <assert.h>
#include <limits>
#include <random>

using F_1 = Doub (*)(Doub);
using F_2 = Doub (*)(Doub, Doub);
using F_3 = Doub (*)(Doub, Doub, Doub);

/*-------------------------------
       Generate a CSV file
-------------------------------*/

void generateCSV(Doub a, Doub b, VecDoub y)
{
    std::ofstream V("V.csv");
    std::ofstream U("U.csv");
    Doub x = (b - a) / Doub(y.size());
    for(int i = 0; i < y.size(); i++)
    {
        V << y[i] << "\n";
        U << a+i*x << "\n";
    }

    V.close();
    U.close();
}



/*-------------------------------
        Richardson Aanlysis
-------------------------------*/

class Richardson
{
private:
    std::vector<double> current;
    std::vector<double> current_error;
    std::vector<double> rich_alpha_k;
    std::vector<double> rich_error;
    std::vector<double> f_computations;
public:
    Richardson(double, double);
    double getError();
    void update(double, double, bool);
    void print();
    ~Richardson();
};

Richardson::Richardson(double curr, double f_comp)
{
    current.push_back(curr);
    current_error.push_back(std::numeric_limits<double>::infinity());
    rich_alpha_k.push_back(std::numeric_limits<double>::infinity());
    rich_error.push_back(std::numeric_limits<double>::infinity());
    f_computations.push_back(f_comp);
}

double Richardson::getError()
{
    return rich_error[rich_error.size() - 1];
}

void Richardson::update(double curr, double f_comp, bool is_order_given=false)
{
    current.push_back(curr);
    f_computations.push_back(f_comp);

    double newest   = current[current.size() - 1];
    double previous = current[current.size() - 2];

    double curr_err  = previous - newest;
    current_error.push_back(curr_err);
    
    if (current.size() <= 2)
    {
        rich_alpha_k.push_back(std::numeric_limits<double>::infinity());
        rich_error.push_back(std::numeric_limits<double>::infinity());
    } 
    else
    {
        double old = current[current.size() - 3];

        double rich_alpk;
        if(is_order_given) rich_alpk = 4.;
        else rich_alpk = (old - previous) / (previous - newest);

        double rich_err  = (newest - previous) /  (rich_alpk - 1);

        rich_alpha_k.push_back(rich_alpk);
        rich_error.push_back(rich_err);
    }
}

void Richardson::print()
{
    std::cout << "i" << setw(15) << "A(hi)" << setw(15) << "A(hi-1)-A(hi)" << setw(15) << "Rich-alpk" << setw(15) << "Rich-error" << setw(15) << "f-computations" << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;

    for (int i = 0; i < current.size(); i++)
    {
        if (i + 1 >= 10)
        {
            std::cout << i + 1 << setw(14) << current[i] << setw(15) << current_error[i] << setw(15) << 
            rich_alpha_k[i] << setw(15) << rich_error[i] << setw(15) << f_computations[i] << std::endl;
        }
        else
        {
            std::cout << i + 1 << setw(15) << current[i] << setw(15) << current_error[i] << setw(15) << 
            rich_alpha_k[i] << setw(15) << rich_error[i] << setw(15) << f_computations[i] << std::endl;
        }
    }
}

Richardson::~Richardson(){}



/*-------------------------------
     Make Vectors Orthonormal
-------------------------------*/

std::vector<VecDoub> GramSchmidt(std::vector<VecDoub> &x)
{
    // Assume the vectors are linearly independent

    std::vector<VecDoub> e = {x[0] / util::norm(x[0])};

    for (int i = 1; i < x.size(); i++)
    {
        VecDoub sum = vecZeros(x[0].size());

        for (int j = 0; j < i; j++)
        {
            sum = sum + util::dot(x[i], e[j]) * e[j];
        }

        VecDoub ei = x[i] - sum;
        ei = ei / util::norm(ei);
        e.push_back(ei);
    }

    return e;
}



/*-------------------------------
        1D Root Finding
-------------------------------*/

void RootLog()
{
    std::cout << "i" << setw(15) << "root" << setw(15) << "droot" << setw(15) << "C" << setw(15) << "e" << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;

}

void RootLogNum(double i, double root, double droot, double C, double e)
{
    if (i + 1 >= 10)
        {
            std::cout << i + 1 << setw(14) << root << setw(15) << droot << setw(15) << 
            C << setw(15) << e << std::endl;
        }
        else
        {
            std::cout << i + 1 << setw(15) << root << setw(15) << droot << setw(15) << 
            C << setw(15) << e << std::endl;
        }
}

Doub Bisection(F_1 F, Doub x_min, Doub x_max, Doub error_threshold)
{
    Doub error = std::numeric_limits<double>::infinity();
    Doub mid_point;
    Doub dummy_mid_point = 0.;
    Doub convergence = 0.5;
    int count = 0;
    RootLog();
    while(abs(error) > error_threshold)
    {
        mid_point = (x_min + x_max) / Doub(2);

        if(F(mid_point) * F(x_max) < 0)
        {
            x_min = mid_point;
        }
        else
        {
            x_max = mid_point;
        }
        Doub dist = mid_point - dummy_mid_point;
        RootLogNum(count, mid_point, dist, convergence, error);
        error = abs(dist);
        dummy_mid_point = mid_point;
        count++;
    }
    return mid_point;
}

Doub Secant(F_1 F, Doub x, Doub error_threshold)
{
    Doub lower_bound = 0.;
    Doub upper_bound = M_PI / Doub(2);
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine engine;
    Doub random = unif(engine);
    Doub order = 1./2. * (1. + sqrt(5.));
    int i = 0;

    RootLog();

    Doub x_curr = random;
    Doub x_next = x_curr - (x_curr - x) / (F(x_curr) - F(x)) * F(x_curr);
    Doub x_prev = x;

    Doub dist_prev = 0.;
    Doub dist_curr = x_curr - x_prev;

    Doub convergence = abs(dist_curr) / pow(abs(dist_prev), order);
    Doub error = convergence * pow(abs(dist_curr), order);

    while(abs(error) > error_threshold)
    {
        RootLogNum(i, x_curr, dist_curr, convergence, error);

        x_prev = x_curr;
        x_curr = x_next;
        x_next = x_curr - (x_curr - x) / (F(x_curr) - F(x)) * F(x_curr);

        dist_prev = dist_curr;
        dist_curr = x_curr - x_prev;

        convergence = abs(dist_curr) / pow(abs(dist_prev), order);
        error = convergence * pow(abs(dist_curr), order);

        i++;
    }

    return x_curr;
}

Doub Newton(F_1 F, F_1 Fd, Doub x, Doub error_threshold)
{
    int i = 0;
    Doub order = 2.;
    Doub x_prev = x;
    x = x_prev - F(x_prev) / Fd(x_prev);
    Doub dist_curr = x - x_prev;
    Doub dist_prev = 1.;
    Doub convergence = abs(dist_curr) / pow(abs(dist_prev), 2);
    Doub error = convergence * pow(abs(dist_curr), 2);
    RootLog();
    while(abs(error) > error_threshold)
    {
        RootLogNum(i, x, dist_curr, convergence, error);
        x_prev = x;
        x = x_prev - F(x_prev) / Fd(x_prev);
        dist_curr = x - x_prev;
        dist_prev = dist_curr;
        convergence = abs(dist_curr) / pow(abs(dist_prev), 2);
        error = convergence * pow(abs(dist_curr), 2);
        i++;
    }
    return x;
}



/*-------------------------------
      Multi-D Root Finding
-------------------------------*/

template <class T>
void NewtonTable(VecDoub_IO &x, Bool &check, T &vecfunc) {
	const Int MAXITS=200;
	const Doub TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
	const Doub TOLX=numeric_limits<Doub>::epsilon();
	Int i,j,its,n=x.size();
	Doub den,f,fold,stpmax,sum,temp,test;
	VecDoub g(n),p(n),xold(n);
	MatDoub fjac(n,n);
	NRfmin<T> fmin(vecfunc);
	NRfdjac<T> fdjac(vecfunc);
	VecDoub &fvec=fmin.fvec;
	f=fmin(x);
	test=0.0;

	//Changed File Here
	std::cout << "i";
	for (int i1 = 0; i1 < x.size(); i1++) std::cout << setw(14) << "x" << i1+1;
	std::cout << setw(15) << "dk";
	std::cout << setw(15) << "C";
	std::cout << setw(15) << "e";
	std::cout << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;

	//

	for (i=0;i<n;i++)
		if (abs(fvec[i]) > test) test=abs(fvec[i]);
	if (test < 0.01*TOLF) {
		check=false;
		return;
	}
	sum=0.0;
	for (i=0;i<n;i++) sum += SQR(x[i]);
	stpmax=STPMX*MAX(sqrt(sum),Doub(n));

	//Changed File Here
	double dk_old = 0;
	//

	for (its=0;its<MAXITS;its++) {
		fjac=fdjac(x,fvec);
		for (i=0;i<n;i++) {
			sum=0.0;
			for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=0;i<n;i++) xold[i]=x[i];
		fold=f;
		for (i=0;i<n;i++) p[i] = -fvec[i];
		LUdcmp alu(fjac);
		alu.solve(p,p);
		lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);

		//Changed file here
		std::cout << its+1;
        if(its+1 > 9) std::cout << setw(14);
        else std::cout << setw(15);
		for (int i2 = 0; i2 < x.size(); i2++) 
        {
            std::cout << x[i2] << setw(15);
        }
		double dk = util::norm(x - xold);
		double convergence = dk / abs(dk_old * dk_old);
		// double convergence = 1000;
		std::cout << dk;
		std::cout << setw(15) << convergence;
		std::cout << setw(15) << convergence * dk * dk;
		std::cout << std::endl;
		dk_old = dk;
		//

		test=0.0;
		for (i=0;i<n;i++)
			if (abs(fvec[i]) > test) test=abs(fvec[i]);
		if (test < TOLF) {
			check=false;
			return;
		}
		if (check) {
			test=0.0;
			den=MAX(f,0.5*n);
			for (i=0;i<n;i++) {
				temp=abs(g[i])*MAX(abs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			check=(test < TOLMIN);
			return;
		}
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(abs(x[i]-xold[i]))/MAX(abs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX)
			return;
	}
	throw("MAXITS exceeded in newt");
}



/*-------------------------------
        Integration
-------------------------------*/

template<typename T>
Doub Midpoint(T F, Doub a, Doub b, Doub error_threshold)
{
    int steps = 2;
    Doub step_size = (b - a) / Doub(steps);
    Doub sum = step_size * F(a + step_size + step_size / Doub(2));
    Richardson rich(sum, steps);
    while(abs(rich.getError()) > error_threshold)
    {
        sum = 0;
        steps *= 2;
        step_size = (b - a) / Doub(steps);
        for(int i = 0; i < steps; i++)
        {
            sum += step_size * F(a + Doub(i) * step_size + step_size / Doub(2));
        }

        rich.update(sum, steps);
    }
    rich.print();
    return sum;
}

template<typename T>
Doub Trapezoidal(T F, Doub a, Doub b, Doub error_threshold)
{
    int steps = 1;
    Doub step_size = (b - a) / Doub(steps);
    Doub sum = step_size / Doub(2) * (F(a) + F(b));
    Richardson rich(sum, steps);
    while(abs(rich.getError()) > error_threshold)
    {
        sum = 0;
        steps *= 2;
        step_size = (b - a) / Doub(steps);
        for (int i = 0; i <= steps; i++)
        {
            sum += step_size * F(a + Doub(i) * step_size);
        }
        sum -= step_size / Doub(2) * (F(a) + F(b));
        rich.update(sum, steps);
    }
    rich.print();
    return sum;
}

template<typename T>
Doub Simpsons(T F, Doub a, Doub b, Doub error_threshold)
{
    int steps = 1;
    Doub step_size = (b - a) / Doub(steps);
    Doub sum = step_size / Doub(3) * (F(a) + F(b));
    Richardson rich(sum, steps);
    while(abs(rich.getError()) > error_threshold)
    {
        sum = 0;
        steps *= 2;
        step_size = (b - a) / Doub(steps);
        for (int i = 1; i < steps; i++)
        {
            if(i % 2 == 1)
            {
                sum += Doub(4) / Doub(3) * step_size * F(a + Doub(i) * step_size); 
            }
            else
            {
                sum += Doub(2) / Doub(3) * step_size * F(a + Doub(i) * step_size);
            }
        }
        sum += step_size / Doub(3) * (F(a) + F(b));
        rich.update(sum, steps);
    }
    rich.print();
    return sum;
}

template<typename T>
Doub DE(T f, Doub a, Doub b, Doub error_threshold){
    Doub c = 1;
    Doub H = 4.3;
    auto F = [=](Doub t){
        Doub dxdt = 0.5 * (b - a) * pow(1. / cosh(c *sinh(t)), 2) * c * cosh(t);
        Doub q = exp(-2 * sinh(abs(t)));
        Doub d = (b - a) * q / (1. + q);
        return t < 0. ? f(a + d)*dxdt : f(b - d)*dxdt;
    };
    return Trapezoidal(F, -H, H, error_threshold);
}



/*-------------------------------
      Initial Value Problem
-------------------------------*/

template<typename T>
VecDoub Euler(T F, Doub start, Doub end, VecDoub init, int N, Doub error_threshold)
{
    VecDoub y;
    Doub t, h;
    
    Richardson rich_x(init[0], 0);
    Richardson rich_y(init[1], 0);
    
    while(abs(rich_x.getError()) > error_threshold || abs(rich_y.getError()) > error_threshold)
    {
        t = start;
        h = (end - start) / Doub(N);
        y = init;
        for(int i = 0; i < N; i++, t+=h)
        {
            y = y + h * F(t, y);
        }

        rich_x.update(y[0], N);
        rich_y.update(y[1], N);
        N *= 2;
    }

    std::cout << "Richardson for x" << std::endl;
    rich_x.print();
    std::cout << std::endl;

    std::cout << "Richardson for y" << std::endl;
    rich_y.print();
    std::cout << std::endl;
    
    return y;
}

template<typename T>
VecDoub MidpointODE(T F, Doub start, Doub end, VecDoub init, int N, Doub error_threshold)
{
    Doub t = start, h = (end - start) / Doub(N);
    VecDoub y = init, k1 = h * F(t, y), k2 = h * F(t + 0.5 * h, y + 0.5 * k1);
    y = y + k2;
    
    Richardson rich(y[2], N);
    
    while(abs(rich.getError()) > error_threshold)
    {
        N *= 2;
        t = start;
        h = (end - start) / Doub(N);
        y = init;
        for(int i = 0; i < N; i++, t+=h)
        {
            k1 = h * F(t, y);
            k2 = h * F(t + 0.5 * h, y + 0.5 * k1);
            y = y + k2;
        }
        rich.update(y[2], N, true);
    }

    std::cout << "Richardson for y(2)" << std::endl;
    rich.print();
    std::cout << std::endl;
    
    return y;
}

template<typename T>
VecDoub RungeKutta4(T F, Doub start, Doub end, VecDoub init, int N, Doub error_threshold)
{
    Doub t, h;
    VecDoub y, k1, k2, k3, k4;
    
    Richardson rich_x(init[0], 0);
    Richardson rich_y(init[1], 0);
    
    while(abs(rich_x.getError()) > error_threshold || abs(rich_y.getError()) > error_threshold)
    {
        t = start;
        h = (end - start) / Doub(N);
        y = init;
        for(int i = 0; i < N; i++, t+=h)
        {
            k1 = h * F(t, y);
            k2 = h * F(t + 0.5 * h, y + 0.5 * k1);
            k3 = h * F(t + 0.5 * h, y + 0.5 * k2);
            k4 = h * F(t + h, y + k3);
            y = y + 1./6. * k1 + 1./3. * k2 + 1./3. * k3 + 1./6. * k4;
        }

        rich_x.update(y[0], N);
        rich_y.update(y[1], N);
        N *= 2;
    }

    std::cout << "Richardson for x" << std::endl;
    rich_x.print();
    std::cout << std::endl;

    std::cout << "Richardson for y" << std::endl;
    rich_y.print();
    std::cout << std::endl;
    
    return y;
}

template<typename T>
VecDoub TrapezoidalODE(T F, Doub start, Doub end, VecDoub init, int N, Doub error_threshold)
{
    VecDoub y, y_next;
    Doub t, h;

    auto phi = [&](VecDoub y_next)
    {
        return y_next - y - h / 2. * (F(t, y) + F(t + h, y_next));
    };

    Richardson rich_x(init[0], 0);
    Richardson rich_y(init[1], 0);

    while(abs(rich_x.getError()) > error_threshold || abs(rich_y.getError()) > error_threshold)
    {
        h = (end - start) / Doub(N);
        y = init;
        t = start;

        for(int i = 0; i < N; i++, t += h)
        {
            y_next = y + h * F(t, y);
            bool check;
            newt(y_next, check, phi);
            assert(!check);
            y = y_next;
        }

        rich_x.update(y[0], N);
        rich_y.update(y[1], N);

        N *= 2;
    }

    std::cout << "Richardson for x" << std::endl;
    rich_x.print();
    std::cout << std::endl;

    std::cout << "Richardson for y" << std::endl;
    rich_y.print();
    std::cout << std::endl;
    
    return y;
}

template<typename T>
VecDoub MidpointPlot(T F, Doub start, Doub end, VecDoub init)
{
    VecDoub y, k1, k2;
    Doub t, h;
    std::vector<Doub> y1;
    std::vector<Doub> y2;

    h = 0.001;
    t = start;
    y = init;

    y1.push_back(y[0]);
    y2.push_back(y[2]);

    int N = (end - start) / h;
    for(int i = 0; i < N; i++, t+=h)
    {
        k1 = h * F(t, y);
        k2 = h * F(t + 0.5 * h, y + 0.5 * k1);
        y = y + k2;
        y1.push_back(y[0]);
        y2.push_back(y[2]);
    }

    assert(y1.size() == y2.size());
    std::ofstream my_file("Trajectory.csv");
    for(int i = 0; i < y1.size(); i++)
    {
        my_file << y1[i] << "," << y2[i] << "\n";
    }
    my_file.close();

    return y;
}



/*-------------------------------
        Finite Difference
-------------------------------*/

VecDoub addYValues(VecDoub &y)
{
    VecDoub new_y(2 * y.size() - 1);
    for (int i = 0; i < y.size() - 1; i++)
    {
        new_y[2 * i] = y[i];
        new_y[2 * i + 1] = (y[i] + y[i+1]) / 2;
    }
    new_y[new_y.size() - 1] = y[y.size() - 1];
    return new_y;
}

struct vec
{
    VecDoub lower;
    VecDoub center;
    VecDoub upper;
};

vec jacobian(Doub &a, Doub &b, VecDoub &y, Doub (*Fy)(Doub, Doub, Doub), Doub (*Fyd)(Doub, Doub, Doub))
{
    Doub h = (b - a) / Doub(y.size() - 1);
    MatDoub j(y.size()-2, y.size()-2);
    j[0][0]                   = 2. + pow(h, 2) * Fy((y[2] - y[0]) / (2. * h), y[1], a + h);
    j[0][1]                   = -1. + h / 2. * Fyd((y[2] - y[0]) / (2. * h), y[1], a + h);
    j[y.size()-3][y.size()-4] = -1. - h / 2. * Fyd((y[y.size()-1] - y[y.size()-3]) / (2. * h), y[y.size()-2], b - h);
    j[y.size()-3][y.size()-3] = 2. + pow(h, 2) * Fy((y[y.size()-1] - y[y.size()-3]) / (2. * h), y[y.size()-2], b - h);
    for (int i = 2; i <= y.size() - 2; i++)
    {
        j[i-1][i-2] = -1. - h / 2. * Fyd((y[i+1] - y[i-1]) / (2. * h), y[i], a + Doub(i) * h);
        j[i-1][i-1] = 2. + pow(h, 2) * Fy((y[i+1] - y[i-1]) / (2. * h), y[i], a + Doub(i) * h);
        j[i-1][i]   = -1. + h / 2. * Fyd((y[i+1] - y[i-1]) / (2. * h), y[i], a + Doub(i) * h);
    }
    VecDoub lower(y.size()-2);
    VecDoub center(y.size()-2);
    VecDoub upper(y.size()-2);
    lower[0] = 0.;
    upper[upper.size() - 1] = 0.;
    for (int row = 1; row < lower.size(); row++)
    {
        lower[row] = j[row][row-1];
    }
    for (int row = 0; row < upper.size() - 1; row++)
    {
        upper[row] = j[row][row+1];
    }
    for (int row = 0; row < center.size(); row++)
    {
        center[row] = j[row][row];
    }
    return {lower, center, upper};
}

VecDoub Phi(Doub &a, Doub &b, VecDoub &y, Doub (*F)(Doub, Doub, Doub))
{
    Doub h = (b - a) / Doub(y.size() - 1);
    VecDoub phi(y.size()-2);
    for (int i = 1; i < y.size() - 1; i++)
    {
        phi[i-1] = -y[i-1] + 2. * y[i] - y[i+1] + pow(h, 2) * F((y[i+1] - y[i-1]) / (2. * h), y[i], a + Doub(i) * h); 
    }
    return phi;
}

VecDoub updateY(VecDoub &y, VecDoub &lower, VecDoub &center, VecDoub &upper, VecDoub &phi)
{
    VecDoub solution(y.size()-2);
    tridag(lower, center, upper, -1 * phi, solution);
    for (int i = 1; i < y.size() - 1; i++)
    {
        y[i] = y[i] + solution[i-1];
    }
    VecDoub new_y = addYValues(y);
    return new_y;
}

VecDoub FiniteDiff(Doub a, Doub b, Doub alpha, Doub beta, F_3 F, F_3 Fy, F_3 Fyd, Doub error_threshold)
{
    VecDoub y = util::vec2(alpha, beta);
    y = addYValues(y);

    Richardson rich(y[y.size() / 2], y.size());

    while(abs(rich.getError()) > error_threshold)
    {
        auto [lower, center, upper] = jacobian(a, b, y, Fy, Fyd);
        VecDoub phi = Phi(a, b, y, F);
        y = updateY(y, lower, center, upper, phi);
        rich.update(y[y.size() / 2], y.size());
    }

    generateCSV(a, b, y);

    rich.print();
    
    return y;
}