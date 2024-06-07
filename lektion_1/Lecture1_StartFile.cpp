#include <iostream>
#include "../source_code/Numerical-Recipes-master/nr3.h"
#include "../source_code/Numerical-Recipes-master/ludcmp.h"
#include "utilities.h"

using namespace std;

int main()
{

	// Exercise 1:
	// Solve A x = b using LU decomposition, and print the result.

	MatDoub A(3, 3);
	A[0][0] = 1.0;
	A[0][1] = 2.0;
	A[0][2] = 3.0;
	A[1][0] = 2.0;
	A[1][1] = -4.0;
	A[1][2] = 6.0;
	A[2][0] = 3.0;
	A[2][1] = -9.0;
	A[2][2] = -3.0;

	VecDoub b(3);
	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;

	util::print(b);

	std::cout << "\nmatrix a\n";
	util::print(A);

	// create LU decomposition
	// allocat x vector
	VecDoub x(3);
	LUdcmp ALU(A);
	ALU.solve(b, x);

	util::print(ALU.lu, "LU");

	// print x
	util::print(x);

	// evaluate x

	// print x

	return 0;
}
