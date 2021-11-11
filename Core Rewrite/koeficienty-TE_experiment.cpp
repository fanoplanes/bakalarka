#define SWITCH 1
#include <cmath>
#include <fstream>
#if SWITCH==1
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <boost/math/constants/constants.hpp>
#endif
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>
#if SWITCH==1
#define PI boost::math::constants::pi<FLOAT_TYPE>()
#define COMPLEX_TYPE boost::multiprecision::number<boost::multiprecision::backends::mpc_complex_backend<500> >
#define FLOAT_TYPE boost::multiprecision::number<boost::multiprecision::backends::mpfr_float_backend<500> >
using namespace boost::multiprecision;
#endif
#if SWITCH==0
#define COMPLEX_TYPE complex<double>
#define FLOAT_TYPE double
#define PI M_PI
#endif
using namespace std;
const double eps_one = 1.;		// one  | vzorka |  air
const COMPLEX_TYPE i (0, 1);

Eigen::Matrix<COMPLEX_TYPE, 2, 2> transfermatrix(FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE);
Eigen::Matrix<COMPLEX_TYPE, 2, 2> intermatrix(FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE);

int main()
{
	const int N = 100;
	const FLOAT_TYPE l_a = 1;
	const FLOAT_TYPE l_b = FLOAT_TYPE(1)/2;
	const FLOAT_TYPE sirkaap = 1;
	const FLOAT_TYPE sirkabp = FLOAT_TYPE(1)/2;
	const FLOAT_TYPE delta = FLOAT_TYPE(1)/100;
	const FLOAT_TYPE theta_delta = FLOAT_TYPE(1)/100;
	const FLOAT_TYPE eps_a = 1;
	const FLOAT_TYPE eps_air = 1;
	const FLOAT_TYPE eps_b = 4;
	const FLOAT_TYPE mi_a = 1;
	const FLOAT_TYPE mi_b = 1;
	const FLOAT_TYPE mi_air = 1;
	const FLOAT_TYPE mi_one = 1;
	const FLOAT_TYPE omega_0 = PI/(2*sqrt(eps_b)*l_b);
	FLOAT_TYPE structure[2*N][2];

	for(int i=0; i<N/2; i++)
	{
		structure[2*i][0] = eps_a;
		structure[2*i][1] = l_a;
	}
	for(int i=0; i<N/2; i++)
	{
		structure[2*i +1][0] = eps_b;
		structure[2*i +1][1] = l_b;
	}

	for(int i=N/2; i<N; i++)
	{
		structure[2*i][0] = eps_b;
		structure[2*i][1] = l_b;
	}
	for(int i=N/2; i<N; i++)
	{
		structure[2*i +1][0] = eps_a;
		structure[2*i +1][1] = l_a;
	}

	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
	Eigen::Matrix<COMPLEX_TYPE, 2, 2> temp;
	Eigen::Matrix<COMPLEX_TYPE, 2, 1> E;

	ofstream fout("Output-TE.dat");
	ofstream eout("Field.dat");

	//for(FLOAT_TYPE omega = delta; omega <= 2.*omega_0; omega+=delta)
	//{
	FLOAT_TYPE omega = 0.99*omega_0;
		//for(FLOAT_TYPE theta=0; theta < M_PI/2.; theta+=theta_delta)
		//{
		FLOAT_TYPE theta = 0;
			E <<	0.5, 0;
			Out = transfermatrix(eps_one, structure[0][0], theta, omega);
			E = Out * E;
			eout << 0 << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;
			for(int i=0; i<2*N-1; i++)
			{
				E << 0.5, 0;
				Out = intermatrix(structure[i][0], structure[i][1], theta, omega) * Out; //
				Out = transfermatrix(structure[i][0], structure[i+1][0], theta, omega) * Out; //
				E = Out * E;
				eout << i << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;
			}

			Out = intermatrix(structure[2*N-1][0], structure[2*N-1][1], theta, omega) * Out; //
			Out = transfermatrix(structure[2*N-1][0], eps_air, theta, omega) * Out; //
			E << 0.5, 0;
			E = Out * E;
			eout << 2*N-1 << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;
			cout << Out.determinant() << endl;

			COMPLEX_TYPE r = -Out(1,0)/Out(1,1);
			COMPLEX_TYPE t = 1./Out(0,0);
			FLOAT_TYPE R = norm(r);
			FLOAT_TYPE T = norm(t);
			if(eps_one == eps_air)
			{
				fout << omega/omega_0 << "\t" << theta*180./PI << "\t" << T << endl;
			}
			else
			{
				fout << omega/omega_0 << "\t" << theta*180./PI << "\t" << 1 - R << endl;
			}
		//}
		//fout << endl;
	//}
	//system("gnuplot -p -c plot-TE.p");
}
Eigen::Matrix<COMPLEX_TYPE, 2, 2> transfermatrix(FLOAT_TYPE permittivity1, FLOAT_TYPE permittivity2, FLOAT_TYPE theta, FLOAT_TYPE omega)
{
	FLOAT_TYPE inter_1 = permittivity1 - eps_one*sin(theta)*sin(theta);
	COMPLEX_TYPE k_1 = omega*sqrt(COMPLEX_TYPE(inter_1, 0));
	FLOAT_TYPE inter_2 = permittivity2 - eps_one*sin(theta)*sin(theta);
	COMPLEX_TYPE k_2 = omega*sqrt(COMPLEX_TYPE(inter_2, 0));
	COMPLEX_TYPE chi_te = k_1/k_2;
	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
	Out <<	1. + chi_te, 1. - chi_te,
			1. - chi_te, 1. + chi_te;
	Out /= FLOAT_TYPE(2);
	return Out;
}

Eigen::Matrix<COMPLEX_TYPE, 2, 2> intermatrix(FLOAT_TYPE permittivity, FLOAT_TYPE width, FLOAT_TYPE theta, FLOAT_TYPE omega)
{
	FLOAT_TYPE inter = permittivity - eps_one*sin(theta)*sin(theta);
	COMPLEX_TYPE k = omega*sqrt(COMPLEX_TYPE(inter, 0));
	COMPLEX_TYPE phi = i*k*width;
	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
	Out <<	exp(phi), 0,
			0, exp(-phi);
	return Out;
}
