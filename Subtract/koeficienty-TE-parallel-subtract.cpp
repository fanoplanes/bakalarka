#define SWITCH 1 //v mode 1 pocita s boost multiprecision, v mode 0 s built-in double presnostou
#include <cmath>
#include <fstream>
#include <omp.h>

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
#define COMPLEX_TYPE boost::multiprecision::number<boost::multiprecision::backends::mpc_complex_backend<100> >
#define FLOAT_TYPE boost::multiprecision::number<boost::multiprecision::backends::mpfr_float_backend<100> >
#define PI boost::math::constants::pi<FLOAT_TYPE>()
using namespace boost::multiprecision;
#endif

#if SWITCH==0
#define COMPLEX_TYPE complex<double>
#define FLOAT_TYPE double
#define PI M_PI
#endif

using namespace std;
const FLOAT_TYPE eps_one = 1.;		// one  | vzorka |  air
const FLOAT_TYPE mi_one = 1.;
const FLOAT_TYPE eps_air = 1.;
const FLOAT_TYPE mi_air = 1.;

const COMPLEX_TYPE i  (0, 1);

Eigen::Matrix<COMPLEX_TYPE, 2, 2> transfermatrix(FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE);
Eigen::Matrix<COMPLEX_TYPE, 2, 2> intermatrix(FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE);
FLOAT_TYPE numeric(FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE, FLOAT_TYPE);

int main()
{
	const int N = 10;
	const FLOAT_TYPE l_a = 1;
	const FLOAT_TYPE l_b = 0.5;
	const FLOAT_TYPE sirkaap = 1.;
	const FLOAT_TYPE sirkabp = 0.5;
	const FLOAT_TYPE delta = 1e-2;
	const FLOAT_TYPE theta_delta = 1e-2;
	const FLOAT_TYPE eps_a = 1.;
	const FLOAT_TYPE eps_b = 4.;
	const FLOAT_TYPE mi_a = 1.;
	const FLOAT_TYPE mi_b = 1.;
	const FLOAT_TYPE omega_0 = PI/(2*sqrt(eps_b)*l_b);
	FLOAT_TYPE eps_eff_parr=0;
	FLOAT_TYPE length;
	FLOAT_TYPE structure[2*N][2];
	unsigned long int thett;

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

	for(int i=0; i<2*N; i++)
	{
		eps_eff_parr+= structure[i][0]*structure[i][1];
	}

	length=0;

	for(int i=0; i<2*N; i++)
	{
		length += structure[i][1];
	}

	eps_eff_parr/=length;

	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;

	ofstream fout("Output-TE.dat");

	for(FLOAT_TYPE omega = delta; omega <= 0.5*omega_0; omega+=delta)
	{
		omp_set_num_threads(omp_get_max_threads());
		#pragma omp parallel for schedule (dynamic) ordered default (shared) private (thett, Out)
		for(thett=0; thett <= (unsigned long int)(PI/(2.*theta_delta)); thett++)
		{
			FLOAT_TYPE theta = thett*theta_delta;
			Out = transfermatrix(eps_one, structure[0][0], theta, omega);
			for(int i=0; i<2*N-1; i++)
			{
				Out = intermatrix(structure[i][0], structure[i][1], theta, omega) * Out; //
				Out = transfermatrix(structure[i][0], structure[i+1][0], theta, omega) * Out; //
			}

			Out = intermatrix(structure[2*N-1][0], structure[2*N-1][1], theta, omega) * Out; //
			Out = transfermatrix(structure[2*N-1][0], eps_air, theta, omega) * Out; //

			COMPLEX_TYPE r = -Out(1,0)/Out(1,1);
			COMPLEX_TYPE t = 1./Out(0,0);
			FLOAT_TYPE R = norm(r);
			FLOAT_TYPE T = norm(t);
			#pragma omp ordered
			{
			if(eps_one == eps_air)
			{
				fout << omega/omega_0 << "\t" << theta*180./PI << "\t" << abs(T-numeric(eps_eff_parr, omega, theta, length)) << endl;
			}
			else
			{
				fout << omega/omega_0 << "\t" << theta*180./PI << "\t" << abs((1 - R) - numeric(eps_eff_parr, omega, theta, length)) << endl;
			}
			}
		}
		fout << endl;
	}
	system("gnuplot -p -c plot-TE.p");
}
Eigen::Matrix<COMPLEX_TYPE, 2, 2> transfermatrix(FLOAT_TYPE permittivity1, FLOAT_TYPE permittivity2, FLOAT_TYPE theta, FLOAT_TYPE omega)
{

	FLOAT_TYPE inter_1 = permittivity1 - eps_one*sin(theta)*sin(theta);
	COMPLEX_TYPE k_1 = omega*sqrt(COMPLEX_TYPE (inter_1, 0));
	FLOAT_TYPE inter_2 = permittivity2 - eps_one*sin(theta)*sin(theta);
	COMPLEX_TYPE k_2 = omega*sqrt(COMPLEX_TYPE (inter_2, 0));
	COMPLEX_TYPE chi_te = k_1/k_2;
	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
	Out <<	1. + chi_te, 1. - chi_te,
			1. - chi_te, 1. + chi_te;
	Out /= 2.;
	return Out;
}

Eigen::Matrix<COMPLEX_TYPE, 2, 2> intermatrix(FLOAT_TYPE permittivity, FLOAT_TYPE width,FLOAT_TYPE theta, FLOAT_TYPE omega)
{
	FLOAT_TYPE inter = permittivity - eps_one*sin(theta)*sin(theta);
	COMPLEX_TYPE k = omega*sqrt(COMPLEX_TYPE (inter, 0));
	COMPLEX_TYPE phi = i*k*width;
	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
	Out <<	exp(phi), 	0,
			0, 		exp(-phi);
	return Out;
}
FLOAT_TYPE numeric(FLOAT_TYPE eps_2, FLOAT_TYPE omega, FLOAT_TYPE theta, FLOAT_TYPE l)
{
	COMPLEX_TYPE k_1z = omega*sqrt(eps_one)*cos(theta);
	FLOAT_TYPE inter_2 = eps_2 - eps_one*sin(theta)*sin(theta);
	FLOAT_TYPE inter_3 = eps_air - eps_one*sin(theta)*sin(theta);
	COMPLEX_TYPE k_2z = omega*sqrt(COMPLEX_TYPE (inter_2, 0.));
	COMPLEX_TYPE k_3z = omega*sqrt(COMPLEX_TYPE(inter_3,0));
	COMPLEX_TYPE r = - ((1. - (k_1z/k_3z))*cos(k_2z*l) - i * ((k_2z/k_3z)-(k_1z/k_2z))*sin(k_2z*l))/((1. + (k_1z/k_3z))*cos(k_2z*l) - i * ((k_2z/k_3z)+(k_1z/k_2z))*sin(k_2z*l));
	return 1 - norm(r);
}
