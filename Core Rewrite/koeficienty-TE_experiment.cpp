#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>

using namespace std;
const long double eps_one = 1.;		// one  | vzorka |  air
const complex<long double> i  (0, 1);

Eigen::Matrix<complex<long double>, 2, 2> transfermatrix(long double, long double, long double, long double);
Eigen::Matrix<complex<long double>, 2, 2> intermatrix(long double, long double, long double, long double);

int main()
{
	const int N = 150;
	const long double l_a = 1;
	const long double l_b = 0.5;
	const long double sirkaap = 1.;
	const long double sirkabp = 0.5;
	const long double delta = 1e-3;
	const long double theta_delta = 1e-3;
	const long double eps_a = 1.;
	const long double eps_air = 1.;
	const long double eps_b = 4.;
	const long double mi_a = 1.;
	const long double mi_b = 1.;
	const long double mi_air = 1.;
	const long double mi_one = 1.;
	const long double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
	long double structure[2*N][2];

	for(int i=0; i<N/3; i++)
	{
		structure[2*i][0] = eps_a;
		structure[2*i][1] = l_a;
	}
	for(int i=0; i<N/3; i++)
	{
		structure[2*i +1][0] = eps_b;
		structure[2*i +1][1] = l_b;
	}

	for(int i=N/3; i<2*N/3; i++)
	{
		structure[2*i][0] = eps_b;
		structure[2*i][1] = l_b;
	}
	for(int i=N/3; i<2*N/3; i++)
	{
		structure[2*i +1][0] = eps_a;
		structure[2*i +1][1] = l_a;
	}

	for(int i=2*N/3; i<N; i++)
	{
		structure[2*i][0] = eps_a;
		structure[2*i][1] = l_a;
	}
	for(int i=2*N/3; i<N; i++)
	{
		structure[2*i +1][0] = eps_b;
		structure[2*i +1][1] = l_b;
	}

	Eigen::Matrix<complex<long double>, 2, 2> Out;
	Eigen::Matrix<complex<long double>, 2, 1> E;

	ofstream fout("Output-TE.dat");
	ofstream eout("Field.dat");

	//for(long double omega = delta; omega <= 2.*omega_0; omega+=delta)
	//{
	long double omega = omega_0;
		//for(long double theta=0; theta < M_PI/2.; theta+=theta_delta)
		//{
		long double theta = 0;
			E <<	1, 0;
			Out = transfermatrix(structure[2*N-1][0], eps_air, theta, omega);

			for(int i=2*N-1; i>0; i--)
			{
				Out *= intermatrix(structure[i][0], structure[i][1], theta, omega);
				Out *= transfermatrix(structure[i-1][0], structure[i][0], theta, omega);
			}

			Out*= intermatrix(structure[0][0], structure[0][1], theta, omega);
			Out *= transfermatrix(eps_one, structure[0][0], theta, omega);

			E = transfermatrix(eps_one, structure[0][0], theta, omega) * E;
			E = intermatrix(structure[0][0], structure[0][1], theta, omega) * E;
			eout << 0 << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;
			for(int i=1; i<2*N; i++)
			{
				E = transfermatrix(structure[i-1][0], structure[i][0], theta, omega) * E;
				E = intermatrix(structure[i][0], structure[i][1], theta, omega) * E;
				eout << i << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;
			}
			//E << 1, 0;
			//E = transfermatrix(structure[2*N-1][0], eps_air, theta, omega).inverse() * E;
			//for(int i=2*N-1; i>=N; i--)
			//{
			//	E = intermatrix(structure[i][0], structure[i][1], theta, omega).inverse() * E;
			//	E = transfermatrix(structure[i-1][0], structure[i][0], theta, omega).inverse() * E;
			//	eout << i << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;
			//}

			complex<long double> r = -Out(1,0)/Out(1,1);
			complex<long double> t = 1.L/Out(0,0);
			long double R = norm(r);
			long double T = norm(t);
			if(eps_one == eps_air)
			{
				fout << omega/omega_0 << "\t" << theta*180./M_PI << "\t" << T << endl;
			}
			else
			{
				fout << omega/omega_0 << "\t" << theta*180./M_PI << "\t" << 1 - R << endl;
			}
		//}
		//fout << endl;
	//}
	//system("gnuplot -p -c plot-TE.p");
}
Eigen::Matrix<complex<long double>, 2, 2> transfermatrix(long double permittivity1, long double permittivity2, long double theta, long double omega)
{
	long double inter_1 = permittivity1 - eps_one*sin(theta)*sin(theta);
	complex<long double> k_1 = omega*sqrt(complex<long double> (inter_1, 0));
	long double inter_2 = permittivity2 - eps_one*sin(theta)*sin(theta);
	complex<long double> k_2 = omega*sqrt(complex<long double> (inter_2, 0));
	complex<long double> chi_te = k_1/k_2;
	Eigen::Matrix<complex<long double>, 2, 2> Out;
	Out <<	1.L + chi_te, 1.L - chi_te,
			1.L - chi_te, 1.L + chi_te;
	Out /= 2.;
	return Out;
}

Eigen::Matrix<complex<long double>, 2, 2> intermatrix(long double permittivity, long double width,long double theta, long double omega)
{
	long double inter = permittivity - eps_one*sin(theta)*sin(theta);
	complex<long double> k = omega*sqrt(complex<long double> (inter, 0));
	complex<long double> phi = i*k*width;
	Eigen::Matrix<complex<long double>, 2, 2> Out;
	Out <<	exp(phi), 0,
			0, exp(-phi);
	return Out;
}
