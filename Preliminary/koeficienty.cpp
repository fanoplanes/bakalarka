#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>

using namespace std;


int main()
{
	const complex<long double> i  (0, 1);
	const int N = 5;
	const double theta = 0;
	const double eps_a = 1;
	const double eps_b = 4;
	const double eps_air = 1;
	const double mi_a = 1;
	const double mi_b = 1;
	const double mi_air = 1;
	const long double l_a = 1;
	const long double l_b = 0.5;
	const double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
	const double delta = 1e-4;
	long double k_a;
	long double k_b;
	long double k_air;
	complex<long double> phi_a = 0;
	complex<long double> phi_b = 0;
	long double chi_te;
	long double chi_te_air_a;
	long double chi_te_air_b;
	const long double sirkaap = 1.;
	const long double sirkabp = 0.8;

	Eigen::Matrix<long double, 2, 2> m_ab;
	Eigen::Matrix<long double, 2, 2> m_ba;
	Eigen::Matrix<complex<long double>, 2, 2> m_a;
	Eigen::Matrix<complex<long double>, 2, 2> m_b;
	Eigen::Matrix<complex<long double>, 2, 2> M;
	Eigen::Matrix<complex<long double>, 2, 2> Out;
	Eigen::Matrix<long double, 2, 2> Air2A;
	Eigen::Matrix<long double, 2, 2> A2Air;
	Eigen::Matrix<complex<long double>, 2, 2> m_ap;
	Eigen::Matrix<complex<long double>, 2, 2> m_bp;

	ofstream fout("Output.dat");

	for(double omega = delta; omega < 2*omega_0; omega+=delta)
	{
		k_air = omega*sqrt(eps_air - sin(theta)*sin(theta));
		k_a = omega*sqrt(eps_a - sin(theta)*sin(theta));
		k_b = omega*sqrt(eps_b - sin(theta)*sin(theta));
		phi_a = i*k_a*l_a;
		phi_b = i*k_b*l_b;
		chi_te = mi_b*k_a/(mi_a*k_b);
		chi_te_air_a = mi_a*k_air/(mi_air*k_a);	//zo vzduchu do a

		m_ab << 1 + chi_te, 1 - chi_te,
				1 - chi_te, 1 + chi_te;
		
		m_ab /= 2;
		m_ba = m_ab.inverse();

		Air2A << 1 + chi_te_air_a, 1 - chi_te_air_a,
				 1 - chi_te_air_a, 1 + chi_te_air_a;

		A2Air = Air2A.inverse();

		m_a << 	exp(phi_a), 0,
				0, exp(-phi_a);
		
		m_b <<  exp(phi_b), 0,
				0, exp(-phi_b);

		m_bp << 	exp(i*k_b*sirkabp), 0,
		     		0, exp(-i*k_b*sirkabp);


		m_ap << 	exp(i*k_a*sirkaap), 0,
		     		0, exp(-i*k_a*sirkaap);

		M = m_a * m_ba * m_b * m_ab;		//|b|a

		//maticu budujeme od konca
		Out = A2Air;					//~

		for(int i=0; i<N; i++)
			{
				Out *= M;	//(|b|a)^N~
			}

		Out *= m_a * m_ba * m_bp * m_ab ;	//|b'|a(|b|a)^N~

		for(int j=0; j<N; j++)
			{
				Out *= M;	//(|b|a)^N|b'|a(|b|a)^N~
			}

		Out *= m_a * Air2A;		//~a(|b|a)^N|b'|a(|b|a)^N~

		fout << omega/omega_0 << "\t" << 1./(norm(Out(0, 0))) << endl;
		//fout << omega/omega_0 << "\t" << abs(Out.determinant()) << endl;
	
	}
	system("xmgrace Output.dat");
}
