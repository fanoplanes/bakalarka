#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>

using namespace std;

int main()
{
	const long double l_a = 15;
	const long double l_b = 0.5;
	const double delta = 1e-3;
	const double theta_delta = 1e-3;
	const double eps_a = 2;
	const double eps_air = 1;
	const double eps_b = 4;
	const double mi_a = 1;
	const double mi_b = 1;
	const double mi_air = 1;
	const double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
	const complex<long double> i  (0, 1);
	long double k_a;
	long double k_air;
	long double chi_te_air_a;
	complex<long double> phi_a;

	Eigen::Matrix<complex<long double>, 2, 2> m_a;
	Eigen::Matrix<complex<long double>, 2, 2> Out;
	Eigen::Matrix<long double, 2, 2> Air2A;
	Eigen::Matrix<long double, 2, 2> A2Air;

	ofstream fout("Output_FP_TE.dat");

      for(double omega = delta; omega < 2.*omega_0; omega+=delta)
      {
            for(double theta = 0; theta < M_PI/2.; theta += theta_delta )
            {
		k_air = omega*sqrt(eps_air - sin(theta)*sin(theta));
		k_a = omega*sqrt(eps_a - sin(theta)*sin(theta));
		phi_a = i*k_a*l_a;
		chi_te_air_a = mi_a*k_air/(mi_air*k_a);

		Air2A << 	1 + chi_te_air_a, 1 - chi_te_air_a,
				1 - chi_te_air_a, 1 + chi_te_air_a;

		A2Air = Air2A.inverse();

		m_a << 	exp(phi_a), 0,
				0, exp(-phi_a);

		Out = A2Air;
		Out *= m_a * Air2A;

		fout << omega/omega_0 << "\t" << theta*180./M_PI << "\t" << 1./(norm(Out(0, 0))) << endl;
		      }
		fout << endl;
	}
	//system("xmgrace Output-TE.dat");
}
