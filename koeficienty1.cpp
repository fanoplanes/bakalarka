#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>

using namespace std;


int main()
{
	const complex<double> i  (0, 1);
	int N = 10;
	double theta = 0;
	double eps_a = 1;
	double eps_b = 4;
	double eps_air = 1;
	double mi_a = 1;
	double mi_b = 1;
	double mi_air = 1;
	double l_a = 1;
	double l_b = 0.5;
	double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
	double delta = 1e-5;
	double k_a;
	double k_b;
	double k_air;
	complex<double> phi_a = 0;
	complex<double> phi_b = 0;
	double chi_te;
	double chi_te_air_a;
	double chi_te_air_b;

	Eigen::Matrix2d m_ab;
	Eigen::Matrix2d m_ba;
	Eigen::Matrix2cd m_a;
	Eigen::Matrix2cd m_b;
	Eigen::Matrix2cd M;
	Eigen::Matrix2cd Out;
	Eigen::Matrix2d Air2A;
	Eigen::Matrix2d A2Air;

	ofstream fout("Output_chybaaa.dat", ios::app);

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

		M = m_ba * m_b * m_ab * m_a;

		Out = Air2A;

		for(int i = 0; i<N; i++)
		{
			Out *= M;
		}

		//Out *= m_a * m_a; //dve vrstvy chyby a

		for(int j = 0; j<N; j++)
		{
			Out *= M;
		}

		Out *= m_a;

		Out *= A2Air;

		fout << omega/omega_0 << "\t" << 1./(norm(Out(0, 0))) << endl;
		//fout << omega/omega_0 << "\t" <<abs(Out.determinant()) << endl;
	
	}
	system("xmgrace Output_chybaaa.dat");
}