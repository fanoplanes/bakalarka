#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>

using namespace std;


int main()
{
	const complex<double> i (0., 1.);
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
	double delta = 1e-3;
	double k_a;
	double k_b;
	double k_air;
	double sirka;
	complex<double> phi_a = 0;
	complex<double> phi_b = 0;
	complex<double> phi_p = 0;
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
	Eigen::Matrix2cd Porucha_a;
	Eigen::Matrix2cd Porucha_b;
	Eigen::Matrix2cd P_b;

	ofstream fout("Output.dat");
	sirka =0.5;

	//for(double omega = delta; omega < 2*omega_0; omega+=delta)
	//{
		double omega = omega_0;
		k_air = omega*sqrt(eps_air - sin(theta)*sin(theta));
		k_a = omega*sqrt(eps_a - sin(theta)*sin(theta));
		k_b = omega*sqrt(eps_b - sin(theta)*sin(theta));
		phi_a = i*k_a*l_a;
		phi_b = i*k_b*l_b;
		chi_te = mi_b*k_a/(mi_a*k_b);
		chi_te_air_a = mi_a*k_air/(mi_air*k_a);	//zo vzduchu do a

		m_ab << 	1 + chi_te, 1 - chi_te,
				1 - chi_te, 1 + chi_te;

		m_ab /= 2;
		m_ba = m_ab.inverse();

		Air2A << 	1 + chi_te_air_a, 1 - chi_te_air_a,
				1 - chi_te_air_a, 1 + chi_te_air_a;

		A2Air = Air2A.inverse();

		m_a << 	exp(phi_a), 0,
				0, exp(-phi_a);

		m_b <<  	exp(phi_b), 0,
				0, exp(-phi_b);

		P_b << 	exp(phi_b*sirka), 0,
				0, exp(-phi_b*sirka);

		//M = m_ba * m_b * m_ab * m_a;
		M = m_a * m_ba * m_b * m_ab;			// |b|a
		Porucha_b  = m_a *m_ba * P_b* m_ab ;	// |b|a

		Out = m_a * Air2A;				// ~a
		cout << phi_a << "\t" << phi_b << endl;
		cout << omega/omega_0 << "\t" << abs(Out.determinant()) << endl;

		for(int a = 0; a<N; a++)
			{
				Out.noalias() = M * Out;				//~a(|b|a)^N
				cout << omega/omega_0 << "\t" << abs(M.determinant()) << endl;
				cout << omega/omega_0 << "\t" << abs(Out.determinant()) << endl;
			}

	//	Out.noalias() = Porucha_b * m_ab * Out; 		//~a(|b|a)|b'

		for(int j = 0; j<N; j++)
			{
			//	Out.noalias() = M * m_a * m_ba * Out;		//~a(|b|a)^N|b'|a(|b|a)^N
				Out.noalias() = M * Out;				//~a(|b|a)^N|b'|a(|b|a)^N
			}

		Out.noalias() = A2Air * Out;				//~a(|b|a)^N|b'|a(|b|a)^N~

		//fout << omega/omega_0 << "\t" << 1./(norm(Out(0, 0))) << endl;
		//fout << omega/omega_0 << "\t" <<abs(Out.determinant()) << endl;
		//cout << omega/omega_0 << "\t" << abs(Out.determinant()) << endl;

	//}
	//system("xmgrace Output.dat");
}
