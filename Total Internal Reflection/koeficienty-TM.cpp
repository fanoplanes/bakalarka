#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>

using namespace std;


int main()
{
	const complex<long double> i (0, 1);
	const int N = 10;
	const long double theta = 4*M_PI/10.;
	const long double eps_a = 1;
	const long double eps_b = 4;
	const long double eps_air = 1;
	const long double eps_one = 4.;
	const long double l_a = 1;
	const long double l_b = 0.5;
	const long double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
	const long double delta = 1e-3;
	const long double theta_delta = 1e-3;
	const long double sirkaap = 1.; 				//šírka poruchy a
	const long double sirkabp = 0.5;				//šírka poruchy b
	complex<long double> k_one;
	complex<long double> k_air;
	complex<long double> k_a;
	complex<long double> k_b;

	Eigen::Matrix<complex<long double>, 2, 2> m_ab;
	Eigen::Matrix<complex<long double>, 2, 2> m_ba;
	Eigen::Matrix<complex<long double>, 2, 2> m_a;
	Eigen::Matrix<complex<long double>, 2, 2> m_b;
	Eigen::Matrix<complex<long double>, 2, 2> M;
	Eigen::Matrix<complex<long double>, 2, 2> Out;
	Eigen::Matrix<complex<long double>, 2, 2> Air2A;
	Eigen::Matrix<complex<long double>, 2, 2> A2Air;
	Eigen::Matrix<complex<long double>, 2, 2> One2A;
	Eigen::Matrix<complex<long double>, 2, 2> m_ap;
	Eigen::Matrix<complex<long double>, 2, 2> m_bp;

	ofstream fout("Output-TM.dat");

	for(long double omega = delta; omega < 2*omega_0; omega+=delta)
	{
		for(long double theta = 0; theta < M_PI/2.; theta += theta_delta)
		{
		k_one = omega*sqrt(eps_one)*cos(theta);
		long double inter_air = eps_air - eps_one*sin(theta)*sin(theta);
		k_air = omega*sqrt(complex<long double> (inter_air, 0));
		long double inter_a = eps_a - eps_one*sin(theta)*sin(theta);
		k_a = omega*sqrt(complex<long double> (inter_a, 0));
		long double inter_b = eps_b - eps_one*sin(theta)*sin(theta);
		k_b = omega*sqrt(complex<long double> (inter_b, 0));
		complex<long double> phi_a = i*k_a*l_a;
		complex<long double> phi_b = i*k_b*l_b;
		complex<long double> chi_tm = eps_b*k_a/(eps_a*k_b);
		complex<long double> chi_tm_air_a = eps_a*k_air/(eps_air*k_a);	//zo vzduchu do a
		complex<long double> chi_tm_one_a = eps_a*k_one/(eps_one*k_a);

		m_ab <<	1.L + chi_tm, 1.L - chi_tm,
				1.L - chi_tm, 1.L + chi_tm;
		m_ab /= 2.L;
		m_ba = m_ab.inverse();

		Air2A <<	1.L + chi_tm_air_a, 1.L - chi_tm_air_a,
				1.L - chi_tm_air_a, 1.L + chi_tm_air_a;

		Air2A /= 2.L;

		A2Air = Air2A.inverse();

		One2A <<	1.L + chi_tm_one_a, 1.L - chi_tm_one_a,
				1.L - chi_tm_one_a, 1.L + chi_tm_one_a;

		m_a <<	exp(phi_a), 0,
				0, exp(-phi_a);

		m_b <<	exp(phi_b), 0,
				0, exp(-phi_b);

		m_bp <<	exp(i*k_b*sirkabp), 0,
				0, exp(-i*k_b*sirkabp);


		m_ap <<	exp(i*k_a*sirkaap), 0,
				0, exp(-i*k_a*sirkaap);

		M = m_a * m_ba * m_b * m_ab;						//|b|a

		//maticu budujem od konca (postupne pridávam do matice Out členy násobením *sprava*)
		Out = A2Air;								//~

		for(int i=0; i<N; i++)
		{
			Out *= M;								//(|b|a)^N~
		}

		//Out *= m_a * m_ba * m_bp * m_ab ;	//porucha b			//|b'|a(|b|a)^N~
		//Out *= m_ap * m_ba * m_b * m_ab;	//porucha a			//|b|a'(|b|a)^N~

		//for(int j=0; j<N; j++)
		//{
		//	Out *= M;								//(|b|a)^N|b'|a(|b|a)^N~
		//}

		Out *= m_a * One2A;							//~a(|b|a)^N|b'|a(|b|a)^N~

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
		//fout << omega/omega_0 << "\t" << abs(Out.determinant()) << endl;
		}
            fout << endl;
	}
	//system("xmgrace Output-TM.dat");
}
