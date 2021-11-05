#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <complex>

using namespace std;

int main()
{
	const int N = 10;
	const long double l_a = 1;
	const long double l_b = 0.5;
	const long double sirkaap = 1.;
	const long double sirkabp = 0.5;
	const long double delta = 1e-3;
	const long double theta_delta = 1e-3;
	const long double eps_a = 1.;			// one  | a | b | a | b | a | b | a |  air
	const long double eps_air = 1.;
	const long double eps_one = 1.;
	const long double eps_b = 4.;
	const long double mi_a = 1.;
	const long double mi_b = 1.;
	const long double mi_air = 1.;
	const long double mi_one = 1.;
	const long double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
	const complex<long double> i	(0, 1);
	complex<long double> k_a;
	complex<long double> k_b;
	complex<long double> k_air;
	long double k_one;
	complex<long double> chi_te;
	complex<long double> chi_te_air_a;
	complex<long double> chi_te_air_b;
	complex<long double> chi_te_one_a;
	complex<long double> chi_te_one_b;
	complex<long double> phi_a;
	complex<long double> phi_b;
	//const double theta = 0*M_PI/24.;

	Eigen::Matrix<complex<long double>, 2, 2> m_ab;
	Eigen::Matrix<complex<long double>, 2, 2> m_ba;
	Eigen::Matrix<complex<long double>, 2, 2> m_a;
	Eigen::Matrix<complex<long double>, 2, 2> m_b;
	Eigen::Matrix<complex<long double>, 2, 2> M;
	Eigen::Matrix<complex<long double>, 2, 2> Out;
	Eigen::Matrix<complex<long double>, 2, 2> Air2A;
	Eigen::Matrix<complex<long double>, 2, 2> A2Air;
	Eigen::Matrix<complex<long double>, 2, 2> Air2B;
	Eigen::Matrix<complex<long double>, 2, 2> B2Air;
	Eigen::Matrix<complex<long double>, 2, 2> One2A;
	Eigen::Matrix<complex<long double>, 2, 2> One2B;
	Eigen::Matrix<complex<long double>, 2, 2> m_ap;
	Eigen::Matrix<complex<long double>, 2, 2> m_bp;
	Eigen::Matrix<complex<long double>, 2, 2> I;

	I << 1,0,0,1;
	string input;
	for(int i=0; i<N; i++)
	{
		input += "ab";
	}

	int input_length = input.length();
	int structure [2*input_length + 1];		// 0 - one2a; 100 - one2b; 10 - a2air; 11 - b2air; 1 - a2b, 2 - b2a; 3 - a; 4 - b; 5 - a2a; 6 - b2b;

	string aa = "aa";
	string ab = "ab";
	string bb = "bb";
	string ba = "ba";
	string a = "a";
	string b = "b";

	string input_0 (1, input[0]);
	if(input_0==a) structure[0]=0;
	if(input_0==b) structure[0]=100;
	string input_last (1, input[input_length-1]);
	if(input_last==a) structure[2*input_length]=10;
	if(input_last==b) structure[2*input_length]=11;

	for(int iter=0; iter<input_length; iter++)
	{
		string input_iter (1, input[iter]);
		if(input_iter==a) structure[2*iter+1]=3;
		if(input_iter==b) structure[2*iter+1]=4;
	}

	for(int iter=0; iter<input_length-1; iter++)
	{
		string part_1 (1,input[iter]);
		string part_2 (1,input[iter+1]);
		string comp = part_1 + part_2;
		if(comp==aa) structure[2*iter+2]=5;
		if(comp==bb) structure[2*iter+2]=6;
		if(comp==ab) structure[2*iter+2]=1;
		if(comp==ba) structure[2*iter+2]=2;
	}

	ofstream fout("Output-TM.dat");

	for(long double omega = delta; omega <= 2*omega_0; omega+=delta)
	{
		for(long double theta=0; theta < M_PI/2.; theta+=theta_delta)
		{
		Out=I;
		k_one = omega*sqrt(eps_one)*cos(theta);
		long double inter_air = eps_air - eps_one*sin(theta)*sin(theta);
		k_air = omega*sqrt(complex<long double> (inter_air, 0));
		long double inter_a = eps_a - eps_one*sin(theta)*sin(theta);
		k_a = omega*sqrt(complex<long double> (inter_a, 0));
		long double inter_b = eps_b - eps_one*sin(theta)*sin(theta);
		k_b = omega*sqrt(complex<long double> (inter_b, 0));
		phi_a = i*k_a*l_a;
		phi_b = i*k_b*l_b;
		chi_te = eps_b*k_a/(eps_a*k_b);
		chi_te_air_a = eps_a*k_air/(eps_air*k_a);
		chi_te_air_b = eps_b*k_air/(eps_air*k_b);
		chi_te_one_a = eps_a*k_one/(eps_one*k_a);
		chi_te_one_b = eps_b*k_one/(eps_one*k_b);

		m_ab <<	1.L + chi_te, 1.L - chi_te,
				1.L - chi_te, 1.L + chi_te;

		m_ab /= 2.;
		m_ba = m_ab.inverse();

		Air2A <<	1.L + chi_te_air_a, 1.L - chi_te_air_a,
				1.L - chi_te_air_a, 1.L + chi_te_air_a;

		Air2A /= 2.L;

		A2Air = Air2A.inverse();

		Air2B <<	1.L + chi_te_air_b, 1.L - chi_te_air_b,
				1.L - chi_te_air_b, 1.L + chi_te_air_b;

		Air2B /= 2.L;

		B2Air = Air2B.inverse();

		One2A <<	1.L + chi_te_one_a, 1.L - chi_te_one_a,
				1.L - chi_te_one_a, 1.L + chi_te_one_a;

		One2A /= 2.L;

		One2B <<	1.L + chi_te_one_b, 1.L - chi_te_one_b,
				1.L - chi_te_one_b, 1.L + chi_te_one_b;

		One2B /= 2.L;

		m_a <<	exp(phi_a), 0,
				0, exp(-phi_a);

		m_b <<	exp(phi_b), 0,
				0, exp(-phi_b);

		m_bp <<	exp(i*k_b*sirkabp), 0,
				0, exp(-i*k_b*sirkabp);


		m_ap <<	exp(i*k_a*sirkaap), 0,
				0, exp(-i*k_a*sirkaap);

		for(int iterat=2*input_length; iterat>=0; iterat--)
		{
			if(structure[iterat]==0) Out *= One2A;
			if(structure[iterat]==1) Out *= m_ab;
			if(structure[iterat]==2) Out *= m_ba;
			if(structure[iterat]==3) Out *= m_a;
			if(structure[iterat]==4) Out *= m_b;
			if(structure[iterat]==5) Out *= I;
			if(structure[iterat]==6) Out *= I;
			if(structure[iterat]==10) Out *= One2A;
			if(structure[iterat]==11) Out *= B2Air;
			if(structure[iterat]==100) Out *= One2B;
		}


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
		//fout <<  << endl;
		//fout << omega/omega_0 << "\t" << abs(Air2A.determinant()) << endl;
		}
		fout << endl;
	}
	system("gnuplot -p -c plot-TM.p");
}
