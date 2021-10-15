#include <fstream>
#include <iostream>
#include <complex>

using namespace std;

int main()
{
	//const double theta = 32*M_PI/128.;
	const long double eps_1 = 16.;
	const long double eps_2 = 4.;
	const long double eps_3 = 1.;
	const long double eps_b = 4.;
	const long double l_b = 0.5;
	const long double l=15;
	const long double delta = 1e-3;
	const long double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
	complex<long double> k_1z;
	complex<long double> k_2z;
	complex<long double> k_3z;
	complex<long double> r1;
	complex<long double> r2;
	complex<long double> t1;
	complex<long double> t2;
	complex<long double> t;
	complex<long double> r;
	const complex<long double> i (0,1);
	long double T;
	long double R;
	const double theta_delta = 1e-3;

	ofstream fout("Output_FP_TE_ANALYT.dat");

	for(long double omega = delta; omega <=0.5*omega_0; omega += delta)
	{
	for(long double theta = 0; theta < M_PI/2.; theta += theta_delta)
		{
		k_1z = omega*sqrt(eps_1)*cos(theta);
		long double inter_2 = eps_2 - eps_1*sin(theta)*sin(theta);
		long double inter_3 = eps_3 - eps_1*sin(theta)*sin(theta);
		k_2z = omega*sqrt(complex<long double> (inter_2, 0.));
		k_3z = omega*sqrt(complex<long double> (inter_3,0));
		r = - ((1.L - (k_1z/k_3z))*cos(k_2z*l) - i * ((k_2z/k_3z)-(k_1z/k_2z))*sin(k_2z*l))/((1.L + (k_1z/k_3z))*cos(k_2z*l) - i * ((k_2z/k_3z)+(k_1z/k_2z))*sin(k_2z*l));
		R = norm(r);
		fout << omega/omega_0 << "\t" << theta*180./M_PI << "\t" << 1 - R << endl;
		}
		fout << endl;
	}
	return 0;
}
