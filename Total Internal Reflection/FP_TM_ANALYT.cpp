#include <fstream>
#include <iostream>
#include <complex>

using namespace std;

int main()
{
      //const double theta = 32*M_PI/128.;
      const long double eps_parr = 2.;
      const long double eps_perp = 4./3.;
      const long double eps_1 = 4.;
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
      complex<long double> r;
      const complex<long double> i (0,1);
      long double R;
      const double theta_delta = 1e-3;

      ofstream fout("Output_FP_TM_ANALYT.dat");

      for(long double omega = delta; omega <=2*omega_0; omega += delta)
      {
      for(long double theta = 0; theta < M_PI/2.; theta += theta_delta)
            {
            k_1z = omega*sqrt(eps_1)*cos(theta);
            long double inter_2 = 1 - sin(theta)*sin(theta)*eps_1/eps_perp;
            long double inter_3 = eps_3 - eps_1*sin(theta)*sin(theta);
            k_2z = sqrt(eps_parr)*omega*sqrt(complex<long double> (inter_2, 0.));
		k_3z = omega*sqrt(complex<long double> (inter_3,0));
		r1 = ((k_2z/eps_parr) - (k_1z/eps_1))/((k_2z/eps_parr) + (k_1z/eps_1));
		r2 = ((k_2z/eps_parr) - (k_3z/eps_3))/((k_2z/eps_parr) + (k_3z/eps_3));
		t1 = (2.L*k_1z/eps_1)/((k_2z/eps_parr) + (k_1z/eps_1));
		r = r1 + (t1*t1*r2*exp(i*2.L*k_2z*l))/(1.L-r1*r2*exp(i*2.L*k_2z*l));
		R = norm(r);
            fout << omega/omega_0 << "\t" << theta*180./M_PI << "\t" << R << endl;
            }
            fout << endl;
      }
      return 0;
}
