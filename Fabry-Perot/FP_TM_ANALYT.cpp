#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

int main()
{
      //const double theta = 32*M_PI/128.;
      const double eps_parr = 2.;
      const double eps_perp = 4./3.;
      const double eps_air = 1.;
      const double eps_b = 4.;
      const double l_b = 0.5;
      const double mi = 1.;
      const double l=15;
      const double delta = 1e-3;
      const double omega_0 = M_PI/(2*sqrt(eps_b)*l_b);
      long double k_1z;
      long double k_2z;
      long double R_0;
      long double r;
      long double T;
      const double theta_delta = 1e-3;

      ofstream fout("Output_FP_TM_ANALYT.dat");

      for(double omega = delta; omega <=2*omega_0; omega += delta)
      {
            for(double theta = 0; theta < M_PI/2.; theta += theta_delta)
            {            k_1z = omega*sqrt(eps_air - sin(theta)*sin(theta));
            k_2z = sqrt(eps_parr)*omega*sqrt(1-sin(theta)*sin(theta)*eps_air/eps_perp);
            r = ((k_2z/eps_parr) - (k_1z/eps_air))/((k_2z/eps_parr)+(k_1z/eps_air));
            R_0 = r*r;
            T = 1/(1+(4*R_0*sin(k_2z*l)*sin(k_2z*l)/((1-R_0)*(1-R_0))));
            fout << omega/omega_0 << "\t" << theta*180./M_PI << "\t" << T << endl;
            }
            fout << endl;
      }
      return 0;
}
