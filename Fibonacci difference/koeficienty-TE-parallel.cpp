#define SWITCH 1
#include "../Common/common.h"

using namespace std;

const FLOAT_TYPE eps_one = 4.; // one  | vzorka |  air
const FLOAT_TYPE mi_one = 1.;

const FLOAT_TYPE eps_a = 1.;
const FLOAT_TYPE mi_a = 1.;

const FLOAT_TYPE eps_b = 5.;
const FLOAT_TYPE mi_b = 1.;

const FLOAT_TYPE eps_air = 4.;
const FLOAT_TYPE mi_air = 1.;

int main() {
  int n = 10; // 377
  int arr0[] = {0};
  int arr1[] = {0, 1};
  vector<bool> vec0;
  vector<bool> vec1;
  vector<bool> vec;
  vec0.assign(arr0, arr0 + 1);
  vec1.assign(arr1, arr1 + 2);
  for (int i = 2; i <= n; i++) {
    vec = vec1;
    vec.insert(vec.end(), vec0.begin(), vec0.end());
    vec0 = vec1;
    vec1 = vec;
  }

  const int N = 144;

  const FLOAT_TYPE l_average = 1.5;

  const FLOAT_TYPE phi = (1 + sqrt(5)) / 2;
  const FLOAT_TYPE nu = 1.0;

  const FLOAT_TYPE l_a = (1 + phi) * l_average / (nu + phi);
  const FLOAT_TYPE l_b = nu * l_a;
  const FLOAT_TYPE l = N * l_average;

  const FLOAT_TYPE delta = 5e-4;
  const FLOAT_TYPE theta_delta = 2e-3;

  FLOAT_TYPE eps_parr = 0;
  FLOAT_TYPE eps_perp = 0;
  FLOAT_TYPE iterator = 0;

  const FLOAT_TYPE omega_0 = PI / (2 * sqrt(eps_b) * l_b);

  FLOAT_TYPE structure[N][2];
  unsigned long int thett;

  for (int i = 0; i < N; i++) {
    structure[i][0] = (i % 2 == 0) ? eps_b : eps_a;
    structure[i][1] = (vec[i] == 0) ? l_b : l_a;
  }

  for (int a = 0; a < N; a++) {
    eps_parr += structure[a][0] * structure[a][1];
    iterator += structure[a][1];
  }

  eps_parr /= iterator;

  for (int a = 0; a < N; a++) {
    eps_perp += structure[a][1] / structure[a][0];
  }

  eps_perp = iterator / eps_perp;

  cout << eps_parr << endl;
  cout << eps_perp << endl;

  Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;

  ofstream fout("Output-TE.dat");
  FLOAT_TYPE start = 0.00 / theta_delta;
  FLOAT_TYPE stop = PI /(2 * theta_delta);

  for (FLOAT_TYPE omega = delta; omega <= 2.5 * omega_0; omega += delta) {
    omp_set_num_threads(omp_get_max_threads());
    #pragma omp parallel for schedule(dynamic) ordered default(shared) private(thett, Out)

    for (thett = (long int)start; thett <= (unsigned long int)(stop); thett++) {

      FLOAT_TYPE theta = thett * theta_delta;
      Out = transfermatrix(eps_one, structure[0][0], theta, omega, eps_one);

      for (int i = 0; i < N - 1; i++) {
        Out = intermatrix(structure[i][0], structure[i][1], theta, omega, eps_one) * Out;
        Out = transfermatrix(structure[i][0], structure[i + 1][0], theta, omega, eps_one) * Out;
      }

      Out = intermatrix(structure[N - 1][0], structure[N - 1][1], theta, omega, eps_one) * Out;
      Out = transfermatrix(structure[N - 1][0], eps_air, theta, omega, eps_one) * Out;

      COMPLEX_TYPE r = -Out(1, 0) / Out(1, 1);
      COMPLEX_TYPE t = 1. / Out(0, 0);
      FLOAT_TYPE R = norm(r);
      FLOAT_TYPE T = norm(t);

      // analytic computation
      COMPLEX_TYPE k_1z_an = omega * sqrt(eps_one) * cos(theta);
      COMPLEX_TYPE inter_2_an = eps_parr - eps_one * sin(theta) * sin(theta);
      COMPLEX_TYPE inter_3_an = eps_air - eps_one * sin(theta) * sin(theta);
      COMPLEX_TYPE k_2z_an = omega * sqrt(inter_2_an);
      COMPLEX_TYPE k_3z_an = omega * sqrt(inter_3_an);
      COMPLEX_TYPE r_an =
          -((1.L - (k_1z_an / k_3z_an)) * cos(k_2z_an * l) -
            i * ((k_2z_an / k_3z_an) - (k_1z_an / k_2z_an)) *
                sin(k_2z_an * l)) /
          ((1.L + (k_1z_an / k_3z_an)) * cos(k_2z_an * l) -
           i * ((k_2z_an / k_3z_an) + (k_1z_an / k_2z_an)) * sin(k_2z_an * l));
      COMPLEX_TYPE R_an = norm(r_an);

      #pragma omp ordered
      {
        if (eps_one == eps_air) {
          fout << omega / omega_0 << "\t" << theta * 180. / PI << "\t"
               << T - (1 - R_an) << endl;
        } else {
          fout << omega / omega_0 << "\t" << theta * 180. / PI << "\t"
               << (1 - R) - (1 - R_an) << endl;
        }
      }
    }
    fout << endl;
  }
  system("gnuplot -p -c plot-TE.p");
}

// Local Variables:
// compile-command: "g++ -march=native -Ofast koeficienty-TE-parallel.cpp -o \
// koeficienty-TE-parallel -lmpc -lmpfr -fopenmp && ./koeficienty-TE-parallel"
// End:
