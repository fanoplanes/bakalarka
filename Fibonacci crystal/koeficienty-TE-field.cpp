#define SWITCH 1
#include "../Common/common.h"
#include <string>

using namespace std;

const FLOAT_TYPE eps_one = 4.; // one  | vzorka |  air
const FLOAT_TYPE mi_one = 1.;

const FLOAT_TYPE eps_a = 1.;
const FLOAT_TYPE mi_a = 1.;

const FLOAT_TYPE eps_b = 5.;
const FLOAT_TYPE mi_b = 1.;

const FLOAT_TYPE eps_air = 4.;
const FLOAT_TYPE mi_air = 1.;

int main(int argc, char *argv[]) {

  FLOAT_TYPE x = 0.0905148 + 0.00458171*(FLOAT_TYPE)argv[1];
  FLOAT_TYPE y = 0.7066*x + 0.9575;
  int num = atoi(argv[1]);

  int n = 10;
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
  const FLOAT_TYPE nu = 0.4;

  const FLOAT_TYPE l_a = (1 + phi) * l_average / (nu + phi);
  const FLOAT_TYPE l_b = nu * l_a;

  const FLOAT_TYPE delta = 1e-2;
  const FLOAT_TYPE theta_delta = 1e-2;

  const FLOAT_TYPE omega_0 = PI / (2 * sqrt(eps_b) * l_b);

  FLOAT_TYPE structure[N][2];
  unsigned long int thett;

  for (int i = 0; i < N; i++) {
    structure[i][0] = (i % 2 == 0) ? eps_a : eps_b;
    structure[i][1] = (vec[i] == 0) ? l_a : l_b;
  }

  Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
  Eigen::Matrix<COMPLEX_TYPE, 2, 1> E;
  Eigen::Matrix<COMPLEX_TYPE, 2, 2> I;

  I << 1, 0, 0, 1;

  string file = "Field";
  string end = ".dat";
  string number = to_string(num);
  string answer = file + number + end;

  ofstream eout(answer);

  FLOAT_TYPE omega = x * omega_0;
  FLOAT_TYPE theta = y;

  E << 1, 0;
  Out = transfermatrix(eps_one, structure[0][0], theta, omega, eps_one);
  E = Out * E;
  eout << 0 << "\t" << sqrt(norm(E(0)) + norm(E(1))) << endl;

  for (int i = 0; i < N - 1; i++) {
    Out = I;
    Out = intermatrix(structure[i][0], structure[i][1], theta, omega, eps_one) * Out; //
    Out = transfermatrix(structure[i][0], structure[i + 1][0], theta, omega, eps_one) *
          Out; //
    E = Out * E;
    eout << i << "\t" << sqrt(norm(E(0)) + norm(E(1))) << endl;
  }

  Out = I;
  Out = intermatrix(structure[N - 1][0], structure[N - 1][1], theta, omega, eps_one) *
        Out;                                                              //
  Out = transfermatrix(structure[N - 1][0], eps_air, theta, omega, eps_one) * Out; //
  E = Out * E;
  eout << N - 1 << "\t" << sqrt(norm(E(0)) + norm(E(1))) << endl;

  system("gnuplot -p -c Field.p");
}

// Local Variables:
// compile-command: "g++ -march=native -Ofast koeficienty-TE-field.cpp -o \
// koeficienty-TE-field -lmpc -lmpfr -fopenmp && ./koeficienty-TE-field"
// End:
