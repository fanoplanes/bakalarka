#include "../Common/common.h"
#include <random>

using namespace std;

const FLOAT_TYPE eps_one = 1.; // one  | vzorka |  air
const FLOAT_TYPE mi_one = 1.;

const FLOAT_TYPE eps_a = 1.;
const FLOAT_TYPE mi_a = 1.;

const FLOAT_TYPE eps_b = 4.;
const FLOAT_TYPE mi_b = 1.;

const FLOAT_TYPE eps_air = 1.;
const FLOAT_TYPE mi_air = 1.;

int main() {
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0,1);

	const int N = 144;

	const FLOAT_TYPE l_average = 1;

	const FLOAT_TYPE phi = (1 + sqrt(5)) / 2;
	const FLOAT_TYPE nu = 0.4;

	const FLOAT_TYPE l_a = (1 + phi) * l_average / (nu + phi);
	const FLOAT_TYPE l_b = nu * l_a;

	const FLOAT_TYPE delta = 1e-2;
	const FLOAT_TYPE theta_delta = 1e-2;

	FLOAT_TYPE eps_parr = 0;
	FLOAT_TYPE eps_perp = 0;
	FLOAT_TYPE iterator = 0;

	const FLOAT_TYPE omega_0 = PI / (2 * sqrt(eps_b) * l_b);

	FLOAT_TYPE structure[N][2];
	unsigned long int thett;

	for (int i = 0; i < N; i++) {
		structure[i][0] = (i % 2 == 0) ? eps_a : eps_b;
		int temp = distribution(generator);
		structure[i][1] = (temp == 0) ? l_a : l_b;
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
	ofstream strout("structure.dat");
	for(int g=0; g<N; g++)
	{
		strout << structure[g][0] << "\t" << structure[g][1] << endl;
	}

	for (FLOAT_TYPE omega = delta; omega <= 2.0 * omega_0; omega += delta) {
		omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic)                                     \
		ordered default(shared) private(thett, Out)
		for (thett = 0; thett <= (unsigned long int)(PI / (2. * theta_delta)); thett++) {
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
#pragma omp ordered
			{
				if (eps_one == eps_air) {
					fout << omega / omega_0 << "\t" << theta * 180. / PI << "\t" << T
							 << endl;
				} else {
					fout << omega / omega_0 << "\t" << theta * 180. / PI << "\t"
							 << (1 - R) << endl;
				}
			}
		}
		fout << endl;
	}
	system("gnuplot -p -c plot-TE.p");
	return 0;
}

// Local Variables:
// compile-command: "g++ -march=native -Ofast koeficienty-TE-parallel.cpp -o \
// koeficienty-TE-parallel -lmpc -lmpfr -fopenmp && ./koeficienty-TE-parallel"
// End:
