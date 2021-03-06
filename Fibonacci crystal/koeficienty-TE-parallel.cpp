#include "../Common/common.h"

using namespace std; // pre túto štruktúru n=144 je kumulovaná hodnota transmitivity 34.9
//pre n=610 je to 16.9

const FLOAT_TYPE eps_one = 4.; // one  | vzorka |  air
const FLOAT_TYPE mi_one = 1.;

const FLOAT_TYPE eps_a = 1.;
const FLOAT_TYPE mi_a = 1.;

const FLOAT_TYPE eps_b = 5.;
const FLOAT_TYPE mi_b = 1.;

const FLOAT_TYPE eps_air = 4.;
const FLOAT_TYPE mi_air = 1.;

int main() {
	FLOAT_TYPE cumulator = 0;
	int n = 12;
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
		structure[i][1] = (vec[i] == 0) ? l_a : l_b;
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

	FLOAT_TYPE beginning_float = asin(sqrt(eps_parr/eps_one))/theta_delta;
	beginning_float = floor(beginning_float);
	unsigned long int beginning = (unsigned long int)beginning_float;

	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;

	ofstream fout("Output-TE.dat");
	ofstream strout("structure.dat");
	for(int g=0; g<N; g++)
	{
		strout << structure[g][0] << "\t" << structure[g][1] << endl;
	}

	for (FLOAT_TYPE omega = delta; omega <= 0.5 * omega_0; omega += delta) {
		omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic) ordered default(shared) private(thett, Out)
		for (thett = beginning; thett <= (unsigned long int)(PI / (2. * theta_delta)); thett++) {
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
					//fout << omega / omega_0 << "\t" << theta * 180. / PI << "\t" << T
					//		 << endl;
					cumulator += T;
				} else {
					//fout << omega / omega_0 << "\t" << theta * 180. / PI << "\t"
					//		 << (1 - R) << endl;
					cumulator += (1-R);
				}
			}
		}
		//fout << endl;
	}
	//system("gnuplot -p -c plot-TE.p");
	cout << cumulator << endl;
	return 0;
}

// Local Variables:
// compile-command: "g++ -march=native -Ofast koeficienty-TE-parallel.cpp -o \
// koeficienty-TE-parallel -lmpc -lmpfr -fopenmp && ./koeficienty-TE-parallel"
// End:
