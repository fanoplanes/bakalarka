#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#define PI M_PI

using namespace std;

const double eps_one = 4.; // one  | vzorka |  air
const double mi_one = 1.;

const double eps_a = 5.;
const double mi_a = 1.;

const double eps_b = 1.;
const double mi_b = 1.;

const double eps_air = 4.;
const double mi_air = 1.;

int main() {
	int n = 15; // 1597
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

	const int N = 1597;

	const double l_average = 1.5;

	const double phi = (1 + sqrt(5)) / 2;
	const double nu = 0.4;

	const double l_a = (1 + phi) * l_average / (nu + phi);
	const double l_b = nu * l_a;

	const double delta = 5e-3;
	const double theta_delta = 5e-3;

	double eps_parr = 0;
	double eps_perp = 0;
	double iterator = 0;

	const double omega_0 = PI / (2 * sqrt(eps_a) * l_b);

	double structure[N][2];
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

	ofstream fout("structure.dat");

	double step = 1e-2;
	int iterant = 0;

	for(int i=0; i<N; i++)
	{
		for(double j=0; j<structure[i][1]; j += step)
		{
			for(int k=0; k<5; k++){
				fout << iterant << "\t" << k << "\t" << structure[i][0] << endl;
			}
			iterant++;
			fout << endl;
		}
	}

	return 0;
}

// Local Variables:
// compile-command: "g++ -march=native -Ofast structure.cpp -o \
// structure && ./structure"
// End:
