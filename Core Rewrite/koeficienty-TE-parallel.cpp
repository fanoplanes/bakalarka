#define SWITCH 1
#include "../Common/common.h"

using namespace std;

const FLOAT_TYPE eps_one = 1.; // one  | vzorka |  air
const FLOAT_TYPE mi_one = 1.;

const FLOAT_TYPE eps_a = 1.;
const FLOAT_TYPE mi_a = 1.;

const FLOAT_TYPE eps_b = 4.;
const FLOAT_TYPE mi_b = 1.;

const FLOAT_TYPE eps_air = 1.;
const FLOAT_TYPE mi_air = 1.;

int main()
{
	const int N = 72;
	const FLOAT_TYPE l_a = 1;
	const FLOAT_TYPE l_b = 1;
	const FLOAT_TYPE delta = 1e-2;
	const FLOAT_TYPE theta_delta = 1e-2;
	const FLOAT_TYPE omega_0 = PI/(2*sqrt(eps_b)*l_b);
	FLOAT_TYPE structure[2*N][2];
	unsigned long int thett;

	FLOAT_TYPE eps_parr=0;
	FLOAT_TYPE eps_perp=0;
	FLOAT_TYPE iterator=0;

	for(int i=0; i<N; i++)
	{
		structure[2*i][0] = eps_a;
		structure[2*i][1] = l_a;
	}
	for(int i=0; i<N; i++)
	{
		structure[2*i +1][0] = eps_b;
		structure[2*i +1][1] = l_b;
	}


	for(int a=0; a<2*N; a++)
	{
		eps_parr += structure[a][0]*structure[a][1];
		iterator += structure[a][1];
	}

	eps_parr /= iterator;

	for(int a=0; a<2*N; a++)
	{
		eps_perp += structure[a][1]/structure[a][0];
	}

	eps_perp = iterator/eps_perp;

	cout << eps_parr << endl;
	cout << eps_perp << endl;

	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;

	ofstream fout("Output-TE.dat");
	ofstream strout("structure.dat");
	for(int g=0; g<2*N; g++)
	{
		strout << structure[g][0] << "\t" << structure[g][1] << endl;
	}

	for(FLOAT_TYPE omega = delta; omega <= 2*omega_0; omega+=delta)
	{
		omp_set_num_threads(omp_get_max_threads());
		#pragma omp parallel for schedule (dynamic) ordered default (shared) private (thett, Out)
		for(thett=0; thett <= (unsigned long int)(PI/(2.*theta_delta)); thett++)
		{
			FLOAT_TYPE theta = thett*theta_delta;
			Out = transfermatrix(eps_one, structure[0][0], theta, omega, eps_one);
			for(int i=0; i<2*N-1; i++)
			{
				Out = intermatrix(structure[i][0], structure[i][1], theta, omega, eps_one) * Out; //
				Out = transfermatrix(structure[i][0], structure[i+1][0], theta, omega, eps_one) * Out; //
			}

			Out = intermatrix(structure[2*N-1][0], structure[2*N-1][1], theta, omega, eps_one) * Out; //
			Out = transfermatrix(structure[2*N-1][0], eps_air, theta, omega, eps_one) * Out; //

			COMPLEX_TYPE r = -Out(1,0)/Out(1,1);
			COMPLEX_TYPE t = 1./Out(0,0);
			FLOAT_TYPE R = norm(r);
			FLOAT_TYPE T = norm(t);
			#pragma omp ordered
			{
			if(eps_one == eps_air)
			{
				fout << omega/omega_0 << "\t" << theta*180./PI << "\t" << T << endl;
			}
			else
			{
				fout << omega/omega_0 << "\t" << theta*180./PI << "\t" << 1 - R << endl;
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
