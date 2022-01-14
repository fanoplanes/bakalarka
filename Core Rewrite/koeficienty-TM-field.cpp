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

int main()
{
	const int N = 10;
	const FLOAT_TYPE l_a = 1;
	const FLOAT_TYPE l_b = 0.5;
	const FLOAT_TYPE delta = 1e-2;
	const FLOAT_TYPE theta_delta = 1e-2;
	const FLOAT_TYPE omega_0 = PI/(2*sqrt(eps_b)*l_b);
	FLOAT_TYPE structure[2*N][2];
	unsigned long int thett;

	for(int i=0; i<N/2; i++)
	{
		structure[2*i][0] = eps_a;
		structure[2*i][1] = l_a;
	}
	for(int i=0; i<N/2; i++)
	{
		structure[2*i +1][0] = eps_b;
		structure[2*i +1][1] = l_b;
	}

	for(int i=N/2; i<N; i++)
	{
		structure[2*i][0] = eps_b;
		structure[2*i][1] = l_b;
	}
	for(int i=N/2; i<N; i++)
	{
		structure[2*i +1][0] = eps_a;
		structure[2*i +1][1] = l_a;
	}


	Eigen::Matrix<COMPLEX_TYPE, 2, 2> Out;
	Eigen::Matrix<COMPLEX_TYPE, 2, 1> E;

	ofstream fout("Output-TE.dat");
	ofstream eout("Field.dat");

	FLOAT_TYPE omega = omega_0;
	FLOAT_TYPE theta = 0;

			E << 1, 0;
			Out = transfermatrix(eps_one, structure[0][0], theta, omega, eps_one);
			E = Out * E;
			eout << 0 << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;

			for(int i=0; i<2*N-1; i++)
			{
				E << 1, 0;
				Out = intermatrix(structure[i][0], structure[i][1], theta, omega, eps_one) * Out; //
				Out = transfermatrix(structure[i][0], structure[i+1][0], theta, omega, eps_one) * Out; //
				E = Out * E;
				eout << i << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;
			}

			E << 1, 0;
			Out = intermatrix(structure[2*N-1][0], structure[2*N-1][1], theta, omega, eps_one) * Out; //
			Out = transfermatrix(structure[2*N-1][0], eps_air, theta, omega, eps_one) * Out; //
			E = Out * E;
			eout << 2*N-1 << "\t" << sqrt(norm(E(0))+norm(E(1))) << endl;

			COMPLEX_TYPE r = -Out(1,0)/Out(1,1);
			COMPLEX_TYPE t = 1./Out(0,0);
			FLOAT_TYPE R = norm(r);
			FLOAT_TYPE T = norm(t);
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
	system("gnuplot -p -c Field.p");
}

// Local Variables:
// compile-command: "g++ -march=native -Ofast koeficienty-TM-field.cpp -o \
// koeficienty-TM-field -lmpc -lmpfr -fopenmp && ./koeficienty-TM-field"
// End:
