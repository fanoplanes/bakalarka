#include <cstdlib>

int main()
{
	system("g++ -march=native -Ofast -Wall koeficienty-TE_experiment.cpp -o koeficienty-TE -lmpfr -lmpc");
	system("./koeficienty-TE");
	return 0;
}
