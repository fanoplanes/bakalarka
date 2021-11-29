#include<iostream>
#include<cmath>

using namespace std;

int main()
{
	int exp=5;
	int L = pow(2, exp);
	bool result[(int)pow(2,exp)];
	result[0]=0;
	for(int i = 0; i<exp; i++)
	{
		for(int j=pow(2,i); j<pow(2, i+1); j++)
		{
			result[j]=!result[j-(int)pow(2,i)];
		}
	}

	for(int i=0; i<pow(2,exp); i++)
	{
		cout << result[i];
	}
	cout << endl;

}
