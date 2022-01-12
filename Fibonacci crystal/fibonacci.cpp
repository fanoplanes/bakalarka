#include<iostream>
#include<vector>

using namespace std;

int main()
{
	int n=6;
	int arr0[] = {0};
	int arr1[] = {0,1};
	vector<bool> vec0;
	vector<bool> vec1;
	vector<bool> vec;
	vec0.assign(arr0, arr0+1);
	vec1.assign(arr1, arr1+2);
	for(int i=2; i<=n; i++)
	{
		vec=vec1;
		vec.insert(vec.end(),vec0.begin(),vec0.end());
		vec0=vec1;
		vec1=vec;
	}
}
