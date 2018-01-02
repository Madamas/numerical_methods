// Маслов КВ-51
// sh(x)
// [-1.8;1.9]
#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <stdio.h>
#include <unistd.h>
#include <cmath>
using std::abs;

double* series1(double eps,double x){
	double sum = x;
	double un = x;
	double prev = x;
	int k = 1;
	while(un >= eps){
	//	printf("%.6e\t%.6e\t%.6e\t%d\n",un,prev,sum,k);
		un = prev*x*x/(2*k*(2*k+1));
		prev = un;
		sum += un;
		k++;
	}
	//printf("%.6e\t%.6e\t%.6e\t%d\n",un,prev,sum,k);
	double *arr = new double[3];
	arr[0] = sum;
	arr[1] = k;
	arr[2] = prev/3.1;
	return arr;
};

double* series2(int k,double x){
	double sum = x;
	double un = x;
	double prev = x;
	for(int i = 1;i < k; i++){
		//printf("%.6e\t%.6e\t%.6e\t%.2e\n",un,prev,sum,x);
		un = prev*x*x/(2*i*(2*i+1));
		prev = un;
		sum += un;
	}
	//printf("%.6e\t%.6e\t%.6e\t%.2e\n",un,prev,sum,x);
	double *arr = new double[3];
	arr[0] = sum;
	arr[1] = k;
	arr[2] = prev/3;
	return arr;
};


int main(){
	double eps = 0.01;
	printf("eps\tn\tdelta\t\tR\n");
	for(int i = 0; i <= 4; i++){
		double* arr = series1(eps,0.05);
		printf("%.0e\t%d\t%.6e\t%.6e\n", eps, (int)arr[1], abs(sinh(0.05)-arr[0]), arr[2]);
		eps *= 0.001;
	}
	
	double h = 0.37;
	printf("\n x\tdelta\t\tR\n");
	for(int i = 0; i<=10;i++){
		double x = -1.8+h*i;
		x = trunc(x*100)/100;
		double* arr = series2(3,x);
		printf("%c%.2f\t%.6e\t%.6e\n", x >= 0 ? ' ' : '\0', x, abs(sinh(x)-arr[0]), arr[2]);
	}
	return 0;
};