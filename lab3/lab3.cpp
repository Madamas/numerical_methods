/*
 * Maslov Vadim, KV-51
 * 23. метод єдиного поділу та метод простої ітерації
 * 10 5 5 16    72
 * 19 38 18 0   76
 * 16 19 53 17  98
 * 10 14 14 4   48
*/
#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <string.h>
using namespace::std;

#define EPS 0.001

double const matrix[4][5]={
	{10,5,5,16,72},
	{19,38,18,0,76},
	{16,19,53,17,98},
	{10,14,14,4,48}
};
double const smatrix[4][5]={
	{-55,-10,2,36,-148},
	{19,38,18,0,76},
	{16,19,53,17,98},
	{25,7,21,61,222}
};

void print(double arr[4][5]){
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 5; j++){
			std::cout << arr[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

int* unifiedSeparation(double arr[4][5]){
	double mrow[5];
	double row[5];
	double buf[4][5];
	double ret[3];
	cout << "START" << endl;
			print(arr);
	for(int k = 0; k <= 3; k++){
		if(arr[k][k] == 0)
			return nullptr;
		else{
			for(int i = 0; i <= 4; i++)
				mrow[i] = arr[k][i] / arr[k][k];
				memcpy(arr[k],mrow,5*sizeof(double));
			for(int i = k + 1; i <= 3; i++){
				for(int j = 0; j < 5; j++)
					if(i != k)
						row[j] = (arr[i][j] - (mrow[j] * arr[i][k]));
					else
						continue;
				memcpy(arr[i],row,5*sizeof(double));
			}
			cout << "AFT" << endl;
			print(arr);
		}
	}
	for (int k = 3; k >= 0; k--){
			if(k == 3){
				ret[k] = arr[3][4]-arr[3][3];
			}else{
				double sum = 0;
				for(int i = k + 1; i <= 3; i++){
					sum += arr[k][i] * ret[i];
				}
				ret[k] = arr[k][4] - sum;
			}
			cout << "X" << k+1 << " = " << ret[k] << endl;
			cout << endl;		
	}
}


void iterations(double arr[4][5], double eps){
	double x[4], buf[4];
	double sum;
	double norm;
	double q = 0;

	for (int i = 0; i < 4; i++){
		sum = 0;
		for (int j = 0; j < 4; j++){
			sum += abs(arr[i][j]);
		}
		if (sum > q)
			q = sum;
	}

	for(int i = 0; i < 4; i++)
		buf[i] = arr[i][4];
	for (int i = 0; i < 4; i++){
		sum = 0;
		for (int j = 0; j < 4; j++)
			sum += arr[i][j] * buf[j];
		buf[i] = arr[i][4] + sum; 
	}

	int k = 2;
	while (true){
		for (int i = 0; i < 4; i++){
			x[i] = arr[i][4]; 
			for (int j = 0; j < 4; j++)
				x[i] += arr[i][j] * buf[j];
		}
		
		norm = 0;
    	for (int i = 0; i < 4; i++) {
        	norm += (x[i] - buf[i])*(x[i] - buf[i]);
    	}
		norm = sqrt(norm);

		if (norm <= eps*(1 - q)/q)
            break;

        k++;
		memcpy(buf, x, 4 * sizeof(double));
	}

	for(int i = 0; i < 4; i++)
		cout << "X[" << i+1 << "]= " << x[i] << endl;
	cout << "For eps = " << EPS << " num of iterations is = " << k << endl;
}

int main(){
double buffer[4][5];

memcpy(buffer,matrix,4*5*sizeof(double));
int* res = unifiedSeparation(buffer);

if(res == nullptr)
	cout << "Couldn't find solution" << endl;


memcpy(buffer, smatrix, 4*5*sizeof(double));

double ii;
	for (int i = 0; i < 4; i++){
		ii = buffer[i][i];
		for (int j = 0; j < 4; j++){
			if (i == j)
				buffer[i][i] = 0;
			else 
				buffer[i][j] = -buffer[i][j]/ii;

		}
		buffer[i][4] = buffer[i][4]/ii;
	}

iterations(buffer,EPS);
return 0;
}
