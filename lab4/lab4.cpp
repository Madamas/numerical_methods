/*
KV-51 Maslov Vadim
23. function - (cos(x) + x*sin(x)/(cos(x))^2)
boundaries - [-1.5;1.5]
primitive - x/cos(x)
Simpson's method
*/
#include <iostream>
#include <cstddef>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <string.h>
using namespace::std;

// #define sec(x) (1/cos(x))
#define A -1.5
#define B 1.5
double sec(double x){
	return (1/cos(x));
}

double f(double x){
	return ((cos(x) + x*sin(x))/(cos(x)*cos(x)));
}

double d4f(double x){
	return ((sec(x)*sec(x))*(x*sin(x)-3*cos(x))+4*x*cos(x)*(16*tan(x)*pow(sec(x),4)+8*pow(tan(x),3)*sec(x)*sec(x))\
          +8*tan(x)*sec(x)*sec(x)*((-2)*sin(x)-x*cos(x))\
          +6*(cos(x)-sin(x))*(2*(pow(sec(x),4))+4*(tan(x)*tan(x))*(sec(x)*sec(x)))+(cos(x)+x*sin(x))*\
          (16*(pow(sec(x),6))+88*(tan(x)*tan(x))*(pow(sec(x),4))+16*(pow(tan(x),4))*(sec(x)*sec(x))));
}
//((sec(x))^2)*(5*cos(x)-x*sin(x))+8*tan(x)*((sec(x))^2)*(4*sin(x)+x*cos(x))+4(-2*sin(x)-x*cos(x))(16*tan(x)*((sec(x))^4)+8*((tan(x))^3)*((sec(x))^2))+6(x*sin(x)-3*cos(x))(2*((sec(x))^4)+4*((tan(x))^2)*((sec(x))^2))+(cos(x)-x*sin(x))(16*((sec(x))^6)+88*((tan(x))^2)((sec(x))^4)+16*((tan(x))^4)*(sec(x))^2)
double F(double x){
	return (x/cos(x));
}

double find_max(){
	double h = 0.00001;
	double iter = B;
	double max = d4f(B);
	while(iter <= A){
		if(fabs(d4f(iter)) > max)
			max = fabs(d4f(iter));
	}
	return max;
}

double simpson(double a,double b,double h_len){
	int n = (b-a)/h_len;
	if(n % 2 == 1) {
        n++;
    }
	double h = (b-a)/n;
	double sum = f(a)+f(b);
	double s = 0;
	for(int i = 1; i <= n-1; i+=2)
		s += f(a + i*h);
	sum += 4*s;
	s = 0;
	for(int i = 2; i <= n-2; i+=2)
		s += f(a + i*h);
	return h*(sum + 2*s)/3;
}

int main(){
	double max = find_max();
	double h,val,In,I2n;
	printf("Newton-Leibniz integral equals %.5e\n", F(B)-F(A));
	printf("\nSimpson's method\n");
	printf("\n   eps\t    h\t     Value\t delta\n");
	for(double eps = 1e-3;eps >= 1e-9; eps *= 1e-3){
		h = 0.00001 * sqrt(sqrt((180*eps)/(B-A)*max));
		val = simpson(A,B,h);
		printf("  %.0e\t %.2e   %.5f\t%.2e\n", eps, h, val, fabs(val - (F(B)-F(A))));
	}
	printf("\nDouble recalculation\n");
	printf("\n   eps\t    h\t      delta\n");
	for(double eps = 1e-3;eps >= 1e-9; eps *= 1e-3){
		int n = (B-A)/sqrt(sqrt(eps));

		In = simpson(A,B,(B-A)/n);
		n *= 2;
		I2n = simpson(A,B,(B-A)/n);
		while(fabs(In - I2n) > 15*eps){
			In = I2n;
			n *= 2;
			I2n = simpson(A,B,(B-A)/n);
		}

		In = I2n;
		printf("  %.0e\t %.2e\t%.2e\n", eps, (B-A)/n, fabs(In - (F(B)-F(A))));
	}
	return 0;
}