/*
 * Maslov Vadim, KV-51
 * 23. sin(2x)+(cos(x))^2-x/2=0, (І,Х)
*/

#include <iostream>
#include <math.h>
#include <unistd.h>

using std::swap;

double func(double x){
	return sin(2*x)+cos(x)*cos(x)-x/2;
}

double dfunc(double x){
	return -2*sin(x)*cos(x)+2*cos(2*x)-0.5;
}

double d2func(double x){
	return 2*sin(x)*sin(x)-4*sin(2*x)-2*cos(x)*cos(x);
}

double phi(double x, double lambda){
	return dfunc(x) > 0.0 ? x - lambda*func(x) : x + lambda*func(x);
}

double iteratable(double a, double b, double eps, double &delta,int &i){
	double m1 = fabs(dfunc(b));
	double M1 = fabs(dfunc(a));
	if(m1 > M1){
		swap(m1,M1);
	}
	double x0, xk = (a+b)/2;
	double lambda = 1.0/M1;
	double q = 1.0 - m1/M1;
	i = 0;
	do{
		x0 = xk;
		xk = phi(x0, lambda);
		i++;
	}while(fabs(xk-x0) > (1-q)/q*eps);
	delta = fabs(xk-x0) * q/(1-q);

	return xk;
}

double chords(double a,double b, double eps, double &delta,int &i){
	double m1 = fabs(dfunc(b));
	double M1 = fabs(dfunc(a));
	if(m1 > M1){
		swap(m1,M1);
	}
	double xk = a;
	double c = b;
	if( d2func(a)*func(a) > 0){
		xk = b;
		c = a;
	}
	i = 0;
	while(fabs(func(xk)/m1) > eps){
		xk -= func(xk)*(xk-c)/(func(xk)-func(c));
		i++;
	};
	delta = fabs(func(xk)/m1);

	return xk;
}

void printIteratable(double a, double b){
	double x;
	double delta;
	int i;
	printf("[%.2f;%.2f], iteratable\neps\tX\t\t\tdelta\n", a, b);
	for (double eps = 1e-2; eps > 1e-14; eps *= 1e-3){
		x = iteratable(a,b,eps,delta,i);
		printf("%.0e\t%.14f\t%.6e\n", eps, x, delta);
	}

}

void printChords(double a, double b){
	double x;
	double delta;
	int i;
	printf("[%.2f;%.2f], chords\neps\tX\t\t\tdelta\n", a, b);
	for (double eps = 1e-2; eps > 1e-14; eps *= 1e-3){
		x = chords(a,b,eps,delta,i);
		printf("%.0e\t%.14f\t%.6e\t %d\n", eps, x, delta, i);
	}
}

void comparationTable(double a, double b){
	double delta;
	int i1,i2;
	printf("[%.2f;%.2f], compare\neps\tIter\tChords\n", a, b);
	for (double eps = 1e-2; eps > 1e-14; eps *= 1e-3){
        iteratable(a, b, eps, delta, i1);
        chords(a, b, eps, delta, i2);
        printf("%.0e\t%d\t%d\n", eps, i1, i2);
    }
}

int main(){

    printIteratable(-1.25, -1.05); printf("\n");
    printChords(-1.25, -1.05); printf("\n");

    printIteratable(-0.74, -0.54); printf("\n");
    printChords(-0.74, -0.54); printf("\n");

    printIteratable(1.23, 1.33); printf("\n");
    printChords(1.23, 1.33); printf("\n");

    comparationTable(-1.25, -1.05);

	return 0;
}
