/*
 *
 * Lab6, "НАБЛИЖЕННЯ ФУНКЦІЙ ЗА ДОПОМОГОЮ ІНТЕРПОЛЯЦІЙНИХ СПЛАЙНІВ"
 * Маслов Вадим, KV-51
 * 23. plot e^((x+5)^1/4)*sin(3x)*cos(x), [1;6.5]
 * (23 % 6 = 5 => 101). Метод ітерацій Зейделя
*/

#include <iostream>
#include <math.h>
#include <cstring>
#include <fstream>

#define A 1
#define B 6.5
#define M 60 
#define N (M+1)

#define eps 1e-3

typedef struct
{
    double f;
    double s;
    double t;
    double ft;
    double x;
} Spline;

double f(double x){
    return (pow(M_E,pow(x+5,1/4))*sin(3*x)*cos(x));
}

void get_matrix(double matrix[][N + 1], double a, double b)
{
    const double h = (b - a)/M;
    double x = a + h;

    memset(matrix, 0, N*(N + 1)*sizeof(double));

    for (int i = 1; i < N - 1; i++, x += h) {
        matrix[i][i - 1] = h; // A
        matrix[i][i] = 4*h; // C
        matrix[i][i + 1] = h; // B
        matrix[i][N] = 6*(f(x - h) - 2*f(x) + f(x + h))/h; // F
    }
    matrix[0][0] = 1.0;
    matrix[N - 1][N - 1] = 1.0;
}

void seidel_iter(double matr[][N + 1], double Xk[N])
{
    double alpha[N][N];
    double beta[N], Xprev[N];
    double q = 0.0;
    double norm = 0.0;

    for (int i = 0; i < N; i++) {
        beta[i] = matr[i][N]/matr[i][i];
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                alpha[i][j] = -matr[i][j]/matr[i][i];
            } else {
                alpha[i][j] = 0;
            }
        }
    }

    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < N; j++) {
            sum += fabs(alpha[i][j]);
        }
        if (sum > q) {
            q = sum;
        }
    }

    for (int i = 0; i < N; i++) {
        Xprev[i] = beta[i];
        for (int j = 0; j < i; j++) {
            Xprev[i] += alpha[i][j]*Xprev[j];
        }
        for (int j = i; j < N; j++) {
            Xprev[i] += alpha[i][j]*matr[j][N];
        }
    }

    while (true) {
        for (int i = 0; i < N; i++) {
            Xk[i] = beta[i];
            for (int j = 0; j < i; j++) {
                Xk[i] += alpha[i][j]*Xk[j];
            }
            for (int j = i; j < N; j++) {
                Xk[i] += alpha[i][j]*Xprev[j];
            }
        }

        norm = 0.0;

        for (int i = 0; i < N; i++) {
            norm += (Xk[i] - Xprev[i])*(Xk[i] - Xprev[i]);
        }
        norm = sqrt(norm);

        if (norm <= eps*(1 - q)/q) {
            break;
        }

        for (int i = 0; i < N; i++) {
            Xprev[i] = Xk[i];
        }
    }
}

void get_splines(double a, double b, const double c[N], Spline splines[N])
{
    const double h = (b - a)/M;
    double x = a;

    splines[0].f = f(x);
    splines[0].x = x;
    for (int i = 1; i < N; i++, x += h) {
        splines[i].f = f(x);
        splines[i].t = c[i];
        splines[i].ft = (c[i] - c[i - 1])/h;
        splines[i].s = h*c[i]/2 - h*h*splines[i].ft/6 + (splines[i].f - splines[i - 1].f)/h;
        splines[i].x = x;
    }
}

double get_spline(Spline spline, double x)
{   
    double ret = spline.f + spline.s*(x - spline.x) + spline.t/2*(x - spline.x)*(x - spline.x) + spline.ft/6*(x - spline.x)*(x - spline.x)*(x - spline.x);
    return ret;
}

int main()
{
    double matrix[N][N + 1], vector[N];
    Spline splines[N];

    get_matrix(matrix, A, B);
    seidel_iter(matrix, vector);
    get_splines(A, B, vector, splines);

    std::ofstream file("table.csv", std::ofstream::out);

    int i = 0;
    const double h = (double)(B - A)/M;
    double y;
    for (double x = A; x <= B; x += h, i++) {
        y = get_spline(splines[i], x);
        file << x << ";" << y << std::endl;
        printf("%.4f, %f\n", x, y);
    }

    file.close();

    std::cout << "Exported table.csv" << std::endl;

    return 0;
}
