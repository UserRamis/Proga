#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;
double a = 1.0;
double b = 1.5;
double T = 3.141592 / 2.0;
int M1 = 100;
int M2 = 100;
int N = 100;
double t = (double)T / N;
double hx = (double)a / M1;
double hy = (double)b / M2;
vector < vector < double>> U(M1 + 1, vector<double>(M2 + 1, 0));
vector < vector < double>> U_half(M1 + 1, vector<double>(M2 + 1, 0));
vector < vector < vector<double>>> f_half(M1 + 1, vector < vector < double>>(M2 + 1, vector<double>(N, 0)));

double u(double x, double y, double t);
double f(double x, double y, double t);
double phi(double x, double y);
double m1(double y, double t);
double m2(double y, double t);
double m3(double x, double t);
double m4(double x, double t);
vector<double> Progoncka(double a, double b, double c, vector<double> f);
int main() {

	double progoncka_A = 1.0;
	double progoncka_B_previous = -2 * (1 + hx * hx / t);
	double progoncka_B_half = -2 * (1 + hy * hy / t);
	double progoncka_C = 1.0;
	double hx2 = hx * hx;
	double hy2 = hy * hy;
	double temp1 = hx * hx / (hy * hy);
	double temp2 = hy * hy / (hx * hx);
	double temp3 = 2 * (1 - hy * hy / t);
	double temp4 = 2 * (1 - hx * hx / t);

	for (int k = 0; k <= M1; k++)
		for (int m = 0; m <= M2; m++)
			for (int n = 0; n < N; n++)
				f_half[k][m][n] = f(k * hx, m * hy, (n + 0.5) * t);
	for (int k = 0; k <= M1; k++)
		for (int m = 0; m <= M2; m++)
			U[k][m] = phi(k * hx, m * hy);
	// double max_error=0;
	// for(int k=0;k<=M1;k++){
	// for(int m=0;m<=M2;m++){
	// max_error=max(max_error,abs(u(k*hx,m*hy,0)-U[k][m]));
	// }
	// }
	// cout«max_error«endl;


	double max_image = 0;
	double max_in_row = 0;
	for (int n = 0; n < N; n++) {
		for (int m = 1; m < M2; m++) {
			U_half[0][m] = (m1(m * hy, (n + 1) * t) + m1(m * hy, n * t)) / 2.0
				- t * ((m1((m - 1) * hy, t * (n + 1)) - 2 * m1(m * hy, (n + 1) * t) + m1((m + 1) * hy, (n + 1) * t))
					- (m1((m - 1) * hy, n * t) - 2 * m1(m * hy, t * n) + m1((m + 1) * hy, t * n))) / (hy * hy * 4.0);
			// if(abs(U_half[0][m]-m1(m*hy,(n+0.5)*t))>max_image)
			// max_image=abs(U_half[0][m]-m1(m*hy,(n+0.5)*t));

			U_half[M1][m] = (m2(m * hy, (n + 1) * t) + m2(m * hy, n * t)) / 2.0
				- t * ((m2((m - 1) * hy, t * (n + 1)) - 2 * m2(m * hy, (n + 1) * t) + m2((m + 1) * hy, (n + 1) * t))
					- (m2((m - 1) * hy, n * t) - 2 * m2(m * hy, t * n) + m2((m + 1) * hy, t * n))) / (hy * hy * 4.0);
			// if(max_image<(abs(U_half[M1][m]-m2(m*hy,(n+0.5)*t))))
			// max_image=abs(U_half[M1][m]-m2(m*hy,(n+0.5)*t));

			vector<double> vec(M1 - 1, 0);
			for (int k = 1; k < M1; k++) {
				vec[k - 1] = temp1 * (-U[k][m - 1] + temp3 * U[k][m] - U[k][m + 1]) - hx2 * f_half[k][m][n];
			}
			vec[0] -= U_half[0][m];
			vec[M1 - 2] -= U_half[M1][m];
			vector<double> z = Progoncka(progoncka_A, progoncka_B_previous, progoncka_C, vec);
			for (int k = 1; k < M1; k++) {
				U_half[k][m] = z[k - 1];
			}
			max_in_row = 0;
			for (int k = 1; k < M1; k++) {
				if (abs(U_half[k][m] - u(k * hx, m * hy, (n + 0.5) * t)) > max_in_row) {
					max_in_row = abs(U_half[k][m] - u(k * hx, m * hy, (n + 0.5) * t));
				}
			}
			//cout«max_in_row«' ';
		}
		//cout«endl;
		for (int m = 0; m <= M2; m++) {
			U[0][m] = m1(m * hy, (n + 1) * t);
			U[M1][m] = m2(m * hy, (n + 1) * t);
		}
		for (int k = 1; k < M1; k++) {
			U[k][0] = m3(hx * k, (n + 1) * t);
			U[k][M2] = m4(hx * k, (n + 1) * t);
			vector<double> vec(M2 - 1, 0);
			for (int m = 1; m < M2; m++) {
				vec[m - 1] = temp2 * (-U_half[k - 1][m] + temp4 * U_half[k][m] - U_half[k + 1][m]) - hy2 * f_half[k][m][n];
			}
			vec[0] -= U[k][0];
			vec[M2 - 2] -= U[k][M2];
			vector<double> z = Progoncka(progoncka_A, progoncka_B_half, progoncka_C, vec);
			for (int m = 1; m < M2; m++) {
				U[k][m] = z[m - 1];
			}
			max_in_row = 0;
			for (int m = 1; m < M2; m++) {
				if (abs(U[k][m] - u(k * hx, hy * m, (n + 1) * t)) > max_in_row) {
					max_in_row = abs(U[k][m] - u(k * hx, hy * m, (n + 1) * t));
				}
			}
			//cout«max_in_row«' ';
		}

		//cout«endl«n«endl;
	}
	double maximum = 0;
	double temp;
	for (int k = 0; k <= M1; k++) {
		for (int m = 0; m <= M2; m++) {
			maximum = max (maximum, abs(u(k * hx, m * hy, N * t) - U[k][m]));
			maximum=maximum / 100000;
		}
	}
	//cout<<max_image<<endl;
	cout<<maximum;

	return 0;
}

double
m1(double y, double t) {
	return y * y + t * t;
}
double m2(double y, double t) {
	return 1 + y * y + t * t;
}
double m3(double x, double t) {
	return x * x + t * t;
}
double m4(double x, double t) {
	return 4 + x * x + t * t;
}
double phi(double x, double y) {
	return x * x + y * y;
}
double f(double x, double y, double t) {
	return 2 * t - 4;
}
vector<double> Progoncka(double a, double b, double c, vector<double>f) {
	int n = f.size();
	vector<double> answer(n, 0);
	vector<double> alpha(n + 1);
	vector<double> betta(n + 1);
	alpha[1] = (double)-c / b;
	betta[1] = (double)f[0] / b;
	double temp;
	for (int i = 1; i < n; i++) {
		temp = (b + alpha[i] * a);
		alpha[i + 1] = (double)-c / temp;
		betta[i + 1] = (double)(f[i] - betta[i] * a) / temp;
	}
	answer[n - 1] = betta[n];
	for (int i = n - 1; i > 0; i--) {
		answer[i - 1] = alpha[i] * answer[i] + betta[i];
	}
	return answer;
}

//double
//m1(double y, double t) {
//	return y * y + t * t;
//}
//double m2(double y, double t) {
//	return 1 + y * y + t * t;
//}
//double m3(double x, double t) {
//	return x * x + t * t;
//}
//double m4(double x, double t) {
//	return 4 + x * x + t * t;
//}
//double phi(double x, double y) {
//	return x * x + y * y;
//}
//double f(double x, double y, double t) {
//	return 2 * t - 4;
//}
double u(double x, double y, double t) {
	return pow(x, 2.71828) + y + cos(2 * t);
}