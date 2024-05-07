#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
struct polinom {
	double cube;
	double two;
	double one;
	double zero;
};

polinom* putData(polinom* eq, double three, double two, double one, double zero) {
	eq->cube = three;
	eq->two = two;
	eq->one = one;
	eq->zero = zero;
	return eq;
}
void renderEq(polinom eq) {
	printf("%fx^3 %fx^2 %fx %f\n", eq.cube, eq.two, eq.one, eq.zero);
}
void renderArr(double* arr, int size) {
	for (int i = 0; i < size; i++) {
		printf("%f  ", arr[i]);
	}
	printf("\n");
}
polinom* ToReduced(polinom primary) {
	polinom* newEq = (polinom*)malloc(sizeof(polinom));
	newEq->cube = 1.0;
	newEq->two = 0.0;
	double a = primary.two/primary.cube;
	double b = primary.one / primary.cube;
	double c = primary.zero / primary.cube;
	newEq->one = (b - pow(a, 2)/3);
	newEq->zero = 2 * pow(a / 3, 3) - a * b / 3 + c;
	return newEq;
}
double* solveLin(polinom pol, int& count) {
	if (pol.one != 0) {
		count = 1;
		double* res = (double*)malloc(sizeof(double));
		res[0] = -pol.zero / pol.one;
		return res;
	}
	else if (pol.zero == 0) {
		count = 4;
		return nullptr;
	}
	else {
		count = 0;
		return nullptr;
	}
}
double* solveSq(polinom pol, int& count) {
	if (pol.two == 0) {
		double* res = solveLin(pol, count);
		return res;
	}
	double a = pol.two;
	double b = pol.one;
	double c = pol.zero;
	double D = pow(b, 2) - 4 * a * c;
	if (D < 0) {
		count = 0;
		return nullptr;
	}
	else if (D == 0) {
		count = 1;
		double* res = (double*)malloc(sizeof(double));
		res[0] = -b / (2 * a);
		return res;
	}
	else {
		count = 2;
		double* res = (double*)malloc(sizeof(double));
		res[0] = (-b + D) / (2 * a);
		res[1] = (-b - D) / (2 * a);
		return res;
	}
}
double* solveEq(polinom pol, int& count) {
	if (pol.cube == 0) {
		double* res = solveSq(pol, count);
		return res;
	}
	double offset = pol.two / (3 * pol.cube);
	polinom* eq = ToReduced(pol);
	double q = eq->zero;
	double p = eq->one;
	double solution;
	double D = pow((q / 2), 2) + pow((p / 3), 3);
	D = round(D * pow(10, 7)) / pow(10, 7);
	if (D > 0) {
		count = 1;
		double* res = (double*)malloc(sizeof(double));
		solution = cbrt(-q / 2 + sqrt(D)) + cbrt(-q / 2 - sqrt(D));
		res[0] = solution - offset;
		return res;
	}
	else if (D==0){
		count = 2; 
		double* res = (double*)malloc(sizeof(double)*2);
		double res1 = 2 * (double)cbrt(-q / 2);
		double res2 = -cbrt(-q / 2);
		res[0] = res1-offset;
		res[1] = res2-offset;
		return res;
	}
	else {
		count = 3;
		double* res = (double*)malloc(sizeof(double) * 3);
		double fi = 0;
		if (q < 0) {
			fi = atan(sqrt(-D) / (q / 2));
		}
		else if (q > 0) {
			fi = atan(sqrt(-D) / (q / 2)) + M_PI;
		}
		else {
			fi = M_PI / 2;
		}
		res[0] = 2 * sqrt(-p / 3) * cos(fi / 3) - offset;
		res[1] = 2 * sqrt(-p / 3) * cos(fi / 3 + 2.0 / 3 * M_PI) - offset;
		res[2] = 2 * sqrt(-p / 3) * cos(fi / 3 + 4.0 / 3 * M_PI) - offset;
		return res;
	}
}
double* get_h(double* h, double* x, int size) {
	for (int i = 0; i < size - 1; i++) {
		h[i] = x[i + 1] - x[i];
	}
	return h;
}
double* progon(double* Y, double* h, double* y, int size) {
	double* a = (double*)malloc(size * sizeof(double));
	double* b = (double*)malloc(size * sizeof(double));
	double* c = (double*)malloc(size * sizeof(double));
	double* f = (double*)malloc(size * sizeof(double));
	for (int i = 1; i < size-1; i++) {
		a[i] = h[i-1] / 6.0;
		b[i] = (h[i-1] + h[i]) / 3.0;
		c[i] = h[i + 1] / 6.0;
		f[i] = (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i-1];
	}
	a[1] = 0.0;
	c[size - 2] = 0.0;
	double* A = (double*)malloc(size * sizeof(double));
	double* B = (double*)malloc(size * sizeof(double));
	A[1] = -c[1] / b[1];
	B[1] = f[1] / b[1];
	for (int i = 2; i < size-1; i++) {
		A[i] = (-c[i] / (b[i] + a[i] * A[i - 1]));
		B[i] = (f[i] - a[i] * B[i - 1]) / (b[i] + a[i] * A[i - 1]);
	}
	Y[size - 1] = 0;
	Y[size - 2] = B[size - 2];
	for (int i = size - 3; i >=1; i--) {
		Y[i] = A[i] * Y[i + 1] + B[i];
	}
	Y[0] = 0;
	free(a);
	free(b);
	free(c);
	free(f);
	return Y;
}
polinom* createSpline(polinom* spline, double* x, double* y, double* h, double* Y, int size) {
	for (int i = 0; i < size - 1; i++) {
		spline[i].cube = (Y[i+1] - Y[i]) / (6 * h[i]);
		spline[i].two = 3 * (x[i + 1] * Y[i] - Y[i+1] * x[i]) / (6 * h[i]);
		spline[i].one = (y[i + 1] - y[i]) / h[i] + \
			(Y[i] * (h[i] * h[i] - 3 * x[i + 1] * x[i + 1]) - Y[i + 1] * (h[i] * h[i] * x[i] - 3 * x[i] * x[i])) / (6 * h[i]);
		spline[i].zero = (y[i] * x[i + 1] - y[i + 1] * x[i]) / h[i] + \
			(Y[i] * (-h[i] * h[i] * x[i + 1] + x[i + 1] * x[i + 1] * x[i + 1]) + Y[i + 1] * (h[i] * h[i] * x[i] - x[i] * x[i] * x[i])) / (6 * h[i]);
	}
	return spline;
}
double min(double a, double b) {
	if (a > b) {
		return b;
	}
	return a;
}
bool solveSpline(polinom* spl1, polinom* spl2, double* h1, double* h2, double begin, double* sol) {
	int posX1, posX2, ind1, ind2, left;
	int count = 0;
	left = posX1 = posX2 = begin;
	ind1 = ind2 = 0;
	polinom* temp=(polinom*)malloc(sizeof(polinom));
	double* res = (double*)malloc(3 * sizeof(double));
	while (h1[ind1] && h2[ind2]) {
		res = solveEq(*putData(temp, spl1->cube - spl2->cube, spl1->two - spl2->two, spl1->one - spl2->one, spl1->zero - spl2->zero), count);
		if (count == 4) {
			printf("Some splines are eq");
			return false;
		}
		if (count > 0) {
			for (int i = 0; i < count;i++) {
				if (res[i] >= left && res[i] <= min(posX1 + h1[ind1], posX2 + h2[ind2])) {
					*sol = res[i];
					printf("The result: %lf", *sol);
					return true;
				}
			}
		}
		if (posX1 + h1[ind1] > posX2 + h2[ind2]) {
			posX2 += h2[ind2];
			ind2++;
			spl2++;
			left = posX2;
		}
		else {
			posX1 += h1[ind1];
			ind1++;
			spl1++;
			left = posX1;
		}
	}
	printf("No roots");
	return false;
}
void main() {
	int size;
	scanf("%d", &size);
	double* x1 = (double*)malloc(size * sizeof(double));
	double* y1 = (double*)malloc(size * sizeof(double));
	double* x2 = (double*)malloc(size * sizeof(double));
	double* y2 = (double*)malloc(size * sizeof(double));
	for (int i = 0; i < size; i++) {
		scanf("%lf %lf", &(x1[i]), &(y1[i]));
	}
	for (int i = 0; i < size; i++) {
		scanf("%lf %lf", &(x2[i]), &(y2[i]));
	}
	double* h1 = (double*)malloc((size) * sizeof(double));
	double* Y1 = (double*)malloc(size * sizeof(double));
	double* h2 = (double*)malloc((size) * sizeof(double));
	double* Y2 = (double*)malloc(size * sizeof(double));
	h1[size] = 0.0;
	h2[size] = 0.0;
	get_h(h1, x1, size);
	get_h(h2, x2, size);
	progon(Y1, h1, y1, size);
	progon(Y2, h2, y2, size);
	polinom* spline1 = (polinom*)malloc(size * sizeof(polinom));
	polinom* spline2 = (polinom*)malloc(size * sizeof(polinom));
	spline1 = createSpline(spline1, x1, y1, h1, Y1, size);
	spline2 = createSpline(spline2, x2, y2, h2, Y2, size);
	for (int i = 0; i < size - 1; i++) {
		renderEq(spline1[i]);
	}
	for (int i = 0; i < size - 1; i++) {
		renderEq(spline2[i]);
	}
	double res;
	bool total = solveSpline(spline1, spline2, h1, h2, x1[0], &res);
}