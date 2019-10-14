#ifndef Stokes5_H
#define Stokes5_H

#include<stdio.h>
#include<math.h>

class Stokes5 {
private:
	double A11, A22, A31, A33, A42, A44, A51, A53, A55,
		B22, B31, B42, B44, B53, B55,
		C0, C2, C4,
		D2, D4,
		E2, E4,
		S, eps, k, kd, kH, c, cs, u_mean, Q, R;

public:
	double wave_length;
	double wave_height;
	double current;
	double depth;
	double gravity;
	double x0;
	double y0;
	void set_stokes5_properties(double wave_length, double wave_height, double current, double depth, double gravity, double x0, double y0);
	double eta(double t, double X);
	double u(double t, double X, double Y);
	double v(double t, double X, double Y);

};

#endif


