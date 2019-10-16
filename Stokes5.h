#ifndef Stokes5_H
#define Stokes5_H

#define PI 3.14159265358979323846
#define G 9.81

class Stokes5 {
private:
	double A11, A22, A31, A33, A42, A44, A51, A53, A55,
		B22, B31, B42, B44, B53, B55,
		C0, C2, C4,
		D2, D4,
		E2, E4,
		S, eps, k, kd, kH, c, cs, u_mean, Q, R;
	void calculate_stokes_coefficients();

public:
	double wave_length;
	double wave_height;
	double current = 0.;
	double depth;
	double gravity = G;
	double x0;
	double y0;
	double z0 = 0.; // Still water line
	double theta = 0.;
	void set_stokes5_properties(double _wave_length, double _wave_height, double _depth, double _x0, double _y0);
	double eta(double t, double X, double Y);
	double u(double t, double X, double Y, double Z);
	double v(double t, double X, double Y, double Z);
	double w(double t, double X, double Y, double Z);

};

#endif


