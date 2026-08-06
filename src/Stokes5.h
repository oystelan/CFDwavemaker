#ifndef Stokes5_H
#define Stokes5_H

#ifndef PI
#define PI 3.14159265358979323846
#endif
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
	bool initialized = false;
	double wave_length;
	double wave_height;
	double wave_period;
	double current = 0.;
	double depth;
	double gravity = G;
	double rho = 1025.0;
	double x0 = 0.;
	double y0 = 0.;
	double z0 = 0.; // Still water line
	double t0 = 0.;
	double theta = 0.;
	void set_stokes5_properties(double _wave_length, double _wave_height);
	double eta(double t, double X, double Y);
	double u(double t, double X, double Y, double Z);
	double v(double t, double X, double Y, double Z);
	double w(double t, double X, double Y, double Z);
	double dp(double t, double X, double Y, double Z);

	// Linear dispersion-relation helpers (use this->depth and this->gravity).
	// wavenumberFromPeriod solves omega^2 = g*k*tanh(k*h) iteratively.
	double wavenumberFromPeriod(double T);
	double wavelengthFromPeriod(double T);
	double wavePeriodFromLength(double L);

};

#endif


