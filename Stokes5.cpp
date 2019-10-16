
#include "Stokes5.h"
#include <math.h>



/* Stokes 5th order waves based on Fenton's paper */
void Stokes5::set_stokes5_properties(double _wave_length, double _wave_height, double _depth, double _x0, double _y0)

/* set minimum required coefficients to calculate stokes wave*/
{
	// compute basic properties, dimensional numbers, and set attributes of stokes5  
	k = 2. * PI / wave_length;
	kd = k * depth;
	eps = k * wave_height / 2.;
	S = 1. / cosh(2. * kd);

	wave_height = _wave_height;
	wave_length = _wave_length;
	depth = _depth;
	kH = k * _wave_height;
	x0 = _x0;
	y0 = _y0;
}

/* Compute coefficients based on wave properties and depth */
void Stokes5::calculate_stokes_coefficients(){

	// compute coefficients
	A11 = 1 / sinh(kd);
	A22 = 3 * pow(S, 2) / (2 * pow(1 - S, 2));
	A31 = (-4 - 20 * S + 10 * pow(S, 2) - 13 * pow(S, 3)) / (8 * sinh(kd) * pow(1 - S, 3));
	A33 = (-2 * pow(S, 2) + 11 * pow(S, 3)) / (8 * sinh(kd) * pow(1 - S, 3));
	A42 = (12 * S - 14 * pow(S, 2) - 264 * pow(S, 3) - 45 * pow(S, 4) - 13 * pow(S, 5)) / (24 * pow(1 - S, 5));
	A44 = (10 * pow(S, 3) - 174 * pow(S, 4) + 291 * pow(S, 5) + 278 * pow(S, 6)) / (48 * (3 + 2 * S) * pow(1 - S, 5));
	A51 = (-1184 + 32 * S + 13232 * pow(S, 2) + 21712 * pow(S, 3) + 20940 * pow(S, 4) + 12554 * pow(S, 5) - 500 * pow(S, 6) - 3341 * pow(S, 7) - 670 * pow(S, 8)) / (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));
	A53 = (4 * S + 105 * pow(S, 2) + 198 * pow(S, 3) - 1376 * pow(S, 4) - 1302 * pow(S, 5) - 117 * pow(S, 6) + 58 * pow(S, 7)) / (32 * sinh(kd) * (3 + 2 * S) * (pow(1 - S, 6)));
	A55 = (-6 * pow(S, 3) + 272 * pow(S, 4) - 1552 * pow(S, 5) + 852 * pow(S, 6) + 2029 * pow(S, 7) + 430 * pow(S, 8)) / (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));

	B22 = (cosh(kd) / sinh(kd)) * (1 + 2 * S) / (2 * (1 - S));
	B31 = -3 * (1 + 3 * S + 3 * pow(S, 2) + 2 * pow(S, 3)) / (8 * pow(1 - S, 3));
	B42 = (cosh(kd) / sinh(kd)) * (6 - 26 * S - 182 * pow(S, 2) - 204 * pow(S, 3) - 25 * pow(S, 4) + 26 * pow(S, 5)) / (6 * (3 + 2 * S) * pow(1 - S, 4));
	B44 = (cosh(kd) / sinh(kd)) * (24 + 92 * S + 122 * pow(S, 2) + 66 * pow(S, 3) + 67 * pow(S, 4) + 34 * pow(S, 5)) / (24 * (3 + 2 * S) * pow(1 - S, 4));
	B53 = 9 * (132 + 17 * S - 2216 * pow(S, 2) - 5897 * pow(S, 3) - 6292 * pow(S, 4) - 2687 * pow(S, 5) + 194 * pow(S, 6) + 467 * pow(S, 7) + 82 * pow(S, 8)) / (128 * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));
	B55 = 5 * (300 + 1579 * S + 3176 * pow(S, 2) + 2949 * pow(S, 3) + 1188 * pow(S, 4) + 675 * pow(S, 5) + 1326 * pow(S, 6) + 827 * pow(S, 7) + 130 * pow(S, 8)) / (384 * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));
	C0 = pow(tanh(kd), 0.5);
	C2 = pow(tanh(kd), 0.5) * (2 + 7 * pow(S, 2)) / (4 * pow(1 - S, 2));
	C4 = pow(tanh(kd), 0.5) * (4 + 32 * S - 116 * pow(S, 2) - 400 * pow(S, 3) - 71 * pow(S, 4) + 146 * pow(S, 5)) / (32 * pow(1 - S, 5));
	D2 = -pow(cosh(kd) / sinh(kd), 0.5) / 2;
	D4 = pow(cosh(kd) / sinh(kd), 0.5) * (2 + 4 * S + pow(S, 2) + 2 * pow(S, 3)) / (8 * pow(1 - S, 3));
	E2 = tanh(kd) * (2 + 2 * S + 5 * pow(S, 2)) / (4 * pow(1 - S, 2));
	E4 = tanh(kd) * (8 + 12 * S - 152 * pow(S, 2) - 308 * pow(S, 3) - 42 * pow(S, 4) + 77 * pow(S, 5)) / (32 * pow(1 - S, 5));

	u_mean = (C0 + pow(eps, 2) * C2 + pow(eps, 4) * C4) / sqrt(k / gravity);
	c = u_mean + current; // wave celerity
	Q = (u_mean) * depth + (D2 * pow(eps, 2) + D4 * pow(eps, 4)) / sqrt(pow(k, 3) / gravity);
	cs = c - (Q) / depth; // cs = Stokes drift (set the current to Q/d (Stokes drift=0) for closed wave tanks
	R = 0.5 * (gravity / k) * pow(C0, 2) + kd + (E2) * pow(eps, 2) + (E4) * pow(eps, 4); // Bernoulli constant

}


double Stokes5::eta(double t, double X, double Y)

/* Given an initialized stokes wave (stokes5), compute the surface elevation eta at time t and coordinate X */
{
	double keta, kx;

	//kx = (k) * (X - (c) * t - (x0)); // C is case sensitive
	kx = (k) * (cos(theta * PI / 180.)*(X - (c)*t - (x0)) + sin(theta * PI / 180.)*(Y - (c)*t - (y0))); // C is case sensitive
	keta = kd + eps * cos(kx) +
		pow(eps, 2) * (B22) * cos(2 * kx) +
		pow(eps, 3) * (B31) * (cos(kx) - cos(3 * kx)) +
		pow(eps, 4) * ((B42) * cos(2 * kx) + (B44) * cos(4 * kx)) +
		pow(eps, 5) * (-((B53) + (B55)) * cos(kx) + (B53) * cos(3 * kx) + (B55) * cos(5 * kx));

	return keta / (k) - depth;
}



double Stokes5::u(double t, double X, double Y, double Z)

/* Given an initialized stokes wave (stokes5), compute the horizontal velocity at time t and coordinate X, Y.  Values above the surface are computed uncritically. */
{

	double kx, ky;

	//kx = (k) * (X - (c) * t - (x0));
	kx = (k) * (cos(theta * PI / 180.) * (X - (c)*t - (x0)) + sin(theta * PI / 180.) * (Y - (c)*t - (y0))); // C is case sensitive
	ky = (k) * (Z + depth - (z0));

	return cos(theta * PI / 180.) * (c - u_mean +
		(C0) * sqrt((gravity) / pow(k, 3)) * (pow(eps, 1) * (1 * k * (A11) * cosh(1 * ky) * cos(1 * kx)) +
			pow(eps, 2) * (2 * k * (A22) * cosh(2 * ky) * cos(2 * kx)) +
			pow(eps, 3) * (1 * k * (A31) * cosh(1 * ky) * cos(1 * kx) + 3 * k * (A33) * cosh(3 * ky) * cos(3 * kx)) +
			pow(eps, 4) * (2 * k * (A42) * cosh(2 * ky) * cos(2 * kx) + 4 * k * (A44) * cosh(4 * ky) * cos(4 * kx)) +
			pow(eps, 5) * (1 * k * (A51) * cosh(1 * ky) * cos(1 * kx) + 3 * k * (A53) * cosh(3 * ky) * cos(3 * kx) +
				5 * k * (A55) * cosh(5 * ky) * cos(5 * kx))));
}




double Stokes5::v(double t, double X, double Y, double Z)

/* Given an initialized stokes wave (stokes5), compute the horizontal velocity at time t and coordinate X, Z.  Values above the surface are computed uncritically. */
{

	double  kx, ky;
	
	//kx = (k) * (X - (c) * t - (x0));
	kx = (k) * (cos(theta * PI / 180.) * (X - (c)*t - (x0)) + sin(theta * PI / 180.) * (Y - (c)*t - (y0))); // C is case sensitive
	ky = (k) * (Z + depth - (z0));

	return sin(theta*PI/180.) * (c - u_mean +
		(C0) * sqrt((gravity) / pow(k, 3)) * (pow(eps, 1) * (1 * k * (A11) * cosh(1 * ky) * cos(1 * kx)) +
			pow(eps, 2) * (2 * k * (A22) * cosh(2 * ky) * cos(2 * kx)) +
			pow(eps, 3) * (1 * k * (A31) * cosh(1 * ky) * cos(1 * kx) + 3 * k * (A33) * cosh(3 * ky) * cos(3 * kx)) +
			pow(eps, 4) * (2 * k * (A42) * cosh(2 * ky) * cos(2 * kx) + 4 * k * (A44) * cosh(4 * ky) * cos(4 * kx)) +
			pow(eps, 5) * (1 * k * (A51) * cosh(1 * ky) * cos(1 * kx) + 3 * k * (A53) * cosh(3 * ky) * cos(3 * kx) +
				5 * k * (A55) * cosh(5 * ky) * cos(5 * kx))));
}

double Stokes5::w(double t, double X, double Y, double Z)

/* Given an initialized stokes wave (stokes5), compute the vertical velocity at time t and coordinate X, Z.  Values above the surface are computed uncritically. */
{

	double kx, ky;
	
	//kx = (k) * (X - (c) * t - (x0));
	kx = (k) * (cos(theta * PI / 180.) * (X - (c)*t - (x0)) + sin(theta * PI / 180.) * (Y - (c)*t - (y0))); // C is case sensitive
	ky = (k) * (Z + depth - (z0));

	return (C0) * sqrt((gravity) / pow(k, 3)) * (pow(eps, 1) * (1 * k * (A11) * sinh(1 * ky) * sin(1 * kx)) +
		pow(eps, 2) * (2 * k * (A22) * sinh(2 * ky) * sin(2 * kx)) +
		pow(eps, 3) * (1 * k * (A31) * sinh(1 * ky) * sin(1 * kx) + 3 * k * (A33) * sinh(3 * ky) * sin(3 * kx)) +
		pow(eps, 4) * (2 * k * (A42) * sinh(2 * ky) * sin(2 * kx) + 4 * k * (A44) * sinh(4 * ky) * sin(4 * kx)) +
		pow(eps, 5) * (1 * k * (A51) * sinh(1 * ky) * sin(1 * kx) + 3 * k * (A53) * sinh(3 * ky) * sin(3 * kx) +
			5 * k * (A55) * sinh(5 * ky) * sin(5 * kx)));
}