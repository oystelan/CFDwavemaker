
#include "Stokes5.h"
#include <math.h>

#define PI 3.14159265358979323846

/* Stokes 5th order waves based on Fenton's paper */



void Stokes5::set_stokes5_properties(double wave_length, double wave_height, double current, double depth, double gravity, double x0, double y0)

/* Compute coefficients based on wave properties and depth */
{

	double eps, S, k, kd;


	// compute basic properties, dimensional numbers, and set attributes of stokes5  
	k = 2. * PI / wave_length;
	kd = k * depth;
	eps = k * wave_height / 2.;
	S = 1. / cosh(2. * kd);

	Stokes5::k = k;
	Stokes5::kd = kd;
	Stokes5::wave_height = wave_height;
	Stokes5::wave_length = wave_length;
	Stokes5::current = current;
	Stokes5::depth = depth;
	Stokes5::gravity = gravity;
	Stokes5::kH = k * wave_height;
	Stokes5::eps = eps;
	Stokes5::S = S;
	Stokes5::x0 = x0;
	Stokes5::y0 = y0;

	// compute coefficients
	Stokes5::A11 = 1 / sinh(kd);
	Stokes5::A22 = 3 * pow(S, 2) / (2 * pow(1 - S, 2));
	Stokes5::A31 = (-4 - 20 * S + 10 * pow(S, 2) - 13 * pow(S, 3)) / (8 * sinh(kd) * pow(1 - S, 3));
	Stokes5::A33 = (-2 * pow(S, 2) + 11 * pow(S, 3)) / (8 * sinh(kd) * pow(1 - S, 3));
	Stokes5::A42 = (12 * S - 14 * pow(S, 2) - 264 * pow(S, 3) - 45 * pow(S, 4) - 13 * pow(S, 5)) / (24 * pow(1 - S, 5));
	Stokes5::A44 = (10 * pow(S, 3) - 174 * pow(S, 4) + 291 * pow(S, 5) + 278 * pow(S, 6)) / (48 * (3 + 2 * S) * pow(1 - S, 5));
	Stokes5::A51 = (-1184 + 32 * S + 13232 * pow(S, 2) + 21712 * pow(S, 3) + 20940 * pow(S, 4) + 12554 * pow(S, 5) - 500 * pow(S, 6) - 3341 * pow(S, 7) - 670 * pow(S, 8)) / (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));
	Stokes5::A53 = (4 * S + 105 * pow(S, 2) + 198 * pow(S, 3) - 1376 * pow(S, 4) - 1302 * pow(S, 5) - 117 * pow(S, 6) + 58 * pow(S, 7)) / (32 * sinh(kd) * (3 + 2 * S) * (pow(1 - S, 6)));
	Stokes5::A55 = (-6 * pow(S, 3) + 272 * pow(S, 4) - 1552 * pow(S, 5) + 852 * pow(S, 6) + 2029 * pow(S, 7) + 430 * pow(S, 8)) / (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));

	Stokes5::B22 = (cosh(kd) / sinh(kd)) * (1 + 2 * S) / (2 * (1 - S));
	Stokes5::B31 = -3 * (1 + 3 * S + 3 * pow(S, 2) + 2 * pow(S, 3)) / (8 * pow(1 - S, 3));
	Stokes5::B42 = (cosh(kd) / sinh(kd)) * (6 - 26 * S - 182 * pow(S, 2) - 204 * pow(S, 3) - 25 * pow(S, 4) + 26 * pow(S, 5)) / (6 * (3 + 2 * S) * pow(1 - S, 4));
	Stokes5::B44 = (cosh(kd) / sinh(kd)) * (24 + 92 * S + 122 * pow(S, 2) + 66 * pow(S, 3) + 67 * pow(S, 4) + 34 * pow(S, 5)) / (24 * (3 + 2 * S) * pow(1 - S, 4));
	Stokes5::B53 = 9 * (132 + 17 * S - 2216 * pow(S, 2) - 5897 * pow(S, 3) - 6292 * pow(S, 4) - 2687 * pow(S, 5) + 194 * pow(S, 6) + 467 * pow(S, 7) + 82 * pow(S, 8)) / (128 * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));
	Stokes5::B55 = 5 * (300 + 1579 * S + 3176 * pow(S, 2) + 2949 * pow(S, 3) + 1188 * pow(S, 4) + 675 * pow(S, 5) + 1326 * pow(S, 6) + 827 * pow(S, 7) + 130 * pow(S, 8)) / (384 * (3 + 2 * S) * (4 + S) * pow(1 - S, 6));
	Stokes5::C0 = pow(tanh(kd), 0.5);
	Stokes5::C2 = pow(tanh(kd), 0.5) * (2 + 7 * pow(S, 2)) / (4 * pow(1 - S, 2));
	Stokes5::C4 = pow(tanh(kd), 0.5) * (4 + 32 * S - 116 * pow(S, 2) - 400 * pow(S, 3) - 71 * pow(S, 4) + 146 * pow(S, 5)) / (32 * pow(1 - S, 5));
	Stokes5::D2 = -pow(cosh(kd) / sinh(kd), 0.5) / 2;
	Stokes5::D4 = pow(cosh(kd) / sinh(kd), 0.5) * (2 + 4 * S + pow(S, 2) + 2 * pow(S, 3)) / (8 * pow(1 - S, 3));
	Stokes5::E2 = tanh(kd) * (2 + 2 * S + 5 * pow(S, 2)) / (4 * pow(1 - S, 2));
	Stokes5::E4 = tanh(kd) * (8 + 12 * S - 152 * pow(S, 2) - 308 * pow(S, 3) - 42 * pow(S, 4) + 77 * pow(S, 5)) / (32 * pow(1 - S, 5));

	Stokes5::u_mean = (Stokes5::C0 + pow(eps, 2) * Stokes5::C2 + pow(eps, 4) * Stokes5::C4) / sqrt(k / gravity);
	Stokes5::c = Stokes5::u_mean + Stokes5::current; // wave celerity
	Stokes5::Q = (Stokes5::u_mean) * depth + (Stokes5::D2 * pow(eps, 2) + Stokes5::D4 * pow(eps, 4)) / sqrt(pow(k, 3) / gravity);
	Stokes5::cs = Stokes5::c - (Stokes5::Q) / depth; // cs = Stokes drift (set the current to Q/d (Stokes drift=0) for closed wave tanks
	Stokes5::R = 0.5 * (gravity / k) * pow(Stokes5::C0, 2) + kd + (Stokes5::E2) * pow(eps, 2) + (Stokes5::E4) * pow(eps, 4); // Bernoulli constant

}


double Stokes5::eta(double t, double X)

/* Given an initialized stokes wave (stokes5), compute the surface elevation eta at time t and coordinate X */
{
	double keta, kx, eps;
	eps = Stokes5::eps;
	kx = (Stokes5::k) * (X - (Stokes5::c) * t - (Stokes5::x0)); // C is case sensitive
	keta = Stokes5::kd + eps * cos(kx) +
		pow(eps, 2) * (Stokes5::B22) * cos(2 * kx) +
		pow(eps, 3) * (Stokes5::B31) * (cos(kx) - cos(3 * kx)) +
		pow(eps, 4) * ((Stokes5::B42) * cos(2 * kx) + (Stokes5::B44) * cos(4 * kx)) +
		pow(eps, 5) * (-((Stokes5::B53) + (Stokes5::B55)) * cos(kx) + (Stokes5::B53) * cos(3 * kx) + (Stokes5::B55) * cos(5 * kx));

	return keta / (Stokes5::k) - Stokes5::depth;
}



double Stokes5::u(double t, double X, double Y)

/* Given an initialized stokes wave (stokes5), compute the horizontal velocity at time t and coordinate X, Y.  Values above the surface are computed uncritically. */
{

	double k, kx, ky, eps;
	eps = Stokes5::eps;
	k = Stokes5::k;
	kx = (Stokes5::k) * (X - (Stokes5::c) * t - (Stokes5::x0));
	ky = (Stokes5::k) * (Y + Stokes5::depth - (Stokes5::y0));

	return Stokes5::c - Stokes5::u_mean +
		(Stokes5::C0) * sqrt((Stokes5::gravity) / pow(k, 3)) * (pow(eps, 1) * (1 * k * (Stokes5::A11) * cosh(1 * ky) * cos(1 * kx)) +
			pow(eps, 2) * (2 * k * (Stokes5::A22) * cosh(2 * ky) * cos(2 * kx)) +
			pow(eps, 3) * (1 * k * (Stokes5::A31) * cosh(1 * ky) * cos(1 * kx) + 3 * k * (Stokes5::A33) * cosh(3 * ky) * cos(3 * kx)) +
			pow(eps, 4) * (2 * k * (Stokes5::A42) * cosh(2 * ky) * cos(2 * kx) + 4 * k * (Stokes5::A44) * cosh(4 * ky) * cos(4 * kx)) +
			pow(eps, 5) * (1 * k * (Stokes5::A51) * cosh(1 * ky) * cos(1 * kx) + 3 * k * (Stokes5::A53) * cosh(3 * ky) * cos(3 * kx) +
				5 * k * (Stokes5::A55) * cosh(5 * ky) * cos(5 * kx)));
}




double Stokes5::v(double t, double X, double Y)

/* Given an initialized stokes wave (stokes5), compute the vertical velocity at time t and coordinate X, Y.  Values above the surface are computed uncritically. */
{

	double k, kx, ky, eps;
	eps = Stokes5::eps;
	k = Stokes5::k;
	kx = (Stokes5::k) * (X - (Stokes5::c) * t - (Stokes5::x0));
	ky = (Stokes5::k) * (Y + Stokes5::depth - (Stokes5::y0));

	return (Stokes5::C0) * sqrt((Stokes5::gravity) / pow(k, 3)) * (pow(eps, 1) * (1 * k * (Stokes5::A11) * sinh(1 * ky) * sin(1 * kx)) +
		pow(eps, 2) * (2 * k * (Stokes5::A22) * sinh(2 * ky) * sin(2 * kx)) +
		pow(eps, 3) * (1 * k * (Stokes5::A31) * sinh(1 * ky) * sin(1 * kx) + 3 * k * (Stokes5::A33) * sinh(3 * ky) * sin(3 * kx)) +
		pow(eps, 4) * (2 * k * (Stokes5::A42) * sinh(2 * ky) * sin(2 * kx) + 4 * k * (Stokes5::A44) * sinh(4 * ky) * sin(4 * kx)) +
		pow(eps, 5) * (1 * k * (Stokes5::A51) * sinh(1 * ky) * sin(1 * kx) + 3 * k * (Stokes5::A53) * sinh(3 * ky) * sin(3 * kx) +
			5 * k * (Stokes5::A55) * sinh(5 * ky) * sin(5 * kx)));
}