#include "Irregular.h"
#include <math.h>

#define PI 3.1415926535897
#define G 9.81
#define RHO 1025.0



double Irregular::u(double t, double xx, double yy, double zz) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
		usum += cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * Irregular::Ampspec[i] * Irregular::D[i] * Irregular::omega[i] * (cosh(Irregular::k[i]
			* (zz + Irregular::depth)) / sinh(Irregular::k[i] * Irregular::depth)) * cos(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) 
				* (xx - Irregular::fpoint[0]) + sin(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
	}

	return usum;

}


/* horizontal velocity U for a sinus Irregular::omegaave */
double Irregular::v(double t, double xx, double yy, double zz) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * tofmax + Irregular::phi[i];
		vsum += sin(thetaA[i] + (mtheta * PI / 180.)) * Ampspec[i] * D[i] * Irregular::omega[i] * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth)) * cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - Irregular::omega[i] * t + phi);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double Irregular::w(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;



	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {

		phi = Irregular::omega[i] * tofmax + Irregular::phi[i];
		wsum += D[i] * Ampspec[i] * Irregular::omega[i] * (sinh(k[i] * (zz + depth)) / sinh(k[i] * depth)) * sin(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - Irregular::omega[i] * t + phi);
	}

	return wsum;

}

double Irregular::dp(double t, double xx, double yy, double zz) {

	double psum = 0.0;
	double phi;


	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * tofmax + Irregular::phi[i];
		psum += Ampspec[i] * D[i] * RHO * G * (cosh(k[i] * (zz + depth)) / cosh(k[i] * depth)) * cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - Irregular::omega[i] * t + phi);
	}

	return psum;

}

/* vertical velocity gradient at z=0 for velocity component U */
double Irregular::phi_dxdz(double t, double xx, double yy) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * tofmax + Irregular::phi[i];
		usum += cos(thetaA[i] + (mtheta * PI / 180.)) * Ampspec[i] * D[i] * Irregular::omega[i] * k[i] * cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - Irregular::omega[i] * t + phi);
	}
	return usum;

}



/* vertical velocity gradient at z=0 for velocity component V */
double Irregular::phi_dydz(double t, double xx, double yy) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * tofmax + Irregular::phi[i];
		vsum += sin(thetaA[i] + (mtheta * PI / 180.)) * Ampspec[i] * D[i] * Irregular::omega[i] * k[i] * cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - Irregular::omega[i] * t + phi);
	}
	return vsum;

}


/* vertical velocity gradient at z=0 for velocity component W */
double Irregular::phi_dzdz(double t, double xx, double yy) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {

		phi = Irregular::omega[i] * tofmax + Irregular::phi[i];
		wsum += D[i] * Ampspec[i] * Irregular::omega[i] * k[i] * (cosh(k[i] * depth) / sinh(k[i] * depth)) * sin(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - Irregular::omega[i] * t + phi);
	}
	return wsum;

}