#include "Wavespectra.h"
#include <cmath>


#ifndef PI
#define PI 3.1415926535897
#endif

void Wavespectra::jonswap3(double* ampl, double* omega, double hs, double tp, double gamma) {

}
void Wavespectra::torsethaugen2004(double* ampl, double* omega, double hs, double tp) {

}
void Wavespectra::torsethaugen1996(double* ampl, double* omega, double hs, double tp) {

}
void Wavespectra::PM(double* ampl, double* omega, double hs, double tp) {

}

void Wavespectra::spreading_uniform(double * D, int nfreq, int ndir){
	// Assign wave spreading based on specified spreading function parameters
	D = new double[nfreq * ndir];
	double dsum = 0.;
	for (int i = 0; i < ndir; i++) {
		D[i] = 1.0;
		dsum += D[i];
	}
	for (int i = 0; i < ndir; i++) {
		D[i] = D[i] / dsum;
	}
	// Repeat Directional distribution nfreq times
	for (int i = 0; i < nfreq; i++) {
		for (int j = 0; j < ndir; j++) {
			D[i * ndir + j] = D[j];
		}
	}
}
void Wavespectra::spreading_cos_theta_n(double* D, int nfreq, int ndir, double s, double mtheta, double* theta) {
	// cos(theta)^n
	D = new double[nfreq * ndir];
	double dsum = 0.;

	for (int i = 0; i < ndir; i++) {
		D[i] = pow(cos((theta[i] - (mtheta * PI / 180.))), s);
		dsum += D[i];
	}
	for (int i = 0; i < ndir; i++) {
		D[i] = D[i] / dsum;
	}
	// Repeat Directional distribution nfreq times
	for (int i = 0; i < nfreq; i++) {
		for (int j = 0; j < ndir; j++) {
			D[i * ndir + j] = D[j];
		}
	}
}

void Wavespectra::spreading_cos_theta05_2s(double* D, int nfreq, int ndir, double s, double mtheta, double* theta) {
	// cos(theta/2)^2s
	D = new double[nfreq * ndir];
	double dsum = 0.;
	
	for (int i = 0; i < ndir; i++) {
		D[i] = pow(cos((theta[i] - (mtheta * PI / 180.)) / 2.), 2.0 * s);
		dsum += D[i];
	}
	for (int i = 0; i < ndir; i++) {
		D[i] = D[i] / dsum;
	}
	// Repeat Directional distribution nfreq times
	for (int i = 0; i < nfreq; i++) {
		for (int j = 0; j < ndir; j++) {
			D[i * ndir + j] = D[j];
		}
	}
}

void Wavespectra::spreading_ewans(double* D, int nfreq, int ndir, double mtheta, double* theta, double* omega, double tp) {
	// Ewans simplified spreading function (frequency dependent spreading)
	// EWANS, Kevin C.Observations of the directional spectrum of fetch - limited waves.Journal of Physical Oceanography, 1998, 28.3: 495 - 512.
	double* dsum2 = new double[nfreq];
	int s;
	double fp = 1. / tp;

	for (int i = 0; i < nfreq; i++) {
		if (((omega[i] / (2. * PI)) / fp) < 1.) {
			s = 15.5 * pow((omega[i] / (2. * PI)) / fp, 9.47);
		}
		else {
			s = 13.1 * pow((omega[i] / (2. * PI)) / fp, -1.94);
		}
		dsum2[i] = 0.;
		for (int j = 0; j < ndir; j++) {
			D[i * ndir + j] = pow(cos((theta[j] - (mtheta * PI / 180.)) / 2.), 2. * s);
			dsum2[i] += D[i * ndir + j];
		}
	}
	for (int i = 0; i < nfreq; i++) {
		for (int j = 0; j < ndir; j++) {
			D[i * ndir + j] = D[i * ndir + j] / dsum2[i];
		}
	}
}

