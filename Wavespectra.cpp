#include "Wavespectra.h"
#include <cmath>


#ifndef PI
#define PI 3.1415926535897
#endif

double Wavespectra::gamma(double f, double Gamma, bool gg) {
	// -------------------------------------------------------------------------
	// Returns a Gamma spectrum for given omega, Hs, Tp, and Gamma
	//
	// Warning: This is not the fully implemented gammaspectrumand is only
	// valied for Mand N = 4.
	// 
	// -------------------------------------------------------------------------


	double N = 4.0;
	double M = 4.0;
	double G_0 = 3.26; // for M = 4 and N = 4;

	double A_gamma = (1. + 1.1 * pow(log(Gamma), 1.19)) / Gamma;
	double sigma = 0.09;
	if (f < 1.) {
		double sigma = 0.07;
	}
	
	double gamma_F = 1;
	if (gg) {
		double gamma_F = pow(Gamma, exp(-(1.0 / (2 * pow(sigma, 2))) * (pow(f - 1, 2))));
	}
	return G_0 * A_gamma * pow(f, -N) * exp(-(N / M) * pow(f, -M)) * gamma_F;

}

void Wavespectra::torsethaugen2004(double* S, double* omega, int nfreq, double hs, double tp) {
	// -------------------------------------------------------------------------
	// Description : Returns the simplified double peaked spectrum.The Algorithm
	// is taken from : Paper NO.2004 - JSC - 193, Knut Torsethaugenand Sverre Haver,
	// "Simplified double peak spectral model for ocean waves" Published 2004
	//
	// S - spectral value S(omega)
	// omega(1xn) - frequency(rad)
	// Hs - Significant Waveheight
	// Tp - Peak period
	//
	//-------------------------------------------------------------------------

	// Empirical Parameters
	// ----------------------
	double af = 6.6;
	double ae = 2.0;
	double au = 25;
	double a10 = 0.7;
	double a1 = 0.5;
	double kg = 35.0;
	double b1 = 2.0;
	double a20 = 0.6;
	double a2 = 0.3;
	double a3 = 6;
	double Tpf = af * pow(hs, (1.0 / 3));
	double Tl = ae * pow(hs, (1.0 / 2));
	double Tu = au;
	double el = (Tpf - tp) / (Tpf - Tl);
	double eu = (tp - Tpf) / (Tu - Tpf);
	double g = 9.81;

	for (int i = 0; i < nfreq; i++) {
		double f = omega[i] / (2. * PI);
		if (tp < Tpf) {
			// WIND DRIVEN SEA
			// -----------------
			double Rw = (1 - a10) * exp(-pow((el / a1), 2)) + a10;
			// 1) Primary Peak
			double Hw1 = Rw * hs;
			double Tpw1 = tp;
			double gamma_w1 = kg * pow((2 * PI / g) * Hw1 / pow(Tpw1, 2), (6.0 / 7));

			double fn1 = f * Tpw1;

			// 2) Secondary peak
			double Hw2 = sqrt(1 - pow(Rw, 2)) * hs;
			double Tpw2 = Tpf + b1;
			double fn2 = f * Tpw2;
			double gamma_w2 = 1;

			S[i] = (1. / (2. * PI)) * (1.0 / 16) * pow(Hw1, 2) * Tpw1 * gamma(fn1, gamma_w1, true) + (1.0 / 16) * pow(Hw2, 2) * Tpw2 * gamma(fn2, gamma_w2, false);
		}
		else {
			// SWELL DRIVEN SEA
			//------------------
			double Rs = (1 - a20) * exp(-pow(eu / a2, 2)) + a20;
			// 1) Primary peak
			double Hs1 = Rs * hs;
			double Tps1 = tp;
			double gamma_s1 = (kg * pow((2 * PI / g) * hs / pow(Tpf, 2), (6.0 / 7))) * (1 + a3 * eu);
			double fn1 = f * Tps1;
			// 2) Secondary Peak
			double Hs2 = sqrt(1 - pow(Rs, 2)) * hs;
			double Tps2 = af * pow(Hs2, (1.0 / 3));
			double fn2 = f * Tps2;
			double gamma_s2 = 1;

			S[i] = (1. / (2. * PI)) * (1.0 / 16) * pow(Hs1, 2) * Tps1 * gamma(fn1, gamma_s1, true) + (1.0 / 16) * pow(Hs2, 2) * Tps2 * gamma(fn2, gamma_s2, false);

		}
	}
}
void Wavespectra::torsethaugen1996(double* S, double* omega, double hs, double tp) {

}

void Wavespectra::PM(double* S, double* omega, int nfreq, double hs, double tp) {
	// The Pierson-Moskowitz spectrum (Ref. DNVGL-RP-C205, Oct 2019)
	double omega_p = (2. * PI) / tp;
	for (int i = 0; i < nfreq; i++) {
		S[i] = (5. / 16.) * pow(hs, 2.) * pow(omega_p, 4.) * pow(omega[i], -5.) * exp((-5. / 4.) * pow(omega[i] / omega_p, -4.));
	}
}
void Wavespectra::jonswap3(double* S, double* omega, int nfreq, double hs, double tp, double gam) {
	// JONSWAP wave spectrum (Ref. DNVGL-RP-C205, oct 2019)
	double omega_p = (2. * PI) / tp;
	for (int i = 0; i < nfreq; i++) {
		double S_pm = (5. / 16.) * pow(hs, 2.) * pow(omega_p, 4.) * pow(omega[i], -5.) * exp((-5. / 4.) * pow(omega[i] / omega_p, -4.));
		
		double A_gamma = 0.2 / (0.065 * pow(gam, 0.803) + 0.135);
		double sigma = 0.09;
		if (omega[i] <= omega_p) {
			sigma = 0.07;
		}
		S[i] = A_gamma * S_pm * pow(gam, exp(-0.5 * pow((omega[i] - omega_p) / (sigma * omega_p), 2)));
	}

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

void Wavespectra::spreading_ewans(double tp, double mtheta) {
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

