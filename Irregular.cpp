#include "Irregular.h"
#include <math.h>
#include <algorithm>

#define PI 3.1415926535897
#define G 9.81
#define RHO 1025.0

/* First order wave elevation */
double Irregular::eta(double t, double xx, double yy) {

	double welev = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {

		phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
		welev += Irregular::Ampspec[i] * Irregular::D[i] * cos(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) 
			+ sin(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
	}


	return welev;

}

/* Second order wave elevation */
double Irregular::eta2(double t, double xx, double yy) {

	//double eta1_t = 0;
	double eta2_t = 0;

	for (int i = 0; i < Irregular::nfreq - 1; i++) {
		for (int j = 0; j < Irregular::ndir; j++) {
			// Second order
			int ci = i * Irregular::ndir + j;
			double phi_i = Irregular::k[ci] * (cos(Irregular::thetaA[ci] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) + sin(Irregular::thetaA[ci] 
				+ (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[ci] * t + Irregular::omega[ci] * Irregular::tofmax + Irregular::phi[ci];
			double Rn = Irregular::k[ci] * tanh(Irregular::k[ci] * Irregular::depth);
			
			for (int m = i + 1; m < std::min(Irregular::nfreq, i + Irregular::bandwidth); m++) {
				int cm = m * Irregular::ndir + j;
				double gamma_nm = cos(Irregular::thetaA[ci] - Irregular::thetaA[cm]);
				double k_nm_plus = sqrt(Irregular::k[ci] * Irregular::k[ci] + Irregular::k[cm] * Irregular::k[cm] + (2. * Irregular::k[ci] * Irregular::k[cm] * gamma_nm));
				double k_nm_minus = sqrt(Irregular::k[ci] * Irregular::k[ci] + Irregular::k[cm] * Irregular::k[cm] - (2. * Irregular::k[ci] * Irregular::k[cm] * gamma_nm));

				double Rm = Irregular::k[cm] * tanh(Irregular::k[cm] * Irregular::depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (Irregular::k[ci] * Irregular::k[ci] - Rn * Rn) + sqrt(Rn) * (Irregular::k[cm] * Irregular::k[cm] - Rm * Rm)) +
					2. * pow(sqrt(Rn) + sqrt(Rm), 2.) * (Irregular::k[ci] * Irregular::k[cm] * gamma_nm - Rn * Rm) / (pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus * tanh(k_nm_plus * Irregular::depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (Irregular::k[ci] * Irregular::k[ci] - Rn * Rn) - sqrt(Rn) * (Irregular::k[cm] * Irregular::k[cm] - Rm * Rm)) + 
					2. * pow(sqrt(Rn) - sqrt(Rm), 2.) * (Irregular::k[ci] * Irregular::k[cm] * gamma_nm + Rn * Rm) / (pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus * tanh(k_nm_minus * Irregular::depth));
				// Catch NaN when two equal frequency components interact
				if (Irregular::omega[ci] == Irregular::omega[cm]) {
					D_nm_minus = 0.;
				}
				double alpha_nm_minus = (((Irregular::omega[ci] / Irregular::omega[cm]) + (Irregular::omega[cm] / Irregular::omega[ci])) + (G * G / (Irregular::omega[ci] * Irregular::omega[cm])) *
					((D_nm_minus - Irregular::k[ci] * Irregular::k[cm] * (gamma_nm + tanh(Irregular::k[ci] * Irregular::depth) * tanh(Irregular::k[cm] * depth))) / (Irregular::omega[ci] * Irregular::omega[cm])));
				double alpha_nm_plus = (((Irregular::omega[ci] / Irregular::omega[cm]) + (Irregular::omega[cm] / Irregular::omega[ci])) + (G * G / (Irregular::omega[ci] * Irregular::omega[cm])) *
					((D_nm_plus - Irregular::k[ci] * Irregular::k[cm] * (gamma_nm - tanh(Irregular::k[ci] * Irregular::depth) * tanh(Irregular::k[cm] * Irregular::depth))) / (Irregular::omega[ci] * Irregular::omega[cm])));

				double phi_m = Irregular::k[cm] * (cos(Irregular::thetaA[cm] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) + sin(Irregular::thetaA[cm] + 
					(Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[cm] * t + Irregular::omega[cm] * tofmax + Irregular::phi[cm];

				eta2_t += ((Irregular::Ampspec[ci] * Irregular::D[ci] * Irregular::Ampspec[cm] * Irregular::D[cm] * Irregular::omega[ci] * Irregular::omega[cm]) / (2. * G)) * (alpha_nm_minus * cos(phi_i - phi_m) + alpha_nm_plus * cos(phi_i + phi_m));



			}
		}
	}

	return eta2_t;

}


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
		phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
		vsum += sin(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * Irregular::Ampspec[i] * Irregular::D[i] * Irregular::omega[i] * (cosh(Irregular::k[i] 
			* (zz + depth)) / sinh(k[i] * Irregular::depth)) * cos(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) 
				+ sin(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double Irregular::w(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
		wsum += Irregular::D[i] * Irregular::Ampspec[i] * Irregular::omega[i] * (sinh(Irregular::k[i] * (zz + Irregular::depth)) / sinh(Irregular::k[i] * Irregular::depth)) 
			* sin(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) + sin(Irregular::thetaA[i] + (Irregular::mtheta 
				* PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
	}

	return wsum;

}

/* Second order horizontal velocity U for a sinus wave */
double uu_2order(double t, double xx, double yy, double zz) {



	//double eta1_t = 0;
	double usum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) + sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) + sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm - Rn * Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus * tanh(k_nm_plus * depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) - sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) - sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm + Rn * Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus * tanh(k_nm_minus * depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				usum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * cos(thetaA[ci] + (mtheta * PI / 180.)) - k[cm] * cos(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * cos(thetaA[ci] + (mtheta * PI / 180.)) + k[cm] * cos(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



			}
		}
	}
	return usum2;
}


/* Second order horizontal velocity V for a sinus wave */
double vv_2order(double t, double xx, double yy, double zz) {



	//double eta1_t = 0;
	double vsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) + sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) + sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm - Rn * Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus * tanh(k_nm_plus * depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) - sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) - sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm + Rn * Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus * tanh(k_nm_minus * depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}
				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				vsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * sin(thetaA[ci] + (mtheta * PI / 180.)) - k[cm] * sin(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * sin(thetaA[ci] + (mtheta * PI / 180.)) + k[cm] * sin(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



			}
		}
	}
	return vsum2;
}

/* Second order vertical velocity component W for a sinus wave */
double ww_2order(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adjusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) + sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) + sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm - Rn * Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus * tanh(k_nm_plus * depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) - sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 2. * pow(sqrt(Rn) - sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm + Rn * Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus * tanh(k_nm_minus * depth));

				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				wsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus * k_nm_minus * sin(phi_i - phi_m) * (sinh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * k_nm_plus * sin(phi_i + phi_m) * (sinh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



			}
		}
	}
	return wsum2;
}

	double Irregular::dp(double t, double xx, double yy, double zz) {

		double psum = 0.0;
		double phi;

		for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
			phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
			psum += Irregular::Ampspec[i] * Irregular::D[i] * RHO * G * (cosh(Irregular::k[i] * (zz + Irregular::depth)) / cosh(Irregular::k[i] * Irregular::depth)) 
				* cos(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) + sin(Irregular::thetaA[i] 
					+ (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
		}

		return psum;

	}

/* vertical velocity gradient at z=0 for velocity component U */
double Irregular::phi_dxdz(double t, double xx, double yy) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
		usum += cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * Irregular::Ampspec[i] * Irregular::D[i] * Irregular::omega[i] * Irregular::k[i] 
			* cos(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) + sin(Irregular::thetaA[i] 
				+ (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
	}
	return usum;

}


/* vertical velocity gradient at z=0 for velocity component V */
double Irregular::phi_dydz(double t, double xx, double yy) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
		vsum += sin(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * Irregular::Ampspec[i] * Irregular::D[i] * Irregular::omega[i] * Irregular::k[i] 
			* cos(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) + sin(Irregular::thetaA[i] 
				+ (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
	}
	return vsum;

}


/* vertical velocity gradient at z=0 for velocity component W */
double Irregular::phi_dzdz(double t, double xx, double yy) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < Irregular::ndir * Irregular::nfreq; i++) {
		phi = Irregular::omega[i] * Irregular::tofmax + Irregular::phi[i];
		wsum += Irregular::D[i] * Irregular::Ampspec[i] * Irregular::omega[i] * Irregular::k[i] * (cosh(Irregular::k[i] * Irregular::depth) / sinh(Irregular::k[i] 
			* Irregular::depth)) * sin(Irregular::k[i] * (cos(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (xx - Irregular::fpoint[0]) 
				+ sin(Irregular::thetaA[i] + (Irregular::mtheta * PI / 180.)) * (yy - Irregular::fpoint[1])) - Irregular::omega[i] * t + phi);
	}
	return wsum;

}