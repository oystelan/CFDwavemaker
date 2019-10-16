#include "Irregular.h"
#include <math.h>
#include <algorithm>

#define PI 3.1415926535897
#define G 9.81
#define RHO 1025.0

/* First order wave elevation */
double Irregular::eta1(double t, double xx, double yy) {

	double welev = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {

		phi = omega[i] * tofmax + phase[i];
		welev += Ampspec[i] * D[i] * cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) 
			+ sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}


	return welev;

}

/* Second order wave elevation */
double Irregular::eta2(double t, double xx, double yy) {

	//double eta1_t = 0;
	double eta2_t = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[ci] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) + sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) +
					2. * pow(sqrt(Rn) + sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm - Rn * Rm) / (pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus * tanh(k_nm_plus * depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) - sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) + 
					2. * pow(sqrt(Rn) - sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm + Rn * Rm) / (pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus * tanh(k_nm_minus * depth));
				// Catch NaN when two equal frequency components interact
				if (omega[ci] == omega[cm]) {
					D_nm_minus = 0.;
				}
				double alpha_nm_minus = (((omega[ci] / omega[cm]) + (omega[cm] / omega[ci])) + (G * G / (omega[ci] * omega[cm])) *
					((D_nm_minus - k[ci] * k[cm] * (gamma_nm + tanh(k[ci] * depth) * tanh(k[cm] * depth))) / (omega[ci] * omega[cm])));
				double alpha_nm_plus = (((omega[ci] / omega[cm]) + (omega[cm] / omega[ci])) + (G * G / (omega[ci] * omega[cm])) *
					((D_nm_plus - k[ci] * k[cm] * (gamma_nm - tanh(k[ci] * depth) * tanh(k[cm] * depth))) / (omega[ci] * omega[cm])));

				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[cm] + 
					(mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				eta2_t += ((Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * omega[ci] * omega[cm]) / (2. * G)) * (alpha_nm_minus * cos(phi_i - phi_m) + alpha_nm_plus * cos(phi_i + phi_m));
			}
		}
	}

	return eta2_t;

}


double Irregular::u1(double t, double xx, double yy, double zz) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		usum += cos(thetaA[i] + (mtheta * PI / 180.)) * Ampspec[i] * D[i] * omega[i] * (cosh(k[i]
			* (zz + depth)) / sinh(k[i] * depth)) * cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) 
				* (xx - fpoint[0]) + sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}

	return usum;

}


/* horizontal velocity U for a sinus Irregular::omegaave */
double Irregular::v1(double t, double xx, double yy, double zz) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		vsum += sin(thetaA[i] + (mtheta * PI / 180.)) * Ampspec[i] * D[i] * omega[i] * (cosh(k[i] 
			* (zz + depth)) / sinh(k[i] * depth)) * cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) 
				+ sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double Irregular::w1(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		wsum += D[i] * Ampspec[i] * omega[i] * (sinh(k[i] * (zz + depth)) / sinh(k[i] * depth)) 
			* sin(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] + (mtheta 
				* PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}

	return wsum;

}

/* Second order horizontal velocity U for a sinus wave */
double Irregular::u2(double t, double xx, double yy, double zz) {



	//double eta1_t = 0;
	double usum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[ci] +
				(mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
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


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) + sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) 
					+ 2. * pow(sqrt(Rn) + sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm - Rn * Rm) / (pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus * tanh(k_nm_plus * depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm)) * (sqrt(Rm) * (k[ci] * k[ci] - Rn * Rn) - sqrt(Rn) * (k[cm] * k[cm] - Rm * Rm)) 
					+ 2. * pow(sqrt(Rn) - sqrt(Rm), 2.) * (k[ci] * k[cm] * gamma_nm + Rn * Rm) / (pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus * tanh(k_nm_minus * depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (omega[ci] - omega[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (omega[ci] + omega[cm]));

				// Catch NaN when two equal frequency components interact
				if (omega[ci] == omega[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				usum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * omega[ci] * omega[cm]) *
					(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * cos(thetaA[ci] + (mtheta * PI / 180.)) - k[cm] * cos(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * cos(thetaA[ci] + (mtheta * PI / 180.)) + k[cm] * cos(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



			}
		}
	}
	return usum2;
}


/* Second order horizontal velocity V for a sinus wave */
double Irregular::v2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double vsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
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


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (omega[ci] - omega[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (omega[ci] + omega[cm]));

				// Catch NaN when two equal frequency components interact
				if (omega[ci] == omega[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}
				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				vsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * omega[ci] * omega[cm]) *
					(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * sin(thetaA[ci] + (mtheta * PI / 180.)) - k[cm] * sin(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * sin(thetaA[ci] + (mtheta * PI / 180.)) + k[cm] * sin(thetaA[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



			}
		}
	}
	return vsum2;
}

/* Second order vertical velocity component W for a sinus wave */
double Irregular::w2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
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

				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (omega[ci] - omega[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (omega[ci] + omega[cm]));

				// Catch NaN when two equal frequency components interact
				if (omega[ci] == omega[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				wsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * omega[ci] * omega[cm]) *
					(beta_nm_minus * k_nm_minus * sin(phi_i - phi_m) * (sinh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * k_nm_plus * sin(phi_i + phi_m) * (sinh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



			}
		}
	}
	return wsum2;
}

	double Irregular::dp1(double t, double xx, double yy, double zz) {

		double psum = 0.0;
		double phi;

		for (int i = 0; i < ndir * nfreq; i++) {
			phi = omega[i] * tofmax + phase[i];
			psum += Ampspec[i] * D[i] * RHO * G * (cosh(k[i] * (zz + depth)) / cosh(k[i] * depth)) 
				* cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] 
					+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
		}

		return psum;

	}

/* vertical velocity gradient at z=0 for velocity component U */
double Irregular::phi1_dxdz(double t, double xx, double yy) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		usum += cos(thetaA[i] + (mtheta * PI / 180.)) * Ampspec[i] * D[i] * omega[i] * k[i] 
			* cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}
	return usum;

}


/* vertical velocity gradient at z=0 for velocity component V */
double Irregular::phi1_dydz(double t, double xx, double yy) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		vsum += sin(thetaA[i] + (mtheta * PI / 180.)) * Ampspec[i] * D[i] * omega[i] * k[i] 
			* cos(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(thetaA[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}
	return vsum;

}


/* vertical velocity gradient at z=0 for velocity component W */
double Irregular::phi1_dzdz(double t, double xx, double yy) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		wsum += D[i] * Ampspec[i] * omega[i] * k[i] * (cosh(k[i] * depth) / sinh(k[i] 
			* depth)) * sin(k[i] * (cos(thetaA[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) 
				+ sin(thetaA[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}
	return wsum;

}

double Irregular::eta(double t, double x, double y) {
	if (pertmet == 2) {
		return (eta1(t, x, y) + eta2(t, x, y));
	}
	else {
		return eta1(t, x, y);
	}
};

double Irregular::u(double t, double x, double y, double z) {
	
	if (pertmet == 2) {
		int z0 = std::min(0., z);
		return (u1(t, x, y, z) + u2(t, x, y, z) + phi1_dxdz(t, x, y) * std::max(0., z));
	}
	else if (pertmet == 1) {
		z = std::min(0., z);
		return (u1(t, x, y, z) + u2(t, x, y, z));
	}
	else {
		return u1(t, x, y, z);
	}
};
double Irregular::v(double t, double x, double y, double z) {
	if (pertmet == 2) {
		int z0 = std::min(0., z);
		return (v1(t, x, y, z) + v2(t, x, y, z) + phi1_dydz(t, x, y) * std::max(0., z));
	}
	else if (pertmet == 1) {
		z = std::min(0., z);
		return (v1(t, x, y, z) + v2(t, x, y, z));
	}
	else {
		return v1(t, x, y, z);
	}
};
double Irregular::w(double t, double x, double y, double z) {
	if (pertmet == 2) {
		int z0 = std::min(0., z);
		return (w1(t, x, y, z) + w2(t, x, y, z) + phi1_dzdz(t, x, y) * std::max(0., z));
	}
	else if (pertmet == 1) {
		z = std::min(0., z);
		return (w1(t, x, y, z) + w2(t, x, y, z));
	}
	else {
		return w1(t, x, y, z);
	}
};