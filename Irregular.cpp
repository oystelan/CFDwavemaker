#include "Irregular.h"
#include <math.h>
#include <algorithm>
#include <numeric>

#define PI 3.1415926535897
#define G 9.81
#define RHO 1025.0


double Irregular::eta(double t, double x, double y) {
	switch(pertmet){
		case 2:
			return (swl + eta1(t, x, y) + eta2(t, x, y));
	
		default:
			//std::cout << eta1(t, x, y) << std::endl;
			return (swl + eta1(t, x, y));
			
	}
};

double Irregular::u(double t, double x, double y, double z) {

	switch (pertmet) {
		case 2:
			{double z0 = std::min(0., z - swl);
			return (u1(t, x, y, z0) + u2(t, x, y, z0) + phi1_dxdz(t, x, y) * std::max(0., z - swl));}
		case 1:
			z = std::min(0., z - swl);
			return (u1(t, x, y, z) + u2(t, x, y, z));
			
		default:
			return u1(t, x, y, z - swl);
			
	}
};
double Irregular::v(double t, double x, double y, double z) {
	switch (pertmet) {
		case 2:
			{double z0 = std::min(0., z - swl);
			return (v1(t, x, y, z0) + v2(t, x, y, z0) + phi1_dydz(t, x, y) * std::max(0., z - swl)); }
			
		case 1:
			z = std::min(0., z - swl);
			return (v1(t, x, y, z) + v2(t, x, y, z));
			
		default:
			return v1(t, x, y, z - swl);
			
	}
};
double Irregular::w(double t, double x, double y, double z) {
	switch (pertmet) {
		case 2:
			{double z0 = std::min(0., z - swl);
			return (w1(t, x, y, z0) + w2(t, x, y, z0) + phi1_dzdz(t, x, y) * std::max(0., z - swl)); }
			
		case 1:
			z = std::min(0., z - swl);
			return (w1(t, x, y, z) + w2(t, x, y, z));
			
		default:
			return w1(t, x, y, z - swl);
			
	}
};

/* velocity profile for horizontal component. For implmentation of sloping bottom*/
double Irregular::profileX(int ind, double x, double y, double z) {
	switch (sloping_bottom) {
	case 0:
		return (cosh(k[ind] * (z + depth)) / cosh(k[ind] * depth));
	case 1:
		return 0.0;
	}
}

double Irregular::profileZ(int ind, double x, double y, double z) {
	switch (sloping_bottom) {
	case 0:
		return (sinh(k[ind] * (z + depth)) / cosh(k[ind] * depth));
	case 1:
		return 0.0;
	}
}

double Irregular::dp(double t, double x, double y, double z) {
	//todo: implement second order pressure component
	return dp1(t, x, y, z - swl);
};

/* First order wave elevation */
double Irregular::eta1(double t, double xx, double yy) {

	double welev = 0.0;
	double phi;
	//std::cout << "ndir" << ndir << std::endl;
	//std::cout << "nfreq" << nfreq << std::endl;
	//std::cout << "mtheta" << mtheta << std::endl;
	
	for (int i = 0; i < ndir * nfreq; i++) {

		phi = omega[i] * tofmax + phase[i];
		welev += A[i] * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) 
			+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
		//std::cout << i << ": A" << A[i] << std::endl;
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
			double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(theta[ci] - theta[cm]);
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

				double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + 
					(mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				eta2_t += ((A[ci] * A[cm] * omega[ci] * omega[cm]) / (2. * G)) * (alpha_nm_minus * cos(phi_i - phi_m) + alpha_nm_plus * cos(phi_i + phi_m));
			}
		}
	}

	return eta2_t;

}


double Irregular::u1(double t, double xx, double yy, double zz) {

	double usum = 0.0;
	double phi;
	//(cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))
	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		usum +=  A[i] * G * k[i] / omega[i] * profileX(i, xx, yy, zz)
			* cos(theta[i] + (mtheta * PI / 180.)) * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.))
				* (xx - fpoint[0]) + sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}

	return usum;

}


/* horizontal velocity U for a sinus Irregular::omegaave */
double Irregular::v1(double t, double xx, double yy, double zz) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		vsum += A[i] * G * k[i] / omega[i] * profileX(i, xx, yy, zz)
			* sin(theta[i] + (mtheta * PI / 180.)) * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0])
				+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double Irregular::w1(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		wsum += A[i] * G * k[i] / omega[i] * profileZ(i, xx, yy, zz)
			* sin(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] + (mtheta 
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
			double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] +
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
				double gamma_nm = cos(theta[ci] - theta[cm]);
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

				double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				usum2 += (A[ci] * A[cm] * omega[ci] * omega[cm]) *
					(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * cos(theta[ci] + (mtheta * PI / 180.)) - k[cm] * cos(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * cos(theta[ci] + (mtheta * PI / 180.)) + k[cm] * cos(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



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
			double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(theta[ci] - theta[cm]);
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
				double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				vsum2 += (A[ci] * A[cm] * omega[ci] * omega[cm]) *
					(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * sin(theta[ci] + (mtheta * PI / 180.)) - k[cm] * sin(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * sin(theta[ci] + (mtheta * PI / 180.)) + k[cm] * sin(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



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
			double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adjusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(theta[ci] - theta[cm]);
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

				double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				wsum2 += (A[ci] * A[cm] * omega[ci] * omega[cm]) *
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
			psum += A[i] * RHO * G * (cosh(k[i] * (zz + depth)) / cosh(k[i] * depth)) 
				* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
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
		usum += cos(theta[i] + (mtheta * PI / 180.)) * A[i] * omega[i] * k[i] 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
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
		vsum += sin(theta[i] + (mtheta * PI / 180.)) * A[i] * omega[i] * k[i] 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
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
		wsum += A[i] * omega[i] * k[i] * (cosh(k[i] * depth) / sinh(k[i] 
			* depth)) * sin(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) 
				+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi);
	}
	return wsum;

}



// First order velocity potential
double Irregular::phi1_pot(double t, double xx, double yy, double zz) {

	double phisum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		phisum += A[i] * (omega[i] / k[i]) * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*sin(k[i] * (cos(theta[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(theta[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - omega[i] * t + phi);
	}

	return phisum;

}


/* Second order component of velocity potential
Warning: Unsure if this is finished...should be checked thoroughly at some point*/
double Irregular::phi2_pot(double t, double xx, double yy, double zz) {
	//double eta1_t = 0;
	double phisum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i * ndir + j;
			double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < std::min(nfreq, i + bandwidth); m++) {
				int cm = m * ndir + j;
				double gamma_nm = cos(theta[ci] - theta[cm]);
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

				double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

				phisum2 += (A[ci] * A[cm] * omega[ci] * omega[cm]) *
					(beta_nm_minus * cos(phi_i - phi_m) * (k[ci] * cos(theta[ci] + (mtheta * PI / 180.)) - k[cm] * cos(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_minus * (zz + depth)) / cosh(k_nm_minus * depth)) +
						beta_nm_plus * cos(phi_i + phi_m) * (k[ci] * cos(theta[ci] + (mtheta * PI / 180.)) + k[cm] * cos(theta[cm] + (mtheta * PI / 180.))) * (cosh(k_nm_plus * (zz + depth)) / cosh(k_nm_plus * depth)));



			}
		}
	}
	return phisum2;
}

double Irregular::sum(double ll[], int nsum) {
	double ss = 0.0;
	for (int i = 0; i < nsum; i++) {
		ss += ll[i];
	}
	return ss;
}

void Irregular::normalize_data() {
	// Normalize and/or amplify the amplitude spectrum if Normalize is switched on
	double sumA = std::accumulate(A.begin(), A.end(), 0.);
	if (normalize) {
		for (int i = 0; i < nfreq*ndir; i++) {
			A[i] = ampl * A[i] / sumA;
		}
	}
	else {
		for (int i = 0; i < nfreq*ndir; i++) {
			A[i] = ampl * A[i];
		}
	}
};


