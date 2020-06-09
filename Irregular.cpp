#include "Irregular.h"
#include <math.h>
#include <algorithm>
#include <numeric>




double Irregular::eta(double t, double x, double y) {
	switch(order){
		case 2:
			return (swl + eta1(t, x, y) + eta2(t, x, y));
	
		default:
			//std::cout << eta1(t, x, y) << std::endl;
			return (swl + eta1(t, x, y));
			
	}
};

double Irregular::u(double t, double x, double y, double z) {

	switch (order) {
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
	switch (order) {
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
	switch (order) {
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
	return 0.0;
}

double Irregular::profileZ(int ind, double x, double y, double z) {
	switch (sloping_bottom) {
	case 0:
		return (sinh(k[ind] * (z + depth)) / cosh(k[ind] * depth));
	case 1:
		return 0.0;
	}
	return 0.0;
}

double Irregular::dp(double t, double x, double y, double z) {
	//todo: implement second order pressure component
	return dp1(t, x, y, z - swl);
};

/* Function for printing wave components (sometimes useful for QA)*/
void Irregular::print()
{
	std::cout << "Irregular Waves" << std::endl;
	std::cout << "===============" << std::endl;
	std::cout << "\n";
	std::cout << "Number of frequencies: " << nfreq << std::endl;
	std::cout << "Number of wave directions " << ndir << std::endl;
	std::cout << "\n";
	std::cout << "Wave components:\n";
	std::cout << "----------------\n";
	std::cout << "|  Omega (rad/s)   |  Amplitude (m)  |  Wave number k  |   phase (rad)   | \n";
	for (int i = 0; i < (nfreq * ndir); i++) {
		std::cout << omega[i] << " " << A[i] << " " << k[i] << " " << phase[i] << std::endl;
	}
	std::cout << std::endl;
}

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

int Irregular::bw_limiter(int i1) {
	int cc = i1;
	while ((omega[cc]-omega[i1] < dw_bandwidth) && (cc < nfreq * ndir)) {	
		cc++;
	}
	return cc;
}

void Irregular::calculate_bwindices() {
	bwlim.clear();
	for (int i = 0; i < nfreq * ndir -1 ; i++) {		
		bwlim.push_back(bw_limiter(i));
		std::cout << "i: " << i << ", bwlim: " << bw_limiter(i) << ", omega: " << omega[i] << ", omegalim: " << omega[bw_limiter(i)-1] << std::endl;
	}
}


/* Second order wave elevation */
double Irregular::eta2(double t, double xx, double yy) {

	//double eta1_t = 0;
	double eta2_t = 0;
	for (int i = 0; i < nfreq*ndir - 1; i++) {
		// Second order
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] 
			+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
			
		eta2_t += 0.5 * A[ci] * A[ci] * k[ci] * cos(phi_i);

		for (int m = i + 1; m < bwlim[i]; m++) {
			int cm = m;	
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

	for (int i = 0; i < nfreq * ndir - 1; i++) {
	
		// Second order
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] +
			(mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
		//cout << "Rn: " << Rn << endl;

		//// Adiusting Bandwidth for 2 order cut-off
		//if (i + 1 + f_bw<nfreqs) {
		//	p = i + 1 + f_bw;
		//}
		//else { p = nfreqs; }
		for (int m = i + 1; m < bwlim[i]; m++) {
			int cm = m;
			
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
	return usum2;
}


/* Second order horizontal velocity V for a sinus wave */
double Irregular::v2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double vsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		
		// Second order
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
		//cout << "Rn: " << Rn << endl;

		//// Adiusting Bandwidth for 2 order cut-off
		//if (i + 1 + f_bw<nfreqs) {
		//	p = i + 1 + f_bw;
		//}
		//else { p = nfreqs; }
		for (int m = i + 1; m < bwlim[i]; m++) {
			int cm = m;
				
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
	return vsum2;
}

/* Second order vertical velocity component W for a sinus wave */
double Irregular::w2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		
		// Second order
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
		
		for (int m = i + 1; m < bwlim[i]; m++) {
			int cm = m;
			
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
			for (int m = i + 1; m < nfreq; m++) {
				int cm = m * ndir + j;
				if (abs(omega[cm] - omega[ci]) <= dw_bandwidth) {
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
}


// trapezoidal integration function
double Irregular::trapz(double x[], double y[], int n) {
	
	double ss = 0.;

	for (int i = 0; i < (n-1); i++) {
		ss += (x[i+1]-x[i])*((y[i+1]+y[i])*0.5);
	}

	return ss;
}

double Irregular::interpolate(std::vector<double>& xData, std::vector<double>& yData, double x, bool extrapolate)
{
	int size = xData.size();

	int i = 0;                                                                  // find left end of interval for interpolation
	if (x >= xData[size - 2])                                                 // special case: beyond right end
	{
		i = size - 2;
	}
	else
	{
		while (x > xData[i + 1]) i++;
	}
	double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1];      // points on either side (unless beyond ends)
	if (!extrapolate)                                                         // if beyond ends of array and not extrapolating
	{
		if (x < xL) yR = yL;
		if (x > xR) yL = yR;
	}

	double dydx = (yR - yL) / (xR - xL);                                    // gradient

	return yL + dydx * (x - xL);                                              // linear interpolation
}

/* Function for calculating phase velocity based integration of the wave spectrum*/
double Irregular::phase_velocity(int opt)
{
	// check 
	if (ndir > 1) {
		std::cerr << "Phase velocity only supported for 2D wave spectra at the moment. Sorry for the inconvenience" << std::endl;
		exit(-1);
	}
	
	double w_t, dw, S;
	double m0 = 0.;
	double m1 = 0.;
	double m2 = 0.;
	
	for (int i = 0; i < (nfreq-1); i++) {
		w_t = (omega[i + 1] + omega[i]) / 2.;
		dw = omega[i + 1] - omega[i];
		S = pow((A[i+1]+A[i])*0.5,2.) / (2. * dw);
		m0 += w_t * S;
		m1 += w_t * S * w_t;
		m2 += w_t * S * w_t * w_t;
	}
	
	double t1 = 2. * PI * (m0 / m1);
	double t2 = 2. * PI * sqrt(m0 / m2);

	double w1 = 2 * PI / t1;
	double w2 = 2 * PI / t2;

	// interpolate wave number
	double k1 = interpolate(omega, k, w1, true);
	double k2 = interpolate(omega, k, w2, true);

	if (opt == 1) {
		return w1 / k1;
	}
	else if(opt == 2) {
		return w2 / k2;
	}
	return 0.;
}

/* Function for calculating the mean wave length based integration of the wave spectrum*/
double Irregular::mean_wave_length(int opt)
{
	// check 
	if (ndir > 1) {
		std::cerr << "mean wave length only supported for 2D wave spectra at the moment. Sorry for the inconvenience" << std::endl;
		exit(-1);
	}

	double w_t, dw, S;
	double m0 = 0.;
	double m1 = 0.;
	double m2 = 0.;

	for (int i = 0; i < (nfreq - 1); i++) {
		w_t = (omega[i + 1] + omega[i]) / 2.;
		dw = omega[i + 1] - omega[i];
		S = pow((A[i + 1] + A[i]) * 0.5, 2.) / (2. * dw);
		m0 += w_t * S;
		m1 += w_t * S * w_t;
		m2 += w_t * S * w_t * w_t;
	}

	double t1 = 2. * PI * (m0 / m1);
	double t2 = 2. * PI * sqrt(m0 / m2);

	double w1 = 2 * PI / t1;
	double w2 = 2 * PI / t2;

	// interpolate wave number
	double l1 = (2. * PI) / interpolate(omega, k, w1, true);
	double l2 = (2. * PI) / interpolate(omega, k, w2, true);

	if (opt == 1) {
		return l1;
	}
	else if (opt == 2) {
		return l2;
	}
	return 0.;
}

/* Function for calculating the mean wave length based integration of the wave spectrum*/
double Irregular::mean_wave_period(int opt)
{
	// check 
	if (ndir > 1) {
		std::cerr << "mean wave period only supported for 2D wave spectra at the moment. Sorry for the inconvenience" << std::endl;
		exit(-1);
	}

	double w_t, dw, S;
	double m0 = 0.;
	double m1 = 0.;
	double m2 = 0.;

	for (int i = 0; i < (nfreq - 1); i++) {
		w_t = (omega[i + 1] + omega[i]) / 2.;
		dw = omega[i + 1] - omega[i];
		S = pow((A[i + 1] + A[i]) * 0.5, 2.) / (2. * dw);
		m0 += w_t * S;
		m1 += w_t * S * w_t;
		m2 += w_t * S * w_t * w_t;
	}

	double t1 = 2. * PI * (m0 / m1);
	double t2 = 2. * PI * sqrt(m0 / m2);

	if (opt == 1) {
		return t1;
	}
	else if (opt == 2) {
		return t2;
	}
	return 0.;
}



