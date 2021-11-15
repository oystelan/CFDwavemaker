#include "Irregular.h"
#include <math.h>
#include <algorithm>
#include <numeric>




double Irregular::eta(double t, double x, double y) {
	switch(order){
		case 2:
			return eta2(t,x,y);
			//return (swl + eta1(t, x, y) + eta2(t, x, y));
	
		default:
			//std::cout << eta1(t, x, y) << std::endl;
			return (swl + eta1(t, x, y));
			
	}
};

double Irregular::u(double t, double x, double y, double z) {

	switch (order) {
		case 2:
			{double z0 = std::min(0., z - swl);
			//std::cout << "ape" << std::endl;
			//return (u1(t, x, y, z0) + u2(t, x, y, z0) + phi1_dxdz(t, x, y) * std::max(0., z - swl));}
			return u2(t, x, y, z0);} 
		case 1:
			{z = std::min(0., z - swl);
			//std::cout << "ape2" << std::endl;
			return (u1(t, x, y, z) + u2(t, x, y, z));} 
			
		default:
			{return u1(t, x, y, z - swl);} 
			
	}
};
double Irregular::v(double t, double x, double y, double z) {
	switch (order) {
		case 2:
			{double z0 = std::min(0., z - swl);
			//return (v1(t, x, y, z0) + v2(t, x, y, z0) + phi1_dydz(t, x, y) * std::max(0., z - swl)); }
			return v2(t, x, y, z0);} 
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
			//return (w1(t, x, y, z0) + w2(t, x, y, z0) + phi1_dzdz(t, x, y) * std::max(0., z - swl)); }
			return w2(t, x, y, z0);} 
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
	for (int i = 0; i < nfreq * ndir ; i++) {		
		bwlim.push_back(bw_limiter(i));
		//std::cout << "i: " << i << ", bwlim: " << bw_limiter(i) << ", omega: " << omega[i] << ", omegalim: " << omega[bw_limiter(i)-1] << std::endl;
	}
}


/* Second order wave elevation */
double Irregular::eta2(double t, double xx, double yy) {

	//double eta1_t = 0;
	double eta2_t = 0;
	for (int i = 0; i < nfreq*ndir; i++) {
		// Second order,
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] 
			+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		
		double Rn = k[ci] * tanh(k[ci] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal
		double gamma_nm = 1;
		double k_nm_plus = sqrt(k[ci] * k[ci] + k[ci] * k[ci] + (2. * k[ci] * k[ci] * gamma_nm));
		double D_nm_plus = ( 8.*Rsn*Rsn * (k[ci] * k[ci] - Rn * Rn)) / (4.* Rsn * Rsn - k_nm_plus * tanh(k_nm_plus * depth)) +
		 (2.* Rsn * (Rsn*(k[ci]*k[ci] - Rn * Rn) + Rsn*(k[ci]*k[ci] - Rn*Rn)))/(4.* Rsn* Rsn - k_nm_plus * tanh(k_nm_plus * depth));
		double alpha_nm_plus = (D_nm_plus-(k[ci]*k[ci]-Rn*Rn))/Rn + 2*Rn;

		eta2_t += 0.25 * A[ci] * A[ci] * alpha_nm_plus * cos(2*phi_i);

		// Off-diagonal (sum over half only due to symmetry)
		for (int m = i + 1; m < nfreq*ndir; m++) {
			int cm = m;

			double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] +
				(mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

			double gamma_nm = cos(theta[ci] - theta[cm]);
			double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
			double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

			double Rm = k[cm] * tanh(k[cm] * depth);
			double Rsm = sqrt(Rm);
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[ci] * k[cm] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci] * k[ci] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[ci] * k[cm] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci]*k[ci] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double alpha_nm_plus = (D_nm_plus - (k[ci]*k[cm]* gamma_nm - Rn*Rm)) / sqrt(Rn*Rm) + (Rn + Rm);

			double alpha_nm_minus = (D_nm_minus - (k[ci]*k[cm]* gamma_nm + Rn*Rm)) / sqrt(Rn*Rm) + (Rn + Rm);

			eta2_t += 0.5* A[ci] * A[cm]  * (alpha_nm_plus * cos(phi_i + phi_m)); // sum term

			eta2_t += 0.5* A[ci] * A[cm]  * (alpha_nm_minus * cos(phi_i - phi_m)); // difference term
            
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
	
	/*double kmaxx = 0.;
	double kmaxy = 0.;
	for (int i = 0; i < nfreq * ndir; i++) {
		if (k[i]*cos(theta[i])  > kmaxx){
			kmaxx = k[i]*cos(theta[i]); 
		}
		if (abs(k[i]*sin(theta[i])) > kmaxy){
			kmaxy = k[i]*sin(theta[i]); 
		}   
	}*/	
	//std::cout << "kmax: " << kmax << std::endl;
	for (int i = 0; i < nfreq * ndir; i++) {
	
		// Second order
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] +
			(mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
		double Rsn = sqrt(Rn);
		
		// Diagonal terms
		double k_nm_plus = 2*k[ci];
		double D_nm_plus = ( 8.*Rsn*Rsn * (k[ci] * k[ci] - Rn * Rn)) / (4.* Rsn * Rsn - k_nm_plus * tanh(k_nm_plus * depth)) +
		 (2.* Rsn * (Rsn*(k[ci]*k[ci] - Rn * Rn) + Rsn*(k[ci]*k[ci] - Rn*Rn)))/(4.* Rsn* Rsn - k_nm_plus * tanh(k_nm_plus * depth));
		
		// Sum term
		//if ((((k[ci]*cos(theta[ci])) + (k[ci]*cos(theta[ci]))) > kmaxx) || (((k[ci]*sin(theta[ci])) + (k[ci]*sin(theta[ci]))) > kmaxy)){ 
			usum2 += 0.25*(A[ci] * G/ omega[ci])*(A[ci] * G / omega[ci])*(2*k[ci]*cos(theta[ci])) * ((D_nm_plus/(2*omega[ci])) * cos(2*phi_i) * cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth)); // sum term only. difference term = 0 on diagonal
		//}	
		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			int cm = m;
			
			double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

			double gamma_nm = cos(theta[ci] - theta[cm]);
			double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
			double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

			double Rm = k[cm] * tanh(k[cm] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[ci] * k[cm] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci] * k[ci] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[ci] * k[cm] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci]*k[ci] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			//if ((k[ci] + k[cm]) <= kmax){
			//if ((((k[ci]*cos(theta[ci])) + (k[cm]*cos(theta[cm]))) > kmaxx) || (((k[ci]*sin(theta[ci])) + (k[cm]*sin(theta[cm]))) > kmaxy)){ 
				// Sum term
				usum2 += 0.5*(A[ci] * G/ omega[ci])*(A[cm] * G / omega[cm])*(k[ci]*cos(theta[ci]+ (mtheta * PI / 180.)) + k[cm]*cos(theta[cm]+ (mtheta * PI / 180.))) * ((D_nm_plus/(omega[ci]+omega[cm])) * cos(phi_i+phi_m) * cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth));	
			//if ((k[ci] - k[cm] ) != 0. ){
			if (k[ci] != k[cm]){ 
			// Diff term
				usum2 += 0.5*(A[ci] * G / omega[ci]) * (A[cm]* G/ omega[cm]) * (k[ci]*cos(theta[ci]+ (mtheta * PI / 180.))-k[cm]*cos(theta[cm]+ (mtheta * PI / 180.))) * ((D_nm_minus/(omega[ci]-omega[cm])) * cos(phi_i-phi_m) * cosh(k_nm_minus * (zz+depth))/cosh(k_nm_minus*depth));

			//std::cout << i << " " << m << " D:"<< D_nm_minus << std::endl;
			}	
		}
	}
	return usum2;
}


/* Second order horizontal velocity V for a sinus wave */
double Irregular::v2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double vsum2 = 0;

	for (int i = 0; i < nfreq * ndir; i++) {
		
		// Second order
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		//double gamma_nm = cos(theta[ci] - theta[ci]);
		//double k_nm_plus = sqrt(k[ci] * k[ci] + k[ci] * k[ci] + (2. * k[ci] * k[ci] * gamma_nm));
		//double D_nm_plus = ( 8.*Rsn*Rsn * (k[ci] * k[ci] - Rn * Rn)) / (4.* Rsn * Rsn - k_nm_plus * tanh(k_nm_plus * depth)) +
		// (2.* Rsn * (Rsn*(k[ci]*k[ci] - Rn * Rn) + Rsn*(k[ci]*k[ci] - Rn*Rn)))/(4.* Rsn* Rsn - k_nm_plus * tanh(k_nm_plus * depth));
		// Diagonal terms
		double k_nm_plus = 2*k[ci];
		double D_nm_plus = ( 8.*Rsn*Rsn * (k[ci] * k[ci] - Rn * Rn)) / (4.* Rsn * Rsn - k_nm_plus * tanh(k_nm_plus * depth)) +
		 (2.* Rsn * (Rsn*(k[ci]*k[ci] - Rn * Rn) + Rsn*(k[ci]*k[ci] - Rn*Rn)))/(4.* Rsn* Rsn - k_nm_plus * tanh(k_nm_plus * depth));
		
		
		vsum2 += 0.25*(A[ci] * G/ omega[ci])*(A[ci] * G / omega[ci])*(2*k[ci]*sin(theta[ci])) * ((D_nm_plus/(2*omega[ci])) * cos(2*phi_i) * cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth)); // sum term only. difference term = 0 on diagonal
		
		for (int m = i + 1; m < bwlim[i]; m++) {
			int cm = m;
			
			double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

			double gamma_nm = cos(theta[ci] - theta[cm]);
			double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
			double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

			double Rm = k[cm] * tanh(k[cm] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[ci] * k[cm] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci] * k[ci] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[ci] * k[cm] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci]*k[ci] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			// Sum term
			vsum2 += 0.5*(A[ci] * G / omega[ci]) * (A[cm] * G / omega[cm])*(k[ci]*sin(theta[ci]+ (mtheta * PI / 180.))+k[cm]*sin(theta[cm]+ (mtheta * PI / 180.))) * ((D_nm_plus/(omega[ci]+omega[cm])) * cos(phi_i+phi_m) * cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth));
			if (k[ci] != k[cm]){   
				// Diff term
				vsum2 += 0.5*(A[ci] * G / omega[ci]) * (A[cm]* G/ omega[cm]) * (k[ci]*sin(theta[ci]+ (mtheta * PI / 180.))-k[cm]*sin(theta[cm]+ (mtheta * PI / 180.))) * ((D_nm_minus/(omega[ci]-omega[cm])) * cos(phi_i-phi_m) * cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth));
			}	
				
		}
		
	}
	return vsum2;
}

/* Second order vertical velocity component W for a sinus wave */
double Irregular::w2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2 = 0;
	double kmax = *max_element(k.begin(), k.end());
	double kmin = *min_element(k.begin(), k.end());
	for (int i = 0; i < nfreq * ndir; i++) {
		
		// Second order
		int ci = i;
		double phi_i = k[ci] * (cos(theta[ci] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[ci] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double k_nm_plus = 2*k[ci];
		double D_nm_plus = ( 8.*Rsn*Rsn * (k[ci] * k[ci] - Rn * Rn)) / (4.* Rsn * Rsn - k_nm_plus * tanh(k_nm_plus * depth)) +
		 (2.* Rsn * (Rsn*(k[ci]*k[ci] - Rn * Rn) + Rsn*(k[ci]*k[ci] - Rn*Rn)))/(4.* Rsn* Rsn - k_nm_plus * tanh(k_nm_plus * depth));
		if (2*k[ci] <= kmax){
			wsum2 += 0.25*(A[ci] * G/ omega[ci])*(A[ci] * G / omega[ci])* k_nm_plus * ((D_nm_plus/(2*omega[ci])) * sin(2*phi_i) * sinh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth)); // sum term only. difference term = 0 on diagonal
		}	
		for (int m = i + 1; m < nfreq * ndir; m++) {
			int cm = m;
			
			double phi_m = k[cm] * (cos(theta[cm] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[cm] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

			double gamma_nm = cos(theta[ci] - theta[cm]);
			double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
			double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

			double Rm = k[cm] * tanh(k[cm] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[ci] * k[cm] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci] * k[ci] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[ci] * k[cm] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci]*k[ci] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			// Sum term
			if ((k[ci] + k[cm]) <= kmax){
				wsum2 += 0.5*(A[ci] * G / omega[ci]) * (A[cm] * G / omega[cm])* k_nm_plus * ((D_nm_plus/(omega[ci]+omega[cm])) * sin(phi_i+phi_m) * sinh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth));
			}	
			if ((k[ci] - k[cm] ) != 0. ){   
				// Diff term
				wsum2 += 0.5*(A[ci] * G / omega[ci]) * (A[cm] * G/ omega[cm]) * k_nm_minus * ((D_nm_minus/(omega[ci]-omega[cm])) * sin(phi_i-phi_m) * sinh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth));
			}

		}
	}
	return wsum2;
}

/* Second order velocity vector (All three components) */
std::vector<double> Irregular::uvw2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double usum2 = 0., vsum2 = 0., wsum2 = 0;

	for (int i = 0; i < nfreq * ndir; i++) {
	
		// Second order parts
		int ci = i;
		double thi = theta[ci]+ (mtheta * PI / 180.);
		double phi_i = k[ci] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[ci] * t + omega[ci] * tofmax + phase[ci];
		double Rn = k[ci] * tanh(k[ci] * depth);
		double Rsn = sqrt(Rn);
		
		// Diagonal terms		
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[ci] * k[ci] - Rn * Rn))/(2.* Rsn* Rsn - k[ci] * tanh(2*k[ci] * depth));		
		// Sum term

		double AA = 0.25*(A[ci] * G/ omega[ci])*(A[ci] * G / omega[ci])*(D_nn_plus/(2*omega[ci]));
		usum2 += AA*(2*k[ci]*cos(thi)) * cos(2*phi_i) * cosh(2*k[ci]*(zz+depth)) / cosh(2*k[ci]*depth);
		vsum2 += AA*(2*k[ci]*sin(thi)) * cos(2*phi_i) * cosh(2*k[ci]*(zz+depth)) / cosh(2*k[ci]*depth); 		
		wsum2 += AA* 2*k[ci] * sin(2*phi_i) * sinh(2*k[ci]*(zz+depth)) / cosh(2*k[ci]*depth); 

		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			int cm = m;
			
			double thm = theta[cm]+ (mtheta * PI / 180.);
			
			double phi_m = k[cm] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[cm] * t + omega[cm] * tofmax + phase[cm];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2. * k[ci] * k[cm] * gamma_nm));
			double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2. * k[ci] * k[cm] * gamma_nm));

			double Rm = k[cm] * tanh(k[cm] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[ci] * k[cm] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci] * k[ci] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[ci] * k[cm] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[cm]*k[cm] - Rm * Rm) + Rsm*(k[ci]*k[ci] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[ci] * G / omega[ci]) * (A[cm] * G / omega[cm]);
			double BB = (D_nm_plus/(omega[ci]+omega[cm]));
			
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			double sinhcosh_plus = sinh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			
			// Sum terms
			usum2 += AA * BB * (k[ci]*cos(thi) + k[cm]*cos(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
			vsum2 += AA * BB * (k[ci]*sin(thi) + k[cm]*sin(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
			wsum2 += AA * BB * k_nm_plus * sin(phi_i+phi_m) * sinhcosh_plus;
	
			if (k[ci] != k[cm]){
				double CC = (D_nm_minus/(omega[ci]-omega[cm])); 
				double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				double sinhcosh_minus = sinh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				
				// Diff terms
				usum2 += AA * CC * (k[ci]*cos(thi)-k[cm]*cos(thm)) * cos(phi_i-phi_m) * coshcosh_minus;
				vsum2 += AA * CC * (k[ci]*sin(thi)-k[cm]*sin(thm)) * cos(phi_i-phi_m) * coshcosh_minus;
				wsum2 += AA * CC * k_nm_minus * sin(phi_i-phi_m) * sinhcosh_minus;
			}	
		}
	}
	
	std::vector<double> res;
	res.push_back(usum2);
	res.push_back(vsum2);
	res.push_back(wsum2);

	return res;
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

/* Function for calculating the mean wave length based integration of the wave spectrum*/
double Irregular::bandwidth_estimator()
{
	// check 
	if (ndir > 1) {
		std::cerr << "Automatic calculation of bandwidth parameter only supported for 2D wave spectra at the moment. Please specify bandiwidth manually. Sorry for the inconvenience" << std::endl;
		exit(-1);
	}

	double w_t, dw, S;
	double m0 = 0.;
	double m1 = 0.;

	for (int i = 0; i < (nfreq - 1); i++) {
		w_t = (omega[i + 1] + omega[i]) / 2.;
		dw = omega[i + 1] - omega[i];
		S = pow((A[i + 1] + A[i]) * 0.5, 2.) / (2. * dw);
		m0 += w_t * S;
		m1 += w_t * S * w_t;
		
	}

	double bw = 0.7 * ( m1 / m0);
	
	return bw;
}

// when called, writes spectral components stored in irregular class to dump file.
void Irregular::dumpSpectralComponents() {
	std::string fpath = "./spectral_components.dat";
	FILE* fp = fopen(fpath.c_str(), "w");
	fprintf(fp, "# spectral wave data output\n");
	fprintf(fp, "# [irregular wave components]\n");
	fprintf(fp, "# nfreq %d\n",nfreq);
	fprintf(fp, "# ndir 0\n");
	fprintf(fp, "# omega [rad/s]    a[m]       k[rad/m]        phi[rad]    theta[rad]\n");
	for (int i = 0; i < nfreq; i++) {
		fprintf(fp, "%12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n", omega[i], A[i], k[i], phase[i], theta[i]);
	}
	fclose(fp);
}



