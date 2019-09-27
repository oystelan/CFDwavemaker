// This program has the soul purpose of providing wave kinematics input to any type
// of CFD program which can be linked up a dynamic link library.
// The program is created and updated by the Oeystein Lande. The link library
// compiles on both linux and windows. The following wave theory/types are currently
// supported:
// linear wave theory (2D and 3D), second order wave theory (2D and 3D) (sharma &
// dean), wave paddle theory (2D only at the moment)
//
//
// Current version: v1.0
// Date: 2017-09-30
// --------------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <ctime>
#include "CFDwavemaker.h"
#include "omp.h"
#include <direct.h>
 
#include "input_file_reader.h"

//#include <fftw3.h>

extern "C" {
#include "Stokes5.h"
}

using namespace std;

double *UX;
double *UY;
double *UZ;
double *UXL;
double *UYL;
double *UZL;
double *ETA;

Stokes5 wave; // declare wave

int initialized;
int initsurf = 0;
int initkin = 0;
double dx, dy, dz, dxl, dyl, dzl;


#define GetCurrentDir _getcwd

/*void testfft() {
	fftw_complex *in, *out;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 64);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 64);

	fftw_free(in); fftw_free(out);

}*/

double sum(double ll[], int nsum) {
	double ss = 0.0;
	for (int i = 0; i < nsum; i++) {
		ss += ll[i];
	}
	return ss;
}

void wait(int seconds)
{
	clock_t endwait;
	endwait = clock() + seconds * CLOCKS_PER_SEC;
	while (clock() < endwait) {}
}


int check_license()
{

	int licensecheck = 1;
  // Date for program to stop working
  int expyear = 2020;
  int expmonth = 12;
  int expday = 31;

  time_t t = time(0);   // get time now

	#if defined(_MSC_VER)

      struct tm now;

      localtime_s(&now, &t);



    	cout << (now.tm_year + 1900) << '-'
    		<< (now.tm_mon + 1) << '-'
    		<< now.tm_mday
    		<< endl;
    	if ((now.tm_year + 1900) > expyear) {
    		licensecheck = 0;
    	}
    	else if ((now.tm_year + 1900) == expyear) {
    		if ((now.tm_mon + 1) > expmonth) {
    			licensecheck = 0;
    		}
    	}
    	else if ((now.tm_year + 1900) == expyear) {
    		if ((now.tm_mon + 1) == expmonth) {
    			if (now.tm_mday > expday) {
    				licensecheck = 0;
    			}
    		}
    	}
  #else
      struct tm * now;
      //time(&t)
      now = localtime(&t);

    	cout << (now->tm_year + 1900) << '-'
    		<< (now->tm_mon + 1) << '-'
    		<< now->tm_mday
    		<< endl;
    	if ((now->tm_year + 1900) > expyear) {
    		licensecheck = 0;
    	}
    	else if ((now->tm_year + 1900) == expyear) {
    		if ((now->tm_mon + 1) > expmonth) {
    			licensecheck = 0;
    		}
    	}
    	else if ((now->tm_year + 1900) == expyear) {
    		if ((now->tm_mon + 1) == expmonth) {
    			if (now->tm_mday > expday) {
    				licensecheck = 0;
    			}
    		}
    	}

  #endif
	return licensecheck;
}


string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}
int main() {
	cout << GetCurrentWorkingDir() << endl;
	read_inputdata_v2();
}






//Define some useful functions
/* Rampfunction */
// NB: Not yet implemented inverse ramp
static double ramp(double x,double xsign, double xstart, double xend) {

	if (xsign > 0.){
		if (x <= xstart) {
			return 1.0;
		}
		else if (x >= xend) {
			return 0.0;
		}
		else {
			return 1. - ((x - xstart) / (xend - xstart));
		}
	}
	else {
		return 1.0;
	}
}


//Define some useful functions
/* Rampfunction */
// NB: Not yet implemented inverse ramp
static double timeramp(double t, double tsign, double tstart, double tend) {

	if (tsign > 0.) {
		if (t <= tstart) {
			return 0.0;
		}
		else if (t >= tend) {
			return 1.0;
		}
		else {
			return ((t - tstart) / (tend - tstart));
		}
	}
	else {
		return 1.0;
	}
}


// Linear interpolation function

int findNearestNeighbourIndex(double value, double *x, int len)
{
	double dist;
	int idx;
	int i;
	idx = -1;
	dist = DBL_MAX;
	for (i = 0; i < len; i++) {
		double newDist = value - x[i];
		if (newDist >= 0 && newDist < dist) {
			dist = newDist;
			idx = i;
		}
	}
	return idx;
}
double interp1(double *x, int x_tam, double *y, double xx)
{
	double dx, dy, slope, intercept, yy;
	int indiceEnVector;


	indiceEnVector = findNearestNeighbourIndex(xx, x, x_tam);
	if (indiceEnVector != -1) {
		dx = x[indiceEnVector + 1] - x[indiceEnVector];
		dy = y[indiceEnVector + 1] - y[indiceEnVector];
		slope = dy / dx;
		intercept = y[indiceEnVector] - x[indiceEnVector] * slope;
		yy = slope * xx + intercept;
	}
	else
		yy = DBL_MAX;

	return yy;
}


//--------------------------------------------------------------------
// -------------------------------------------------------------------
//
//  WAVEFUNCTIONS
//
//



/* wave elevation for a sinus wave */
double waveelev(double t, double xx, double yy) {

	double welev = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {

		phi = w[i] * tofmax+phas[i];
		welev += Ampspec[i] * D[i] *cos(k[i] * (cos(thetaA[i]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);

		//welev += A*cos(w[i] * t - k[i] * (cos(thetaA[i])*(xx - fpoint[0]) + sin(thetaA[i])*(yy - fpoint[1])) + phi);
	}


	return welev;

}


double waveelev_2order(double t, double xx, double yy) {



	//double eta1_t = 0;
	double eta2_t = 0;

	for (int i = 0; i < nfreq-1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));
				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
				}
				double alpha_nm_minus = (((w[ci] / w[cm]) + (w[cm] / w[ci])) + (G*G / (w[ci] * w[cm]))*
					((D_nm_minus - k[ci] * k[cm] * (gamma_nm + tanh(k[ci] * depth)*tanh(k[cm] * depth))) / (w[ci] * w[cm])));
				double alpha_nm_plus = (((w[ci] / w[cm]) + (w[cm] / w[ci])) + (G*G / (w[ci] * w[cm]))*
					((D_nm_plus - k[ci] * k[cm] * (gamma_nm - tanh(k[ci] * depth)*tanh(k[cm] * depth))) / (w[ci] * w[cm])));

				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				eta2_t += ((Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm])/(2.*G)) * (alpha_nm_minus*cos(phi_i - phi_m) + alpha_nm_plus*cos(phi_i + phi_m));



			}
		}
	}

	return eta2_t;

}


/* Second order horizontal velocity U for a sinus wave */
double uu_2order(double t, double xx, double yy, double zz) {



	//double eta1_t = 0;
	double usum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				usum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*cos(phi_i - phi_m)*(k[ci]*cos(thetaA[ci]+ (mtheta*PI / 180.))-k[cm]*cos(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth)) +
					beta_nm_plus*cos(phi_i + phi_m)*(k[ci] * cos(thetaA[ci]+ (mtheta*PI / 180.)) + k[cm] * cos(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



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
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}
				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				vsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*cos(phi_i - phi_m)*(k[ci] * sin(thetaA[ci]+ (mtheta*PI / 180.)) - k[cm] * sin(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_minus*(zz + depth)) / cosh(k_nm_minus*depth)) +
					beta_nm_plus*cos(phi_i + phi_m)*(k[ci] * sin(thetaA[ci]+ (mtheta*PI / 180.)) + k[cm] * sin(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



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
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adjusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));

				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				wsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*k_nm_minus*sin(phi_i - phi_m)*(sinh(k_nm_minus*(zz + depth)) / cosh(k_nm_minus*depth)) +
					beta_nm_plus*k_nm_plus*sin(phi_i + phi_m)*(sinh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



			}
		}
	}
	return wsum2;
}

/* Second order component of velocity potential 
Warning: Unsure if this is finished...should be checked thoroughly at some point*/
double phi_pot_2order(double t, double xx, double yy, double zz) {
	//double eta1_t = 0;
	double phisum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				phisum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*cos(phi_i - phi_m)*(k[ci] * cos(thetaA[ci] + (mtheta*PI / 180.)) - k[cm] * cos(thetaA[cm] + (mtheta*PI / 180.)))*(cosh(k_nm_minus*(zz + depth)) / cosh(k_nm_minus*depth)) +
						beta_nm_plus*cos(phi_i + phi_m)*(k[ci] * cos(thetaA[ci] + (mtheta*PI / 180.)) + k[cm] * cos(thetaA[cm] + (mtheta*PI / 180.)))*(cosh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



			}
		}
	}
	return phisum2;
}

// First order velocity potential
double phi_pot(double t, double xx, double yy, double zz) {

	double phisum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		phisum += Ampspec[i] * D[i] * (w[i] / k[i]) * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*sin(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return phisum;

}

double uu(double t, double xx, double yy, double zz) {

	double usum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		usum += cos(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return usum;

}


/* horizontal velocity U for a sinus wave */
double vv(double t, double xx, double yy, double zz) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		vsum += sin(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double ww(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;



	for (int i = 0; i< ndir*nfreq; i++) {

		phi = w[i] * tofmax + phas[i];
		wsum += D[i] * Ampspec[i] * w[i] * (sinh(k[i] * (zz + depth)) / sinh(k[i] * depth))*sin(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return wsum;

}

double pp(double t, double xx, double yy, double zz) {

	double psum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		psum += Ampspec[i] * D[i] * RHO * G * (cosh(k[i] * (zz + depth)) / cosh(k[i] * depth))*cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return psum;

}

/* vertical velocity gradient at z=0 for velocity component U */
double phi_dxdz(double t, double xx, double yy) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		usum += cos(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * k[i] * cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}
	return usum;

}



/* vertical velocity gradient at z=0 for velocity component V */
double phi_dydz(double t, double xx, double yy) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		vsum += sin(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * k[i] *cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}
	return vsum;

}


/* vertical velocity gradient at z=0 for velocity component W */
double phi_dzdz(double t, double xx, double yy) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {

		phi = w[i] * tofmax + phas[i];
		wsum += D[i] * Ampspec[i] * w[i] * k[i] * (cosh(k[i] * depth) / sinh(k[i] * depth))*sin(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}
	return wsum;

}


/* Horizontal velocity taken directly from the timeseries*/
double u_piston(double t) {
	double ux = interp1(PD_time, n_timesteps, PD_velo, t);

	return ux+alpha_u*ux;
}

/* Wave elevation taken directly from piston timeseries*/
double wave_elev_piston(double t) {
	return alpha_z*interp1(PD_time, n_timesteps, PD_eta, t);
}


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void initialize_kinematics(double tpt) {
	// Allocating memory for storage of surface elevation and velocities
	UX = new double[NX*NY*NZ];
	UY = new double[NX*NY*NZ];
	UZ = new double[NX*NY*NZ];

	UXL = new double[NXL*NYL*NZL];
	UYL = new double[NXL*NYL*NZL];
	UZL = new double[NXL*NYL*NZL];

	cout << "Memory allocation successful for storage of kinematics." << endl;

	dx = (domainsize[1] - domainsize[0]) / double(NX-1);
	dy = (domainsize[3] - domainsize[2]) / double(NY-1);
	dz = (domainsize[6] - domainsize[5]) / double(NZ-1);
	
	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel // start parallell initialization
	{
		#pragma omp master
		cout << "Number of available threads: " << omp_get_num_threads() << endl;

		double xpt, ypt, zpt;
		double eta_temp;

		// Main grid
		#pragma omp for
		for (int i = 0; i < NX; i++) {
			xpt = domainsize[0] + dx*i;
			for (int j = 0; j < NY; j++) {
				ypt = domainsize[2] + dy*j;
				eta_temp = waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt);

				double Ux0 = uu(tpt, xpt, ypt, 0.0) + uu_2order(tpt, xpt, ypt, 0.0);
				double Uy0 = vv(tpt, xpt, ypt, 0.0) + vv_2order(tpt, xpt, ypt, 0.0);
				double Uz0 = ww(tpt, xpt, ypt, 0.0) + ww_2order(tpt, xpt, ypt, 0.0);

				double PHI_dxdz = phi_dxdz(tpt, xpt, ypt);
				double PHI_dydz = phi_dydz(tpt, xpt, ypt);
				double PHI_dzdz = phi_dzdz(tpt, xpt, ypt);

				for (int m = 0; m < NZ; m++) {
					zpt = domainsize[5] + dz*m;
					if (zpt > (eta_temp + dz)) {
						UX[i*NY*NZ + j*NZ + m] = 0.0;
						UY[i*NY*NZ + j*NZ + m] = 0.0;
						UZ[i*NY*NZ + j*NZ + m] = 0.0;
					}
					else if (zpt > 0.) {
						UX[i*NY*NZ + j*NZ + m] = Ux0 + PHI_dxdz*zpt;
						UY[i*NY*NZ + j*NZ + m] = Uy0 + PHI_dydz*zpt;
						UZ[i*NY*NZ + j*NZ + m] = Uz0 + PHI_dzdz*zpt;
					}
					else {
						UX[i*NY*NZ + j*NZ + m] = uu(tpt, xpt, ypt, zpt) + uu_2order(tpt, xpt, ypt, zpt);
						UY[i*NY*NZ + j*NZ + m] = vv(tpt, xpt, ypt, zpt) + vv_2order(tpt, xpt, ypt, zpt);
						UZ[i*NY*NZ + j*NZ + m] = ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt);
					}
					/*UX[i*NY*NZ + j*NZ + m] = uu(tpt, xpt, ypt, zpt);
					UY[i*NY*NZ + j*NZ + m] = vv(tpt, xpt, ypt, zpt);
					UZ[i*NY*NZ + j*NZ + m] = ww(tpt, xpt, ypt, zpt);*/
				}
			}
		}
	} // End parallel initialization

	cout << "Generation of upper domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	cout << "Initialization time: " << dd << " seconds." << endl;

	dd = omp_get_wtime();
	#pragma omp parallel // start parallel initialization
	{
		double xpt, ypt, zpt;
		// Secondary grid (coarse res at depth)
		dxl = (domainsize[1] - domainsize[0]) / double(NXL - 1);
		dyl = (domainsize[3] - domainsize[2]) / double(NYL - 1);
		dzl = (domainsize[5] - domainsize[4]) / double(NZL - 1);
		#pragma omp for
		for (int i = 0; i < NXL; i++) {
			xpt = domainsize[0] + dxl*i;
			for (int j = 0; j < NYL; j++) {
				ypt = domainsize[2] + dyl*j;
				for (int m = 0; m < NZL; m++) {
					zpt = domainsize[4] + dzl*m;
					UXL[i*NYL*NZL + j*NZL + m] = uu(tpt, xpt, ypt, zpt) + uu_2order(tpt, xpt, ypt, zpt);
					UYL[i*NYL*NZL + j*NZL + m] = vv(tpt, xpt, ypt, zpt) + vv_2order(tpt, xpt, ypt, zpt);
					UZL[i*NYL*NZL + j*NZL + m] = ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt);
				}
			}
		}
	} // End parallel
	cout << "Generation of lower domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	cout << "Initialization time: " << dd << " seconds." << endl;
	cout << "Interpolation can commence..." << endl;
	initkin = 1;
}


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void initialize_surface_elevation(double tpt) {

	// Allocating memory for storage of surface elevation and velocities
	ETA = new double[NX*NY];

	cout << "Memory allocation successful for Surface elevation storage." << endl;

	dx = (domainsize[1] - domainsize[0]) / double(NX - 1);
	dy = (domainsize[3] - domainsize[2]) / double(NY - 1);
	
	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel
	{
		double xpt, ypt;
		// Main grid
		#pragma omp for
		for (int i = 0; i < NX; i++) {
			xpt = domainsize[0] + dx*i;
			for (int j = 0; j < NY; j++) {
				ypt = domainsize[2] + dy*j;
				ETA[i*NY + j] = waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt);
			}
		}
	}
	dd = omp_get_wtime()-dd;

	cout << "Surface Elevation generated successfully. ";
	cout << "Initialization time: " << dd << " seconds." << endl;
	initsurf = 1;
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double trilinear_interpolation(double *VAR, double xpt, double ypt, double zpt) {
	double nxp = min(double(NX), max(0., (xpt - domainsize[0]) / dx));
	double nyp = min(double(NY), max(0., (ypt - domainsize[2]) / dy));
	double nzp = min(double(NZ), max(0., (zpt - domainsize[5]) / dz));

	double C000 = VAR[int(floor(nxp)*NY*NZ + floor(nyp)*NZ + floor(nzp))];
	double C001 = VAR[int(floor(nxp)*NY*NZ + floor(nyp)*NZ + ceil(nzp))];
	double C010 = VAR[int(floor(nxp)*NY*NZ + ceil(nyp)*NZ + floor(nzp))];
	double C011 = VAR[int(floor(nxp)*NY*NZ + ceil(nyp)*NZ + ceil(nzp))];
	double C100 = VAR[int(ceil(nxp)*NY*NZ + floor(nyp)*NZ + floor(nzp))];
	double C101 = VAR[int(ceil(nxp)*NY*NZ + floor(nyp)*NZ + ceil(nzp))];
	double C110 = VAR[int(ceil(nxp)*NY*NZ + ceil(nyp)*NZ + floor(nzp))];
	double C111 = VAR[int(ceil(nxp)*NY*NZ + ceil(nyp)*NZ + ceil(nzp))];
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double zd = nzp - floor(nzp);

	double C00 = C000*(1. - xd) + C100*xd;
	double C01 = C001*(1. - xd) + C101*xd;
	double C10 = C010*(1. - xd) + C110*xd;
	double C11 = C011*(1. - xd) + C111*xd;

	double C0 = C00*(1. - yd) + C10*yd;
	double C1 = C01*(1. - yd) + C11*yd;

	return C0*(1. - zd) + C1*zd;
}
/* Function for trilinear interpolation on a cartesian evenly spaced mesh on the lower part of the domain*/
double trilinear_interpolationL(double *VAR, double xpt, double ypt, double zpt) {
	
	double nxp = min(double(NXL),max(0.,(xpt - domainsize[0]) / dxl));
	double nyp = min(double(NYL), max(0., (ypt - domainsize[2]) / dyl));
	double nzp = min(double(NZL), max(0., (zpt - domainsize[4]) / dzl));

	double C000 = VAR[int(floor(nxp)*NYL*NZL + floor(nyp)*NZL + floor(nzp))];
	double C001 = VAR[int(floor(nxp)*NYL*NZL + floor(nyp)*NZL + ceil(nzp))];
	double C010 = VAR[int(floor(nxp)*NYL*NZL + ceil(nyp)*NZL + floor(nzp))];
	double C011 = VAR[int(floor(nxp)*NYL*NZL + ceil(nyp)*NZL + ceil(nzp))];
	double C100 = VAR[int(ceil(nxp)*NYL*NZL + floor(nyp)*NZL + floor(nzp))];
	double C101 = VAR[int(ceil(nxp)*NYL*NZL + floor(nyp)*NZL + ceil(nzp))];
	double C110 = VAR[int(ceil(nxp)*NYL*NZL + ceil(nyp)*NZL + floor(nzp))];
	double C111 = VAR[int(ceil(nxp)*NYL*NZL + ceil(nyp)*NZL + ceil(nzp))];
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double zd = nzp - floor(nzp);

	double C00 = C000*(1. - xd) + C100*xd;
	double C01 = C001*(1. - xd) + C101*xd;
	double C10 = C010*(1. - xd) + C110*xd;
	double C11 = C011*(1. - xd) + C111*xd;

	double C0 = C00*(1. - yd) + C10*yd;
	double C1 = C01*(1. - yd) + C11*yd;

	return C0*(1. - zd) + C1*zd;
}

/* bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid*/
double bilinear_interpolation(double *VAR, double xpt, double ypt) {
	
	double nxp = min(double(NX), max(0.,(xpt - domainsize[0]) / dx));
	double nyp = min(double(NY), max(0.,(ypt - domainsize[2]) / dy));

	double C00 = VAR[int(floor(nxp)*NY + floor(nyp))];
	double C01 = VAR[int(floor(nxp)*NY + ceil(nyp))];
	double C10 = VAR[int(ceil(nxp)*NY + floor(nyp))];
	double C11 = VAR[int(ceil(nxp)*NY + ceil(nyp))];
	
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);

	double C0 = C00*(1. - xd) + C10*xd;
	double C1 = C01*(1. - xd) + C11*xd;

	return C0*(1. - yd) + C1*yd;
}



//extern "C" {

// For comflow 3.9 pointer are sent. For later releases, the entire variables are sent

//EXPORT double VelocityX(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
double wave_VeloX(double xpt, double ypt, double zpt, double tpt)
{

	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
	// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*uu(tpt, xpt, ypt, zpt)*timeramp(tpt,rampswitch,0.,ramp_time);
	// Linear wave theory, constant profile used above free surface
	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*uu(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
	// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (uu(tpt, xpt, ypt, zpt) + uu_2order(tpt, xpt, ypt, zpt) )*timeramp(tpt, rampswitch, 0., ramp_time);
	// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * ( uu(tpt, xpt, ypt, z) + uu_2order(tpt, xpt, ypt, z) )*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (uu(tpt, xpt, ypt, z) + uu_2order(tpt, xpt, ypt, z) )*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((uu(tpt, xpt, ypt, z) + uu_2order(tpt, xpt, ypt, z)) + phi_dxdz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return u_piston(tpt);
	case 7:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*Stokes5_u(&wave, tpt, xpt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
	case 8:
		if (initkin == 0) {
			cout << "Generating kinematics for interpolation:" << endl;
			initialize_kinematics(0.0);
		}
		/*
		// Check if coordinates are within bounding box
		if (xpt < domainsize[0] || xpt > domainsize[1]) {
			cerr << "xpt: " << xpt << " out of bounds! Please extend interpolation box boundaries in x-direction" << endl;
		}
		else if (ypt < domainsize[2] || ypt > domainsize[3]) {
			cerr << "ypt: " << ypt << " out of bounds! Please extend interpolation box boundaries in y-direction" << endl;
		}
		else if (zpt < domainsize[4] || zpt > domainsize[6]) {
			cerr << "zpt: " << zpt << " out of bounds! Please extend interpolation box boundaries in z-direction" << endl;
		}
		*/
		if (zpt < domainsize[5]) {
			return trilinear_interpolationL(UXL, xpt, ypt, zpt);
		}
		else {
			return trilinear_interpolation(UX, xpt, ypt, zpt);
		}
	default:
		return 0.0;
	}


}

//
//EXPORT double VelocityY(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
double wave_VeloY(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*vv(tpt, xpt, ypt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface
	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*vv(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * ( vv(tpt, xpt, ypt, zpt) + vv_2order(tpt, xpt, ypt, zpt) )*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * ( vv(tpt, xpt, ypt, z) + vv_2order(tpt, xpt, ypt, z) )*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (vv(tpt, xpt, ypt, z) + vv_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((vv(tpt, xpt, ypt, z) + vv_2order(tpt, xpt, ypt, z)) + phi_dydz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return 0.0;
	case 7:
		return 0.0; // todo: implement directionality of stokes5th wave
	case 8:
		if (zpt < domainsize[5]) {
			return trilinear_interpolationL(UYL, xpt, ypt, zpt);
		}
		else {
			return trilinear_interpolation(UY, xpt, ypt, zpt);
		}
	default:
		return 0.0;
	}
}



//
//EXPORT double VelocityZ(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
double wave_VeloZ(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*ww(tpt, xpt, ypt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface
	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*ww(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt))*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*(ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z)) + phi_dzdz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return 0.0;
	case 7:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*Stokes5_v(&wave, tpt, xpt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
	case 8:
		if (zpt < domainsize[5]) {
			return trilinear_interpolationL(UZL, xpt, ypt, zpt);
		}
		else {
			return trilinear_interpolation(UZ, xpt, ypt, zpt);
		}
	default:
		return 0.0;
	}
}

//
double wave_DynPres(double xpt, double ypt, double zpt, double tpt)
//EXPORT double DynamicPressure(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*pp(tpt, xpt, ypt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface

	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*pp(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
	/*
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt))*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*(ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z)) + phi_dzdz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return 0.0;
		*/
	default:
		//cout << "WARNING: HIGHER ORDER PRESSURE CALCULATION NOT YET SUPPORTED." << endl;
		return 0.0;
	}
}

//
//EXPORT double SurfaceElevation(int& ii, int& jj, double& xpt, double& ypt, double& tpt)
double wave_SurfElev(double xpt, double ypt, double tpt)
{

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		//return waveelev(tpt, xpt, ypt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*waveelev(tpt, xpt, ypt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface
	case 2:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*waveelev(tpt, xpt, ypt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt))*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 6:
		return wave_elev_piston(tpt);
	case 7:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*Stokes5_eta(&wave, tpt, xpt)*timeramp(tpt, rampswitch, 0., ramp_time);
	case 8:
		if (initsurf == 0) {
			cout << "Initializing surface elevation storage:" << endl;
			initialize_surface_elevation(0.0);
			return bilinear_interpolation(ETA, xpt, ypt);
		}
		else {
			//cout << "asking for surface elevation..." << endl;
			//cout << xpt << " " << ypt << " " << tpt << endl;
			return bilinear_interpolation(ETA, xpt, ypt);
		}
	default:
		return 0.0;
	}
}



double wave_VFrac(double xpt, double ypt, double zpt, double tpt, double delta_cell)

	// Function which calculates volume fraction
	// assumes still water level z = 0
{
	double wwelev = wave_SurfElev(xpt,ypt,tpt);

	//cout << wwelev << endl;
	if (wwelev < zpt - (delta_cell / 2.)) {
		return 0.0;
	}
	else if (wwelev > zpt + (delta_cell / 2.)) {
		return 1.0;
	}
	else {
		// Calculate volume fraction for the given cell with size delta_cell and position zpt
		return (wwelev - (zpt - (delta_cell / 2.))) / delta_cell;
	}
	//return 0;
}




//EXPORT int Init(double& tmin_in, double& tmax_in)
int wave_Initialize()
{
	cout << "---------------------------------------" << endl;
	cout << "CFD WAVEMAKER v.1.07" << endl;
	cout << "---------------------------------------" << endl;
	
	// Check if license has expired
	if (check_license() == 1){
		cout << "License checks out...carry on..." << endl;
	}
	else {
		cout << "License for CFDwavemaker has expired. Please contact Oeystein Lande to update the program." << endl << endl;
		cout << "This program will auto-distruct in \n5..." << endl;
		wait(1);
		cout << "4..." << endl;
		wait(1);
		cout << "3..." << endl;
		wait(1);
		cout << "2..." << endl;
		wait(1);
		cout << "1..." << endl;
		wait(1);
		cout << "Bang!" << endl;

		return -1;
	}
	//for (int i = 0; i < nfreq; i++) {
	//	k[i] = pow(2. * pi * f[i], 2.) / 9.81;
	//}
	int i = read_inputdata();

	return 0;
}


int wave_Cleanup()
{
	if (wavetype == 1) {
		delete[] w, Ampspec, D, k, thetaA, phas;
	}
	else if (wavetype == 2) {
		delete[] w, Ampspec, k, thetaA, phas;
	}
	else if (wavetype == 3) {
		delete[] w, Ampspec, k, thetaA, phas;
		delete[] ETA, UX, UY, UZ, index;
	}
	else if (wavetype == 4) {
		delete[] PD_time, PD_ampl, PD_velo, PD_eta;
	}
	return 0;
}
//}



// external functions used by COMFLOW
double VelocityX(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	//cout << xpt << " " << ypt << " " << zpt << " " << wave_VeloX(xpt, ypt, zpt, time) << endl;
	return wave_VeloX(xpt, ypt, zpt, time);
}
double VelocityY(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return wave_VeloY(xpt, ypt, zpt, time);
}
double VelocityZ(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return wave_VeloZ(xpt, ypt, zpt, time);
}
double DynamicPressure(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return wave_DynPres(xpt, ypt, zpt, time);
}
double SurfaceElevation(int i, int j, double xpt, double ypt, double time) {
	return wave_SurfElev(xpt, ypt, time);
}
double VolumeFraction(double xpt, double ypt, double zpt, double time, double delta_xyz) {
	return wave_VFrac(xpt, ypt, zpt, time, delta_xyz);
}
double LiquidFillRatio(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return 0.0;
}
int Init(void* fptr) {
	return wave_Initialize();
}
int Prepare(double time, int mode) {
	return 0;
}
int Cleanup() {
	return wave_Cleanup();
}

