#ifndef Irregular_H
#define Irregular_H

#include <iostream>
#include "Wavespectra.h"


class Irregular {
private:
	double phi1_pot(double t, double xx, double yy, double zz);
	double phi2_pot(double t, double xx, double yy, double zz);

	double sum(double ll[], int nsum);

	Wavespectra wavespec;

public:
	Irregular() {
		ampl = 1.;
		normalize = 0;
		mtheta = 0.;
		extmet = 0;
		pertmet = 0;
		bandwidth = 1000;
		tofmax = 0.;
		fpoint[0] = 0.;
		fpoint[1] = 0.;
	};

	~Irregular() {
		std::cout << "Irregular class destroyed." << std::endl;
		delete[] omega, A, k, theta, phase;
	};

	// Variables
	int nfreq, ndir, extmet, bandwidth, normalize;
	int pertmet = 0;
	double ampl, depth, mtheta, tofmax, fpoint[2];
	double swl; // still water level

	// Declaration of pointers where data will be stored
	double* omega;
	double* A;
	double* k;
	double* theta;
	double* phase;

	// First order
	double eta1(double t, double xx, double yy); // wave elevation
	double u1(double t, double xx, double yy, double zz); // velocity component x
	double v1(double t, double xx, double yy, double zz); // velocity component y
	double w1(double t, double xx, double yy, double zz); // velocity component z
	double dp1(double t, double xx, double yy, double zz); // dynamic pressure component	

	// second order
	double eta2(double t, double xx, double yy);
	double u2(double t, double xx, double yy, double zz);
	double v2(double t, double xx, double yy, double zz);
	double w2(double t, double xx, double yy, double zz);

	// gradients
	double phi1_dxdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component U */
	double phi1_dydz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component V */
	double phi1_dzdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component W */


	// main kinematics functions
	double eta(double t, double x, double y);
	double u(double t, double x, double y, double z); 
	double v(double t, double x, double y, double z);
	double w(double t, double x, double y, double z);
	double dp(double t, double x, double y, double z);

	void normalize_data();

	void allocate_arrays(int array_length) {
		// allocate memory for storage of spectral component data
		omega = new double[array_length];
		k = new double[array_length];
		phase = new double[array_length];
		A = new double[array_length];
		theta = new double[array_length];
	}	
};

#endif