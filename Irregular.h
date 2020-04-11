#ifndef Irregular_H
#define Irregular_H

#include <iostream>
#include "Wavespectra.h"
#include <vector> 

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
		extrapolation_met = 0;
		order = 2;
		bandwidth = 1000;
		tofmax = 0.;
		fpoint[0] = 0.;
		fpoint[1] = 0.;
		swl = 0.;
	};

	~Irregular() {
		std::cout << "Irregular class destroyed." << std::endl;
	};

	// Variables
	int nfreq, ndir, extrapolation_met, bandwidth, normalize;
	int order = 2; // 1 = linear airy wave theory; 2= second order wave theory
	int sloping_bottom = 0;
	double ampl, depth, mtheta, tofmax, fpoint[2];
	double swl; // still water level

	// Declaration of vectors to store spectral data components
	std::vector<double> omega;
	std::vector<double> A;
	std::vector<double> k;
	std::vector<double> theta;
	std::vector<double> phase;

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
	double profileX(int ind, double x, double y, double z);
	double profileZ(int ind, double x, double y, double z);
	double dp(double t, double x, double y, double z);

	void normalize_data();
};

#endif