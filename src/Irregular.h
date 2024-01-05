#ifndef Irregular_H
#define Irregular_H

#include <iostream>
#include <vector>
#include <iomanip>
#ifndef PI
#define PI 3.1415926535897
#endif

class Irregular {
private:
	double phi1_pot(double t, double xx, double yy, double zz);
	double phi2_pot(double t, double xx, double yy, double zz);

	double sum(double ll[], int nsum);

public:
	Irregular() {
		ampl = 1.;
		normalize = false;
		mtheta = 0.;
		extrapolation_met = 0;
		order = 0;
		tofmax = 0.;
		fpoint[0] = 0.;
		fpoint[1] = 0.;
		swl = 0.;
		depth = 300.;
		ndir = 1;
		g = 9.81;
		rho = 1025.0;

	};

	~Irregular() {
		//std::cout << "Irregular class destroyed." << std::endl;
	};
	// bools
	int initialized = false;

	// Variables
	int nfreq, ndir;
	bool normalize, sloping_bottom;
	int extrapolation_met = 0; // 0 = exponential all components. 1 
	int order; // 0 = linear airy wave theory; 2= second order wave theory
	double ampl, depth, mtheta, tofmax, fpoint[2], rho, g;
	double swl = 0.; // still water level
	double dw_bandwidth = 100.; // default, frequencies bandwidth will interact
	double dw_cutoff = 100000.;

	// Declaration of vectors to store spectral data components
	std::vector<double> omega;
	std::vector<double> A;
	std::vector<double> k;
	std::vector<double> theta;
	std::vector<double> phase;
	std::vector<int> bwlim;

	void print();

	// First order
	double eta1(double t, double xx, double yy); // wave elevation
	int bw_limiter(int i1); // calculates the index of the maximum bandwidth frequency
	void calculate_bwindices();
	double u1(double t, double xx, double yy, double zz); // velocity component x
	double v1(double t, double xx, double yy, double zz); // velocity component y
	double w1(double t, double xx, double yy, double zz); // velocity component z
	double dp1(double t, double xx, double yy, double zz); // dynamic pressure component	

	// second order
	double eta2(double t, double xx, double yy);
	double u2(double t, double xx, double yy, double zz);
	double v2(double t, double xx, double yy, double zz);
	double w2(double t, double xx, double yy, double zz);
	std::vector<double> uvw2(double t, double xx, double yy, double zz);

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
	double trapz(double x[], double y[], int n);
	double interpolate(std::vector<double>& xData, std::vector<double>& yData, double x, bool extrapolate);
	
	// Calculation of various spectral properties
	double phase_speed(int opt);
	double mean_wave_length(int opt);
	double mean_wave_period(int opt);
	double bandwidth_estimator();

	void dumpSpectralComponents();
};

#endif
