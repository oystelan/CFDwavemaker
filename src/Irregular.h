#ifndef Irregular_H
#define Irregular_H

#include <iostream>
#include <vector>
#include <iomanip>
#ifndef SMALL
#define SMALL 1.e-12
#endif
#ifndef HUGE
#define HUGE 1.e12
#endif
#ifndef PI
#define PI 3.1415926535897
#endif

struct UVWP {
    double u;
    double v;
    double w;
    double p;
};

struct PairCacheEntry {
    int i;
    int m;

    bool has_plus;
    bool has_minus;

    double gamma_nm;

    double k_nm_plus;
    double k_nm_minus;

    double D_nm_plus;
    double D_nm_minus;

    double alpha_nm_plus;
    double alpha_nm_minus;

    // Smit et al. (2017) s-coordinate interaction coefficients
    double H_nm_plus;
    double H_nm_minus;
};

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
		order = 2;
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
	int vertical_theory = 0; // 0 = conventional Sharma & Dean, 1 = Smit s-coordinate
	double ampl, depth, mtheta, tofmax, fpoint[2], rho, g;
	double swl = 0.; // still water level
	double dw_bandwidth = 100.; // default, frequencies bandwidth will interact
	double dw_cutoff = 100000.;
	bool individual_ramp = false; // rampup time for each frequency component

	// Declaration of vectors to store spectral data components
	std::vector<double> omega;
	std::vector<double> A;
	std::vector<double> k;
	std::vector<double> theta;
	std::vector<double> phase;
	std::vector<double> freq_ramptime;
	std::vector<double> freq_starttime;
	std::vector<double> freq_duration;
	std::vector<int> bwlim;

	// vectors for storing initial parts of second order computation that is independent of time and space
	std::vector<double> eta2_diag;
	std::vector<double> eta2_sum;
	std::vector<double> eta2_diff;
	std::vector<double> uvw2_diag;
	std::vector<double> uvw2_sum;
	std::vector<double> uvw2_diff;

	// New uvw2 caches
    std::vector<double> cth_cache;
    std::vector<double> sth_cache;
    std::vector<double> R_cache;
    std::vector<double> Rs_cache;
    std::vector<double> ag_over_om_cache;
    std::vector<double> cosh2kd_cache;
    std::vector<double> phase0_cache;
	double kmax_cache = 0.0;

	std::vector<PairCacheEntry> pair_cache;

    void rebuild_component_cache();
	void rebuild_second_order_cache();
	void rebuild_pair_cache();

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
	double dp2(double t, double xx, double yy, double zz);
	UVWP uvw2(double t, double xx, double yy, double zz);

	// gradients
	double phi1_dxdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component U */
	double phi1_dydz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component V */
	double phi1_dzdz(double t, double xx, double yy); /* vertical velocity gradient at z=0 for velocity component W */
	double dp1_dz(double t, double xx, double yy); /* vertical pressure gradient at z=0*/

	// main kinematics functions
	double eta(double t, double x, double y);
	double u(double t, double x, double y, double z); 
	double v(double t, double x, double y, double z);
	double w(double t, double x, double y, double z);
	double profileX(int ind, double x, double y, double z);
	double profileZ(int ind, double x, double y, double z);
	// Smit s-coordinate stretched profiles: argument becomes k*h*(h+z)/(h+eta)
	double profileX_s(int ind, double z, double eta);
	double profileZ_s(int ind, double z, double eta);
	double profileX_s2(double k_nm, double z, double eta); // for combined wavenumber
	double profileZ_s2(double k_nm, double z, double eta);
	// Smit s-coordinate kinematics
	double eta_smit(double t, double x, double y);   // derived from Smit's Eq. 10 (per-pair C from H)
	double u_smit(double t, double x, double y, double z);
	double v_smit(double t, double x, double y, double z);
	double w_smit(double t, double x, double y, double z);
	double dp_smit(double t, double x, double y, double z);
	double dp(double t, double x, double y, double z);

	void normalize_data();
	double trapz(double x[], double y[], int n);
	double interpolate(std::vector<double>& xData, std::vector<double>& yData, double x, bool extrapolate);
	
	// Calculation of various spectral properties
	double phase_speed(int opt);
	std::vector<double> spectral_moments();
	double mean_wave_length(int opt);
	double mean_wave_period(int opt);
	double bandwidth_estimator();

	// ramp function for each individual frequency component
	double freq_ramp(double t, double ramptime, double tstart, double duration);

	void dumpSpectralComponents();
	void print_spectral_parameters();

	// new functions for faster computation of second order terms
	void second_order_init();
	double eta2_boost(double t, double xx, double yy);
	double u2_boost(double t, double xx, double yy, double zz);
	double v2_boost(double t, double xx, double yy, double zz);
	double w2_boost(double t, double xx, double yy, double zz);
};

#endif
