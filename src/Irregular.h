#ifndef Irregular_H
#define Irregular_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <complex>
#include <mutex>
#include <atomic>
#include <memory>
#include "seabed.h"
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

    // Sharma & Dean truncation flag: ki+km <= kmax. The S&D evaluation
    // paths (eta2, uvw2, factored) drop sum interactions beyond kmax; the
    // smit velocity/pressure paths must NOT (their first-order stretching
    // and chain content cannot be truncated per-pair, so consistency
    // requires the full pair set).
    bool kmax_plus;

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
	bool normalize;
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

	// Factorized second-order elevation for long-crested components on an
	// equal-dk wavenumber grid (waveinput keyword "equal_dk 1"): all pair
	// interactions collapse into wavenumber bins j*dk, so the O(N^2) pair
	// work is spent once per time level (independent of evaluation points)
	// and each eta2() call costs O(n_bins).
	bool equal_dk = false;              // set by the input file parser
	// new 4-column table with equal_dk: the k column is primary and
	// omega is derived from the dispersion relation (post-parse)
	bool k_first = false;
	// legacy 5-column table: k was supplied in the file (recomputed;
	// used only to decide whether a mismatch warning makes sense)
	bool k_column_read = false;
	std::atomic<bool> fact_checked{ false };
	bool fact_ok = false;
	double fact_dk = 0.0;
	int fact_nbins = 0;
	std::vector<int> fact_n;            // k[i] = fact_n[i] * fact_dk
	std::mutex fact_mutex;
	void factored_init();

	// Per-time-level coefficient sets. Callers (e.g. the spline lsgrid)
	// evaluate SEVERAL time levels concurrently, so coefficients are kept
	// in a small cache of immutable entries handed out as shared pointers:
	// an entry is fully built under fact_mutex before it becomes visible
	// and is never mutated afterwards (except adding the independent uvw
	// part under the same mutex).
	struct FactCoeffs {
		double t = 0.0;
		bool has_eta2 = false;
		bool has_uvw = false;
		std::vector<double> ramp;                    // per-component ramp
		std::vector<std::complex<double>> eta2_bins;
		std::vector<std::complex<double>> Qk;        // u,v: k-weighted
		std::vector<std::complex<double>> Qw;        // w: sin-type
		std::vector<std::complex<double>> Qp;        // p: omega-weighted
	};
	std::vector<std::shared_ptr<FactCoeffs>> fact_cache;
	std::shared_ptr<FactCoeffs> fact_get(double t, bool need_uvw);
	void fact_build_eta2(FactCoeffs& c, double t);
	void fact_build_uvw(FactCoeffs& c, double t);
	double eta2_factored_eval(const FactCoeffs& c, double xx, double yy);
	UVWP uvw2_factored_eval(const FactCoeffs& c, double xx, double yy,
	                        double zz);

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
	double eta2_full(double t, double xx, double yy); // no kmax truncation (smit path)
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

	// Sloping seabed (variable bathymetry). When seabed.active(), the
	// kinematics dispatch to the slope-corrected first-order evaluation
	// below: per point, the local depth h(x,y) and bed slope grad h come
	// from the Seabed probe, each component's wavenumber is re-solved
	// from the local dispersion relation, and the Airy vertical
	// eigenfunction is replaced by (see seabed_development/):
	//   F_j(z) = [cosh k(z+h) - i b sinh k(z+h)]
	//            / [cosh kh - i b sinh kh],   b = khat_j . grad h,
	// which satisfies the sloping-bottom kinematic condition
	// w + grad h . u_H = 0 at z = -h exactly (pointwise identity).
	Seabed seabed;
	// local wavenumber from om^2 = g k tanh(k h) (Newton)
	double dispersion_k(double om, double h) const;
	// slope-corrected first-order kinematics; zz measured from swl
	UVWP uvwp1_slope(double t, double xx, double yy, double zz);
	// first-order elevation with LOCAL wavenumbers in the phase (the
	// dynamic surface BC of the slope-corrected potential; F(0) = 1
	// exactly, so no slope term remains at z = 0). Must be used instead
	// of eta1() whenever the seabed is active, so the surface and the
	// kinematics/eta2_slope share the same local phase function.
	double eta1_slope(double t, double xx, double yy);

	// Depth-local Sharma & Dean second order (variable bathymetry):
	// the same interaction coefficients as rebuild_pair_cache()/uvw2(),
	// but built for the LOCAL depth h(x, y) with local wavenumbers.
	// The interaction limit (Johannessen, APOR_736: sum interactions
	// with the high-k spectral tail misbehave for wide-banded spectra,
	// e.g. Torsethaugen's w^-4 decay) is applied to CROSS sum
	// interactions: k_i + k_m <= kmax (max local component k). The
	// DIAGONAL self-interaction (each component's own Stokes second
	// harmonic - bounded and physical) is kept in both the elevation
	// and the kinematics, unlike the constant-depth S&D path which
	// drops it from eta2 only; this keeps the exported surface and the
	// kinematics mutually consistent. A per-thread coefficient set is
	// cached on h, so evaluations down a water column (same x, y)
	// rebuild nothing.
	struct LocalPairSet {
		double h = -1.0;               // depth this set was built for
		const void* owner = nullptr;   // Irregular instance it belongs to
		std::vector<double> kl, Rl, Rsl, c2kd, cth, sth;
		std::vector<double> diagD;     // diagonal D_nn+
		struct P {
			int i, m;
			bool has_minus;
			bool kmax_p;               // ki + km <= local kmax
			double k_p, k_mn, D_p, D_mn, a_p, a_mn;
		};
		std::vector<P> pairs;
	};
	void build_local_pair_set(double h, LocalPairSet& s) const;
	double eta2_slope(double t, double xx, double yy);
	UVWP uvwp2_slope(double t, double xx, double yy, double zz);

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
	// First-order-only dynamic pressure with the smit stretched profile:
	// the dp1 part of dp_smit(), O(nfreq) per call (no O(nfreq^2) pair
	// sums). Used to serve p_nh in lsgrid mode (wavetype 5), where the
	// spline grid stores no pressure field and the full second-order
	// evaluation is too expensive per boundary query.
	double dp1_smit(double t, double x, double y, double z);
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
