#ifndef LAYER_STRETCH_H
#define LAYER_STRETCH_H

/* ---------------------------------------------------------------------------
   Vertical sigma-layer stretching for the lsgrid, and the automatic choice of
   its parameter.

   The mapping is exponential: the computational coordinate u in [0,1] (0 =
   seabed, 1 = surface) maps to sigma in [-1,0] via

       sigma(u) = -(exp(beta*(1-u)) - 1) / (exp(beta) - 1).

   Adjacent layer thicknesses form a geometric progression with ratio
   exp(beta/(nl-1)), so the bottom-to-surface thickness ratio actually realised
   by nl layers is exp(beta*(nl-2)/(nl-1)). exp(beta) is only the nl -> inf
   limit of that, and overstates the grading noticeably at the small nl this
   grid targets (nl=4, beta=3 gives 7.4, not 20.1). beta = 0 is uniform.

   optimal_beta() chooses beta by minimising the error that the vertical cubic
   Hermite reconstruction (lsGridSpline::spline_interp_velo) makes on the sea
   state's representative velocity profile -- the spectral RMS profile

       sigma_u^2(z) = 1/2 sum_j V_j(z)^2,
       V_j(z)       = A_j g k_j / omega_j * cosh(k_j(z+h)) / cosh(k_j h),

   using the same per-component vertical factor as the direct evaluation
   (Irregular::u1 / profileX). The overall wave amplitude cancels out of the
   argmin, so the answer depends only on the spectral shape; the sum over
   components handles bimodal seas (swell and wind sea decay very differently
   with depth) without any effective-wavenumber fit.

   Why a direct search rather than a formula: a prototype study
   (sigma_developments/layer_optimization) compared this against closed forms
   beta* = ke*h/4, against classical 4th-derivative knot equidistribution, and
   against richer knot-density families. The closed forms were unreliable (a
   factor ~6 spread over the sea states tested); the richer families overfitted
   the RMS proxy, driving its error down without improving -- sometimes while
   worsening -- the error on actual crest profiles. Minimising directly within
   this one-parameter family was the most robust option and costs microseconds.

   Header-only and free of any library dependency, so it can be exercised
   standalone (see sigma_developments/layer_optimization/verify_cpp.cpp).
   --------------------------------------------------------------------------- */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace layer_stretch {

// Upper bound of the automatic search; a grading of exp(12) is already far
// beyond anything a real spectrum calls for.
const double BETA_MAX = 12.0;

inline double u2sigma(double u, double beta) {
	if (std::abs(beta) < 1.0e-9) return u - 1.0;
	return -std::expm1(beta * (1.0 - u)) / std::expm1(beta);
}

inline double sigma2u(double sigma, double beta) {
	if (std::abs(beta) < 1.0e-9) return sigma + 1.0;
	return 1.0 - std::log1p(-sigma * std::expm1(beta)) / beta;
}

inline double dsigma_du(double u, double beta) {
	if (std::abs(beta) < 1.0e-9) return 1.0;
	return beta * std::exp(beta * (1.0 - u)) / std::expm1(beta);
}

// Bottom-to-surface layer thickness ratio realised by nl layers.
inline double thickness_ratio(double beta, int nl) {
	if (nl < 3) return 1.0;
	return std::exp(beta * double(nl - 2) / double(nl - 1));
}

/* cosh(k*x)/cosh(k*h) and sinh(k*x)/cosh(k*h) for 0 <= x <= h, written so that
   neither numerator nor denominator can overflow at large k*h. */
inline void profile_ratios(double k, double x, double h,
	double& c_ratio, double& s_ratio) {
	const double e = std::exp(k * (x - h));      // <= 1
	const double a = std::exp(-2.0 * k * x);     // <= 1
	const double b = std::exp(-2.0 * k * h);     // <= 1
	c_ratio = e * (1.0 + a) / (1.0 + b);
	s_ratio = e * (1.0 - a) / (1.0 + b);
}

/* sigma_u(z) and d sigma_u/dz for a component set. c[j] is the squared
   surface velocity amplitude (A_j g k_j / omega_j)^2 of component j. */
inline void rms_profile(const std::vector<double>& c, const std::vector<double>& k,
	double h, double z, double& s, double& dsdz) {
	double p = 0.0, dp = 0.0;
	for (std::size_t j = 0; j < k.size(); j++) {
		double cr, sr;
		profile_ratios(k[j], z + h, h, cr, sr);
		p += c[j] * cr * cr;
		dp += c[j] * 2.0 * k[j] * cr * sr;
	}
	s = std::sqrt(0.5 * p);
	dsdz = (s > 0.0) ? 0.25 * dp / s : 0.0;
}

/* Relative L2 error of the cubic Hermite reconstruction of the profile for one
   candidate beta. zq/fq hold the reference profile on a fine column and den is
   sum(fq^2); they do not depend on beta, so the caller builds them once. */
inline double hermite_cost(double beta, int nl,
	const std::vector<double>& c, const std::vector<double>& k, double h,
	const std::vector<double>& zq, const std::vector<double>& fq, double den) {

	// nodal values and nodal d/du, exactly as the grid stores them
	std::vector<double> f(nl), m(nl);
	const double du = 1.0 / double(nl - 1);
	for (int i = 0; i < nl; i++) {
		const double u = du * i;
		const double z = u2sigma(u, beta) * h;          // at eta = 0
		double s, dsdz;
		rms_profile(c, k, h, z, s, dsdz);
		f[i] = s;
		m[i] = dsdz * dsigma_du(u, beta) * h;           // chain rule to d/du
	}

	double num = 0.0;
	for (std::size_t q = 0; q < zq.size(); q++) {
		double u = sigma2u(zq[q] / h, beta);
		u = std::min(std::max(u, 0.0), 1.0);
		int i = int(u / du);
		if (i > nl - 2) i = nl - 2;
		if (i < 0) i = 0;
		const double t = (u - i * du) / du;
		// same Hermite form as spline_interp_velo
		const double d = f[i + 1] - f[i];
		const double a = du * m[i] - d;
		const double b = -du * m[i + 1] + d;
		const double fi = f[i] * (1.0 - t) + f[i + 1] * t
			+ t * (1.0 - t) * (a * (1.0 - t) + b * t);
		const double e = fi - fq[q];
		num += e * e;
	}
	return (den > 0.0) ? std::sqrt(num / den) : 0.0;
}

/* Beta minimising hermite_cost for the given component set and layer count.
   Returns a negative value if the inputs cannot define a profile, in which
   case the caller should keep whatever beta it already had. If err is
   non-null it receives the relative L2 error at the optimum. */
inline double optimal_beta(int nl, const std::vector<double>& c,
	const std::vector<double>& k, double h, double* err = 0, int nq = 257) {

	if (nl < 3 || h <= 0.0 || k.empty() || c.size() != k.size()) return -1.0;

	std::vector<double> zq(nq), fq(nq);
	double den = 0.0;
	for (int q = 0; q < nq; q++) {
		zq[q] = -h + h * double(q) / double(nq - 1);
		double dsdz;
		rms_profile(c, k, h, zq[q], fq[q], dsdz);
		den += fq[q] * fq[q];
	}
	if (!(den > 0.0)) return -1.0;

	// coarse scan, then successive bracket refinement. The cost is smooth in
	// beta but not always unimodal -- a bimodal spectrum can put two comparable
	// basins in [0, BETA_MAX] -- so the first pass is a full sweep at spacing
	// BETA_MAX/120 and only the refinement is local.
	double lo = 0.0, hi = BETA_MAX, best = 0.0, best_cost = 1e300;
	int n = 121;
	for (int pass = 0; pass < 4; pass++) {
		const double blo = lo, bhi = hi;
		int ibest = 0;
		best_cost = 1e300;
		for (int i = 0; i < n; i++) {
			const double beta = blo + (bhi - blo) * double(i) / double(n - 1);
			const double cost = hermite_cost(beta, nl, c, k, h, zq, fq, den);
			if (cost < best_cost) { best_cost = cost; best = beta; ibest = i; }
		}
		const double step = (bhi - blo) / double(n - 1);
		lo = std::max(blo, blo + step * (ibest - 1));
		hi = std::min(bhi, blo + step * (ibest + 1));
		n = 21;
	}
	if (err) *err = best_cost;
	return best;
}

/* Relative L2 error of the nl-layer Hermite reconstruction at a FIXED beta --
   the companion to optimal_beta() for callers that keep a prescribed grading
   (used by the automatic layer-count search). Negative if the inputs cannot
   define a profile. */
inline double profile_error(double beta, int nl, const std::vector<double>& c,
	const std::vector<double>& k, double h, int nq = 257) {

	if (nl < 3 || h <= 0.0 || k.empty() || c.size() != k.size()) return -1.0;

	std::vector<double> zq(nq), fq(nq);
	double den = 0.0;
	for (int q = 0; q < nq; q++) {
		zq[q] = -h + h * double(q) / double(nq - 1);
		double dsdz;
		rms_profile(c, k, h, zq[q], fq[q], dsdz);
		den += fq[q] * fq[q];
	}
	if (!(den > 0.0)) return -1.0;
	return hermite_cost(beta, nl, c, k, h, zq, fq, den);
}

}  // namespace layer_stretch

#endif
