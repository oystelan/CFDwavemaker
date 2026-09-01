#ifndef SIGMA_FFT_H
#define SIGMA_FFT_H

#include <complex>
#include <vector>

class Irregular;

/* Layer-by-layer FFT fill of the smit s-coordinate second-order
   kinematics on an equal-dk wave-vector lattice.

   Supported component grids (waveinput "equal_dk 1" +
   "vertical_theory smit", swl = 0):
     - 1D: long-crested, propagation along +x, k_n = n*dkx
     - 2D: short-crested lattice (kx, ky) = (i*dkx, j*dky), i > 0,
       j signed (directions folded into the lattice; includes mtheta)

   The FFT evaluates one full spatial period Lx = 2*pi/dkx (and
   Ly = 2*pi/dky in 2D) at the node spacing Lx/nfftx, anchored at the
   lsgrid origin. The lsgrid stores a WINDOW of those nodes: its grid
   spacing must equal the tile node spacing (validated in setup), and
   its node count may be anything up to one full period + 1 (the wrap
   point). A full-period grid is exactly periodic (queries outside
   wrap); a smaller window behaves like a regular bounded lsgrid.

   Why it is exact (see sigma_developments/fft_prototype/README.md and
   the sigma note): on a sigma layer every vertical factor is a layer
   constant (the bin magnitude |k12| = hypot(jx*dkx, jy*dky) depends
   only on the bin); bound-wave pair sums collapse into lattice bins;
   the particular solution is the pointwise product
   phi_P = g(1+s)*eta1*G; chain and Jacobian terms are pointwise
   products of IFFT'd fields. One O(N^2) pair pass per time level plus
   a handful of (2D) FFTs per layer. */
class SigmaFFTFill {
public:
	bool enabled = false;
	bool twodim = false;
	int nfftx = 0, nffty = 1;   // grid intervals per period (pow2; 1 in 1D)
	double Lx = 0.0, Ly = 0.0;
	double dkx = 0.0, dky = 0.0;

	// Dry-run lattice detection (for the domain lock, before allocation).
	// Returns false if the input is not FFT-fillable. needx/needy are the
	// minimum FFT sizes (2*(2*jmax)+2 rounded to pow2 happens in caller).
	static bool probe(Irregular& irr, double swl, bool& twodim_out,
	                  double& dkx_out, double& dky_out,
	                  int& jxmax_out, int& jymax_out);

	// Check preconditions against the actual grid and build tables.
	// (dxg, dyg): lsgrid node spacings; each must equal L/nfft for a
	// power-of-two nfft large enough to hold the signed second-order
	// bins (2*(2*jmax)+2), otherwise setup fails and the caller falls
	// back to the pointwise fill. dyg is ignored for 1D lattices.
	// (x0, y0): physical coordinates of grid node (0,0) -- absolute;
	// fpoint offsets are handled internally.
	// svals: sigma value of each layer.
	// nxg/nyg: node count of the requested lsgrid window, used only to
	// report how much of the FFT lattice it covers (0 = whole period).
	bool setup(Irregular& irr, double dxg, double dyg,
	           double x0, double y0,
	           const std::vector<double>& svals, double swl,
	           int nxg = 0, int nyg = 0);

	// Fill one time level. All arrays are sized npts = nfftx*nffty with
	// row-major layout [ix*nffty + iy]; u/v/w are [layer][npts].
	// p_out (optional): the second-order DYNAMIC pressure (rho included,
	// Bernoulli term included) on the same [layer][npts] layout -- the
	// FFT-field form of Irregular::dp_smit(). Costs ~1/3 more transforms
	// per layer, so it is computed only when requested.
	void fill(Irregular& irr, double t,
	          std::vector<double>& eta_out,
	          std::vector<double>& u_out,
	          std::vector<double>& v_out,
	          std::vector<double>& w_out,
	          std::vector<double>* p_out = nullptr);

private:
	struct PairT {
		int i, m;
		int jpx, jpy, jmx, jmy;     // sum / difference lattice bins
		double alpha_p, alpha_m;    // S&D elevation coefficients
		double H_p, H_m;            // smit potential coefficients
		double om_p, om_m;
	};
	std::vector<PairT> pairs;
	std::vector<double> diag_alpha, diag_H;
	std::vector<int> bx, by;            // component lattice indices
	int jxmax2 = 0, jymax2 = 0;         // max |bin| at second order
	double x0_ = 0.0, y0_ = 0.0;
	double depth_ = 0.0;
	std::vector<double> svals_;
	// per-layer vertical factors: [layer][component] and [layer][bin grid]
	std::vector<std::vector<double>> Ch1, Sh1, ChB, ShB;

	int npts() const { return nfftx * nffty; }
	int bin_index(int jx, int jy) const;  // signed bins -> FFT index
	void ifft2(std::vector<std::complex<double>>& a) const;
	void eval1(const std::vector<std::complex<double>>& c,
	           std::vector<double>& out,
	           std::vector<std::complex<double>>& work) const;
	void evalB(const std::vector<std::complex<double>>& cbins,
	           std::vector<double>& out,
	           std::vector<std::complex<double>>& work) const;
};

#endif
