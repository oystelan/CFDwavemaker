#include "sigma_fft.h"
#include "mpi_utils.h"
#include "Irregular.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

// In-place iterative radix-2 complex DFT with e^{+i 2 pi j m / N} kernel
// (unnormalized inverse transform). Self-contained: no external FFT
// dependency (SWD keeps its optional FFTW linkage independently).
void fft_pow2_inverse(std::complex<double>* a, size_t n)
{
	if (n < 2) return;
	for (size_t i = 1, j = 0; i < n; ++i) {
		size_t bit = n >> 1;
		for (; j & bit; bit >>= 1) j ^= bit;
		j ^= bit;
		if (i < j) std::swap(a[i], a[j]);
	}
	for (size_t len = 2; len <= n; len <<= 1) {
		const double ang = 2.0 * PI / double(len);
		const std::complex<double> wlen(std::cos(ang), std::sin(ang));
		for (size_t i = 0; i < n; i += len) {
			std::complex<double> w(1.0, 0.0);
			for (size_t k = 0; k < len / 2; ++k) {
				const std::complex<double> u = a[i + k];
				const std::complex<double> v = a[i + k + len / 2] * w;
				a[i + k] = u + v;
				a[i + k + len / 2] = u - v;
				w *= wlen;
			}
		}
	}
}

// Fit an equal-spacing lattice v_i = n_i * dk (n_i integer, possibly
// signed). Returns false if no consistent fit exists.
bool fit_lattice(const std::vector<double>& v, double& dk,
                 std::vector<int>& n)
{
	const int N = (int)v.size();
	n.assign(N, 0);
	double vmax = 0.0;
	for (double x : v) vmax = std::max(vmax, std::abs(x));
	if (vmax < 1e-14) { dk = 0.0; return true; }   // all zero

	// guess: smallest |v| clearly above the file-rounding noise floor
	// (input files carry rounded k/theta, so a nominally-zero component
	// of the wave vector can come out as ~1e-9 - it must be treated as
	// zero, NOT as the lattice spacing), and smallest positive gap
	// between sorted distinct |v|
	const double zero_tol = 1e-6 * vmax;
	std::vector<double> mags;
	for (double x : v) if (std::abs(x) > zero_tol) mags.push_back(std::abs(x));
	if (mags.empty()) { dk = 0.0; return true; }
	std::sort(mags.begin(), mags.end());
	double guess = mags[0];
	for (size_t i = 1; i < mags.size(); ++i) {
		const double gap = mags[i] - mags[i - 1];
		if (gap > zero_tol && gap < guess) guess = gap;
	}

	dk = guess;
	for (int pass = 0; pass < 3; ++pass) {
		double skn = 0.0, snn = 0.0;
		for (int i = 0; i < N; ++i) {
			n[i] = (std::abs(v[i]) <= zero_tol)
				? 0 : (int)std::lround(v[i] / dk);
			skn += v[i] * n[i];
			snn += double(n[i]) * n[i];
		}
		if (snn <= 0.0) return false;
		dk = skn / snn;
	}
	if (dk < 1e-4 * vmax) return false;   // degenerate spacing
	for (int i = 0; i < N; ++i) {
		n[i] = (std::abs(v[i]) <= zero_tol)
			? 0 : (int)std::lround(v[i] / dk);
		if (std::abs(v[i] - n[i] * dk) > 1e-5 * vmax) return false;
	}
	return true;
}

} // namespace

int SigmaFFTFill::bin_index(int jx, int jy) const
{
	const int ix = ((jx % nfftx) + nfftx) % nfftx;
	const int iy = ((jy % nffty) + nffty) % nffty;
	return ix * nffty + iy;
}

void SigmaFFTFill::ifft2(std::vector<std::complex<double>>& a) const
{
	/* The row/column 1-D transforms are independent -> OpenMP-parallel.
	   When ifft2 is called from inside an already-parallel region (the
	   per-layer loop in fill()), nested parallelism is disabled by
	   default and these become one-thread teams -- correct, just inert. */
	// rows (y-direction, contiguous)
	if (nffty > 1) {
#pragma omp parallel for
		for (int ix = 0; ix < nfftx; ++ix)
			fft_pow2_inverse(&a[(size_t)ix * nffty], nffty);
	}
	// columns (x-direction, stride nffty)
#pragma omp parallel
	{
		std::vector<std::complex<double>> col(nfftx);
#pragma omp for
		for (int iy = 0; iy < nffty; ++iy) {
			for (int ix = 0; ix < nfftx; ++ix) col[ix] = a[(size_t)ix * nffty + iy];
			fft_pow2_inverse(col.data(), nfftx);
			for (int ix = 0; ix < nfftx; ++ix) a[(size_t)ix * nffty + iy] = col[ix];
		}
	}
}

void SigmaFFTFill::eval1(const std::vector<std::complex<double>>& c,
                         std::vector<double>& out,
                         std::vector<std::complex<double>>& work) const
{
	std::fill(work.begin(), work.end(), std::complex<double>(0.0, 0.0));
	for (size_t nn = 0; nn < c.size(); ++nn)
		work[bin_index(bx[nn], by[nn])] += c[nn];
	ifft2(work);
	out.resize(npts());
	for (int m = 0; m < npts(); ++m) out[m] = work[m].real();
}

void SigmaFFTFill::evalB(const std::vector<std::complex<double>>& cbins,
                         std::vector<double>& out,
                         std::vector<std::complex<double>>& work) const
{
	work.assign(cbins.begin(), cbins.end());
	ifft2(work);
	out.resize(npts());
	for (int m = 0; m < npts(); ++m) out[m] = work[m].real();
}

bool SigmaFFTFill::probe(Irregular& irr, double swl, bool& twodim_out,
                         double& dkx_out, double& dky_out,
                         int& jxmax_out, int& jymax_out)
{
	if (!irr.equal_dk || irr.vertical_theory != 1) return false;
	if (std::abs(swl) > 1e-12) return false;

	const int N = irr.nfreq * irr.ndir;
	if ((int)irr.cth_cache.size() != N) irr.rebuild_second_order_cache();

	std::vector<double> kx(N), ky(N);
	for (int i = 0; i < N; ++i) {
		kx[i] = irr.k[i] * irr.cth_cache[i];
		ky[i] = irr.k[i] * irr.sth_cache[i];
	}

	double dkx, dky;
	std::vector<int> nx_, ny_;
	if (!fit_lattice(kx, dkx, nx_)) return false;
	if (!fit_lattice(ky, dky, ny_)) return false;
	if (dkx == 0.0) return false;   // waves purely along y not supported

	twodim_out = (dky != 0.0);
	dkx_out = dkx;
	dky_out = dky;
	jxmax_out = 0;
	jymax_out = 0;
	for (int i = 0; i < N; ++i) {
		jxmax_out = std::max(jxmax_out, std::abs(nx_[i]));
		jymax_out = std::max(jymax_out, std::abs(ny_[i]));
	}
	// pathological lattices (and int overflow in size computations)
	if (jxmax_out > 2048 || jymax_out > 2048) return false;
	return true;
}

bool SigmaFFTFill::setup(Irregular& irr, double dxg, double dyg,
                         double x0, double y0,
                         const std::vector<double>& svals, double swl,
                         int nxg, int nyg)
{
	enabled = false;
	bool twod;
	int jxmax, jymax;
	if (!probe(irr, swl, twod, dkx, dky, jxmax, jymax)) return false;

	const int N = irr.nfreq * irr.ndir;
	std::vector<double> kx(N), ky(N);
	for (int i = 0; i < N; ++i) {
		kx[i] = irr.k[i] * irr.cth_cache[i];
		ky[i] = irr.k[i] * irr.sth_cache[i];
	}
	fit_lattice(kx, dkx, bx);
	fit_lattice(ky, dky, by);

	twodim = twod;
	jxmax2 = 2 * jxmax;
	jymax2 = 2 * jymax;
	Lx = 2.0 * PI / dkx;
	Ly = twodim ? 2.0 * PI / dky : 0.0;

	// FFT sizes from the grid spacing: nfft = L/dx must be a power of
	// two large enough that the signed second-order bins are distinct
	// frequencies. The grid then samples FFT nodes exactly (any window
	// of them; the input reader snaps the lsgrid accordingly).
	if (dxg <= 0.0) return false;
	nfftx = (int)std::lround(Lx / dxg);
	if (std::fabs(Lx / dxg - nfftx) > 1e-6 * (Lx / dxg)) return false;
	if (nfftx < 2 * jxmax2 + 2) return false;
	if (nfftx & (nfftx - 1)) return false;
	if (twodim) {
		if (dyg <= 0.0) return false;
		nffty = (int)std::lround(Ly / dyg);
		if (std::fabs(Ly / dyg - nffty) > 1e-6 * (Ly / dyg)) return false;
		if (nffty < 2 * jymax2 + 2) return false;
		if (nffty & (nffty - 1)) return false;
	} else {
		nffty = 1;
	}

	x0_ = x0 - irr.fpoint[0];
	y0_ = y0 - irr.fpoint[1];
	depth_ = irr.depth;
	svals_ = svals;

	// pair tables (alpha for the untruncated eta2, H for the bound
	// potential; pair_cache carries them for the full pair set)
	pairs.clear();
	for (const auto& p : irr.pair_cache) {
		PairT q{};
		q.i = p.i; q.m = p.m;
		q.jpx = bx[p.i] + bx[p.m];
		q.jpy = by[p.i] + by[p.m];
		q.jmx = bx[p.m] - bx[p.i];
		q.jmy = by[p.m] - by[p.i];
		q.alpha_p = p.has_plus ? p.alpha_nm_plus : 0.0;
		q.H_p = p.has_plus ? p.H_nm_plus : 0.0;
		q.alpha_m = p.has_minus ? p.alpha_nm_minus : 0.0;
		q.H_m = p.has_minus ? p.H_nm_minus : 0.0;
		q.om_p = irr.omega[p.i] + irr.omega[p.m];
		q.om_m = irr.omega[p.i] - irr.omega[p.m];
		pairs.push_back(q);
	}
	diag_alpha.resize(N);
	diag_H.resize(N);
	for (int i = 0; i < N; ++i) {
		const double ki = irr.k[i], oi = irr.omega[i];
		const double Rni = irr.R_cache[i];
		const double D_nn = (6.0 * Rni * (ki * ki - Rni * Rni))
			/ (2.0 * Rni - ki * std::tanh(2.0 * ki * irr.depth));
		diag_alpha[i] = (D_nn - (ki * ki - Rni * Rni)) / Rni + 2.0 * Rni;
		diag_H[i] = -irr.g / (2.0 * oi) * (diag_alpha[i] / 2.0
			- (3.0 * oi * oi) / (2.0 * irr.g)
			+ irr.g * ki * ki / (2.0 * oi * oi));
	}

	// per-layer vertical factors (time-independent). ChB/ShB are full
	// FFT-index grids, populated only where second-order bins can live.
	const int nl = (int)svals_.size();
	Ch1.assign(nl, std::vector<double>(N));
	Sh1.assign(nl, std::vector<double>(N));
	ChB.assign(nl, std::vector<double>(npts(), 0.0));
	ShB.assign(nl, std::vector<double>(npts(), 0.0));
	for (int m = 0; m < nl; ++m) {
		const double arg = depth_ * (1.0 + svals_[m]);   // layer constant
		for (int i = 0; i < N; ++i) {
			const double ki = irr.k[i];
			Ch1[m][i] = std::cosh(ki * arg) / std::cosh(ki * depth_);
			Sh1[m][i] = std::sinh(ki * arg) / std::cosh(ki * depth_);
		}
		for (int ix = 0; ix < nfftx; ++ix) {
			const int jx = (ix <= nfftx / 2) ? ix : ix - nfftx;
			if (std::abs(jx) > jxmax2) continue;
			for (int iy = 0; iy < nffty; ++iy) {
				const int jy = (iy <= nffty / 2) ? iy : iy - nffty;
				if (std::abs(jy) > jymax2) continue;
				const double kj = std::hypot(jx * dkx, jy * dky);
				// stable ratio form (avoids cosh overflow at large kj*d)
				const double e1 = std::exp(kj * (arg - depth_));
				const double e2 = std::exp(-kj * (arg + depth_));
				const double e3 = std::exp(-2.0 * kj * depth_);
				ChB[m][ix * nffty + iy] = (e1 + e2) / (1.0 + e3);
				ShB[m][ix * nffty + iy] = (e1 - e2) / (1.0 + e3);
			}
		}
	}

	enabled = true;
	// Nodes of the periodic FFT lattice that the requested lsgrid window
	// actually samples. The window is snapped to the lattice, and node
	// nx-1 wraps onto node 0, so a window wider than one period re-visits
	// nodes rather than adding new ones -- worth reporting, because a
	// window that covers only a fraction of the period is the usual sign
	// of bounds that do not match the CFD domain.
	const int nwx = (nxg > 0) ? std::min(nxg, nfftx) : nfftx;
	const int nwy = twodim ? ((nyg > 0) ? std::min(nyg, nffty) : nffty) : 1;
	cfdwm_mpi::out()
	    << "sigma-fft: FFT layer fill enabled ("
	    << (twodim ? "2D short-crested lattice, " : "1D long-crested, ")
	    << nfftx << "x" << nffty << " grid, Lx=" << Lx
	    << (twodim ? " m, Ly=" : " m")
	    << (twodim ? std::to_string(Ly) + " m" : std::string())
	    << ", " << nl << " layers)." << std::endl;
	// Does the node spacing actually resolve the spectrum? The lattice size
	// is set by the Nyquist requirement for the second-order sum-frequency
	// bins, which is NOT the same as resolving them for interpolation: the
	// stored values are exact collocation values, but the cubic Hermite
	// between nodes needs several points per wavelength. Report it rather
	// than silently choosing a finer (and much more expensive) grid.
	{
		double kmax = 0.0;
		for (int i = 0; i < N; ++i) kmax = std::max(kmax, irr.k[i]);
		if (kmax > 0.0) {
			const double dxn = Lx / nfftx;
			const double L1 = 2.0 * PI / kmax;          // shortest 1st order
			cfdwm_mpi::out()
			    << "sigma-fft: node spacing " << dxn << " m resolves the "
			       "shortest first-order wave (" << L1 << " m) with "
			    << L1 / dxn << " points, the shortest second-order sum ("
			    << 0.5 * L1 << " m) with " << 0.5 * L1 / dxn << ". Below "
			       "about 6 points the Hermite interpolation between nodes "
			       "starts to smear that component; raise lsgrid nx/ny to "
			       "refine (cost grows as nx*ny)." << std::endl;
		}
	}
	cfdwm_mpi::out()
	    << "sigma-fft: " << (long)nwx * nwy << " of " << (long)nfftx * nffty
	    << " FFT nodes lie inside the lsgrid bounds ["
	    << x0 << ", " << x0 + (nxg > 1 ? (nxg - 1) * dxg : 0.) << "] x ["
	    << y0 << ", " << y0 + (twodim && nyg > 1 ? (nyg - 1) * dyg : 0.)
	    << "] (" << nwx << "x" << nwy << " of " << nfftx << "x" << nffty
	    << "). Depth is the constant " << depth_
	    << " m: the FFT path never reads the [seabed] bathymetry."
	    << std::endl;
	return true;
}

void SigmaFFTFill::fill(Irregular& irr, double t,
                        std::vector<double>& eta_out,
                        std::vector<double>& u_out,
                        std::vector<double>& v_out,
                        std::vector<double>& w_out,
                        std::vector<double>* p_out)
{
	const int N = irr.nfreq * irr.ndir;
	const int nl = (int)svals_.size();
	const int NP = npts();
	const double g = irr.g;
	const bool want_p = (p_out != nullptr);
	const std::complex<double> I(0.0, 1.0);

	// component phasors at grid node (0,0); er includes the ramp
	std::vector<std::complex<double>> b(N), er(N);
	std::vector<double> kxv(N), kyv(N);
	for (int n = 0; n < N; ++n) {
		kxv[n] = irr.k[n] * irr.cth_cache[n];
		kyv[n] = irr.k[n] * irr.sth_cache[n];
		const double ph = kxv[n] * x0_ + kyv[n] * y0_
			- irr.omega[n] * t + irr.phase0_cache[n];
		b[n] = irr.A[n] * std::exp(I * ph);
		const double fr = irr.freq_ramp(t, irr.freq_ramptime[n],
			irr.freq_starttime[n], irr.freq_duration[n]);
		er[n] = b[n] * fr;
	}

	// second-order bin coefficients: one O(N^2) pass per time level.
	// OpenMP: each thread accumulates into private bin arrays, merged
	// under a critical section (bin write collisions between pairs make
	// direct sharing impossible). MPI: each rank sums its block of the
	// pair list (rank 0 also the diagonal), Allreduce combines -- the
	// summation order differs from serial, so results match to roundoff
	// (~1e-16 relative), not bit-exactly.
	// Cpt (pressure only): the time derivative of the BOUND second-order
	// potential at fixed s, i.e. the +H*om12*Ch12*cos(phi12) rows of
	// dp_smit() -- same bins as Cub/Cvb/Cws, frequency factor instead of
	// wavenumber factor.
	std::vector<std::complex<double>> Ceta2(NP), Cub(NP), Cvb(NP), Cws(NP);
	std::vector<std::complex<double>> Cpt(want_p ? NP : 0);
	long p0 = 0, p1 = (long)pairs.size();
	bool add_diag = true;
#if defined(MPI_enable)
	if (cfdwm_mpi::active() && cfdwm_mpi::size() > 1) {
		int a0, a1;
		cfdwm_mpi::range((int)pairs.size(), a0, a1);
		p0 = a0; p1 = a1;
		add_diag = (cfdwm_mpi::rank() == 0);
	}
#endif
#pragma omp parallel
	{
		std::vector<std::complex<double>> lE(NP), lU(NP), lV(NP), lW(NP);
		std::vector<std::complex<double>> lP(want_p ? NP : 0);
#pragma omp for nowait
		for (long qq = p0; qq < p1; ++qq) {
			const PairT& q = pairs[qq];
			const std::complex<double> Ekin = er[q.i] * er[q.m];
			const std::complex<double> Eeta = er[q.i] * b[q.m]; // ramp[i] only
			const std::complex<double> EkinM = er[q.m] * std::conj(er[q.i]);
			const std::complex<double> EetaM = b[q.m] * std::conj(er[q.i]);
			if (q.alpha_p != 0.0 || q.H_p != 0.0) {
				const int bi = bin_index(q.jpx, q.jpy);
				const double kabs = std::hypot(q.jpx * dkx, q.jpy * dky);
				lE[bi] += 0.5 * q.alpha_p * Eeta;
				lU[bi] += -q.H_p * (q.jpx * dkx) * Ekin;
				lV[bi] += -q.H_p * (q.jpy * dky) * Ekin;
				lW[bi] += -q.H_p * (kabs * depth_) * (-I) * Ekin;
				if (want_p) lP[bi] += q.H_p * q.om_p * Ekin;
			}
			if (q.alpha_m != 0.0 || q.H_m != 0.0) {
				const int bi = bin_index(q.jmx, q.jmy);
				const double kabs = std::hypot(q.jmx * dkx, q.jmy * dky);
				lE[bi] += 0.5 * q.alpha_m * EetaM;
				// kx12 = kx_i - kx_m = -(jmx*dkx), likewise for y
				lU[bi] += q.H_m * (q.jmx * dkx) * EkinM;
				lV[bi] += q.H_m * (q.jmy * dky) * EkinM;
				lW[bi] += -q.H_m * (kabs * depth_) * (+I) * EkinM;
				// cos is even: the (phi_m - phi_i) phasor carries
				// cos(phi_i - phi_m) with om_m = omega_i - omega_m
				if (want_p) lP[bi] += q.H_m * q.om_m * EkinM;
			}
		}
#pragma omp critical (sigma_fft_bin_merge)
		{
			for (int m = 0; m < NP; ++m) {
				Ceta2[m] += lE[m];
				Cub[m]  += lU[m];
				Cvb[m]  += lV[m];
				Cws[m]  += lW[m];
			}
			if (want_p)
				for (int m = 0; m < NP; ++m) Cpt[m] += lP[m];
		}
	}
	if (add_diag)
		for (int i = 0; i < N; ++i) {
			const int bi = bin_index(2 * bx[i], 2 * by[i]);
			const std::complex<double> E2 = er[i] * er[i];
			Ceta2[bi] += 0.25 * diag_alpha[i] * (er[i] * b[i]);
			Cub[bi] += -0.5 * diag_H[i] * (2.0 * kxv[i]) * E2;
			Cvb[bi] += -0.5 * diag_H[i] * (2.0 * kyv[i]) * E2;
			Cws[bi] += -0.5 * diag_H[i] * (2.0 * irr.k[i] * depth_) * (-I) * E2;
			if (want_p) Cpt[bi] += 0.5 * diag_H[i] * (2.0 * irr.omega[i]) * E2;
		}
#if defined(MPI_enable)
	if (cfdwm_mpi::active() && cfdwm_mpi::size() > 1) {
		MPI_Allreduce(MPI_IN_PLACE, reinterpret_cast<double*>(Ceta2.data()),
		              2 * NP, MPI_DOUBLE, MPI_SUM, cfdwm_mpi::comm());
		MPI_Allreduce(MPI_IN_PLACE, reinterpret_cast<double*>(Cub.data()),
		              2 * NP, MPI_DOUBLE, MPI_SUM, cfdwm_mpi::comm());
		MPI_Allreduce(MPI_IN_PLACE, reinterpret_cast<double*>(Cvb.data()),
		              2 * NP, MPI_DOUBLE, MPI_SUM, cfdwm_mpi::comm());
		MPI_Allreduce(MPI_IN_PLACE, reinterpret_cast<double*>(Cws.data()),
		              2 * NP, MPI_DOUBLE, MPI_SUM, cfdwm_mpi::comm());
		if (want_p)
			MPI_Allreduce(MPI_IN_PLACE, reinterpret_cast<double*>(Cpt.data()),
			              2 * NP, MPI_DOUBLE, MPI_SUM, cfdwm_mpi::comm());
	}
#endif

	std::vector<std::complex<double>> work(NP), c(N);

	// x-fields independent of layer
	std::vector<double> eta1, eta2v, eta_x, eta_y, eta1_t;
	eval1(er, eta1, work);
	evalB(Ceta2, eta2v, work);
	for (int n = 0; n < N; ++n) c[n] = I * kxv[n] * er[n];
	eval1(c, eta_x, work);
	for (int n = 0; n < N; ++n) c[n] = I * kyv[n] * er[n];
	eval1(c, eta_y, work);
	if (want_p) {
		// d(eta1)/dt = sum A*omega*sin(phi)  (pressure only)
		for (int n = 0; n < N; ++n) c[n] = irr.omega[n] * (-I) * er[n];
		eval1(c, eta1_t, work);
	}

	eta_out.resize(NP);
	for (int m = 0; m < NP; ++m) eta_out[m] = eta1[m] + eta2v[m];

	u_out.resize((size_t)nl * NP);
	v_out.resize((size_t)nl * NP);
	w_out.resize((size_t)nl * NP);
	if (want_p) p_out->resize((size_t)nl * NP);

	// Per-layer transforms: layers are independent -> OpenMP-parallel with
	// per-layer private buffers. Under MPI each rank computes only the
	// layers it owns (round-robin) and the rows are broadcast afterwards,
	// so the redundant per-rank full fill is gone.
	int mpi_np = 1, mpi_rk = 0;
#if defined(MPI_enable)
	if (cfdwm_mpi::active()) {
		mpi_np = cfdwm_mpi::size();
		mpi_rk = cfdwm_mpi::rank();
	}
#endif
#pragma omp parallel for schedule(dynamic)
	for (int lay = 0; lay < nl; ++lay) {
		if (mpi_np > 1 && lay % mpi_np != mpi_rk) continue;   // MPI: not ours
		std::vector<std::complex<double>> wk(NP), cc(N), cb(NP);
		std::vector<double> u1, v1, phi1_s, Gf, Gf_s, Gx, Gy, ub, vb, p2s;
		std::vector<double> p1f, Gt, p2t;
		const double sp1 = 1.0 + svals_[lay];
		const double g1s = g * sp1;
		const std::vector<double>& C1v = Ch1[lay];
		const std::vector<double>& S1v = Sh1[lay];

		for (int n = 0; n < N; ++n)
			cc[n] = (g * irr.k[n] / irr.omega[n]) * C1v[n]
				* irr.cth_cache[n] * er[n];
		eval1(cc, u1, wk);
		for (int n = 0; n < N; ++n)
			cc[n] = (g * irr.k[n] / irr.omega[n]) * C1v[n]
				* irr.sth_cache[n] * er[n];
		eval1(cc, v1, wk);
		for (int n = 0; n < N; ++n)
			cc[n] = (g / irr.omega[n]) * irr.k[n] * depth_ * S1v[n]
				* (-I) * er[n];
		eval1(cc, phi1_s, wk);
		for (int n = 0; n < N; ++n)
			cc[n] = (irr.k[n] * S1v[n] / irr.omega[n]) * (-I) * er[n];
		eval1(cc, Gf, wk);
		for (int n = 0; n < N; ++n)
			cc[n] = (irr.k[n] / irr.omega[n])
				* (irr.k[n] * depth_ * C1v[n]) * (-I) * er[n];
		eval1(cc, Gf_s, wk);
		for (int n = 0; n < N; ++n)
			cc[n] = (irr.k[n] * S1v[n] / irr.omega[n]) * kxv[n] * er[n];
		eval1(cc, Gx, wk);
		for (int n = 0; n < N; ++n)
			cc[n] = (irr.k[n] * S1v[n] / irr.omega[n]) * kyv[n] * er[n];
		eval1(cc, Gy, wk);

		for (int j = 0; j < NP; ++j) cb[j] = Cub[j] * ChB[lay][j];
		evalB(cb, ub, wk);
		for (int j = 0; j < NP; ++j) cb[j] = Cvb[j] * ChB[lay][j];
		evalB(cb, vb, wk);
		for (int j = 0; j < NP; ++j) cb[j] = Cws[j] * ShB[lay][j];
		evalB(cb, p2s, wk);

		if (want_p) {
			// first-order dynamic pressure / rho: g * A * Ch * cos(phi)
			for (int n = 0; n < N; ++n) cc[n] = g * C1v[n] * er[n];
			eval1(cc, p1f, wk);
			// dG/dt = -sum A * k * Sh * cos(phi)
			for (int n = 0; n < N; ++n) cc[n] = -irr.k[n] * S1v[n] * er[n];
			eval1(cc, Gt, wk);
			// bound-potential time derivative at fixed s
			for (int j = 0; j < NP; ++j) cb[j] = Cpt[j] * ChB[lay][j];
			evalB(cb, p2t, wk);
		}

		double* urow = &u_out[(size_t)lay * NP];
		double* vrow = &v_out[(size_t)lay * NP];
		double* wrow = &w_out[(size_t)lay * NP];
		double* prow = want_p ? &(*p_out)[(size_t)lay * NP] : nullptr;
		for (int m = 0; m < NP; ++m) {
			const double eta = eta_out[m];
			const double Dcol = depth_ + eta;
			const double phi2_s = p2s[m] + g * eta1[m] * Gf[m]
				+ g1s * eta1[m] * Gf_s[m];
			const double phis = phi1_s[m] + phi2_s;
			const double cf = -(sp1 / Dcol);
			urow[m] = u1[m] + ub[m]
				+ g1s * (eta_x[m] * Gf[m] + eta1[m] * Gx[m])
				+ cf * eta_x[m] * phis;
			vrow[m] = v1[m] + vb[m]
				+ g1s * (eta_y[m] * Gf[m] + eta1[m] * Gy[m])
				+ cf * eta_y[m] * phis;
			wrow[m] = phis / Dcol;
			if (want_p) {
				/* dp_smit() in field form:
				     p = rho*[ p1 - dphi2/dt|_s
				               + (1+s)/D * eta1_t * dphi/ds
				               - 0.5*|u1|^2 ],
				   dphi2/dt|_s = p2t (bound) + g(1+s)*d(eta1*G)/dt
				   (the particular solution phi_P = g(1+s)*eta1*G
				   differentiates as a pointwise product). The Bernoulli
				   term uses the FIRST-order stretched velocities, and
				   w1 = dphi1/ds / D, exactly as in dp_smit(). */
				const double dphi2dt = p2t[m]
					+ g1s * (eta1_t[m] * Gf[m] + eta1[m] * Gt[m]);
				const double w1 = phi1_s[m] / Dcol;
				prow[m] = irr.rho * (p1f[m] - dphi2dt
					+ (sp1 / Dcol) * eta1_t[m] * phis
					- 0.5 * (u1[m] * u1[m] + v1[m] * v1[m] + w1 * w1));
			}
		}
	}
#if defined(MPI_enable)
	if (mpi_np > 1)
		for (int lay = 0; lay < nl; ++lay) {
			const int owner = lay % mpi_np;
			MPI_Bcast(&u_out[(size_t)lay * NP], NP, MPI_DOUBLE, owner, cfdwm_mpi::comm());
			MPI_Bcast(&v_out[(size_t)lay * NP], NP, MPI_DOUBLE, owner, cfdwm_mpi::comm());
			MPI_Bcast(&w_out[(size_t)lay * NP], NP, MPI_DOUBLE, owner, cfdwm_mpi::comm());
			if (want_p)
				MPI_Bcast(&(*p_out)[(size_t)lay * NP], NP, MPI_DOUBLE, owner, cfdwm_mpi::comm());
		}
#endif
}
