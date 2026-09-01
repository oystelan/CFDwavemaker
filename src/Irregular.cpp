#include "Irregular.h"
#include <math.h>
#include <algorithm>
#include <numeric>
#include <omp.h>



double Irregular::eta(double t, double x, double y) {
	// variable bathymetry: slope-corrected first order plus the
	// depth-local Sharma & Dean second order (full pair set)
	if (seabed.active()) {
		if (order >= 2)
			return (swl + eta1_slope(t, x, y) + eta2_slope(t, x, y));
		return (swl + eta1_slope(t, x, y));
	}
	// smit s-coordinate: export the same (untruncated) second-order
	// surface that the smit kinematics stretch and clamp against, so the
	// CFD-visible surface and the kinematics are consistent
	if (vertical_theory == 1 && order >= 1)
		return (swl + eta1(t, x, y) + eta_smit(t, x, y));
	switch(order){
		case 2:
			//return eta2(t,x,y);
			return (swl + eta1(t, x, y) + eta2(t, x, y));
		case 3:
			//return eta2(t,x,y);
			return (swl + eta1(t, x, y) + eta2_boost(t, x, y));
		default:
			//std::cout << eta1(t, x, y) << std::endl;
			return (swl + eta1(t, x, y));
			
	}
};

double Irregular::u(double t, double x, double y, double z) {
	if (seabed.active()) {
		double u1v = uvwp1_slope(t, x, y, z - swl).u;
		if (order >= 2) u1v += uvwp2_slope(t, x, y, z - swl).u;
		return u1v;
	}
	if (vertical_theory == 1) return u_smit(t, x, y, z);

	switch (order) {
		case 3:
			{double z0 = std::min(0., z - swl);
			//return u2(t,x,y,z0);}
			//std::cout << "ape" << std::endl;
			return (u1(t, x, y, z0) + u2_boost(t, x, y, z0) + phi1_dxdz(t, x, y) * std::max(0., z - swl));}
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[0];} 
		case 2:
			{double z0 = std::min(0., z - swl);
			//return u2(t,x,y,z0);}
			//std::cout << "ape" << std::endl;
			return (u1(t, x, y, z0) + u2(t, x, y, z0) + phi1_dxdz(t, x, y) * std::max(0., z - swl));}
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[0];} 
		case 1:
			{z = std::min(0., z - swl);
			//std::cout << "ape2" << std::endl;
			return (u1(t, x, y, z) + u2(t, x, y, z));} 
			
		default:
			{return u1(t, x, y, z - swl);} 
			
	}
};
double Irregular::v(double t, double x, double y, double z) {
	if (seabed.active()) {
		double v1v = uvwp1_slope(t, x, y, z - swl).v;
		if (order >= 2) v1v += uvwp2_slope(t, x, y, z - swl).v;
		return v1v;
	}
	if (vertical_theory == 1) return v_smit(t, x, y, z);
	switch (order) {
		case 3:
			{double z0 = std::min(0., z - swl);
			//return v2(t,x,y,z0);}
			return (v1(t, x, y, z0) + v2_boost(t, x, y, z0) + phi1_dydz(t, x, y) * std::max(0., z - swl)); }
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[1];} 
		case 2:
			{double z0 = std::min(0., z - swl);
			//return v2(t,x,y,z0);}
			return (v1(t, x, y, z0) + v2(t, x, y, z0) + phi1_dydz(t, x, y) * std::max(0., z - swl)); }
			//std::vector<double> U = uvw2(t, x, y, z);
			//return U[1];}  
		case 1:
			z = std::min(0., z - swl);
			return (v1(t, x, y, z) + v2(t, x, y, z));
			
		default:
			return v1(t, x, y, z - swl);
			
	}
};
double Irregular::w(double t, double x, double y, double z) {
	if (seabed.active()) {
		double w1v = uvwp1_slope(t, x, y, z - swl).w;
		if (order >= 2) w1v += uvwp2_slope(t, x, y, z - swl).w;
		return w1v;
	}
	if (vertical_theory == 1) return w_smit(t, x, y, z);
	switch (order) {
		case 3:
			{double z0 = std::min(0., z - swl);
			return (w1(t, x, y, z0) + w2_boost(t, x, y, z0) + phi1_dzdz(t, x, y) * std::max(0., z - swl)); }
			//return w2(t,x,y,z0);}  
		case 2:
			{double z0 = std::min(0., z - swl);
			return (w1(t, x, y, z0) + w2(t, x, y, z0) + phi1_dzdz(t, x, y) * std::max(0., z - swl)); }
			//return w2(t,x,y,z0);}   
		case 1:
			z = std::min(0., z - swl);
			return (w1(t, x, y, z) + w2(t, x, y, z));
			
		default:
			return w1(t, x, y, z - swl);
			
	}
};

/* velocity profile for horizontal component */
double Irregular::profileX(int ind, double x, double y, double z) {
	return (cosh(k[ind] * (z + depth)) / cosh(k[ind] * depth));
}

double Irregular::profileZ(int ind, double x, double y, double z) {
	return (sinh(k[ind] * (z + depth)) / cosh(k[ind] * depth));
}

/* Local wavenumber from the dispersion relation om^2 = g k tanh(k h)
   (Newton iteration from the deep-water start k = om^2/g). */
double Irregular::dispersion_k(double om, double h) const {
	double kk = om * om / g;
	for (int it = 0; it < 50; it++) {
		double th = tanh(std::min(kk * h, 50.));
		double f = g * kk * th - om * om;
		double df = g * th + g * kk * h * (1. - th * th);
		double dk = f / df;
		kk -= dk;
		if (std::fabs(dk) < 1e-14 * kk) break;
	}
	return kk;
}

/* Slope-corrected first-order kinematics over variable bathymetry
   (theory: seabed_development/sloping_seabed_notes.tex; validation:
   seabed_development/rb_residual/). Per evaluation point the depth and
   bed slope are probed from the seabed surface, each component's
   wavenumber is re-solved from the local dispersion relation, and the
   complex vertical eigenfunction

     F(z) = [cosh k(z+h) - i b sinh k(z+h)] / [cosh kh - i b sinh kh]

   with b = khat . grad h replaces the Airy profile. All horizontal
   derivatives use frozen coefficients (grad phi ~ i k phi). This makes
   the bottom kinematic BC w + grad h . u_H = 0 exact at z = -h, while
   the free-surface conditions keep an O(b) defect that vanishes in
   deep water (see FINDINGS_seabed.md). zz is measured from the swl. */
namespace {
// per-thread local-wavenumber cache: a water column (fixed x, y, hence
// fixed h) evaluates uvwp1_slope once per layer, order term and time
// level -- the Newton dispersion solves only depend on (owner, h, N)
struct FirstOrderLocal {
	const void* owner = nullptr;
	double h = 0.;
	std::vector<double> kl;
};
thread_local FirstOrderLocal first_order_cache;
} // namespace

UVWP Irregular::uvwp1_slope(double t, double xx, double yy, double zz) {
	double h = seabed.depth(xx, yy);
	double hx, hy;
	seabed.gradient(xx, yy, hx, hy);

	const int N = ndir * nfreq;
	FirstOrderLocal& f = first_order_cache;
	if (f.owner != this || f.h != h || (int)f.kl.size() != N) {
		f.owner = this; f.h = h; f.kl.resize(N);
		for (int i = 0; i < N; i++)
			f.kl[i] = dispersion_k(omega[i], h);
	}

	UVWP res;
	res.u = 0.; res.v = 0.; res.w = 0.; res.p = 0.;
	const std::complex<double> I(0., 1.);
	for (int i = 0; i < N; i++) {
		double thi = theta[i] + (mtheta * PI / 180.);
		double cth = cos(thi), sth = sin(thi);
		double kl = f.kl[i];
		double beta = cth * hx + sth * hy; // khat . grad h
		double arg = kl * (cth * (xx - fpoint[0]) + sth * (yy - fpoint[1]))
			- omega[i] * t + omega[i] * tofmax + phase[i];
		double kzh = kl * std::max(zz + h, 0.); // clamp below the bed
		double kh = kl * h;
		double C = cosh(kzh) / cosh(kh);
		double S = sinh(kzh) / cosh(kh);
		std::complex<double> den(1., -beta * tanh(kh));
		std::complex<double> F = (C - I * beta * S) / den;
		std::complex<double> Fp = (S - I * beta * C) / den; // dF/dz / k
		std::complex<double> E(cos(arg), sin(arg));
		double ramp = freq_ramp(t, freq_ramptime[i], freq_starttime[i],
			freq_duration[i]);
		double fac = A[i] * g / omega[i] * ramp;
		double FE = std::real(F * E);
		res.u += fac * kl * cth * FE;
		res.v += fac * kl * sth * FE;
		res.w += fac * kl * std::imag(Fp * E);
		res.p += rho * g * A[i] * ramp * FE;
	}
	return res;
}

/* First-order elevation over variable bathymetry. The slope-corrected
   potential gives eta1 = -(1/g) dphi/dt at z = 0 with F(0) = 1
   identically, i.e. the plain component sum but with the LOCAL
   wavenumbers in the phase - the same phase function as uvwp1_slope
   and eta2_slope, keeping surface and kinematics in phase on a slope. */
double Irregular::eta1_slope(double t, double xx, double yy) {
	double h = seabed.depth(xx, yy);

	const int N = ndir * nfreq;
	FirstOrderLocal& f = first_order_cache;
	if (f.owner != this || f.h != h || (int)f.kl.size() != N) {
		f.owner = this; f.h = h; f.kl.resize(N);
		for (int i = 0; i < N; i++)
			f.kl[i] = dispersion_k(omega[i], h);
	}

	double welev = 0.;
	for (int i = 0; i < N; i++) {
		double thi = theta[i] + (mtheta * PI / 180.);
		double arg = f.kl[i] * (cos(thi) * (xx - fpoint[0])
			+ sin(thi) * (yy - fpoint[1]))
			- omega[i] * t + omega[i] * tofmax + phase[i];
		welev += A[i] * cos(arg)
			* freq_ramp(t, freq_ramptime[i], freq_starttime[i],
				freq_duration[i]);
	}
	return welev;
}

/* Build the depth-local Sharma & Dean interaction coefficients for
   depth h: local wavenumbers from the local dispersion relation, then
   the same D/alpha formulas as rebuild_pair_cache(). Bandwidth/cutoff
   gates (dw_bandwidth, dw_cutoff) are honoured; cross sum interactions
   carry the interaction-limit flag kmax_p (see Irregular.h). */
void Irregular::build_local_pair_set(double h, LocalPairSet& s) const {
	const int N = nfreq * ndir;
	const double mt = mtheta * PI / 180.;

	s.h = h;
	s.owner = this;
	s.kl.resize(N); s.Rl.resize(N); s.Rsl.resize(N); s.c2kd.resize(N);
	s.cth.resize(N); s.sth.resize(N);
	s.diagD.resize(N);
	s.pairs.clear();
	s.pairs.reserve(N * (N - 1) / 2);

	for (int i = 0; i < N; i++) {
		const double kl = dispersion_k(omega[i], h);
		s.kl[i] = kl;
		s.Rl[i] = kl * tanh(kl * h);
		s.Rsl[i] = sqrt(s.Rl[i]);
		s.c2kd[i] = cosh(2. * kl * h);
		s.cth[i] = cos(theta[i] + mt);
		s.sth[i] = sin(theta[i] + mt);
		s.diagD[i] = (6. * s.Rl[i] * (kl * kl - s.Rl[i] * s.Rl[i]))
			/ (2. * s.Rl[i] - kl * tanh(2. * kl * h));
	}
	double kmax_l = 0.;
	for (int i = 0; i < N; i++) kmax_l = std::max(kmax_l, s.kl[i]);

	for (int i = 0; i < N; i++) {
		if (omega[i] > dw_cutoff) continue;
		const double ki = s.kl[i], Rni = s.Rl[i], Rsni = s.Rsl[i];
		for (int m = i + 1; m < N; m++) {
			if (omega[m] > dw_cutoff) continue;
			if (std::abs(omega[i] - omega[m]) > dw_bandwidth) continue;

			const double km = s.kl[m], Rnm = s.Rl[m], Rsm = s.Rsl[m];
			const double gamma_nm = s.cth[i] * s.cth[m] + s.sth[i] * s.sth[m];
			const double ki2 = ki * ki, km2 = km * km, kikm = ki * km;

			LocalPairSet::P p{};
			p.i = i; p.m = m;
			p.kmax_p = ((ki + km) <= kmax_l); // interaction limit
			p.k_p = sqrt(ki2 + km2 + 2. * kikm * gamma_nm);
			p.k_mn = sqrt(ki2 + km2 - 2. * kikm * gamma_nm);

			const double sp = Rsni + Rsm, sm2 = Rsni - Rsm;
			{
				const double R2p = sp * sp;
				const double den = R2p - p.k_p * tanh(p.k_p * h);
				p.D_p = (2. * R2p * (kikm * gamma_nm - Rni * Rnm)) / den
					+ (sp * (Rsni * (km2 - Rnm * Rnm)
						+ Rsm * (ki2 - Rni * Rni))) / den;
				p.a_p = (p.D_p - (kikm * gamma_nm - Rni * Rnm))
					/ sqrt(Rni * Rnm) + (Rni + Rnm);
			}
			p.has_minus = (ki != km);
			if (p.has_minus) {
				const double R2m = sm2 * sm2;
				const double den = R2m - p.k_mn * tanh(p.k_mn * h);
				p.D_mn = (2. * R2m * (kikm * gamma_nm + Rni * Rnm)) / den
					+ (sm2 * (-Rsni * (km2 - Rnm * Rnm)
						+ Rsm * (ki2 - Rni * Rni))) / den;
				p.a_mn = (p.D_mn - (kikm * gamma_nm + Rni * Rnm))
					/ sqrt(Rni * Rnm) + (Rni + Rnm);
			}
			s.pairs.push_back(p);
		}
	}
}

namespace {
// per-thread coefficient cache keyed on the local depth: a water
// column (fixed x, y) reuses the set across layers and time levels
thread_local Irregular::LocalPairSet local_pair_cache;
}

/* Depth-local second-order surface elevation (full pair set). */
double Irregular::eta2_slope(double t, double xx, double yy) {
	const int N = nfreq * ndir;
	const double h = seabed.depth(xx, yy);
	LocalPairSet& s = local_pair_cache;
	if (s.owner != this || s.h != h || (int)s.kl.size() != N)
		build_local_pair_set(h, s);

	const double dxp = xx - fpoint[0], dyp = yy - fpoint[1];
	std::vector<double> cphi(N), sphi(N), rmp(N);
	for (int i = 0; i < N; i++) {
		const double ph = s.kl[i] * (s.cth[i] * dxp + s.sth[i] * dyp)
			- omega[i] * t + omega[i] * tofmax + phase[i];
		cphi[i] = cos(ph); sphi[i] = sin(ph);
		rmp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i],
			freq_duration[i]);
	}

	double e2 = 0.;
	for (int i = 0; i < N; i++) {
		if (omega[i] > dw_cutoff) continue;
		const double alpha_nn = (s.diagD[i]
			- (s.kl[i] * s.kl[i] - s.Rl[i] * s.Rl[i])) / s.Rl[i]
			+ 2. * s.Rl[i];
		const double cos2 = cphi[i] * cphi[i] - sphi[i] * sphi[i];
		e2 += 0.25 * A[i] * A[i] * alpha_nn * cos2 * rmp[i];
	}
	for (size_t q = 0; q < s.pairs.size(); q++) {
		const LocalPairSet::P& p = s.pairs[q];
		const double amp = 0.5 * A[p.i] * A[p.m] * rmp[p.i];
		if (p.kmax_p)
			e2 += amp * p.a_p
				* (cphi[p.i] * cphi[p.m] - sphi[p.i] * sphi[p.m]);
		if (p.has_minus)
			e2 += amp * p.a_mn
				* (cphi[p.i] * cphi[p.m] + sphi[p.i] * sphi[p.m]);
	}
	return e2;
}

/* Depth-local second-order kinematics (full pair set). Mirrors uvw2()
   with local coefficients; the Bernoulli term uses the slope-corrected
   first-order velocities (phi = phi1_mod + phi2_SD). The returned p
   INCLUDES the rho factor (unlike uvw2). zz measured from swl. */
UVWP Irregular::uvwp2_slope(double t, double xx, double yy, double zz) {
	const int N = nfreq * ndir;
	const double h = seabed.depth(xx, yy);
	LocalPairSet& s = local_pair_cache;
	if (s.owner != this || s.h != h || (int)s.kl.size() != N)
		build_local_pair_set(h, s);

	const double dxp = xx - fpoint[0], dyp = yy - fpoint[1];
	const double zpd = std::max(zz + h, 0.);
	std::vector<double> cphi(N), sphi(N), rmp(N);
	for (int i = 0; i < N; i++) {
		const double ph = s.kl[i] * (s.cth[i] * dxp + s.sth[i] * dyp)
			- omega[i] * t + omega[i] * tofmax + phase[i];
		cphi[i] = cos(ph); sphi[i] = sin(ph);
		rmp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i],
			freq_duration[i]);
	}

	double us = 0., vs = 0., ws = 0., ps = 0.;
	for (int i = 0; i < N; i++) {
		if (omega[i] > dw_cutoff) continue;
		const double ki = s.kl[i], oi = omega[i];
		const double agi = A[i] * g / oi;
		const double AA = 0.25 * agi * agi * (s.diagD[i] / (2. * oi));
		const double cos2 = cphi[i] * cphi[i] - sphi[i] * sphi[i];
		const double sin2 = 2. * sphi[i] * cphi[i];
		const double cc = cosh(2. * ki * zpd) / s.c2kd[i];
		const double sc = sinh(2. * ki * zpd) / s.c2kd[i];
		us += AA * (2. * ki * s.cth[i]) * cos2 * cc * rmp[i];
		vs += AA * (2. * ki * s.sth[i]) * cos2 * cc * rmp[i];
		ws += AA * (2. * ki) * sin2 * sc * rmp[i];
		ps += AA * (2. * oi) * cos2 * cc * rmp[i];
	}
	for (size_t q = 0; q < s.pairs.size(); q++) {
		const LocalPairSet::P& p = s.pairs[q];
		const int i = p.i, m = p.m;
		const double ki = s.kl[i], km = s.kl[m];
		const double oi = omega[i], om = omega[m];
		const double AA = 0.5 * (A[i] * g / oi) * (A[m] * g / om);
		const double ri = rmp[i];
		if (p.kmax_p) {
			const double BB = p.D_p / (oi + om);
			const double cp = cphi[i] * cphi[m] - sphi[i] * sphi[m];
			const double sp2 = sphi[i] * cphi[m] + cphi[i] * sphi[m];
			const double den = cosh(p.k_p * h);
			const double cc = cosh(p.k_p * zpd) / den;
			const double sc = sinh(p.k_p * zpd) / den;
			us += AA * BB * (ki * s.cth[i] + km * s.cth[m]) * cp * cc * ri;
			vs += AA * BB * (ki * s.sth[i] + km * s.sth[m]) * cp * cc * ri;
			ws += AA * BB * p.k_p * sp2 * sc * ri;
			ps += AA * BB * (oi + om) * cp * cc * ri;
		}
		if (p.has_minus) {
			const double CC = p.D_mn / (oi - om);
			const double cm = cphi[i] * cphi[m] + sphi[i] * sphi[m];
			const double sm2 = sphi[i] * cphi[m] - cphi[i] * sphi[m];
			const double den = cosh(p.k_mn * h);
			const double cc = cosh(p.k_mn * zpd) / den;
			const double sc = sinh(p.k_mn * zpd) / den;
			us += AA * CC * (ki * s.cth[i] - km * s.cth[m]) * cm * cc * ri;
			vs += AA * CC * (ki * s.sth[i] - km * s.sth[m]) * cm * cc * ri;
			ws += AA * CC * p.k_mn * sm2 * sc * ri;
			ps += AA * CC * (oi - om) * cm * cc * ri;
		}
	}

	// quadratic Bernoulli term from the slope-corrected first order
	{
		UVWP s1 = uvwp1_slope(t, xx, yy, zz);
		double u1v = s1.u, v1v = s1.v, w1v = s1.w;
		ps -= 0.5 * (u1v * u1v + v1v * v1v + w1v * w1v);
	}

	UVWP res;
	res.u = us; res.v = vs; res.w = ws;
	res.p = rho * ps;
	return res;
}

/* Smit et al. (2017) s-coordinate stretched profiles.
   The vertical argument k*(z+h) is replaced by k*h*(h+z)/(h+eta),
   which maps the instantaneous water column [−h, eta] onto the
   reference column [−h, 0]. This is equivalent to Wheeler stretching. */
double Irregular::profileX_s(int ind, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return cosh(k[ind] * sf) / cosh(k[ind] * depth);
}

double Irregular::profileZ_s(int ind, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return sinh(k[ind] * sf) / cosh(k[ind] * depth);
}

double Irregular::profileX_s2(double k_nm, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return cosh(k_nm * sf) / cosh(k_nm * depth);
}

double Irregular::profileZ_s2(double k_nm, double z, double eta) {
	double sf = depth * (depth + z) / (depth + eta);
	return sinh(k_nm * sf) / cosh(k_nm * depth);
}

double Irregular::dp(double t, double x, double y, double z) {
	if (seabed.active()) {
		double p1v = uvwp1_slope(t, x, y, z - swl).p;
		if (order >= 2) p1v += uvwp2_slope(t, x, y, z - swl).p;
		return p1v;
	}
	if (vertical_theory == 1) return dp_smit(t, x, y, z);
	switch (order) {
		case 2:
			{double z0 = std::min(0., z - swl);
			return (dp1(t, x, y, z0) + dp2(t, x, y, z0) + dp1_dz(t, x, y) * std::max(0., z - swl)); }
			//return w2(t,x,y,z0);}   
		case 1:
			z = std::min(0., z - swl);
			return (dp1(t, x, y, z) + dp2(t, x, y, z));
			
		default:
			return dp1(t, x, y, z - swl);
			
	}
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
			+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi) *
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
		//std::cout << i << ": A" << A[i] << std::endl;
	}


	return welev;

}

double Irregular::freq_ramp(double t, double ramptime, double tstart, double duration) {
	double rampdown = std::max(std::min(1. - ((t - tstart- duration + ramptime) / (ramptime + SMALL)),1.), 0.);
	double rampup = std::max(std::min(((t - tstart) / (ramptime + SMALL)), 1.), 0.);
	return rampup*rampdown;
}

int Irregular::bw_limiter(int i1) {
	int cc = i1;
	// bounds check FIRST: with the old order omega[cc] was read one past
	// the end before cc < N stopped the loop
	while ((cc < nfreq * ndir) && (omega[cc] - omega[i1] < dw_bandwidth)) {
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


/* Detect whether the components form a long-crested equal-dk wavenumber
   grid; if so, precompute the bin mapping for factorized evaluation. */
void Irregular::factored_init()
{
    // NOTE: fact_checked is stored by the caller AFTER this returns (inside
    // fact_mutex); storing it here/early would let concurrent threads bypass
    // the mutex and read half-built state without memory ordering.
    fact_ok = false;
    const int N = nfreq * ndir;
    if (N < 2) return;

    // long-crested: all components propagate in the same direction
    for (int i = 1; i < N; ++i) {
        if (std::abs(theta[i] - theta[0]) > 1.e-12) {
            std::cout << "equal_dk: components are not long-crested, "
                      << "falling back to direct second-order evaluation."
                      << std::endl;
            return;
        }
    }

    // initial dk estimate: smallest spacing (k is sorted ascending)
    double dk0 = HUGE;
    for (int i = 1; i < N; ++i) {
        dk0 = std::min(dk0, std::abs(k[i] - k[i - 1]));
    }
    if (dk0 <= SMALL) return;

    // integer bin indices, then refine dk by least squares through origin
    fact_n.resize(N);
    double skn = 0.0, snn = 0.0;
    for (int i = 0; i < N; ++i) {
        fact_n[i] = (int)std::lround(k[i] / dk0);
        if (fact_n[i] < 1) { fact_n.clear(); return; }
        skn += k[i] * fact_n[i];
        snn += (double)fact_n[i] * (double)fact_n[i];
    }
    const double dk = skn / snn;

    int nmax = 0;
    for (int i = 0; i < N; ++i) {
        if (std::abs(fact_n[i] * dk - k[i]) > 1.e-5 * std::max(1.0, k[i])) {
            std::cout << "equal_dk: wavenumbers are not on an equal-dk grid, "
                      << "falling back to direct second-order evaluation."
                      << std::endl;
            fact_n.clear();
            return;
        }
        nmax = std::max(nmax, fact_n[i]);
    }

    fact_dk = dk;
    fact_nbins = 2 * nmax + 1;

    fact_ok = true;
    std::cout << "Second order: factorized equal-dk evaluation enabled "
              << "(dk = " << dk << " 1/m, " << fact_nbins << " bins)."
              << std::endl;
}

/* Fetch (or build) the immutable coefficient set for time level t.
   Multiple time levels can be requested concurrently (the spline lsgrid
   fills four knots in one parallel region), so entries live in a small
   cache; an entry is completed under fact_mutex before it is returned and
   never mutated afterwards, except for adding the independent uvw part
   (also under the mutex). */
std::shared_ptr<Irregular::FactCoeffs> Irregular::fact_get(double t,
                                                           bool need_uvw)
{
    std::lock_guard<std::mutex> lock(fact_mutex);
    for (auto& e : fact_cache) {
        if (e->t == t) {
            if (need_uvw && !e->has_uvw) fact_build_uvw(*e, t);
            if (!need_uvw && !e->has_eta2) fact_build_eta2(*e, t);
            return e;
        }
    }
    auto e = std::make_shared<FactCoeffs>();
    e->t = t;
    if (need_uvw) fact_build_uvw(*e, t);
    else fact_build_eta2(*e, t);
    if ((int)fact_cache.size() >= 8) {
        fact_cache.erase(fact_cache.begin());
    }
    fact_cache.push_back(e);
    return e;
}

/* Accumulate the complex wavenumber-bin coefficients C_j(t) such that
   eta2(x) = sum_j Re{ C_j * exp(i*j*dk*X) }, X along the wave direction.
   O(N^2) once per time level; mirrors the direct eta2() term by term. */
void Irregular::fact_build_eta2(FactCoeffs& c, double t)
{
    const int N = nfreq * ndir;
    c.eta2_bins.assign(fact_nbins, std::complex<double>(0.0, 0.0));
    if (c.ramp.empty()) {
        c.ramp.resize(N);
        for (int i = 0; i < N; ++i) {
            c.ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i],
                                  freq_duration[i]);
        }
    }

    std::vector<std::complex<double>> E(N);
    for (int i = 0; i < N; ++i) {
        const double psi = -omega[i] * t + phase0_cache[i];
        E[i] = std::complex<double>(std::cos(psi), std::sin(psi));
    }

    // diagonal terms -> bin 2*n_i
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) continue;

        const double ki = k[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];

        const double D_nn_plus =
            (6.0 * Rsni * Rsni * (ki * ki - Rni * Rni)) /
            (2.0 * Rsni * Rsni - ki * std::tanh(2.0 * ki * depth));
        const double alpha_nn_plus =
            (D_nn_plus - (ki * ki - Rni * Rni)) / Rni + 2.0 * Rni;

        if (2.0 * ki <= kmax_cache) {
            c.eta2_bins[2 * fact_n[i]] +=
                (0.25 * A[i] * A[i] * alpha_nn_plus * c.ramp[i]) * E[i] * E[i];
        }
    }

    // pair terms -> bins n_i + n_m (sum) and n_m - n_i (difference)
    for (const auto& pair : pair_cache) {
        const int i = pair.i;
        const int m = pair.m;
        const double amp = 0.5 * A[i] * A[m] * c.ramp[i];

        if (pair.has_plus && pair.kmax_plus) {
            c.eta2_bins[fact_n[i] + fact_n[m]] +=
                (amp * pair.alpha_nm_plus) * E[i] * E[m];
        }
        if (pair.has_minus) {
            // cos(phi_i - phi_m) with k_m > k_i: bin n_m - n_i, conjugated
            c.eta2_bins[fact_n[m] - fact_n[i]] +=
                (amp * pair.alpha_nm_minus) * E[m] * std::conj(E[i]);
        }
    }
    c.has_eta2 = true;
}

/* Evaluate the collapsed Fourier series at one point: O(n_bins). */
double Irregular::eta2_factored_eval(const FactCoeffs& c, double xx,
                                     double yy)
{
    const double X = cth_cache[0] * (xx - fpoint[0])
                   + sth_cache[0] * (yy - fpoint[1]);
    const std::complex<double> rot(std::cos(fact_dk * X),
                                   std::sin(fact_dk * X));
    std::complex<double> Ej(1.0, 0.0);
    double s = 0.0;
    for (int j = 0; j < fact_nbins; ++j) {
        s += c.eta2_bins[j].real() * Ej.real()
           - c.eta2_bins[j].imag() * Ej.imag();
        Ej *= rot;
    }
    return s;
}

double Irregular::eta2(double t, double xx, double yy)
{
    const int N = nfreq * ndir;
    const double dxp = xx - fpoint[0];
    const double dyp = yy - fpoint[1];

    // factorized path for long-crested equal-dk grids: the O(N^2) pair sum
    // runs once per time level; every point evaluation is O(n_bins).
    // The cache rebuild must be serialized here: the factored update writes
    // through bin indices derived from pair_cache, so a concurrent lazy
    // rebuild (benign for the read-only direct path) would corrupt memory.
    if (equal_dk) {
        if (!fact_checked.load()) {
            std::lock_guard<std::mutex> lock(fact_mutex);
            if (!fact_checked.load()) {
                if ((int)cth_cache.size() != N) {
                    rebuild_second_order_cache();
                }
                factored_init();
                fact_checked.store(true);   // release: init fully visible
            }
        }
        if (fact_ok) {
            auto c = fact_get(t, false);
            return eta2_factored_eval(*c, xx, yy);
        }
    }

    if ((int)cth_cache.size() != N) {
        rebuild_second_order_cache();
    }

    double eta2_total = 0.0;

    std::vector<double> cphi(N), sphi(N), ramp(N);

    for (int i = 0; i < N; ++i) {
        const double ph =
            k[i] * (cth_cache[i] * dxp + sth_cache[i] * dyp)
            - omega[i] * t
            + phase0_cache[i];

        cphi[i] = std::cos(ph);
        sphi[i] = std::sin(ph);
        ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
    }

    // diagonal terms
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) {
            continue;
        }

        const double ki = k[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];

        const double D_nn_plus =
            (6.0 * Rsni * Rsni * (ki * ki - Rni * Rni)) /
            (2.0 * Rsni * Rsni - ki * std::tanh(2.0 * ki * depth));

        const double alpha_nn_plus =
            (D_nn_plus - (ki * ki - Rni * Rni)) / Rni + 2.0 * Rni;

        if (2.0 * ki <= kmax_cache) {
            const double cos_2phi_i = cphi[i] * cphi[i] - sphi[i] * sphi[i];
            eta2_total += 0.25 * A[i] * A[i] * alpha_nn_plus * cos_2phi_i * ramp[i];
        }
    }

    // off-diagonal terms from pair cache
    for (const auto& pair : pair_cache) {
        const int i = pair.i;
        const int m = pair.m;

        const double amp = 0.5 * A[i] * A[m];
        const double ri = ramp[i];  // preserves your original weighting convention

        if (pair.has_plus && pair.kmax_plus) {
            const double cos_phi_plus =
                cphi[i] * cphi[m] - sphi[i] * sphi[m];
            eta2_total += amp * pair.alpha_nm_plus * cos_phi_plus * ri;
        }

        if (pair.has_minus) {
            const double cos_phi_minus =
                cphi[i] * cphi[m] + sphi[i] * sphi[m];
            eta2_total += amp * pair.alpha_nm_minus * cos_phi_minus * ri;
        }
    }

    return eta2_total;
}

/* Untruncated second-order elevation for the smit s-coordinate path.
   Same Sharma & Dean interaction coefficients as eta2(), but WITHOUT the
   kmax truncation: the sigma-theory free surface must carry the full
   bound harmonics for the dynamic BC p(eta) = rho*g*eta to close (the
   kmax convention drops ALL self-interactions of the upper half of the
   spectrum - for a single component it removes the entire Stokes second
   harmonic). Mean set-down (i=m difference terms) remains excluded,
   matching the S&D convention. */
double Irregular::eta2_full(double t, double xx, double yy)
{
    const int N = nfreq * ndir;
    const double dxp = xx - fpoint[0];
    const double dyp = yy - fpoint[1];

    if ((int)cth_cache.size() != N) {
        rebuild_second_order_cache();
    }

    double eta2_total = 0.0;

    std::vector<double> cphi(N), sphi(N), ramp(N);

    for (int i = 0; i < N; ++i) {
        const double ph =
            k[i] * (cth_cache[i] * dxp + sth_cache[i] * dyp)
            - omega[i] * t
            + phase0_cache[i];

        cphi[i] = std::cos(ph);
        sphi[i] = std::sin(ph);
        ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
    }

    // diagonal terms (no 2k <= kmax gate)
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) {
            continue;
        }

        const double ki = k[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];

        const double D_nn_plus =
            (6.0 * Rsni * Rsni * (ki * ki - Rni * Rni)) /
            (2.0 * Rsni * Rsni - ki * std::tanh(2.0 * ki * depth));

        const double alpha_nn_plus =
            (D_nn_plus - (ki * ki - Rni * Rni)) / Rni + 2.0 * Rni;

        const double cos_2phi_i = cphi[i] * cphi[i] - sphi[i] * sphi[i];
        eta2_total += 0.25 * A[i] * A[i] * alpha_nn_plus * cos_2phi_i * ramp[i];
    }

    // off-diagonal terms (no kmax gate on the sum interaction)
    for (const auto& pair : pair_cache) {
        const int i = pair.i;
        const int m = pair.m;

        const double amp = 0.5 * A[i] * A[m];
        const double ri = ramp[i];

        if (pair.has_plus) {
            const double cos_phi_plus =
                cphi[i] * cphi[m] - sphi[i] * sphi[m];
            eta2_total += amp * pair.alpha_nm_plus * cos_phi_plus * ri;
        }

        if (pair.has_minus) {
            const double cos_phi_minus =
                cphi[i] * cphi[m] + sphi[i] * sphi[m];
            eta2_total += amp * pair.alpha_nm_minus * cos_phi_minus * ri;
        }
    }

    return eta2_total;
}

/* Function which precomputes all second order terms which are not a function of space and time. Speeds up second order calculation substantially */
void Irregular::second_order_init() {
	// we start by allocating memory
	eta2_diag.resize(nfreq*ndir);
	eta2_sum.resize(nfreq*ndir*nfreq*ndir);
	eta2_diff.resize(nfreq*ndir*nfreq*ndir);
	uvw2_diag.resize(nfreq*ndir);
	uvw2_sum.resize(nfreq*ndir*nfreq*ndir);
	uvw2_diff.resize(nfreq*ndir*nfreq*ndir);
	double kmax = *max_element(k.begin(), k.end());
	
	for (int i = 0; i < nfreq*ndir; i++) {
		//double eta2_t = 0;
		if (omega[i] <= dw_cutoff) {
			double thi = theta[i] + (mtheta * PI / 180.);
			//double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
			double Rn = k[i] * tanh(k[i] * depth);
			double Rsn = sqrt(Rn);

			// Diagonal terms
			double D_nn_plus = (6. * Rsn * Rsn * (k[i] * k[i] - Rn * Rn)) / (2. * Rsn * Rsn - k[i] * tanh(2. * k[i] * depth));
			double alpha_nn_plus = (D_nn_plus - (k[i] * k[i] - Rn * Rn)) / Rn + 2. * Rn;

			if (2 * k[i] <= kmax) {
				eta2_diag[i] = 0.25 * A[i] * A[i] * alpha_nn_plus;// * cos(2 * phi_i);
				uvw2_diag[i] = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]))*2*k[i];
			}
			



			// Off-diagonal (sum over half only due to symmetry)
			for (int m = i + 1; m < nfreq * ndir; m++) {
				
				if ((omega[m] <= dw_cutoff) && abs(omega[i] - omega[m]) <= dw_bandwidth) {
					double thm = theta[m] + (mtheta * PI / 180.);
					//double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

					double gamma_nm = cos(thi - thm);
					double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
					double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

					double Rm = k[m] * tanh(k[m] * depth);
					double Rsm = sqrt(Rm);
					double Rsnm2plus = pow(Rsn + Rsm, 2.);
					double Rsnm2minus = pow(Rsn - Rsm, 2.);


					if ((k[i] + k[m]) <= kmax) {
						double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) +
							((Rsn + Rsm) * (Rsn * (k[m] * k[m] - Rm * Rm) + Rsm * (k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));
						double alpha_nm_plus = (D_nm_plus - (k[i] * k[m] * gamma_nm - Rn * Rm)) / sqrt(Rn * Rm) + (Rn + Rm);
						eta2_sum[i*(nfreq * ndir)+ m] = 0.5 * A[i] * A[m] * alpha_nm_plus;// * cos(phi_i + phi_m);
						uvw2_sum[i*(nfreq * ndir)+ m] = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]) * (D_nm_plus/(omega[i]+omega[m]));
						
					}
					
					if (k[i] != k[m]) {
						double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) 
							+ ((Rsn - Rsm) * (-Rsn * (k[m] * k[m] - Rm * Rm) + Rsm * (k[i] * k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));
						double alpha_nm_minus = (D_nm_minus - (k[i] * k[m] * gamma_nm + Rn * Rm)) / sqrt(Rn * Rm) + (Rn + Rm);
						eta2_diff[i*(nfreq * ndir)+ m] = 0.5 * A[i] * A[m] * alpha_nm_minus;// * cos(phi_i - phi_m); // difference term
						uvw2_diff[i*(nfreq * ndir)+ m] = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]) * (D_nm_minus/(omega[i]-omega[m]));
						
					}
					
				}

			}
		}
	}
}

/* Second order wave elevation */
double Irregular::eta2_boost(double t, double xx, double yy) {
	double eta2_ti = 0.;
	#pragma omp parallel for reduction(+:eta2_ti)
	for (int i = 0; i < nfreq*ndir; i++) {
		double eta2_t = 0;
		
		double thi = theta[i] + (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		
		eta2_t += eta2_diag[i]*cos(2 * phi_i);

		// Off-diagonal (sum over half only due to symmetry)
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m] + (mtheta * PI / 180.);
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			eta2_t += eta2_sum[i*(nfreq * ndir)+ m] * cos(phi_i + phi_m); // sum term			
			eta2_t += eta2_diff[i*(nfreq * ndir)+ m] * cos(phi_i - phi_m); // difference term
		}
		eta2_ti += eta2_t*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return eta2_ti;
}


double Irregular::u1(double t, double xx, double yy, double zz) {

	double usum = 0.0;
	double phi;
	//(cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))
	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		usum +=  A[i] * g * k[i] / omega[i] * profileX(i, xx, yy, zz)
			* cos(theta[i] + (mtheta * PI / 180.)) * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.))
				* (xx - fpoint[0]) + sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return usum;

}


/* horizontal velocity U for a sinus Irregular::omegaave */
double Irregular::v1(double t, double xx, double yy, double zz) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		vsum += A[i] * g * k[i] / omega[i] * profileX(i, xx, yy, zz)
			* sin(theta[i] + (mtheta * PI / 180.)) * cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0])
				+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double Irregular::w1(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		wsum += A[i] * g * k[i] / omega[i] * profileZ(i, xx, yy, zz)
			* sin(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] + (mtheta 
				* PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
			freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return wsum;

}

/* Second order horizontal velocity U for a sinus wave */
double Irregular::u2(double t, double xx, double yy, double zz) {

	double usum2_t = 0;
	
	for (int i = 0; i < nfreq * ndir; i++) {
		double usum2 = 0;
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));		
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		usum2 += AA*(2*k[i]*cos(thi)) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth);

		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);

			// Sum terms
			usum2 += AA * BB * (k[i]*cos(thi) + k[m]*cos(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				
				// Diff terms
				usum2 += AA * CC * (k[i]*cos(thi)-k[m]*cos(thm)) * cos(phi_i-phi_m) * coshcosh_minus;
			}	
		}
		usum2_t += usum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return usum2_t;
}

/* Second order horizontal velocity U for a sinus wave */
double Irregular::u2_boost(double t, double xx, double yy, double zz) {

	double usum2_t = 0;
	#pragma omp parallel for reduction(+:usum2_t)
	for (int i = 0; i < nfreq * ndir; i++) {
		double usum2 = 0;
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];

		// Diagonal terms
		usum2 += uvw2_diag[i] * cos(thi) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth);

		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];
			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);

			// Sum terms
			usum2 +=  uvw2_sum[i*(nfreq * ndir)+ m] * (k[i]*cos(thi) + k[m]*cos(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
	
			// Diff terms
			usum2 += uvw2_diff[i*(nfreq * ndir)+ m] * (k[i]*cos(thi) - k[m]*cos(thm)) * cos(phi_i - phi_m) * coshcosh_minus;
		}
		usum2_t += usum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return usum2_t;
}


/* Second order horizontal velocity V for a sinus wave */
double Irregular::v2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double vsum2_t = 0;

	for (int i = 0; i < nfreq * ndir; i++) {
		double vsum2 = 0;		
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));		
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		vsum2 += AA*(2*k[i]*sin(thi)) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) 
			+ ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			
			// Sum terms
			vsum2 += AA * BB * (k[i]*sin(thi) + k[m]*sin(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				
				// Diff terms
				vsum2 += AA * CC * (k[i]*sin(thi)-k[m]*sin(thm)) * cos(phi_i-phi_m) * coshcosh_minus;
			}					
		}
		vsum2_t += vsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);	
	}
	return vsum2_t;
}

/* Second order horizontal velocity V for a sinus wave */
double Irregular::v2_boost(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double vsum2_t = 0;
	#pragma omp parallel for reduction(+:vsum2_t)
	for (int i = 0; i < nfreq * ndir; i++) {
		double vsum2 = 0;		
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
	
		// Diagonal terms
		vsum2 += uvw2_diag[i] * sin(thi) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
			
			vsum2 += uvw2_sum[i*(nfreq * ndir)+ m] * (k[i]*sin(thi) + k[m]*sin(thm)) * cos(phi_i + phi_m) * coshcosh_plus;
			vsum2 += uvw2_diff[i*(nfreq * ndir)+ m] * (k[i]*sin(thi) - k[m]*sin(thm)) * cos(phi_i - phi_m) * coshcosh_minus;					
		}
		vsum2_t += vsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);	
	}
	return vsum2_t;
}

/* Second order vertical velocity component W for a sinus wave */
double Irregular::w2(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2_t = 0;
	//double kmax = *max_element(k.begin(), k.end());
	//double kmin = *min_element(k.begin(), k.end());
	for (int i = 0; i < nfreq * ndir; i++) {
		double wsum2 = 0;
		// Second order
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));		
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		wsum2 += AA* 2*k[i] * sin(2*phi_i) * sinh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			
			double sinhcosh_plus = sinh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);

			// Sum term
			wsum2 += AA * BB * k_nm_plus * sin(phi_i+phi_m) * sinhcosh_plus;	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double sinhcosh_minus = sinh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);				
				// Diff terms
				wsum2 += AA * CC * k_nm_minus * sin(phi_i-phi_m) * sinhcosh_minus;
			}
		}
		wsum2_t += wsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return wsum2_t;
}

/* Second order vertical velocity component W for a sinus wave */
double Irregular::w2_boost(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2_t = 0;
	#pragma omp parallel for reduction(+:wsum2_t)
	for (int i = 0; i < nfreq * ndir; i++) {
		double wsum2 = 0;
		// Second order
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		

		// Diagonal terms
		
		wsum2 += uvw2_diag[i] * sin(2*phi_i) * sinh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); 

		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double sinhcosh_plus = sinh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);
			double sinhcosh_minus = sinh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);	

			// Sum term
			wsum2 += uvw2_sum[i*(nfreq * ndir)+ m] * k_nm_plus * sin(phi_i+phi_m) * sinhcosh_plus;	
				
			// Diff terms
			wsum2 += uvw2_diff[i*(nfreq * ndir)+ m] * k_nm_minus * sin(phi_i-phi_m) * sinhcosh_minus;
			
		}
		wsum2_t += wsum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return wsum2_t;
}


/* Second order dynamic pressure */
double Irregular::dp2(double t, double xx, double yy, double zz) {

	double psum2_t = 0;
	
	for (int i = 0; i < nfreq * ndir; i++) {
		double psum2 = 0;
	
		double thi = theta[i]+ (mtheta * PI / 180.);
		double phi_i = k[i] * (cos(thi) * (xx - fpoint[0]) + sin(thi) * (yy - fpoint[1])) - omega[i] * t + omega[i] * tofmax + phase[i];
		double Rn = k[i] * tanh(k[i] * depth);
		double Rsn = sqrt(Rn);

		// Diagonal terms
		double D_nn_plus = ( 6.*Rsn*Rsn * (k[i] * k[i] - Rn * Rn))/(2.* Rsn* Rsn - k[i] * tanh(2*k[i] * depth));
		double AA = 0.25*(A[i] * g/ omega[i])*(A[i] * g / omega[i])*(D_nn_plus/(2*omega[i]));
		psum2 += AA*(2*omega[i]) * cos(2*phi_i) * cosh(2*k[i]*(zz+depth)) / cosh(2*k[i]*depth); // -d(phi2)/dt, diagonal

		// Off-diagonals
		for (int m = i + 1; m < nfreq * ndir; m++) {
			double thm = theta[m]+ (mtheta * PI / 180.);
			
			double phi_m = k[m] * (cos(thm) * (xx - fpoint[0]) + sin(thm) * (yy - fpoint[1])) - omega[m] * t + omega[m] * tofmax + phase[m];

			double gamma_nm = cos(thi - thm);
			double k_nm_plus = sqrt(k[i] * k[i] + k[m] * k[m] + (2. * k[i] * k[m] * gamma_nm));
			double k_nm_minus = sqrt(k[i] * k[i] + k[m] * k[m] - (2. * k[i] * k[m] * gamma_nm));

			double Rm = k[m] * tanh(k[m] * depth);

			double Rsm = sqrt(Rm);
			
			double Rsnm2plus = pow(Rsn + Rsm, 2.);
			double Rsnm2minus = pow(Rsn - Rsm, 2.);

			double D_nm_plus = (2. * Rsnm2plus * (k[i] * k[m] * gamma_nm - Rn * Rm)) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth)) + 
			((Rsn + Rsm) * (Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i] * k[i] - Rn * Rn))) / (Rsnm2plus - k_nm_plus * tanh(k_nm_plus * depth));

			double D_nm_minus = (2. * Rsnm2minus * (k[i] * k[m] * gamma_nm + Rn * Rm)) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth)) + ((Rsn - Rsm) * (-Rsn*(k[m]*k[m] - Rm * Rm) + Rsm*(k[i]*k[i] - Rn * Rn))) / (Rsnm2minus - k_nm_minus * tanh(k_nm_minus * depth));

			double AA = 0.5*(A[i] * g / omega[i]) * (A[m] * g / omega[m]);
			double BB = (D_nm_plus/(omega[i]+omega[m]));
			double coshcosh_plus = cosh(k_nm_plus*(zz+depth))/cosh(k_nm_plus*depth);

			// Sum terms
			psum2 += AA * BB * (omega[i] + omega[m]) * cos(phi_i + phi_m) * coshcosh_plus; // -d(phi2)/dt, sum
	
			if (k[i] != k[m]){
				double CC = (D_nm_minus/(omega[i]-omega[m])); 
				double coshcosh_minus = cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth);
				
				// Diff terms
				psum2 += AA * CC * (omega[i]-omega[m]) * cos(phi_i-phi_m) * coshcosh_minus;
			}	
		}
		psum2_t += psum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	// Quadratic Bernoulli term -0.5*|grad(phi1)|^2, evaluated exactly from
	// the first-order velocity sums (their per-component ramps included).
	// This replaces earlier per-pair profile products that lacked both the
	// component amplitudes and the phase dependence.
	{
		double u1v = u1(t, xx, yy, zz);
		double v1v = v1(t, xx, yy, zz);
		double w1v = w1(t, xx, yy, zz);
		psum2_t -= 0.5 * (u1v * u1v + v1v * v1v + w1v * w1v);
	}

	// dp2 terms above were assembled without the density factor; dp1 includes rho,
	// so multiply here to keep the second-order dynamic pressure dimensionally
	// consistent (dynamic pressure = rho * [potential + kinetic terms]).
	return rho * psum2_t;
}

void Irregular::rebuild_component_cache()
{
    const int N = nfreq * ndir;
    const double mt = mtheta * PI / 180.0;

    cth_cache.resize(N);
    sth_cache.resize(N);
    R_cache.resize(N);
    Rs_cache.resize(N);
    ag_over_om_cache.resize(N);
    cosh2kd_cache.resize(N);
    phase0_cache.resize(N);

    kmax_cache = 0.0;

    for (int i = 0; i < N; ++i) {
        const double th = theta[i] + mt;

        cth_cache[i] = std::cos(th);
        sth_cache[i] = std::sin(th);

        const double R = k[i] * std::tanh(k[i] * depth);
        R_cache[i] = R;
        Rs_cache[i] = std::sqrt(R);

        ag_over_om_cache[i] = A[i] * g / omega[i];
        cosh2kd_cache[i] = std::cosh(2.0 * k[i] * depth);
        phase0_cache[i] = omega[i] * tofmax + phase[i];

        if (k[i] > kmax_cache) {
            kmax_cache = k[i];
        }
    }
}

void Irregular::rebuild_pair_cache()
{
    const int N = nfreq * ndir;
    pair_cache.clear();
    // Count the surviving pairs first and reserve exactly: a blind
    // N*(N-1)/2 reserve commits gigabytes for large component sets
    // (e.g. dense equal-dk lattices) before any filtering happens.
    size_t npairs = 0;
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) continue;
        for (int m = i + 1; m < N; ++m) {
            if (omega[m] > dw_cutoff) continue;
            if (std::abs(omega[i] - omega[m]) > dw_bandwidth) continue;
            ++npairs;
        }
    }
    if (npairs > (size_t)20000000)
        std::cerr << "WARNING: second-order pair cache has " << npairs
            << " entries (" << N << " components, ~"
            << (npairs * sizeof(PairCacheEntry)) / (1024 * 1024 * 1024)
            << " GB per process). The second-order machinery is O(N^2) in "
            "the component count - reduce the number of components (coarser "
            "lattice / narrower band) or tighten 'bandwidth'." << std::endl;
    pair_cache.reserve(npairs);

    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) {
            continue;
        }

        const double ki = k[i];
        const double oi = omega[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];

        for (int m = i + 1; m < N; ++m) {
            if (omega[m] > dw_cutoff) {
                continue;
            }
            if (std::abs(oi - omega[m]) > dw_bandwidth) {
                continue;
            }

            const double km = k[m];
            const double om = omega[m];
            const double Rnm = R_cache[m];
            const double Rsm = Rs_cache[m];

            const double gamma_nm =
                cth_cache[i] * cth_cache[m] +
                sth_cache[i] * sth_cache[m];

            const double ki2 = ki * ki;
            const double km2 = km * km;
            const double kikm = ki * km;

            const double k_nm_plus =
                std::sqrt(ki2 + km2 + 2.0 * kikm * gamma_nm);

            const double k_nm_minus =
                std::sqrt(ki2 + km2 - 2.0 * kikm * gamma_nm);

            const double sp = Rsni + Rsm;
            const double sm = Rsni - Rsm;
            const double Rsnm2plus = sp * sp;
            const double Rsnm2minus = sm * sm;

            PairCacheEntry entry{};
            entry.i = i;
            entry.m = m;
            entry.gamma_nm = gamma_nm;
            entry.k_nm_plus = k_nm_plus;
            entry.k_nm_minus = k_nm_minus;
            entry.has_plus = false;
            entry.has_minus = false;
            entry.kmax_plus = false;
            entry.D_nm_plus = 0.0;
            entry.D_nm_minus = 0.0;
            entry.alpha_nm_plus = 0.0;
            entry.alpha_nm_minus = 0.0;
            entry.H_nm_plus = 0.0;
            entry.H_nm_minus = 0.0;

            {
                const double denom_plus =
                    Rsnm2plus - k_nm_plus * std::tanh(k_nm_plus * depth);

                const double D_nm_plus =
                    (2.0 * Rsnm2plus * (kikm * gamma_nm - Rni * Rnm)) / denom_plus +
                    (sp * (Rsni * (km2 - Rnm * Rnm) + Rsm * (ki2 - Rni * Rni))) / denom_plus;

                const double alpha_nm_plus =
                    (D_nm_plus - (kikm * gamma_nm - Rni * Rnm)) / std::sqrt(Rni * Rnm)
                    + (Rni + Rnm);

                entry.has_plus = true;
                entry.kmax_plus = ((ki + km) <= kmax_cache);
                entry.D_nm_plus = D_nm_plus;
                entry.alpha_nm_plus = alpha_nm_plus;

                // Smit et al. (2017) H coefficient for sum interaction
                // Derived by inverting Eq. 10: C = (iω₁₂/g)H + ... with C = α/2
                // H = i·H_real, where H_real = -g/ω₁₂ · [α/2 - (ω₁²+ω₂²+ω₁ω₂)/(2g) + g·k₁·k₂/(2ω₁ω₂)]
                {
                    const double om12 = oi + om;
                    const double k1dk2 = kikm * gamma_nm;
                    entry.H_nm_plus = -g / om12 * (
                        alpha_nm_plus / 2.0
                        - (oi * oi + om * om + oi * om) / (2.0 * g)
                        + g * k1dk2 / (2.0 * oi * om));
                }
            }

            if (ki != km) {
                const double denom_minus =
                    Rsnm2minus - k_nm_minus * std::tanh(k_nm_minus * depth);

                const double D_nm_minus =
                    (2.0 * Rsnm2minus * (kikm * gamma_nm + Rni * Rnm)) / denom_minus +
                    (sm * (-Rsni * (km2 - Rnm * Rnm) + Rsm * (ki2 - Rni * Rni))) / denom_minus;

                const double alpha_nm_minus =
                    (D_nm_minus - (kikm * gamma_nm + Rni * Rnm)) / std::sqrt(Rni * Rnm)
                    + (Rni + Rnm);

                entry.has_minus = true;
                entry.D_nm_minus = D_nm_minus;
                entry.alpha_nm_minus = alpha_nm_minus;

                // Smit H coefficient for difference interaction
                // Derived by inverting Eq. 10 with ω₂→−ω₂, k₂→−k₂
                {
                    const double om12 = oi - om;
                    const double k1dk2 = kikm * gamma_nm;
                    if (std::abs(om12) > SMALL) {
                        entry.H_nm_minus = -g / om12 * (
                            alpha_nm_minus / 2.0
                            - (oi * oi + om * om - oi * om) / (2.0 * g)
                            + g * k1dk2 / (2.0 * oi * om));
                    }
                }
            }

            if (entry.has_plus || entry.has_minus) {
                pair_cache.push_back(entry);
            }
        }
    }
}

void Irregular::rebuild_second_order_cache()
{
    rebuild_component_cache();
    rebuild_pair_cache();
}

/* Accumulate wavenumber-bin coefficients for the oscillatory part of the
   second-order kinematics. Mirrors uvw2() term by term:
     u(x,z) = cth0 * sum_j Re{Qk_j E_j(X)} * CR_j(z)
     v(x,z) = sth0 * sum_j Re{Qk_j E_j(X)} * CR_j(z)
     w(x,z) =        sum_j Re{Qw_j E_j(X)} * SR_j(z)
     p(x,z) =        sum_j Re{Qp_j E_j(X)} * CR_j(z)  + O(N) DC terms
   with E_j = exp(i*j*dk*X), CR_j = cosh(j*dk*(z+d))/cosh(j*dk*d) and
   SR_j the sinh analogue. O(N^2) once per time level. */
void Irregular::fact_build_uvw(FactCoeffs& c, double t)
{
    const int N = nfreq * ndir;
    const std::complex<double> I(0.0, 1.0);

    c.Qk.assign(fact_nbins, std::complex<double>(0.0, 0.0));
    c.Qw.assign(fact_nbins, std::complex<double>(0.0, 0.0));
    c.Qp.assign(fact_nbins, std::complex<double>(0.0, 0.0));
    if (c.ramp.empty()) {
        c.ramp.resize(N);
        for (int i = 0; i < N; ++i) {
            c.ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i],
                                  freq_duration[i]);
        }
    }

    std::vector<std::complex<double>> E(N);
    for (int i = 0; i < N; ++i) {
        const double psi = -omega[i] * t + phase0_cache[i];
        E[i] = std::complex<double>(std::cos(psi), std::sin(psi));
    }

    // diagonal terms -> bin 2*n_i (note: no kmax condition, as in uvw2)
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) continue;

        const double ki = k[i];
        const double oi = omega[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];
        const double agi = ag_over_om_cache[i];

        const double D_nn_plus =
            (6.0 * Rsni * Rsni * (ki * ki - Rni * Rni)) /
            (2.0 * Rsni * Rsni - ki * std::tanh(2.0 * ki * depth));
        const double AA_diag = 0.25 * agi * agi * (D_nn_plus / (2.0 * oi));

        const std::complex<double> P = (AA_diag * c.ramp[i]) * E[i] * E[i];
        const int j = 2 * fact_n[i];
        c.Qk[j] += (2.0 * ki) * P;
        c.Qw[j] += -I * (2.0 * ki) * P;
        c.Qp[j] += (2.0 * oi) * P;
    }

    // pair terms
    for (const auto& pair : pair_cache) {
        const int i = pair.i;
        const int m = pair.m;
        const double ki = k[i];
        const double km = k[m];
        const double oi = omega[i];
        const double om = omega[m];
        const double AA = 0.5 * ag_over_om_cache[i] * ag_over_om_cache[m]
                        * c.ramp[i];

        if (pair.has_plus && pair.kmax_plus) {
            const double BB = pair.D_nm_plus / (oi + om);
            const std::complex<double> P = (AA * BB) * E[i] * E[m];
            const int j = fact_n[i] + fact_n[m];
            c.Qk[j] += (ki + km) * P;
            c.Qw[j] += -I * pair.k_nm_plus * P;
            c.Qp[j] += (oi + om) * P;
        }
        if (pair.has_minus) {
            // cos(phi_i - phi_m): bin n_m - n_i with conjugated phasor
            const double CC = pair.D_nm_minus / (oi - om);
            const std::complex<double> P = (AA * CC) * E[m] * std::conj(E[i]);
            const int j = fact_n[m] - fact_n[i];
            c.Qk[j] += (ki - km) * P;
            c.Qw[j] += I * pair.k_nm_minus * P;
            c.Qp[j] += (oi - om) * P;
        }
    }
    c.has_uvw = true;
}

/* Evaluate factorized second-order kinematics at one point.
   O(n_bins) for the oscillatory part + O(N) for the phase-free quadratic
   pressure terms (which are products of first-order profile sums). */
UVWP Irregular::uvw2_factored_eval(const FactCoeffs& c, double xx,
                                   double yy, double zz)
{
    const int N = nfreq * ndir;
    const double zpd = zz + depth;
    const double gg = g * g;
    const double X = cth_cache[0] * (xx - fpoint[0])
                   + sth_cache[0] * (yy - fpoint[1]);

    // per-bin recurrences: phase rotation and stable vertical ratios
    //   CR_j = cosh(a*zpd)/cosh(a*d) = (e^{a(zpd-d)} + e^{-a(zpd+d)}) /
    //                                  (1 + e^{-2ad}),  a = j*dk
    const std::complex<double> rot(std::cos(fact_dk * X),
                                   std::sin(fact_dk * X));
    const double g1 = std::exp(fact_dk * (zpd - depth));
    const double g2 = std::exp(-fact_dk * (zpd + depth));
    const double g3 = std::exp(-2.0 * fact_dk * depth);

    std::complex<double> Ej(1.0, 0.0);
    double e1 = 1.0, e2 = 1.0, e3 = 1.0;

    double u = 0.0, w = 0.0, p = 0.0;
    for (int j = 0; j < fact_nbins; ++j) {
        const double CR = (e1 + e2) / (1.0 + e3);
        const double SR = (e1 - e2) / (1.0 + e3);
        const double reQk = c.Qk[j].real() * Ej.real()
                          - c.Qk[j].imag() * Ej.imag();
        const double reQw = c.Qw[j].real() * Ej.real()
                          - c.Qw[j].imag() * Ej.imag();
        const double reQp = c.Qp[j].real() * Ej.real()
                          - c.Qp[j].imag() * Ej.imag();
        u += reQk * CR;
        w += reQw * SR;
        p += reQp * CR;
        Ej *= rot;
        e1 *= g1; e2 *= g2; e3 *= g3;
    }

    // quadratic Bernoulli term -0.5*|grad(phi1)|^2 from the first-order
    // velocity sums (O(N), exact; per-component ramps included in u1/v1/w1)
    const double cth0 = cth_cache[0];
    const double sth0 = sth_cache[0];
    {
        // the caller's time is the entry's time level
        const double t = c.t;
        const double u1v = u1(t, xx, yy, zz);
        const double v1v = v1(t, xx, yy, zz);
        const double w1v = w1(t, xx, yy, zz);
        p -= 0.5 * (u1v * u1v + v1v * v1v + w1v * w1v);
    }

    return { cth0 * u, sth0 * u, w, p };
}

UVWP Irregular::uvw2(double t, double xx, double yy, double zz)
{
    const int N = nfreq * ndir;
    const double dxp = xx - fpoint[0];
    const double dyp = yy - fpoint[1];
    const double zpd = zz + depth;
    const double gg = g * g;

    // factorized path for long-crested equal-dk grids; cache rebuild is
    // serialized under the mutex (see eta2 for the rationale)
    if (equal_dk) {
        if (!fact_checked.load()) {
            std::lock_guard<std::mutex> lock(fact_mutex);
            if (!fact_checked.load()) {
                if ((int)cth_cache.size() != N) {
                    rebuild_second_order_cache();
                }
                factored_init();
                fact_checked.store(true);   // release: init fully visible
            }
        }
        if (fact_ok) {
            auto c = fact_get(t, true);
            return uvw2_factored_eval(*c, xx, yy, zz);
        }
    }

    if ((int)cth_cache.size() != N) {
        rebuild_second_order_cache();
    }

    double usum2_t = 0.0;
    double vsum2_t = 0.0;
    double wsum2_t = 0.0;
    double psum2_t = 0.0;

    std::vector<double> cphi(N), sphi(N), ramp(N);

    for (int i = 0; i < N; ++i) {
        const double ph =
            k[i] * (cth_cache[i] * dxp + sth_cache[i] * dyp)
            - omega[i] * t
            + phase0_cache[i];

        cphi[i] = std::cos(ph);
        sphi[i] = std::sin(ph);
        ramp[i] = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
    }

    // diagonal terms
    for (int i = 0; i < N; ++i) {
        if (omega[i] > dw_cutoff) {
            continue;
        }

        const double ki = k[i];
        const double oi = omega[i];
        const double cti = cth_cache[i];
        const double sti = sth_cache[i];
        const double Rni = R_cache[i];
        const double Rsni = Rs_cache[i];
        const double agi = ag_over_om_cache[i];

        const double D_nn_plus =
            (6.0 * Rsni * Rsni * (ki * ki - Rni * Rni)) /
            (2.0 * Rsni * Rsni - ki * std::tanh(2.0 * ki * depth));

        const double AA_diag = 0.25 * agi * agi * (D_nn_plus / (2.0 * oi));

        const double cos_2phi_i = cphi[i] * cphi[i] - sphi[i] * sphi[i];
        const double sin_2phi_i = 2.0 * sphi[i] * cphi[i];

        const double coshcosh_diag = std::cosh(2.0 * ki * zpd) / cosh2kd_cache[i];
        const double sinhcosh_diag = std::sinh(2.0 * ki * zpd) / cosh2kd_cache[i];

        double u = 0.0, v = 0.0, w = 0.0, p = 0.0;

        u += AA_diag * (2.0 * ki * cti) * cos_2phi_i * coshcosh_diag;
        v += AA_diag * (2.0 * ki * sti) * cos_2phi_i * coshcosh_diag;
        w += AA_diag * (2.0 * ki) * sin_2phi_i * sinhcosh_diag;
        p += AA_diag * (2.0 * oi) * cos_2phi_i * coshcosh_diag;

        const double ri = ramp[i];
        usum2_t += u * ri;
        vsum2_t += v * ri;
        wsum2_t += w * ri;
        psum2_t += p * ri;
    }

    // pair terms
    for (const auto& pair : pair_cache) {
        const int i = pair.i;
        const int m = pair.m;

        const double ki = k[i];
        const double km = k[m];
        const double oi = omega[i];
        const double om = omega[m];
        const double cti = cth_cache[i];
        const double sti = sth_cache[i];
        const double ctm = cth_cache[m];
        const double stm = sth_cache[m];
        const double ri = ramp[i];

        const double AA = 0.5 * ag_over_om_cache[i] * ag_over_om_cache[m];

        if (pair.has_plus && pair.kmax_plus) {
            const double BB = pair.D_nm_plus / (oi + om);
            const double cos_phi_plus = cphi[i] * cphi[m] - sphi[i] * sphi[m];
            const double sin_phi_plus = sphi[i] * cphi[m] + cphi[i] * sphi[m];

            const double denom = std::cosh(pair.k_nm_plus * depth);
            const double coshcosh_plus = std::cosh(pair.k_nm_plus * zpd) / denom;
            const double sinhcosh_plus = std::sinh(pair.k_nm_plus * zpd) / denom;

            usum2_t += AA * BB * (ki * cti + km * ctm) * cos_phi_plus * coshcosh_plus * ri;
            vsum2_t += AA * BB * (ki * sti + km * stm) * cos_phi_plus * coshcosh_plus * ri;
            wsum2_t += AA * BB * pair.k_nm_plus * sin_phi_plus * sinhcosh_plus * ri;
            psum2_t += AA * BB * (oi + om) * cos_phi_plus * coshcosh_plus * ri;
        }

        if (pair.has_minus) {
            const double CC = pair.D_nm_minus / (oi - om);
            const double cos_phi_minus = cphi[i] * cphi[m] + sphi[i] * sphi[m];
            const double sin_phi_minus = sphi[i] * cphi[m] - cphi[i] * sphi[m];

            const double denom = std::cosh(pair.k_nm_minus * depth);
            const double coshcosh_minus = std::cosh(pair.k_nm_minus * zpd) / denom;
            const double sinhcosh_minus = std::sinh(pair.k_nm_minus * zpd) / denom;

            usum2_t += AA * CC * (ki * cti - km * ctm) * cos_phi_minus * coshcosh_minus * ri;
            vsum2_t += AA * CC * (ki * sti - km * stm) * cos_phi_minus * coshcosh_minus * ri;
            wsum2_t += AA * CC * pair.k_nm_minus * sin_phi_minus * sinhcosh_minus * ri;
            psum2_t += AA * CC * (oi - om) * cos_phi_minus * coshcosh_minus * ri;
        }
    }

    // quadratic Bernoulli term -0.5*|grad(phi1)|^2 from the first-order
    // velocity sums (per-component ramps included in u1/v1/w1)
    {
        const double u1v = u1(t, xx, yy, zz);
        const double v1v = v1(t, xx, yy, zz);
        const double w1v = w1(t, xx, yy, zz);
        psum2_t -= 0.5 * (u1v * u1v + v1v * v1v + w1v * w1v);
    }

    return {usum2_t, vsum2_t, wsum2_t, psum2_t};
}




double Irregular::dp1(double t, double xx, double yy, double zz) {

	double psum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		psum += A[i] * rho * g * (cosh(k[i] * (zz + depth)) / cosh(k[i] * depth)) 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
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
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return usum;

}

/* vertical velocity gradient at z=0 for velocity component U */
double Irregular::dp1_dz(double t, double xx, double yy) {

	double psum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		psum +=  A[i] * rho * g * k[i] 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return psum;

}


/* vertical velocity gradient at z=0 for velocity component V */
double Irregular::phi1_dydz(double t, double xx, double yy) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i < ndir * nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		vsum += sin(theta[i] + (mtheta * PI / 180.)) * A[i] * omega[i] * k[i] 
			* cos(k[i] * (cos(theta[i] + (mtheta * PI / 180.)) * (xx - fpoint[0]) + sin(theta[i] 
				+ (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
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
				+ sin(theta[i] + (mtheta * PI / 180.)) * (yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return wsum;

}

// First order velocity potential
double Irregular::phi1_pot(double t, double xx, double yy, double zz) {

	double phisum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = omega[i] * tofmax + phase[i];
		phisum += A[i] * (omega[i] / k[i]) * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*sin(k[i] 
		* (cos(theta[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(theta[i] + (mtheta*PI / 180.))*
		(yy - fpoint[1])) - omega[i] * t + phi)*
				freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}

	return phisum;

}


/* Second order component of velocity potential
Warning: Unsure if this is finished...should be checked thoroughly at some point*/
double Irregular::phi2_pot(double t, double xx, double yy, double zz) {
	//double eta1_t = 0;
	double phisum2_t = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		double phisum2 = 0;
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
		phisum2_t += phisum2*freq_ramp(t,freq_ramptime[i],freq_starttime[i], freq_duration[i]);
	}
	return phisum2_t;
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
			phase[i] = 0.;
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

std::vector<double>  Irregular::spectral_moments(){
	// Computes the spectral moments (m0, m1, m2) from spectral compnents
	if (ndir > 1) {
		std::cerr << "WARNING: mean wave length only supported for 1D wave spectra at the moment. Sorry for the inconvenience" << std::endl;
	}

	double w_t, dw, S;
	double m0 = 0.;
	double m1 = 0.;
	double m2 = 0.;

	// not very robust
	// for (int i = 0; i < (nfreq - 1); i++) {
	// 	w_t = (omega[i + 1] + omega[i])/2.;
	// 	dw = omega[i + 1] - omega[i];
	// 	S = pow((A[i + 1] + A[i]) * 0.5, 2.) / (2. * dw);
	// 	m0 += dw * S;
	// 	m1 += dw * S * w_t;
	// 	m2 += dw * S * w_t * w_t;
	// }

	for (int i = 0; i < (nfreq - 1); i++) {
		m0 += A[i]*A[i]/2.;
		m1 += omega[i]*A[i]*A[i]/2.;
		m2 += omega[i]*omega[i]*A[i]*A[i]/2.;
	}

	std::vector<double> res;
	res.push_back(m0);
	res.push_back(m1);
	res.push_back(m2);

	return res;
}

/* Function for calculating phase speed based integration of the wave spectrum*/
double Irregular::phase_speed(int opt)
{
	std::vector<double> moment = spectral_moments();
	
	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);

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
	std::vector<double> moment = spectral_moments();
	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);

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
	std::vector<double> moment = spectral_moments();

	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);

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
	std::vector<double> moment = spectral_moments();

	double bw = 0.7 * ( moment[1] / moment[0]);
	
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

void Irregular::print_spectral_parameters() {
	// Print some spectral parameters which can be used to validate the input spectral components
	// Spectral moments 0, 1, 2
	// Tm, Tz, Hs,
	if (A.empty() || omega.empty()) {
		std::cerr << "WARNING: no spectral components available -- skipping spectral parameter report. "
		             "Check that the [wave input irregular] section was parsed (correct header and frequency count)." << std::endl;
		return;
	}
	std::vector<double> moment = spectral_moments();
	double t1 = 2. * PI * (moment[0] / moment[1]);
	double t2 = 2. * PI * sqrt(moment[0] / moment[2]);
	double hs = 4*sqrt(moment[0]);
	std::vector<double>::iterator it = std::max_element(A.begin(), A.end());
    int index = static_cast<int>(std::distance(A.begin(), it));
	double tp = 2*PI/omega[index];
	

	std::cout << "----------------------------------------------------------------" << std::endl;
	std::cout << "Spectral parameters computed from specified spectral components:" << std::endl << std::endl;
	std::cout << "Spectral moments m0: " << moment[0] << ", m1: " << moment[1] << ", m2: " << moment[2] << std::endl;
	std::cout << "Significant Wave height Hs (m):" << hs << std::endl;
	std::cout << "Mean period Tz (sec, from m0 and m2):" << t2 << std::endl;
	std::cout << "Mean period Tm:(sec, from m0 and m1)" << t1 << std::endl;
	std::cout << "Spectral peak period (sec):" << tp << std::endl;
	std::cout << "Mean phase speed (m/s, based on Tm):" << phase_speed(1) << std::endl;
	std::cout << "----------------------------------------------------------------" << std::endl;



}


/* ==========================================================================
   Smit et al. (2017) s-coordinate second-order kinematics.

   The velocity field is (Eq. 11):
     u = Σ η̂₁ U₁¹ e₁  +  Σ η̂₁η̂₂ U¹² e₁e₂

   First-order uses Wheeler-stretched profiles: Ch^z = cosh[k·h·(h+z)/(h+η)] / cosh(kh)
   Second-order has the bound wave (H coefficient + stretched Ch_{1+2}) plus
   a coordinate-transformation correction term g(h+z)/(2(h+η)) · (...)
   ========================================================================== */

/* Second-order surface elevation for the smit s-coordinate path.

   Historical note: this function originally implemented Smit's Eq. 9 for
   H (with typo corrections) as an independent cross-check route. That
   investigation (sigma_developments/) established that Eq. 9 is
   structurally wrong - the correct H comes from inverting Eq. 10 with
   C = alpha/2 (Sharma & Dean), which is what the kinematics use. The
   correct sigma-consistent surface is therefore simply the S&D
   second-order elevation WITHOUT the kmax truncation (eta2_full), the
   same surface the smit kinematics stretch and clamp against. */
double Irregular::eta_smit(double t, double xx, double yy)
{
	return eta2_full(t, xx, yy);
}


double Irregular::u_smit(double t, double xx, double yy, double zz) {

	// Step 1: compute surface elevation for the stretch factor
	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2_full(t, xx, yy);

	// Step 2: clamp to water column
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	// Precompute phases and first-order stretched profiles
	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// Step 3: first-order velocity with stretched profiles (Wheeler stretching)
	double u1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		u1_val += A[i] * g * k[i] / omega[i] * profileX_s(i, zz, eta_val)
			* cth_cache[i] * cos(phi) * fr;
	}

	// Step 4: second-order contributions (skip if order < 1, giving pure Wheeler stretching)
	double u2_val = 0.0;
	if (order >= 1) {

	// Accumulate ∂_s φ² so the Eulerian chain term can map the FULL potential
	// (φ¹+φ²). The φ² part is the necessary O(ε³) term: −(1+s)/D·η_x·∂_s φ².
	double dphi2ds = 0.0;

	// Diagonal terms (self-interaction: component i with itself)
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;

		double k2 = 2.0 * ki; // combined wavenumber magnitude
		double om2 = 2.0 * oi; // combined frequency
		// Compute diagonal alpha (Sharma & Dean) then invert Eq. 10 for H
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));

		// Bound wave: -kx₁₂ · H · Ch(2k) · cos(2φ)  [vertically-Lagrangian: ∂_x* at constant s]
		double kx2 = 2.0 * ki * cth_cache[i];
		double Ch2 = profileX_s2(k2, zz, eta_val);
		u2_val += 0.5 * AA * H_diag * kx2 * (-cos(2.0 * phi_i)) * Ch2;

		// P-derivative ∂ₓP (particular-solution velocity) — coefficient validated
		// against exact Stokes at all depths (papers/sigma_to_euler_derivation.py).
		double sp1 = (depth + zz) / (depth + eta_val);              // (1+s)
		double g1s = g * sp1;                                       // g(1+s)
		double Sh_i = profileZ_s(i, zz, eta_val);
		u2_val += 0.5 * AA * kx2 * g1s * (ki * Sh_i / oi) * cos(2.0 * phi_i);

		// ∂_s φ² (diagonal): bound  φ=−0.5·AA·H·Ch2·sin2φ,  ∂_sCh2 = k2·h·Sh2
		//                    P-pot  φ= 0.5·AA·g(1+s)(k·Sh/ω)·sin2φ
		double Ch_i = profileX_s(i, zz, eta_val);
		double Sh2  = profileZ_s2(k2, zz, eta_val);
		dphi2ds += -0.5 * AA * H_diag * (k2 * depth * Sh2) * sin(2.0 * phi_i);
		dphi2ds +=  0.5 * AA * g * (ki * Sh_i / oi
			+ sp1 * ki * (ki * depth * Ch_i) / oi) * sin(2.0 * phi_i);
	}

	// Off-diagonal terms (pair interactions)
	for (const auto& p : pair_cache) {
		const int i = p.i;
		const int m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		double sp1 = (depth + zz) / (depth + eta_val);                // (1+s)
		double g1s = g * sp1;                                         // g(1+s)
		double Sh_i = profileZ_s(i, zz, eta_val);
		double Sh_m = profileZ_s(m, zz, eta_val);
		double Ch_i = profileX_s(i, zz, eta_val);
		double Ch_m = profileX_s(m, zz, eta_val);

		// Sum interaction (ω₁ + ω₂): bound wave + P-derivative ∂ₓP
		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double kx12 = ki * cth_cache[i] + km * cth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_plus, zz, eta_val);
			u2_val += AA * p.H_nm_plus * kx12 * (-cos(phi_plus)) * Ch12;
			// P-derivative ∂ₓP: factor 0.5 like the diagonal (the paper's P is
			// already the total per unordered pair; eps-scaling budget gives
			// xP = 1/2 exactly, see sigma_developments/FINDINGS_eps_scaling.md)
			u2_val += 0.5 * AA * kx12 * g1s * (ki * Sh_i / oi + km * Sh_m / om) * cos(phi_plus);

			// ∂_s φ² (sum): bound  φ=−AA·H·Ch12·sinφ₊,  P-pot φ=0.5·AA·g(1+s)(...)·sinφ₊
			double Sh12 = profileZ_s2(p.k_nm_plus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_plus * (p.k_nm_plus * depth * Sh12) * sin(phi_plus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi + km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi + km * (km * depth * Ch_m) / om)) * sin(phi_plus);
		}

		// Difference interaction (ω₁ − ω₂): k₂→−k₂, ω₂→−ω₂
		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double kx12 = ki * cth_cache[i] - km * cth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_minus, zz, eta_val);
			u2_val += AA * p.H_nm_minus * kx12 * (-cos(phi_minus)) * Ch12;
			// P-derivative ∂ₓP with the same 0.5 factor as the sum interaction
			u2_val += 0.5 * AA * kx12 * g1s * (ki * Sh_i / oi - km * Sh_m / om) * cos(phi_minus);

			// ∂_s φ² (difference): k₂→−k₂, ω₂→−ω₂
			double Sh12 = profileZ_s2(p.k_nm_minus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_minus * (p.k_nm_minus * depth * Sh12) * sin(phi_minus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi - km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi - km * (km * depth * Ch_m) / om)) * sin(phi_minus);
		}
	}

	// Lagrangian→Eulerian chain term:  u += -(1+s)/D · η_x · φ_s
	// (field product of first-order surface slope η_x and ∂φ₁/∂s — automatically
	//  generates all sum/difference cross terms. Derived & verified symbolically in
	//  papers/sigma_to_euler_derivation.py: for one component this equals u1_quad,
	//  the missing half of the Eulerian 2nd-order velocity from φ₁.)
	{
		double eta_x = 0.0, phi_s = 0.0;
		for (int n = 0; n < N; n++) {
			double phi_n = k[n] * (cth_cache[n] * (xx - fpoint[0]) + sth_cache[n] * (yy - fpoint[1]))
				- omega[n] * t + phase0_cache[n];
			double fr = freq_ramp(t, freq_ramptime[n], freq_starttime[n], freq_duration[n]);
			double sphi = sin(phi_n);
			eta_x += -A[n] * k[n] * cth_cache[n] * sphi * fr;                       // ∂η₁/∂x
			phi_s += (g * A[n] / omega[n]) * k[n] * depth
				* profileZ_s(n, zz, eta_val) * sphi * fr;                           // ∂φ₁/∂s
		}
		phi_s += dphi2ds;   // include ∂_s φ²: chain maps the FULL potential (Eulerian, O(ε³) consistent)
		double oneplus_s_over_D = (depth + zz) / ((depth + eta_val) * (depth + eta_val));
		u2_val += -oneplus_s_over_D * eta_x * phi_s;
	}
	} // end if (order >= 1)

	return u1_val + u2_val;
}


double Irregular::v_smit(double t, double xx, double yy, double zz) {

	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2_full(t, xx, yy);
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// First order with stretched profiles
	double v1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		v1_val += A[i] * g * k[i] / omega[i] * profileX_s(i, zz, eta_val)
			* sth_cache[i] * cos(phi) * fr;
	}

	// Second order — mirrors u_smit with y-components (∂φ/∂y mapping:
	// v = ∂ᵧφ̃ − (1+s)/D·η_y·∂ₛφ̃)
	double v2_val = 0.0;
	if (order >= 1) {

	double dphi2ds = 0.0;

	// Diagonal terms
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;
		double k2 = 2.0 * ki, om2 = 2.0 * oi;
		// Diagonal H from inverting Eq. 10 with Sharma & Dean alpha (NOT the buggy Eq. 9)
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));
		double ky2 = 2.0 * ki * sth_cache[i];
		double Ch2 = profileX_s2(k2, zz, eta_val);
		v2_val += 0.5 * AA * H_diag * ky2 * (-cos(2.0 * phi_i)) * Ch2;

		// P-derivative ∂ᵧP (same 0.5 convention as u_smit)
		double sp1 = (depth + zz) / (depth + eta_val);
		double g1s = g * sp1;
		double Sh_i = profileZ_s(i, zz, eta_val);
		v2_val += 0.5 * AA * ky2 * g1s * (ki * Sh_i / oi) * cos(2.0 * phi_i);

		double Ch_i = profileX_s(i, zz, eta_val);
		double Sh2  = profileZ_s2(k2, zz, eta_val);
		dphi2ds += -0.5 * AA * H_diag * (k2 * depth * Sh2) * sin(2.0 * phi_i);
		dphi2ds +=  0.5 * AA * g * (ki * Sh_i / oi
			+ sp1 * ki * (ki * depth * Ch_i) / oi) * sin(2.0 * phi_i);
	}

	// Off-diagonal terms
	for (const auto& p : pair_cache) {
		const int i = p.i, m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		double sp1 = (depth + zz) / (depth + eta_val);
		double g1s = g * sp1;
		double Sh_i = profileZ_s(i, zz, eta_val);
		double Sh_m = profileZ_s(m, zz, eta_val);
		double Ch_i = profileX_s(i, zz, eta_val);
		double Ch_m = profileX_s(m, zz, eta_val);

		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double ky12 = ki * sth_cache[i] + km * sth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_plus, zz, eta_val);
			v2_val += AA * p.H_nm_plus * ky12 * (-cos(phi_plus)) * Ch12;
			v2_val += 0.5 * AA * ky12 * g1s * (ki * Sh_i / oi + km * Sh_m / om) * cos(phi_plus);

			double Sh12 = profileZ_s2(p.k_nm_plus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_plus * (p.k_nm_plus * depth * Sh12) * sin(phi_plus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi + km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi + km * (km * depth * Ch_m) / om)) * sin(phi_plus);
		}

		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double ky12 = ki * sth_cache[i] - km * sth_cache[m];
			double Ch12 = profileX_s2(p.k_nm_minus, zz, eta_val);
			v2_val += AA * p.H_nm_minus * ky12 * (-cos(phi_minus)) * Ch12;
			v2_val += 0.5 * AA * ky12 * g1s * (ki * Sh_i / oi - km * Sh_m / om) * cos(phi_minus);

			double Sh12 = profileZ_s2(p.k_nm_minus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_minus * (p.k_nm_minus * depth * Sh12) * sin(phi_minus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi - km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi - km * (km * depth * Ch_m) / om)) * sin(phi_minus);
		}
	}

	// Lagrangian→Eulerian chain term:  v += -(1+s)/D · η_y · φ_s
	{
		double eta_y = 0.0, phi_s = 0.0;
		for (int n = 0; n < N; n++) {
			double phi_n = k[n] * (cth_cache[n] * (xx - fpoint[0]) + sth_cache[n] * (yy - fpoint[1]))
				- omega[n] * t + phase0_cache[n];
			double fr = freq_ramp(t, freq_ramptime[n], freq_starttime[n], freq_duration[n]);
			double sphi = sin(phi_n);
			eta_y += -A[n] * k[n] * sth_cache[n] * sphi * fr;                   // ∂η₁/∂y
			phi_s += (g * A[n] / omega[n]) * k[n] * depth
				* profileZ_s(n, zz, eta_val) * sphi * fr;                       // ∂φ₁/∂s
		}
		phi_s += dphi2ds;
		double oneplus_s_over_D = (depth + zz) / ((depth + eta_val) * (depth + eta_val));
		v2_val += -oneplus_s_over_D * eta_y * phi_s;
	}
	} // end if (order >= 1)

	return v1_val + v2_val;
}


double Irregular::w_smit(double t, double xx, double yy, double zz) {

	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2_full(t, xx, yy);
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// Exact vertical mapping:  w = ∂φ/∂z|ₓ = ∂ₛφ̃ / (h+η)  with φ̃ = φ¹+φ²
	// in σ-coordinates. The 1/(h+η) Jacobian and the ∂ₛ of the particular
	// solution P are required for Eulerian consistency (eps-scaling budget,
	// sigma_developments/FINDINGS_eps_scaling.md).

	// ∂φ¹/∂s
	double phi_s = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		phi_s += (g * A[i] / omega[i]) * k[i] * depth
			* profileZ_s(i, zz, eta_val) * sin(phi) * fr;
	}

	// ∂φ²/∂s (bound + P parts, diagonal and pairs) — same accumulation as u_smit
	double dphi2ds = 0.0;
	if (order >= 1) {

	double sp1 = (depth + zz) / (depth + eta_val);              // (1+s)

	// Diagonal terms
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;
		double k2 = 2.0 * ki, om2 = 2.0 * oi;
		// Diagonal H from inverting Eq. 10 with Sharma & Dean alpha (NOT the buggy Eq. 9)
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));

		double Ch_i = profileX_s(i, zz, eta_val);
		double Sh_i = profileZ_s(i, zz, eta_val);
		double Sh2  = profileZ_s2(k2, zz, eta_val);
		dphi2ds += -0.5 * AA * H_diag * (k2 * depth * Sh2) * sin(2.0 * phi_i);
		dphi2ds +=  0.5 * AA * g * (ki * Sh_i / oi
			+ sp1 * ki * (ki * depth * Ch_i) / oi) * sin(2.0 * phi_i);
	}

	// Off-diagonal terms
	for (const auto& p : pair_cache) {
		const int i = p.i, m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		double Sh_i = profileZ_s(i, zz, eta_val);
		double Sh_m = profileZ_s(m, zz, eta_val);
		double Ch_i = profileX_s(i, zz, eta_val);
		double Ch_m = profileX_s(m, zz, eta_val);

		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double Sh12 = profileZ_s2(p.k_nm_plus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_plus * (p.k_nm_plus * depth * Sh12) * sin(phi_plus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi + km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi + km * (km * depth * Ch_m) / om)) * sin(phi_plus);
		}

		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double Sh12 = profileZ_s2(p.k_nm_minus, zz, eta_val);
			dphi2ds += -AA * p.H_nm_minus * (p.k_nm_minus * depth * Sh12) * sin(phi_minus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi - km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi - km * (km * depth * Ch_m) / om)) * sin(phi_minus);
		}
	}
	} // end if (order >= 1)

	return (phi_s + dphi2ds) / (depth + eta_val);
}


// First-order dynamic pressure, smit stretched profile: the dp1 block of
// dp_smit() with the stretching referenced to the FIRST-order elevation
// (consistent truncation, and eta2_full would cost the O(N^2) pair sums
// this fast path exists to avoid).
double Irregular::dp1_smit(double t, double xx, double yy, double zz) {

	double eta_val = eta1(t, xx, yy);
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	double dp1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		dp1_val += A[i] * rho * g * profileX_s(i, zz, eta_val) * cos(phi) * fr;
	}
	return dp1_val;
}

double Irregular::dp_smit(double t, double xx, double yy, double zz) {

	double eta_val = eta1(t, xx, yy);
	if (order >= 1) eta_val += eta2_full(t, xx, yy);
	if (zz > eta_val + swl) return 0.0;
	zz = std::max(-depth, zz - swl);

	const int N = ndir * nfreq;
	if ((int)cth_cache.size() != N) rebuild_second_order_cache();

	// First order pressure with stretched profiles
	double dp1_val = 0.0;
	for (int i = 0; i < N; i++) {
		double phi = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		dp1_val += A[i] * rho * g * profileX_s(i, zz, eta_val) * cos(phi) * fr;
	}

	// Second order pressure via the exact σ-mapping
	//   p = −ρ[∂ₜφ̃ − (1+s)/D·η_t·∂ₛφ̃] − ½ρ|∇φ¹|²
	// (∂ₜ at fixed z picks up the s_t = −(1+s)/D·η_t Jacobian term; the
	// Bernoulli part is the exact field product of the first-order
	// stretched velocities, matching the corrected S&D dp2 formulation.)
	double dp2_val = 0.0;
	if (order >= 1) {

	double sp1 = (depth + zz) / (depth + eta_val);              // (1+s)
	double g1s = g * sp1;                                       // g(1+s)

	double dphi2dt = 0.0;   // ∂φ²/∂t (σ-frame, at fixed s)
	double dphi2ds = 0.0;   // ∂φ²/∂s

	// Diagonal terms
	for (int i = 0; i < N; i++) {
		double phi_i = k[i] * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- omega[i] * t + phase0_cache[i];
		double fr = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double ki = k[i], oi = omega[i];
		double AA = A[i] * A[i] * fr * fr;
		double k2 = 2.0 * ki, om2 = 2.0 * oi;
		// Diagonal H from inverting Eq. 10 with Sharma & Dean alpha (NOT the buggy Eq. 9)
		double Rni = R_cache[i], Rsni = Rs_cache[i];
		double D_nn = (6.0*Rsni*Rsni*(ki*ki - Rni*Rni)) / (2.0*Rsni*Rsni - ki*std::tanh(2.0*ki*depth));
		double alpha_nn = (D_nn - (ki*ki - Rni*Rni)) / Rni + 2.0*Rni;
		double H_diag = -g / om2 * (alpha_nn / 2.0
			- (3.0 * oi * oi) / (2.0 * g)
			+ g * ki * ki / (2.0 * oi * oi));

		double Ch2 = profileX_s2(k2, zz, eta_val);
		double Ch_i = profileX_s(i, zz, eta_val);
		double Sh_i = profileZ_s(i, zz, eta_val);
		double Sh2  = profileZ_s2(k2, zz, eta_val);

		// φ²(diag) = −0.5·AA·H·Ch₂·sin2φ + 0.5·AA·g(1+s)(k·Sh/ω)·sin2φ
		dphi2dt += 0.5 * AA * H_diag * om2 * Ch2 * cos(2.0 * phi_i);
		dphi2dt += -0.5 * AA * g1s * (ki * Sh_i / oi) * om2 * cos(2.0 * phi_i);
		dphi2ds += -0.5 * AA * H_diag * (k2 * depth * Sh2) * sin(2.0 * phi_i);
		dphi2ds +=  0.5 * AA * g * (ki * Sh_i / oi
			+ sp1 * ki * (ki * depth * Ch_i) / oi) * sin(2.0 * phi_i);
	}

	// Off-diagonal terms
	for (const auto& p : pair_cache) {
		const int i = p.i, m = p.m;
		const double Ai = A[i], Am = A[m];
		const double ki = k[i], km = k[m];
		const double oi = omega[i], om = omega[m];

		double phi_i = ki * (cth_cache[i] * (xx - fpoint[0]) + sth_cache[i] * (yy - fpoint[1]))
			- oi * t + phase0_cache[i];
		double phi_m = km * (cth_cache[m] * (xx - fpoint[0]) + sth_cache[m] * (yy - fpoint[1]))
			- om * t + phase0_cache[m];
		double fr_i = freq_ramp(t, freq_ramptime[i], freq_starttime[i], freq_duration[i]);
		double fr_m = freq_ramp(t, freq_ramptime[m], freq_starttime[m], freq_duration[m]);
		double AA = Ai * Am * fr_i * fr_m;

		double Ch_i = profileX_s(i, zz, eta_val);
		double Ch_m = profileX_s(m, zz, eta_val);
		double Sh_i = profileZ_s(i, zz, eta_val);
		double Sh_m = profileZ_s(m, zz, eta_val);

		if (p.has_plus) {
			double phi_plus = phi_i + phi_m;
			double om12 = oi + om;
			double Ch12 = profileX_s2(p.k_nm_plus, zz, eta_val);
			double Sh12 = profileZ_s2(p.k_nm_plus, zz, eta_val);

			// φ²(sum) = −AA·H·Ch₁₂·sinφ₊ + 0.5·AA·g(1+s)(k₁Sh₁/ω₁+k₂Sh₂/ω₂)·sinφ₊
			dphi2dt += AA * p.H_nm_plus * om12 * Ch12 * cos(phi_plus);
			dphi2dt += -0.5 * AA * g1s * (ki * Sh_i / oi + km * Sh_m / om) * om12 * cos(phi_plus);
			dphi2ds += -AA * p.H_nm_plus * (p.k_nm_plus * depth * Sh12) * sin(phi_plus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi + km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi + km * (km * depth * Ch_m) / om)) * sin(phi_plus);
		}

		if (p.has_minus) {
			double phi_minus = phi_i - phi_m;
			double om12 = oi - om;
			double Ch12 = profileX_s2(p.k_nm_minus, zz, eta_val);
			double Sh12 = profileZ_s2(p.k_nm_minus, zz, eta_val);

			// difference: k₂→−k₂, ω₂→−ω₂
			dphi2dt += AA * p.H_nm_minus * om12 * Ch12 * cos(phi_minus);
			dphi2dt += -0.5 * AA * g1s * (ki * Sh_i / oi - km * Sh_m / om) * om12 * cos(phi_minus);
			dphi2ds += -AA * p.H_nm_minus * (p.k_nm_minus * depth * Sh12) * sin(phi_minus);
			dphi2ds +=  0.5 * AA * g * ((ki * Sh_i / oi - km * Sh_m / om)
				+ sp1 * (ki * (ki * depth * Ch_i) / oi - km * (km * depth * Ch_m) / om)) * sin(phi_minus);
		}
	}

	// First-order η_t, ∂ₛφ¹ and velocities (stretched profiles)
	double eta_t = 0.0, phi_s1 = 0.0, u1v = 0.0, v1v = 0.0;
	for (int n = 0; n < N; n++) {
		double phi_n = k[n] * (cth_cache[n] * (xx - fpoint[0]) + sth_cache[n] * (yy - fpoint[1]))
			- omega[n] * t + phase0_cache[n];
		double fr = freq_ramp(t, freq_ramptime[n], freq_starttime[n], freq_duration[n]);
		double cphi = cos(phi_n), sphi = sin(phi_n);
		eta_t += A[n] * omega[n] * sphi * fr;                              // ∂η₁/∂t
		phi_s1 += (g * A[n] / omega[n]) * k[n] * depth
			* profileZ_s(n, zz, eta_val) * sphi * fr;                      // ∂φ₁/∂s
		double uampl = A[n] * g * k[n] / omega[n] * profileX_s(n, zz, eta_val) * cphi * fr;
		u1v += uampl * cth_cache[n];
		v1v += uampl * sth_cache[n];
	}
	double w1v = phi_s1 / (depth + eta_val);

	double oneplus_s_over_D = (depth + zz) / ((depth + eta_val) * (depth + eta_val));
	dp2_val = rho * (-dphi2dt
		+ oneplus_s_over_D * eta_t * (phi_s1 + dphi2ds)
		- 0.5 * (u1v * u1v + v1v * v1v + w1v * w1v));
	} // end if (order >= 1)

	return dp1_val + dp2_val;
}

