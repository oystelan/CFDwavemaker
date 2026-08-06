#include "FentonStream.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Dense linear solve  A x = b  by Gaussian elimination with partial pivoting.
// 1-indexed (A[1..n][1..n], b[1..n]); replaces b with the solution. Returns
// false if the matrix is singular. (Stand-in for LINPACK dgefa/dgesl.)
// ---------------------------------------------------------------------------
bool FentonStream::gauss_solve(std::vector<std::vector<double>>& A, std::vector<double>& b, int n)
{
	for (int col = 1; col <= n; ++col) {
		// pivot
		int piv = col;
		double best = std::fabs(A[col][col]);
		for (int r = col + 1; r <= n; ++r) {
			if (std::fabs(A[r][col]) > best) { best = std::fabs(A[r][col]); piv = r; }
		}
		if (best < 1e-300) return false;
		if (piv != col) { std::swap(A[piv], A[col]); std::swap(b[piv], b[col]); }

		const double dinv = 1.0 / A[col][col];
		for (int r = col + 1; r <= n; ++r) {
			double f = A[r][col] * dinv;
			if (f != 0.0) {
				for (int cc = col; cc <= n; ++cc) A[r][cc] -= f * A[col][cc];
				b[r] -= f * b[col];
			}
		}
	}
	// back substitution
	for (int r = n; r >= 1; --r) {
		double s = b[r];
		for (int cc = r + 1; cc <= n; ++cc) s -= A[r][cc] * b[cc];
		b[r] = s / A[r][r];
	}
	return true;
}

// ---------------------------------------------------------------------------
// Initial linear solution (Fenton 1988 subroutine INIT, Table 1).
// z is 1-indexed, length 2N+11 (index 0 unused).
// ---------------------------------------------------------------------------
void FentonStream::init(std::vector<double>& z)
{
	const int n = N;
	if (finite) {
		if (!case_wavelength) {           // period specified: Eckart + one Newton step -> kd
			double a = 4.0 * PI * PI * height / hoverd;
			double b = a / std::sqrt(std::tanh(a));
			double t = std::tanh(b);
			z[1] = b + (a - b * t) / (t + b * (1.0 - t * t));
		} else {                          // wavelength specified: kd = 2 pi d / lambda
			z[1] = 2.0 * PI * height / hoverd;
		}
		z[2] = z[1] * hoverd;             // kH = kd (H/d)
		z[4] = std::sqrt(std::tanh(z[1]));// c (k/g)^1/2 = sqrt(tanh kd)
	} else {
		z[1] = -1.0;                      // dummy kd
		z[4] = 1.0;
		z[2] = case_wavelength ? (2.0 * PI * height) : (4.0 * PI * PI * height);
	}
	z[3] = 2.0 * PI / z[4];               // tau (gk)^1/2
	if (current_euler) { z[5] = value * std::sqrt(z[2]); z[6] = 0.0; }
	else               { z[6] = value * std::sqrt(z[2]); z[5] = 0.0; }
	z[7] = z[4];                          // ubar (k/g)^1/2  (initial = c)
	z[8] = 0.0;                           // q
	z[9] = 0.5 * z[7] * z[7];             // rk/g
	cosa[0] = 1.0; sina[0] = 0.0;
	z[10] = 0.5 * z[2];                   // k eta_0 (crest)
	for (int i = 1; i <= n; ++i) {
		cosa[i]     = std::cos(i * PI / n);
		cosa[i + n] = std::cos((i + n) * PI / n);
		sina[i]     = std::sin(i * PI / n);
		sina[i + n] = std::sin((i + n) * PI / n);
		z[n + i + 10] = 0.0;              // B_i = 0
		z[i + 10]     = 0.5 * z[2] * cosa[i]; // k eta_m = 0.5 kH cos(m pi/N)
	}
	z[n + 11] = 0.5 * z[2] / z[7];        // B_1 = 0.5 kH / (ubar (k/g)^1/2)
}

// ---------------------------------------------------------------------------
// The 2N+10 nonlinear residuals f_i(z) (Fenton 1988 subroutine EQNS).
// ---------------------------------------------------------------------------
void FentonStream::eqns(const std::vector<double>& z, std::vector<double>& rhs)
{
	const int n = N;

	rhs[1] = finite ? (z[2] - z[1] * hoverd) : (z[1] + 1.0);
	rhs[2] = case_wavelength ? (z[2] - 2.0 * PI * height)
	                         : (z[2] - height * z[3] * z[3]);
	rhs[3] = z[4] * z[3] - 2.0 * PI;
	rhs[4] = z[5] + z[7] - z[4];
	rhs[5] = z[6] + z[7] - z[4];
	if (finite) rhs[5] -= z[8] / z[1];

	// coeff_j = B_j / cosh(j kd)  (finite);  deep branch uses z directly
	std::vector<double> coeff(n + 1, 0.0);
	if (finite) for (int i = 1; i <= n; ++i) coeff[i] = z[n + i + 10] / std::cosh(i * z[1]);

	const int it = current_euler ? 5 : 6;
	rhs[6] = z[it] - value * std::sqrt(z[2]);

	rhs[7] = z[10] + z[n + 10];
	for (int i = 1; i <= n - 1; ++i) rhs[7] += z[10 + i] + z[10 + i];
	rhs[8] = z[10] - z[n + 10] - z[2];

	for (int m = 0; m <= n; ++m) {
		double psi = 0.0, uu = 0.0, vv = 0.0;
		if (finite) {
			for (int j = 1; j <= n; ++j) {
				int nm = (m * j) % (2 * n);
				double e = std::exp(j * (z[1] + z[10 + m]));
				double s = 0.5 * (e - 1.0 / e);   // sinh(j(kd + k eta_m))
				double cc = 0.5 * (e + 1.0 / e);  // cosh(j(kd + k eta_m))
				psi += coeff[j] * s * cosa[nm];
				uu  += j * coeff[j] * cc * cosa[nm];
				vv  += j * coeff[j] * s * sina[nm];
			}
		} else {
			for (int j = 1; j <= n; ++j) {
				int nm = (m * j) % (2 * n);
				double e = std::exp(j * z[10 + m]);
				psi += z[n + j + 10] * e * cosa[nm];
				uu  += j * z[n + j + 10] * e * cosa[nm];
				vv  += j * z[n + j + 10] * e * sina[nm];
			}
		}
		// kinematic free-surface BC:  psi - ubar*k*eta_m - q = 0
		rhs[9 + m] = psi - z[7] * z[10 + m] - z[8];
		// dynamic (Bernoulli) free-surface BC
		double uw = -z[7] + uu;                   // (U - c)(k/g)^1/2 at surface
		rhs[n + 10 + m] = 0.5 * (uw * uw + vv * vv) + z[10 + m] - z[9];
	}
}

// ---------------------------------------------------------------------------
// Configure and solve.
// ---------------------------------------------------------------------------
void FentonStream::set_stream_properties(double _wave_height, double _wave_length, double _wave_period,
	bool specify_length, int current_type, int n_fourier, int n_steps)
{
	wave_height = _wave_height;
	N = (n_fourier > 0) ? n_fourier : 16;
	nstep = (n_steps > 0) ? n_steps : 5;
	case_wavelength = specify_length;
	current_euler = (current_type == 0);
	finite = (depth > 0.0 && depth < 1e8);

	if (wave_height <= 0.0) throw std::invalid_argument("FentonStream: wave_height must be > 0.");
	if (finite && depth <= 0.0) throw std::invalid_argument("FentonStream: depth must be > 0.");

	// dimensionless target ratios for the solver
	if (specify_length) {
		if (_wave_length <= 0.0) throw std::invalid_argument("FentonStream: wave_length must be > 0.");
		wave_length = _wave_length;
		height_full = wave_height / wave_length;                 // H / lambda
	} else {
		if (_wave_period <= 0.0) throw std::invalid_argument("FentonStream: wave_period must be > 0.");
		wave_period = _wave_period;
		height_full = wave_height / (gravity * _wave_period * _wave_period); // H / (g T^2)
	}
	hoverd_full = finite ? (wave_height / depth) : 0.0;
	value = current / std::sqrt(gravity * wave_height);          // c_E or c_S nondim (full H)

	const int num = 2 * N + 10;
	cosa.assign(2 * N + 1, 0.0);
	sina.assign(2 * N + 1, 0.0);
	std::vector<double> z(num + 1, 0.0), sol1(num + 1, 0.0), sol2(num + 1, 0.0);

	const double dhc = height_full / nstep;
	const double dho = hoverd_full / nstep;
	const double crit = 1e-3;

	for (int ns = 1; ns <= nstep; ++ns) {
		height = ns * dhc;
		hoverd = ns * dho;
		if (ns == 1) init(z);
		else for (int i = 1; i <= num; ++i) z[i] = 2.0 * sol2[i] - sol1[i];

		const int maxit = (nstep == 1) ? 20 : 9;
		bool converged = false;
		for (int iter = 0; iter < maxit; ++iter) {
			std::vector<double> rhs1(num + 1, 0.0), rhs2(num + 1, 0.0);
			eqns(z, rhs1);

			std::vector<std::vector<double>> A(num + 1, std::vector<double>(num + 1, 0.0));
			std::vector<double> b(num + 1, 0.0);
			for (int i = 1; i <= num; ++i) {
				double hh = (std::fabs(z[i]) > 1e-4) ? 0.01 * z[i] : 1e-5;
				double zi = z[i];
				z[i] = zi + hh;
				eqns(z, rhs2);
				z[i] = zi;
				for (int j = 1; j <= num; ++j) A[j][i] = (rhs2[j] - rhs1[j]) / hh;
			}
			for (int j = 1; j <= num; ++j) b[j] = -rhs1[j];

			if (!gauss_solve(A, b, num))
				throw std::runtime_error("FentonStream: singular Jacobian in Newton solve.");

			double sum = 0.0;
			for (int i = 1; i <= num; ++i) { z[i] += b[i]; sum += std::fabs(b[i]); }
			double criter = (ns == nstep) ? 0.01 * crit : crit;
			if (sum < criter) { converged = true; break; }
		}
		if (!converged && ns == nstep)
			std::cerr << "FentonStream: warning - Newton iteration did not fully converge.\n";

		for (int i = 1; i <= num; ++i) { sol1[i] = sol2[i]; sol2[i] = z[i]; }
	}

	// --- extract dimensional results ---
	k = finite ? (z[1] / depth) : (z[2] / wave_height);
	c = z[4] * std::sqrt(gravity / k);
	ubar = z[7] * std::sqrt(gravity / k);
	wave_length = 2.0 * PI / k;
	wave_period = z[3] / std::sqrt(gravity * k);
	// Bernoulli constant: z[9] = r k/g,  R = r + g d
	Rc = z[9] * gravity / k + (finite ? gravity * depth : 0.0);

	B.assign(N + 1, 0.0);
	cosh_jkd.assign(N + 1, 1.0);
	for (int j = 1; j <= N; ++j) {
		B[j] = z[N + 10 + j];
		cosh_jkd[j] = finite ? std::cosh(j * k * depth) : 1.0;
	}

	// surface-elevation Fourier coefficients Y_j (Fenton 1988 subroutine OUTPUT)
	Ycoef.assign(N + 1, 0.0);
	for (int j = 1; j <= N; ++j) {
		double sum = 0.5 * (z[10] + z[N + 10] * ((j % 2 == 0) ? 1.0 : -1.0));
		for (int m = 1; m <= N - 1; ++m) sum += z[10 + m] * cosa[(m * j) % (2 * N)];
		Ycoef[j] = 2.0 * sum / N;
	}

	initialized = true;
}

// phase argument kx = k (X - c t), including direction and reference point.
static inline double phase_kx(double k, double theta, double x0, double y0, double c,
                              double t0, double t, double X, double Y)
{
	const double th = theta * PI / 180.0;
	return k * (std::cos(th) * (X - x0) + std::sin(th) * (Y - y0) - c * (t - t0));
}

double FentonStream::eta(double t, double X, double Y)
{
	const double kx = phase_kx(k, theta, x0, y0, c, t0, t, X, Y);
	double keta = 0.5 * Ycoef[N] * std::cos(N * kx);
	for (int j = 1; j <= N - 1; ++j) keta += Ycoef[j] * std::cos(j * kx);
	return keta / k + z0;
}

// Fourier sums of the horizontal (su) and vertical (sw) velocity profiles at a point.
// Finite depth: cosh/sinh(jk(d+y))/cosh(jkd);  deep water: exp(jk y) for both.
void FentonStream::velocity_sums(double kx, double Z, double& su, double& sw) const
{
	su = 0.0; sw = 0.0;
	if (finite) {
		const double ky = k * (Z + depth - z0);      // = k(d + y_Fenton)
		for (int j = 1; j <= N; ++j) {
			double bj = B[j] / cosh_jkd[j];
			su += j * bj * std::cosh(j * ky) * std::cos(j * kx);
			sw += j * bj * std::sinh(j * ky) * std::sin(j * kx);
		}
	} else {
		const double ky = k * (Z - z0);              // y from SWL
		for (int j = 1; j <= N; ++j) {
			double e = B[j] * std::exp(j * ky);
			su += j * e * std::cos(j * kx);
			sw += j * e * std::sin(j * kx);
		}
	}
}

double FentonStream::u(double t, double X, double Y, double Z)
{
	const double kx = phase_kx(k, theta, x0, y0, c, t0, t, X, Y);
	double su, sw; velocity_sums(kx, Z, su, sw);
	double U = c - ubar + std::sqrt(gravity / k) * su;
	return std::cos(theta * PI / 180.0) * U;
}

double FentonStream::v(double t, double X, double Y, double Z)
{
	const double kx = phase_kx(k, theta, x0, y0, c, t0, t, X, Y);
	double su, sw; velocity_sums(kx, Z, su, sw);
	double U = c - ubar + std::sqrt(gravity / k) * su;
	return std::sin(theta * PI / 180.0) * U;
}

double FentonStream::w(double t, double X, double Y, double Z)
{
	const double kx = phase_kx(k, theta, x0, y0, c, t0, t, X, Y);
	double su, sw; velocity_sums(kx, Z, su, sw);
	return std::sqrt(gravity / k) * sw;
}

double FentonStream::dp(double t, double X, double Y, double Z)
{
	// Non-hydrostatic pressure via steady Bernoulli in the wave frame (cf. Stokes5::dp):
	//   p_nh = rho ( R - g(eta + d - z0) - 1/2 q'^2 ),   q'^2 = (U - c)^2 + w^2.
	const double kx = phase_kx(k, theta, x0, y0, c, t0, t, X, Y);
	double su, sw; velocity_sums(kx, Z, su, sw);
	const double sg = std::sqrt(gravity / k);
	const double u_wave = -ubar + sg * su;   // U - c  (wave-frame horizontal)
	const double w_val  = sg * sw;
	const double q_sq = u_wave * u_wave + w_val * w_val;
	const double eta_val = eta(t, X, Y);
	return rho * (Rc - gravity * (eta_val + depth - z0) - 0.5 * q_sq);
}
