#ifndef FentonStream_H
#define FentonStream_H

#include <vector>

#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef G
#define G 9.81
#endif

// Fenton (1988) stream-function (Rienecker & Fenton) steady wave theory.
//   J.D. Fenton, "The numerical solution of steady water wave problems",
//   Computers & Geosciences 14(3), 357-368, 1988.
//
// Standalone class analogous to Stokes5: given wave height + (length OR period),
// water depth and (optionally) a mean current, it solves the nonlinear system for
// the Fourier coefficients B_j and returns eta, u, v, w and non-hydrostatic pressure.
class FentonStream {
private:
	// --- solver configuration ---
	int N = 16;             // number of Fourier components
	int nstep = 5;          // wave-height continuation steps
	bool finite = true;     // finite depth (true) or deep water (false)
	bool case_wavelength = true;   // input is wavelength (true) or period (false)
	bool current_euler = true;     // current given as Eulerian mean (true) or Stokes/mass-transport (false)

	// stepped (working) and full target values used inside eqns()
	double height = 0.0;    // H/lambda  (wavelength case) or H/(g T^2) (period case)
	double hoverd = 0.0;    // H/d
	double height_full = 0.0;
	double hoverd_full = 0.0;
	double value = 0.0;     // c_E/(gH)^1/2  or  c_S/(gH)^1/2

	// --- solved quantities ---
	std::vector<double> B;         // B[1..N] Fourier coefficients (dimensionless)
	std::vector<double> Ycoef;     // Y[1..N] surface-elevation Fourier coefficients (=k*eta coeffs)
	std::vector<double> cosh_jkd;  // cosh(j k d), j=1..N  (=1 for deep water branch)
	std::vector<double> cosa, sina;// cos/sin(i pi/N), i=0..2N (collocation)

	double ubar = 0.0;      // mean fluid speed in the wave frame  (dimensional)
	double Rc = 0.0;        // Bernoulli constant R (dimensional)

	void init(std::vector<double>& z);
	void eqns(const std::vector<double>& z, std::vector<double>& rhs);
	void velocity_sums(double kx, double Z, double& su, double& sw) const;
	static bool gauss_solve(std::vector<std::vector<double>>& A, std::vector<double>& b, int n);

public:
	bool initialized = false;
	double wave_length = 0.0;
	double wave_height = 0.0;
	double wave_period = 0.0;
	double current = 0.0;   // mean current (m/s); interpreted per current_euler
	double depth = 0.0;
	double gravity = G;
	double rho = 1025.0;
	double x0 = 0.0;
	double y0 = 0.0;
	double z0 = 0.0;        // still water line
	double t0 = 0.0;
	double theta = 0.0;     // propagation direction [deg]

	double k = 0.0;         // wavenumber
	double c = 0.0;         // wave celerity

	// Configure and solve. Specify EITHER wave_length OR wave_period (the other is
	// computed). current_type: 0 = Eulerian mean current, 1 = Stokes (mass-transport).
	// n_fourier: number of Fourier modes; n_steps: wave-height continuation steps.
	void set_stream_properties(double _wave_height, double _wave_length, double _wave_period,
		bool specify_length, int current_type, int n_fourier, int n_steps);

	double eta(double t, double X, double Y);
	double u(double t, double X, double Y, double Z);
	double v(double t, double X, double Y, double Z);
	double w(double t, double X, double Y, double Z);
	double dp(double t, double X, double Y, double Z);

	// diagnostics
	int    num_components() const { return N; }
	double coefficient(int j) const { return (j >= 1 && j <= N) ? B[j] : 0.0; }
	double mean_velocity() const { return ubar; }
	double bernoulli_R() const { return Rc; }
};

#endif
