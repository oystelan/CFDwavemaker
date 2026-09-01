#ifndef PI
#define PI 3.1415926535897
#endif

#ifndef SGridSpline_H
#define SGridSpline_H
#include <vector>
#include "Irregular.h"
#if SWD_enable
#include "SpectralWaveData.h"
#endif

#define largeval 1.E12

#include "sigma_fft.h"
#include "layer_stretch.h"

class lsGridSpline {
public:
	lsGridSpline() {
		// Set default values
		stretch_beta = 3.0;
		stretch_auto = true;
		nl = 10;
		dt = 0.5;
		swl = 0.;
	};
	~lsGridSpline() {
		if (dump_vtk) {
			write_vtk(true);
		}
		// NB: was `delete[] UX, UY, UZ, ETA, IGNORE;` -- the comma operator
		// only ever freed UX
		delete[] UX; delete[] UY; delete[] UZ; delete[] ETA;
		delete[] ETAdx; delete[] ETAdy; delete[] ETAdt;
		delete[] Udt; delete[] Udx; delete[] Udy; delete[] Uds;
		delete[] Vdt; delete[] Vdx; delete[] Vdy; delete[] Vds;
		delete[] Wdt; delete[] Wdx; delete[] Wdy; delete[] Wds;
		delete[] P; delete[] Pdt; delete[] Pdx; delete[] Pdy; delete[] Pds;
		delete[] IGNORE;
	};

	// Primary field variables
	double* UX = nullptr, * UY = nullptr, * UZ = nullptr; // 4D dimensions, (time, x, y, z)
	double* ETA = nullptr;// Surface grid (3D) (time, x, y)
	// Secondary field variables (derivatives)
	double* ETAdx = nullptr,* ETAdy = nullptr,* ETAdt = nullptr;
	double* Udt = nullptr, * Udx = nullptr,* Udy = nullptr,* Uds = nullptr;
	double* Vdt = nullptr, * Vdx = nullptr,* Vdy = nullptr,* Vds = nullptr;
	double* Wdt = nullptr, * Wdx = nullptr,* Wdy = nullptr,* Wds = nullptr;
	// Dynamic pressure field (rho included; -rho*dphi/dt - 0.5*rho*|u1|^2,
	// i.e. what Irregular::dp() returns), stored only when store_pressure
	// is set: the non-hydrostatic incident pressure a generating-absorbing
	// boundary condition subtracts is p_nh = P - rho*g*eta, and without a
	// stored field wave_nhPres() can only serve a first-order fallback.
	double* P = nullptr;
	double* Pdt = nullptr, * Pdx = nullptr, * Pdy = nullptr, * Pds = nullptr;
	// Whether P is allocated and filled. Wired from [lsgrid]
	// 'store_pressure' (default: on for border-mode grids -- the GABC use
	// case -- off for full-domain initialization grids).
	bool store_pressure = false;

	int* IGNORE; // matrix with cells to ignore when updating surface elevation
	int nx, ny, nl;

	// ---- storage compaction -------------------------------------------
	// The field arrays are indexed by STORAGE node index, not by the dense
	// i*ny + j. Until compact() runs, slot is the identity and nnodes ==
	// nx*ny, i.e. exactly the old dense layout. After compaction the nodes
	// outside the selected borders are not stored at all (slot holds -1)
	// and the arrays shrink to the retained strips plus the halo that the
	// interpolation and gradient stencils need.
	//
	// ignore_at_init decides WHEN: with it set, compaction happens before
	// the first fill, so the large arrays are never allocated; without it,
	// the whole grid is filled at t = 0 as before (a solver initialising
	// its full domain needs that) and compaction happens on the first
	// time advance, once the t = 0 field has served its purpose.
	// The LAST node of every time slab is a write-only scratch node that
	// dropped nodes resolve to, so a stray access can never run off the
	// end of an array: nnodes counts the stored nodes PLUS that scratch.
	std::vector<int> slot;      // nx*ny -> storage index, or -1 if dropped
	std::vector<long> col_off;  // nx+1: first storage index of column i
	long nnodes = 0;            // slab stride = kept nodes + 1 scratch
	long nkept = 0;             // nodes actually carrying data
	long trash = 0;             // index of the scratch node (== nnodes - 1)
	bool compacted = false;
	inline long nid(int i, int j) const {
		const long s = slot[(long)i * ny + j];
		return s >= 0 ? s : trash;
	}
	inline bool stored(int i, int j) const { return slot[(long)i * ny + j] >= 0; }
	// true when the 4-neighbourhood of (i, j) is stored, i.e. when its
	// finite-difference gradients are meaningful
	inline bool grad_ok(int i, int j) const {
		return stored(i, j)
			&& stored(i > 0 ? i - 1 : i, j) && stored(i < nx - 1 ? i + 1 : i, j)
			&& stored(i, j > 0 ? j - 1 : j) && stored(i, j < ny - 1 ? j + 1 : j);
	}
	// Set by the fill entry points: initialize_kinematics() fills the whole
	// grid, the _with_ignore variant and update() skip the ignore box.
	bool fill_all = true;
	// Which nodes a fill must produce. Once compacted this is "every node
	// that exists" -- note that includes the halo ring INSIDE the ignore
	// box, which must be filled or the gradients of the retained nodes
	// next to it are wrong.
	inline bool fill_node(int i, int j) const {
		if (compacted) return stored(i, j);
		return fill_all || !IGNORE[(long)i * ny + j];
	}
	void build_identity_slot();
	// Drop the nodes the ignore box covers. Safe both before the first
	// fill (nothing to preserve) and after it (retained data is carried
	// over). No-op when there is no ignore box.
	void compact();

	// MPI column decomposition (work-balanced by active-cell count). Computed once
	// after the IGNORE mask is set; reused by initialize_kinematics and update.
	int mpi_i0 = 0, mpi_i1 = 0;      // this rank's column range [mpi_i0, mpi_i1)
	std::vector<int> mpi_bounds;     // partition boundaries for all ranks (size+1)
	void compute_mpi_partition();
	double domain[4] = {};
	// Border mode: store only a strip along the named borders. {left,
	// right, bottom, top}. Declaring WHICH borders need kinematics is
	// resolution-independent, unlike an ignore rectangle in metres.
	int  border[4] = {0, 0, 0, 0};
	bool border_mode = false;
	// Nodes retained inward from a selected border, before the stencil
	// halo is added. Two covers the interpolation cell the client's border
	// queries land in (the GABC samples the boundary face and, for the
	// outgoing-residual decomposition, the first interior cell centre --
	// well under one lsgrid cell in any sane configuration). compact()
	// then adds HALO for the gradient stencils.
	static const int BORDER_NODES = 2;
	// Client-declared inward query reach [m]; 0 = use BORDER_NODES alone.
	double border_margin = 0.0;
	bool ignore_at_init = false;
	int tstep = 0;

	bool disable_checkbounds = true;
	bool init_only = false;

	std::string vtk_directory_path = "./";
	std::string vtk_prefix = "kinematics";
	std::string vtk_timelabel = "TimeValue";
	bool dump_vtk = false;

	double water_depth;
	// Variable bathymetry: when set (and active), every sigma-mapping
	// uses the local depth h(x, y) from the seabed surface instead of
	// the constant water_depth, so the layer distribution follows the
	// bed. Set by CFDwavemaker when a [seabed] section is present.
	const Seabed* seabed = nullptr;
	double local_depth(double xpt, double ypt) const {
		return (seabed && seabed->active())
			? seabed->depth(xpt, ypt) : water_depth;
	}
	double swl;
	double dx, dy, ds;
	double t0, dt;
	// Vertical layer stretching: exponential mapping with a single
	// parameter. Adjacent layer thicknesses form a geometric progression with
	// ratio exp(stretch_beta/(nl-1)), so the bottom-to-surface thickness ratio
	// actually realised is exp(stretch_beta*(nl-2)/(nl-1)) -- see
	// layer_thickness_ratio(). stretch_beta = 0 gives uniform layers.
	double stretch_beta;
	// When set, stretch_beta is chosen automatically at initialization to
	// minimise the vertical interpolation error for the loaded sea state
	// (resolve_auto_stretch). Cleared once resolved.
	bool stretch_auto;
	// When set, the layer count itself is chosen at initialization: the
	// smallest nl whose vertical Hermite reconstruction of the sea-state RMS
	// velocity profile meets nl_tol (re-optimising the grading per candidate
	// under stretch_auto, else evaluating the prescribed beta). Cleared once
	// resolved; the storage is reallocated if the count changes.
	bool nl_auto = false;
	double nl_tol = 2e-4;

	double bxmin, bxmax, bymin, bymax;

	double slayer(int layerno);
	double clayer(int layerno);

	// sGrid transform functions
	double z2s(double z, double wave_elev, double depth);
	double s2z(double s, double wave_elev, double depth);

	// exponential stretching: u in [0,1] (0 = seabed, 1 = surface)
	// <-> sigma in [-1,0]
	double u2sigma(double u);
	double sigma2u(double sigma);

	double sigmaS2cart(double s, double wave_elev, double depth);
	double cart2sigmaS(double zpt, double wave_elev, double depth);

	void square_vals(double* C, double* data, int nxp, int nyp, int tid);
	void cube_vals(double* C, double* data, int nxp, int nyp, int nlp, int tid);
	double spline_interp_velo(double* U, double* Udt, double* Udx, double* Udy, double* Uds, int nxp, int nyp, int nsp0, int nsp1, double xd, double yd, double sd0, double sd1, double td, int tid1, int tid2);

	
	void write_vtk(bool endtime);
	void allocate();
	void reallocate();   // free + allocate after a grid resize

	// Sigma-FFT fast fill (equal-dk + smit): fills a whole time level
	// with FFTs instead of pointwise sums. periodic_x wraps queries
	// outside the domain by the spatial period.
	SigmaFFTFill sfft;
	bool periodic_x = false;
	bool periodic_y = false;   // 2D short-crested lattice mode
	// Bottom-to-surface layer thickness ratio actually realised by nl layers.
	double layer_thickness_ratio() const;
	// Picks stretch_beta from the loaded sea state when stretch_auto is set.
	// Must run before anything that consumes the mapping (sigma-FFT setup).
	void resolve_auto_stretch(Irregular& irregular);

	bool try_sigma_fft_setup(Irregular& irregular);
	void sigma_fft_fill_level(Irregular& irregular, double t, int tt);

	// Second order theory
	void initialize_kinematics(Irregular& irregular);
	void initialize_kinematics_with_ignore(Irregular& irregular);
	void update(Irregular& irregular, double t_target);

#if SWD_enable
	// SWD as kinematics source: sample the swd field onto the grid.
	// NOTE: the SWD library is stateful and supports only ONE open object per
	// process (fixed Fortran file unit in the bundled build), so all SWD
	// sampling below is strictly serial within a process. MPI ranks are
	// separate processes (each with its own swd object), so the column
	// partition still parallelizes the fill across ranks.
	void initialize_kinematics(SpectralWaveData* swd);
	void initialize_kinematics_with_ignore(SpectralWaveData* swd);
	void update(SpectralWaveData* swd, double t_target);
#endif

	void set_ignore();
	// gradient update functions
	void update_gradient_dt(double* data, double* graddt);
	void update_gradient_eta_dt(double* data, double* graddt);
	void update_gradient_eta_dxdy(double* eta, double* gradx,double* grady);
	void update_gradient_dxdydz(double* data, double* gradx, double* grady, double* gradz);

	// Grid interpolation functions. Returns {eta, u, v, w}; with
	// with_pressure set (requires store_pressure) a 5th element is
	// appended: the interpolated dynamic pressure P.
	std::vector<double> get_kinematics_at_point(double tpt, double xpt, double ypt, double zpt, double h, bool with_pressure = false);


	// routines for extracting kinematics at walls
	bool CheckTime(double tpt);
	bool CheckBounds();

	void update_bounds(double xpt, double ypt);
	
	// export_functions to validate kinematics
	void export_vts(FILE* fp, bool last);

};

#endif