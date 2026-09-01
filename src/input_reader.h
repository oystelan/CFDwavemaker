#ifndef CFDWM_INPUT_READER_H
#define CFDWM_INPUT_READER_H

// waveinput.dat reader for CFDwavemaker.
//
// Two-phase design:
//   Phase 1 (parse):   every [section] is read into plain storage.
//                      Nothing is initialized and no section reads
//                      state written by another section, so the order
//                      in which sections appear in the file does not
//                      matter.
//   Phase 2 (resolve): the wave type is derived from [wave type] plus
//                      the presence of [lsgrid], all cross-section
//                      consistency checks run, and the wave objects
//                      are initialized in a fixed, correct sequence.
//
// Robustness rules enforced here (and unit-tested by
// tests/t010_input_reader/):
//   - [section] order is free;
//   - keyword order within a section is free (component tables must
//     stay contiguous blocks, frequency rows before direction rows);
//   - unknown [sections] and unknown keywords raise a warning instead
//     of being silently ignored;
//   - malformed numbers are a hard error instead of silently
//     becoming 0.0;
//   - a duplicated section is a hard error.

#include <string>

#include "Utils.h"     // Ramp, SwdSpectralProps (SWD_enable)

// Version string: normally injected by the makefiles (-DCFDWM_VERSION=...).
#ifndef CFDWM_VERSION
#define CFDWM_VERSION "(unversioned build)"
#endif

#ifndef largeval
#define largeval 1.E12
#endif
#ifndef HUGE
#define HUGE 1.e12
#endif

// Storage class of input data (phase-1 target; owned by
// CFDwavemaker.cpp as the global `inputdata`).
class CFDwavemakerInputdata {
public:
	double depth = 0.;
	bool depth_read = false;
	double x_pos = 0., y_pos = 0., tofmax = 0.;
	double mtheta = 0.;
	double swl = 0.;
	bool bw_auto_calc = false;
	int wavetype = 0;
	bool property_read = false;
	double rho = 1.;
	double gravity = 9.81;
	double ampl = 1.;
	// SWD parameters
	double nsumx = -1, nsumy = -1, impl = 0, ipol = 0, norder = 0;
	bool dc_bias = false;
	std::string swdFileName;

	// sloping water related parameters
	double slope_angle = 0.; // tilt angle of the slope
	double slope_direction = 0.; // direction of water slope (degrees)
	double slope_rotmat[9] = {}; // store rotation matrix (flattened)
	double slope_radius_limiter = HUGE;
	// lsgrid related data
	double lsgrid_domain[4] = {};
	// The grid extent must come from somewhere: either a 'bounds' line here
	// or wave_SetLsgridBounds() before wave_Initialize(). Neither is not a
	// usable default -- a zero-extent grid stores nothing.
	bool   lsgrid_domain_given = false;
	int lsgrid_nx = 5;
	int lsgrid_ny = 0;
	// Whether nx/ny were given explicitly. When they are not, and the
	// sigma-FFT applies, the lattice is taken from the component set alone
	// (the generator decides how far the spectral tail reaches, and the
	// lattice follows). The other fill paths have no lattice to fall back
	// on, so there nx/ny are required.
	// Whether vertical_theory was stated explicitly. An equal-dk lattice
	// implies smit -- the sigma-FFT fills the grid with the smit
	// sigma-transform, so anything else would be inconsistent with the
	// stored field -- so it is defaulted rather than demanded.
	bool vertical_theory_set = false;
	bool lsgrid_nx_set = false;
	bool lsgrid_ny_set = false;
	int lsgrid_nl = 4;
	// nl omitted or 'nl auto': pick the smallest layer count whose vertical
	// profile reconstruction meets lsgrid_nl_tol (resolved at initialization
	// from the loaded sea state; needs the discrete component set, so the
	// SWD backend keeps the default instead).
	bool lsgrid_nl_set = false;
	double lsgrid_nl_tol = 2e-4;   // 'nl_tol': target relative L2 error
	double lsgrid_t0 = 0.;
	double lsgrid_dt = 0.1;
	bool lsgrid_dt_set = false;
	// Points per shortest wavelength / period used when nx, ny or dt are
	// omitted. Cubic Hermite converges at 4th order: ~1.7% velocity error at
	// 4 points, 0.5% at 6, 0.14% at 8.
	double lsgrid_pts_per_wave = 6.0;
	// exponential layer stretching; 0 = uniform layers. The bottom-to-surface
	// thickness ratio realised by nl layers is exp(beta*(nl-2)/(nl-1)).
	double lsgrid_stretch_beta = 3.0;
	// 'stretch_beta auto': pick beta at initialization from the loaded sea
	// state instead of using the value above. Default, since the fixed legacy
	// beta=3.0 is a poor choice in shallow water.
	bool lsgrid_stretch_auto = true;
	// [seabed] variable bathymetry from a triangulated surface file
	std::string seabed_file;
	double seabed_zshift = 0.;
	double seabed_xshift = 0.;
	double seabed_yshift = 0.;
	double seabed_min_depth = 0.1;
	// false = geometry-only mode ([seabed] kinematics 0): the surface is
	// loaded for wave_Seabed*() queries but the kinematics assume the
	// constant depth from [general input data] (allows the sigma-FFT path)
	bool seabed_kinematics = true;
	bool seabed_read = false;
	// Border mode: instead of describing the region NOT to store (the
	// complement, in absolute metres, which has to be re-derived whenever
	// the client's grid changes), name the borders that actually need
	// kinematics. The library then keeps a strip along each of them and
	// drops the interior. {left, right, bottom, top}, matching the basilisk
	// boundary order used by the linked solver.
	int  lsgrid_border[4] = {0, 0, 0, 0};
	bool lsgrid_border_mode = false;
	// Set programmatically through wave_SetLsgridBounds/Borders() before
	// wave_Initialize(); applied after the file is parsed, so a linked
	// solver that knows its own domain wins over the file.
	bool   lsgrid_domain_override = false;
	double lsgrid_domain_set[4] = {0., 0., 0., 0.};
	// How far INWARD from a border the client will ask for kinematics
	// (metres). The GABC samples the boundary face and the first interior
	// cell centre, so this is Delta_solver/2 -- which the library cannot
	// know. Sizes the stored strip; 0 = use the built-in minimum.
	double lsgrid_border_margin = 0.0;
	bool   lsgrid_border_override = false;
	int    lsgrid_border_set[4] = {0, 0, 0, 0};
	bool lsgrid_ignore_at_init = 0;
	bool lsgrid_init_only = 0;
	// 'store_pressure': store the dynamic pressure as a 4th grid field so
	// wave_nhPres()/wave_DynPres() are served by spline interpolation at
	// full second order (the GABC's incident p_nh). -1 = auto: on for
	// border-mode grids (the GABC use case), off for full-domain
	// initialization grids; 0/1 = explicit off/on.
	int lsgrid_store_pressure = -1;

	// Nonhydrostatic pressure options
	// When false (default): p_nh = -rho*(dphi/dt + g*eta)  [linearized, suitable for non-hydrostatic CFD solvers]
	// When true:            p_nh = -rho*(dphi/dt + 0.5*|u|^2 + g*eta)  [full Bernoulli form]
	bool nhpres_include_kinetic = false;

	// Regular wave (wavetype 21) theory selection under [wave input regular]:
	//   0 = Stokes 5th (default), 1 = Fenton stream function.
	int  regular_theory = 0;
	int  wave_order = 16;                 // Fourier order N for stream function
	bool regular_specify_length = true;   // true: user gave wave_length; false: wave_period
	bool regular_length_given = false;
	bool regular_period_given = false;

	// ---- reader bookkeeping (phase 1 -> phase 2) ----
	// base type from [wave type]; the final `wavetype` additionally
	// folds in the presence of [lsgrid] (1 -> 5, 31 -> 34)
	int wavetype_base = 0;
	bool wavetype_read = false;
	bool general_read = false;
	bool refpoint_read = false;
	bool ramps_read = false;
	bool irregular_read = false;
	bool regular_read = false;
	bool piston_read = false;
	bool swd_input_read = false;
	bool vtk_input_read = false;
	bool lsgrid_read = false;
	bool lsgrid2vtk_read = false;
	bool sloping_read = false;

#if defined(SWD_enable)
	// Spectral properties derived from the SWD file at init (long-crested shapes only)
	SwdSpectralProps swd_spectral;
#endif

	CFDwavemakerInputdata() {
	}

	~CFDwavemakerInputdata() {
	}
};

// Locate waveinput.dat (same search paths as before: ./, ./constant/,
// ../input_files/, ../), strip # comments and blank lines, check the
// @v version tag, then run both reader phases. Returns 0 on success,
// -1 when no input file could be opened. Fatal format errors keep the
// historical behavior and terminate (via cfdwm_mpi::abort, so an MPI
// job is taken down as a whole).
int cfdwm_read_input();

#endif // CFDWM_INPUT_READER_H
