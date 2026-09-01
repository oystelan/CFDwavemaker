//=================================
// include guard
#ifndef __CFDWAVEMAKER_H_INCLUDED__
#define __CFDWAVEMAKER_H_INCLUDED__

// Export settings for dynamic linking using extern C on linux and windows
#if defined(_MSC_VER)
//  Microsoft 
#define EXPORT __declspec(dllexport)
#elif defined(__GNUC__)
//  GCC
#define EXPORT __attribute__((visibility("default")))
#else
//  do nothing and hope for the best?
#define EXPORT
#pragma warning Unknown dynamic link import/export semantics.
#endif



// Functions available in CFDwavemaker.cpp
#ifdef __cplusplus
extern "C" {
#endif
	int CFDwavemaker_is_initialized();
	// standard callable functions for static linking
	int wave_Initialize(); //function needs to be called at startup (reading wave input file)
	int wave_Cleanup(); // Needs to be called before program end to clear variables stored in memory
	double wave_VeloX(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_VeloY(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_VeloZ(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_DynPres(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_nhPres(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time} - nonhydrostatic pressure. With an lsgrid carrying a pressure field ([lsgrid] store_pressure, auto-on for border-mode grids) this is spline-interpolated at full second order; without one, irregular+lsgrid falls back to first-order direct summation
	double wave_SurfElev(double, double, double); // input variables are {xpoint,ypoint,time}
	double wave_VFrac(double, double, double, double, double);// input variables are {xpoint,ypoint,zpoint,time, delta_cellsize}
	double wave_Seabed(double, double); // input variables are {xpoint,ypoint}; min-depth clamped bed level (what the kinematics see)
	double wave_SeabedRaw(double, double); // input variables are {xpoint,ypoint}; unclamped bed level, positive above swl on land (for solvers importing the bathymetry)
	int wave_SetLsgridBounds(double, double, double, double); // input {xmin,xmax,ymin,ymax}: set the lsgrid extent from the linked solver instead of waveinput.dat. MUST be called BEFORE wave_Initialize(); returns 0 (and does nothing) if called later
	int wave_SetBorderMargin(double); // input: how far INWARD from a border the caller will request kinematics [m]. The stored strip is then sized in metres instead of grid nodes, so a coarse CFD grid over a finely resolved wave lattice cannot outrun it. For a Basilisk multilayer GABC this is Delta/2. MUST be called BEFORE wave_Initialize()
	int wave_SetLsgridBorders(int, int, int, int); // input {left,right,bottom,top} as 0/1: which domain borders need kinematics. Only a strip along those is stored, the interior is dropped. Resolution-independent, unlike ignore_subdomain. MUST be called BEFORE wave_Initialize(). Equivalent to border_left/right/bottom/top in the [lsgrid] section
	int wave_SeabedActive(void); // returns 1 if a [seabed] surface (variable bathymetry) is loaded
	void wave_SeabedRelease(void); // free the [seabed] mesh after init-time sampling (geometry-only mode); no-op/refused when the kinematics use it; queries after release abort
	double wave_LocalDepth(double, double); // input {xpoint,ypoint}: local still-water depth [m], positive down. The [seabed] surface whenever one is loaded - whatever generates the kinematics, geometry-only mode included - else the constant 'depth' from [general input data]. Unclamped: negative on land
	int wave_SeabedBounds(double*, double*, double*, double*); // outputs {xmin,xmax,ymin,ymax} of the [seabed] footprint bbox; returns 1 if a surface is loaded, else 0 and outputs untouched. Necessary but not sufficient - the mesh need not fill its bbox
	int wave_SeabedContains(double, double); // input {xpoint,ypoint}: 1 if the point is covered by a triangle of the [seabed] surface, 0 if outside (or no surface). Exact test, no clamping and no warning - use it to validate a solver domain up front
	long wave_SeabedOutsideCount(void); // number of depth queries served by out-of-footprint clamping so far. Non-zero = part of the sampled region is a flat zero-slope skirt, not real bathymetry (a common cause of boundary blow-up)
	double wave_OptimalLayerRatio(int); // input: number of solver layers. Optimal geometric layer grading from the sea-state RMS velocity profile, returned as the top-layer-to-uniform thickness ratio in (0,1] (Basilisk geometric_beta / 'cell_aspect_ratio' convention); 1.0 = uniform, negative on failure. Evaluated at the reference ('depth') depth
	double wave_OptimalLayerRatioAtDepth(int, double); // input {nlayers, depth}: same grading at an ARBITRARY depth, with the components carried there by conserved frequency + Green's law. Table-interpolated (built on first call for a given nlayers), so it is cheap enough to call per grid column; prime it with one call before entering a parallel region
	double wave_OptimalLayerRatioAt(int, double, double); // input {nlayers, xpoint, ypoint}: wave_OptimalLayerRatioAtDepth() at wave_LocalDepth(x,y)
	double* wave_Kinematics(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_piston_position(double); // function for extracting paddle position in case of wave piston
	
	// sgrid updater function
	// Checks the timestep and, if necessary, updates the sgrid kinematics.
	// MPI note: in the MPI build this is a COLLECTIVE call for the lsgrid
	// wave types (5, 34) - every rank in the communicator must call it
	// with the same time before querying kinematics, or the ranks
	// deadlock in the grid assembly (MPI_Allgatherv).
	void wave_force_update(double);

	// MPI (optional): set the communicator used to distribute the lsGridSpline
	// second-order kinematics computation. Takes a Fortran integer communicator
	// handle (MPI_Comm_c2f in C, or the native handle in Fortran); defaults to
	// MPI_COMM_WORLD. Must be called before wave_Initialize(). No-op unless the
	// library was built with -DMPI_enable=1.
	EXPORT void wave_set_mpi_comm(int fortran_comm);
	
	// update probes
	EXPORT void update_probes(double); // writes wave kinematics data to probe files for a given time step as argument.

	// some helpful function in case of irregular waves
	double wave_phase_speed(int); // returns the phase velocity based on spectral mean wave period t1 for opt=1, and spectral zero crossing period t2 for opt=2
	double wave_mean_period(int); // calculate mean wave period
	double wave_mean_length(int); // calculate mean wave length

	double wave_water_depth(); // returns the water depth used to wave calculation
	double wave_PeakPeriod(void); // spectral peak period Tp [s] (largest-amplitude component); 0 if not irregular
	double wave_PhaseVelocity(double, double, double); // input {xpoint, ypoint, T0}: linear phase speed at the LOCAL depth ([seabed]-aware, geometry-only mode included) for period T0 (T0 <= 0 -> spectral mean period T1, then Tp); 0 on dry/land

	// external functions used by COMFLOW
	EXPORT double VelocityX(int, int, int, double, double, double, double);
	EXPORT double VelocityY(int, int, int, double, double, double, double); 
	EXPORT double VelocityZ(int, int, int, double, double, double, double);
	EXPORT double DynamicPressure(int, int, int, double, double, double, double);
	EXPORT double SurfaceElevation(int, int, double, double, double);
	EXPORT double VolumeFraction(double, double, double, double, double);
	EXPORT int Init(void*);
	EXPORT int Prepare(double, int);
	EXPORT int Cleanup();
#ifdef __cplusplus
}
#endif
#endif