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
	double wave_SurfElev(double, double, double); // input variables are {xpoint,ypoint,time}
	double wave_VFrac(double, double, double, double, double);// input variables are {xpoint,ypoint,zpoint,time, delta_cellsize}
	double wave_Seabed(double, double); // input variables are {xpoint,ypoint}
	double* wave_Kinematics(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}

	// sgrid updater function
	void wave_force_update(double); // function checking timestep and updating if neccessary the sgrid kinematics
	
	// update probes
	EXPORT void update_probes(double); // writes wave kinematics data to probe files for a given time step as argument.

	// some helpful function in case of irregular waves
	double wave_phase_velocity(int); // returns the phase velocity based on spectral mean wave period t1 for opt=1, and spectral zero crossing period t2 for opt=2
	double wave_mean_period(int); // calculate mean wave period
	double wave_mean_length(int); // calculate mean wave length

	double wave_water_depth(); // returns the water depth used to wave calculation

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