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
	// standard callable functions for static linking
	int wave_Initialize(); //function needs to be called at startup (reading wave input file)
	int wave_Cleanup(); // Needs to be called before program end to clear variables stored in memory
	double wave_VeloX(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_VeloY(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_VeloZ(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_DynPres(double, double, double, double); // input variables are {xpoint,ypoint,zpoint,time}
	double wave_SurfElev(double, double, double); // input variables are {xpoint,ypoint,time}
	double wave_VFrac(double, double, double, double, double);// input variables are {xpoint,ypoint,zpoint,time, delta_cellsize}
	

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