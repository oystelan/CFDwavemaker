#pragma once

#include <string>
#include <vector>
#include "Irregular.h"
#include "Utils.h"
#include "lsgrid.h"
#include "Wavemaker.h"
#include "Stokes5.h"

#if defined(SWD_enable)
#include "SpectralWaveData.h"
#endif

typedef struct {
	double x;
	double y;
	double z;
} probe;


class Probes 
{
public:

	std::string ts_path, ts_filename;
	int ts_nprobes;
	double dt, t0, update_time;
	std::vector<probe> coords;

	Probes() {
		// set default values
		ts_path = "./";
		ts_filename = "probe";
		ts_nprobes = 0;
		update_time = 0.;
		t0 = 0.;
		dt = 0.5;
	}
	void init_probes();
	bool checkTime(double tpt);
	void write(double tpt, Irregular& irregular, Ramp& ramp);
	void write(double tpt, lsGrid& sgrid, Irregular& irregular, Ramp& ramp);
	void write(double tpt, Wavemaker& wavemaker, Ramp& ramp);
	void write(double tpt, Stokes5& stokes5, Ramp& ramp);
#if defined(SWD_enable)
	void write(double tpt, SpectralWaveData* swd, Ramp& ramp);
	void write(double tpt,  lsGrid& sgrid, SpectralWaveData* swd, Ramp& ramp);
#endif

};

