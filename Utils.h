#ifndef Utils_H
#define Utils_H

#include "Irregular.h"
#include <sstream>

// Ramp class contains various ramp functions
class Ramp {
private:
	double ramp1d(double x, double xstart, double xend, bool inv);
public:
	int timeramp, xramp, yramp; // logical operators for switching on and of rampfunctions
	double time_rampup_start, time_rampup_end, time_rampdown_start, time_rampdown_end;
	double x_rampup_start, x_rampup_end, x_rampdown_start, x_rampdown_end;
	double y_rampup_start, y_rampup_end, y_rampdown_start, y_rampdown_end;
	bool ramp_init = false;
	bool ramp_init_time_up = false;
	bool ramp_init_time_down = false;
	bool ramp_init_x_up = false;
	bool ramp_init_x_down = false;
	bool ramp_init_y_up = false;
	bool ramp_init_y_down = false;

	double ramp(double time, double x, double y);
	
};

int dirExists(const char* const path);

void createDirectory(std::string sPath);

#endif
