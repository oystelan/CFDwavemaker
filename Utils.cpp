
#include "Utils.h"
#include <algorithm>
#include <cmath>

// -------------------------------------------------------------------------------------------------
// ramp class function
// -------------------------------------------------------------------------------------------------

//Define some useful functions
/* Rampfunction */
// NB: Not yet implemented inverse ramp
double Ramp::ramp1d(double x, double xstart, double xend, bool inv) {

	if (inv){
		return std::max(std::min(1. - ((x - xstart) / (xend - xstart)),1.), 0.);
	}
	else {
		return std::max(std::min(((x - xstart) / (xend - xstart)), 1.), 0.);
	}
}

bool comp(double a, double b)
{
	return (a < b);
}

// Todo: implement ramps
double Ramp::ramp(double t, double x, double y) {
	if (ramp_init) {
		double ramps[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
		if (ramp_init_time_up) {
			ramps[0] = ramp1d(t, time_rampup_start, time_rampup_end, false);
		}
		if (ramp_init_time_down) {
			ramps[1] = ramp1d(t, time_rampdown_start, time_rampdown_end, true);
		}
		if (ramp_init_x_up) {
			ramps[2] = ramp1d(x, x_rampup_start, x_rampup_end, false);
		}
		if (ramp_init_x_down) {
			ramps[3] = ramp1d(x, x_rampdown_start, x_rampdown_end, true);
		}
		if (ramp_init_y_up) {
			ramps[4] = ramp1d(y, y_rampup_start, y_rampup_end, false);
		}
		if (ramp_init_y_down) {
			ramps[5] = ramp1d(y, y_rampdown_start, y_rampdown_end, true);
		}
		return *std::min_element(std::begin(ramps), std::end(ramps));
	}
	else {
		return 1.0;
	}
}
