
#include "Utils.h"
#include <algorithm>
#include <cmath>
#if defined(_WIN32)
#include < direct.h >
#endif

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


/******************************************************************************
 * Checks to see if a directory exists. Note: This method only checks the
 * existence of the full path AND if path leaf is a dir.
 *
 * @return  >0 if dir exists AND is a dir,
 *           0 if dir does not exist OR exists but not a dir,
 *          <0 if an error occurred (errno is also set)
 *****************************************************************************/
int dirExists(const char* const path)
{
	struct stat info;

	int statRC = stat(path, &info);
	if (statRC != 0)
	{
		if (errno == ENOENT) { return 0; } // something along the path does not exist
		if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
		return -1;
	}

	return (info.st_mode & S_IFDIR) ? 1 : 0;
}

void createDirectory(std::string sPath) {


	int nError = 0;
#if defined(_WIN32)

	nError = _mkdir(sPath.c_str()); // can be used on Windows
#else
	mode_t nMode = 0733; // UNIX style permissions
	nError = mkdir(sPath.c_str(), nMode); // can be used on non-Windows
#endif
	if (nError != 0) {
		// handle your error here
	}
}
