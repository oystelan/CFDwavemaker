// Tiny probe for the input-reader tests: initialize from ./waveinput.dat
// and print a deterministic kinematics fingerprint. run_test.py compares
// the fingerprint across section-order permutations of the same file.
//
// Build:
//   g++ -std=c++11 -O2 -fopenmp probe.cpp -I../../src \
//       ../../builds/linux64/libCFDwavemaker_openmp.a -lgfortran -lfftw3 \
//       -o probe
#include <cstdio>
#include "CFDwavemaker.h"

int main() {
	if (wave_Initialize() != 0) {
		printf("PROBE-FAIL: wave_Initialize\n");
		return 1;
	}
	for (int i = 0; i < 5; i++) {
		double x = 11.0 + 17.0 * i, y = 3.0 + 5.0 * i, t = 1.0 + 2.3 * i;
		wave_force_update(t);
		double eta = wave_SurfElev(x, y, t);
		double z = eta - 2.0 - 1.5 * i;
		printf("PROBE %.12e %.12e %.12e %.12e\n", eta,
			wave_VeloX(x, y, z, t), wave_VeloY(x, y, z, t),
			wave_VeloZ(x, y, z, t));
	}
	printf("PROBE depth %.12e\n", wave_water_depth());
	wave_Cleanup();
	return 0;
}
