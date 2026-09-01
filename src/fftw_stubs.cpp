// Stubs replacing the FFTW-3 entry points referenced by the SWD library's
// optional grid-FFT backend (swd_fft_fftw3). CFDwavemaker only uses the
// pointwise SWD API, which never reaches these — the stubs exist so the
// library links without FFTW (GPL), keeping the distributed binaries
// MIT-clean. If a future code path does reach them, fail loudly instead of
// returning garbage.
#include <cstdio>
#include <cstdlib>

namespace {
void swd_fft_unavailable() {
	std::fprintf(stderr,
	    "CFDwavemaker: SWD grid-FFT API (elev_fft) is not available in this "
	    "build (compiled without FFTW to keep the library MIT-clean). "
	    "Use the pointwise SWD API.\n");
	std::exit(1);
}
}

extern "C" {
void dfftw_plan_dft_c2r_2d_(void) { swd_fft_unavailable(); }
void dfftw_execute_dft_c2r_(void) { swd_fft_unavailable(); }
void dfftw_destroy_plan_(void)    { swd_fft_unavailable(); }
void fftw_cleanup(void)           { /* harmless no-op */ }
}
