## CFDwavemaker 3.1.0

New in this release:

- **Stream-function (Fenton) regular wave theory** (`theory_type stream`) — accurate for steep and shallow-water waves where Stokes theory breaks down.
- **VTK multilayer support**: per-layer thickness (`h`) and non-hydrostatic pressure (`phi`) from `.vts` output, with PPR or linear vertical reconstruction (`vertical_interp`). Faster VTK interpolation.
- **Optional MPI build** (`make mpi`): distribute the second-order s-grid kinematics across MPI ranks (`wave_set_mpi_comm`); results are bit-identical to the serial build.
- **SWD + lsgrid**: kinematics from a Spectral-Wave-Data file (e.g. pregenerated HOSM output) can again be sampled onto the spline interpolation grid, so CFD-side queries are served by thread-safe spline interpolation instead of serialized direct swd access. Validated against direct swd evaluation.
- Consistent `[wave reference point]` handling on the VTK path.
- **New user manual** (PDF attached below); includes a Python `ctypes` linking guide and a ready-to-use wrapper (`examples/python/cfdwavemaker.py`).
- All bundled examples updated to the current input-file format and verified to run out of the box (case2 SWD+lsgrid, case4/case5 VTK with bundled `.vts` data).
- Ramp input simplified: ramp lines are `<keyword> <start> <end>` (the `enable` column has been removed).

**Citing**: Lande, Ø., & Helmers, J. B. (2022). *CFDwavemaker: An Open-Source Library for Efficient Generation of Higher Order Wave Kinematics.* OMAE2022, V007T08A002. ASME.

### Prebuilt Linux binaries (x64)

| File | Description |
|---|---|
| `libCFDwavemaker_swd.a` / `.so` | Default build: OpenMP + SWD bundled. Link with `-lgfortran -lfftw3 -lm -fopenmp`. |
| `libCFDwavemaker_mpi.a` / `.so` | MPI build (no OpenMP): distributes lsgrid kinematics across ranks; SWD bundled. Link with `mpicxx` and `-lgfortran -lfftw3 -lm`. |
| `CFDwavemaker.h` | The C header declaring the extern "C" API. |

The static `.a` archives bundle the Spectral-Wave-Data objects, so they can be copied into your project and linked directly. Requires `libgfortran` and `libfftw3` on the system.
