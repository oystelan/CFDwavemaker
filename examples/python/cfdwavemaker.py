"""
cfdwavemaker.py
================
A thin ``ctypes`` wrapper around the CFDwavemaker shared library.

It mirrors the C API declared in ``src/CFDwavemaker.h`` (the ``extern "C"``
block) and exposes it as a small Python class with correct ``restype`` /
``argtypes`` set for every function, so you get proper double/int marshalling
and can call the library from plain Python or NumPy loops.

Usage
-----
    from cfdwavemaker import CFDwavemaker

    # 'waveinput.dat' is read from the current working directory, so point the
    # wrapper at the folder that contains it (it chdirs for you).
    wm = CFDwavemaker("../../builds/linux64/libCFDwavemaker_openmp.so",
                      workdir="my_case_folder")
    wm.initialize()

    u = wm.velocity(x=0.0, y=0.0, z=-1.0, t=5.0)     # (ux, uy, uz)
    eta = wm.surface_elevation(x=0.0, y=0.0, t=5.0)
    p = wm.dynamic_pressure(0.0, 0.0, -1.0, 5.0)

    wm.cleanup()

The library reads its configuration from ``waveinput.dat`` in the current
working directory when ``initialize()`` (``wave_Initialize``) is called.
"""

import os
from ctypes import (
    CDLL, POINTER, c_double, c_int, c_void_p,
)

# Convenience aliases for the two common signatures.
_XYZT = [c_double, c_double, c_double, c_double]   # (x, y, z, t)
_XYT = [c_double, c_double, c_double]              # (x, y, t)


class CFDwavemaker:
    """ctypes binding for the CFDwavemaker shared library.

    Parameters
    ----------
    libpath : str
        Path to the shared library (e.g. ``libCFDwavemaker_openmp.so`` on
        Linux, ``CFDwavemaker.dll`` on Windows).
    workdir : str, optional
        If given, the process changes into this directory before use so that
        ``wave_Initialize`` finds ``waveinput.dat`` there.
    """

    def __init__(self, libpath, workdir=None):
        # Resolve the library path *before* changing directory, so a relative
        # libpath is interpreted from the caller's current directory.
        libpath = os.path.abspath(libpath)
        if workdir is not None:
            os.chdir(workdir)
        self.lib = CDLL(libpath)
        self._bind()

    def _bind(self):
        """Declare restype/argtypes for every exported function."""
        lib = self.lib

        # --- lifecycle -------------------------------------------------
        lib.wave_Initialize.restype = c_int
        lib.wave_Initialize.argtypes = []
        lib.wave_Cleanup.restype = c_int
        lib.wave_Cleanup.argtypes = []
        lib.CFDwavemaker_is_initialized.restype = c_int
        lib.CFDwavemaker_is_initialized.argtypes = []

        # --- point queries (x, y, z, t) --------------------------------
        for name in ("wave_VeloX", "wave_VeloY", "wave_VeloZ",
                     "wave_DynPres", "wave_nhPres"):
            fn = getattr(lib, name)
            fn.restype = c_double
            fn.argtypes = _XYZT

        # --- surface elevation (x, y, t) -------------------------------
        lib.wave_SurfElev.restype = c_double
        lib.wave_SurfElev.argtypes = _XYT

        # --- volume fraction (x, y, z, t, delta_cellsize) --------------
        lib.wave_VFrac.restype = c_double
        lib.wave_VFrac.argtypes = [c_double] * 5

        # --- seabed (x, y) ---------------------------------------------
        lib.wave_Seabed.restype = c_double
        lib.wave_Seabed.argtypes = [c_double, c_double]

        # --- bulk kinematics: returns pointer to doubles ---------------
        # Layout is wave-type dependent (see kinematics() docstring).
        lib.wave_Kinematics.restype = POINTER(c_double)
        lib.wave_Kinematics.argtypes = _XYZT

        # --- piston paddle position (t) --------------------------------
        lib.wave_piston_position.restype = c_double
        lib.wave_piston_position.argtypes = [c_double]

        # --- grid updater / probes (t) ---------------------------------
        lib.wave_force_update.restype = None
        lib.wave_force_update.argtypes = [c_double]
        lib.update_probes.restype = None
        lib.update_probes.argtypes = [c_double]

        # --- irregular-wave helpers (opt) ------------------------------
        for name in ("wave_phase_speed", "wave_mean_period", "wave_mean_length"):
            fn = getattr(lib, name)
            fn.restype = c_double
            fn.argtypes = [c_int]

        lib.wave_water_depth.restype = c_double
        lib.wave_water_depth.argtypes = []

        # --- MPI communicator (optional; no-op unless built -DMPI_enable)
        lib.wave_set_mpi_comm.restype = None
        lib.wave_set_mpi_comm.argtypes = [c_int]

    # ------------------------------------------------------------------
    # lifecycle
    # ------------------------------------------------------------------
    def initialize(self):
        """Read ``waveinput.dat`` (from CWD) and set up the wave field."""
        return self.lib.wave_Initialize()

    def cleanup(self):
        """Free memory held by the library. Call before program exit."""
        return self.lib.wave_Cleanup()

    def is_initialized(self):
        return bool(self.lib.CFDwavemaker_is_initialized())

    def set_mpi_comm(self, fortran_comm):
        """Set the MPI communicator (Fortran handle) before initialize().

        No-op unless the library was built with ``-DMPI_enable=1``.
        """
        self.lib.wave_set_mpi_comm(int(fortran_comm))

    # ------------------------------------------------------------------
    # point queries
    # ------------------------------------------------------------------
    def velocity_x(self, x, y, z, t):
        return self.lib.wave_VeloX(x, y, z, t)

    def velocity_y(self, x, y, z, t):
        return self.lib.wave_VeloY(x, y, z, t)

    def velocity_z(self, x, y, z, t):
        return self.lib.wave_VeloZ(x, y, z, t)

    def velocity(self, x, y, z, t):
        """Return the velocity vector ``(ux, uy, uz)`` at a point."""
        return (self.lib.wave_VeloX(x, y, z, t),
                self.lib.wave_VeloY(x, y, z, t),
                self.lib.wave_VeloZ(x, y, z, t))

    def dynamic_pressure(self, x, y, z, t):
        return self.lib.wave_DynPres(x, y, z, t)

    def nonhydrostatic_pressure(self, x, y, z, t):
        return self.lib.wave_nhPres(x, y, z, t)

    def surface_elevation(self, x, y, t):
        return self.lib.wave_SurfElev(x, y, t)

    def volume_fraction(self, x, y, z, t, delta_cellsize):
        return self.lib.wave_VFrac(x, y, z, t, delta_cellsize)

    def seabed(self, x, y):
        return self.lib.wave_Seabed(x, y)

    def piston_position(self, t):
        return self.lib.wave_piston_position(t)

    def kinematics(self, x, y, z, t, n=4):
        """Bulk kinematics query (``wave_Kinematics``).

        Returns the first ``n`` doubles of the library's internal result
        array. The layout is wave-type dependent: for the lsGridSpline path
        it is ``[eta, ux, uy, uz]`` (n=4); the VTK path returns
        ``[eta, seabed, ux, uy, uz]`` (n=5). Prefer the individual
        ``velocity``/``surface_elevation`` accessors unless you specifically
        need this batched call.
        """
        ptr = self.lib.wave_Kinematics(x, y, z, t)
        return [ptr[i] for i in range(n)]

    # ------------------------------------------------------------------
    # grid / probes / irregular helpers
    # ------------------------------------------------------------------
    def force_update(self, t):
        """Advance the s-grid (lsGridSpline) kinematics to time ``t``.

        Call this once per time step, on every rank when running under MPI
        (it is a collective operation there).
        """
        self.lib.wave_force_update(t)

    def update_probes(self, t):
        self.lib.update_probes(t)

    def phase_speed(self, opt=1):
        """Phase velocity: opt=1 uses mean period T1, opt=2 uses T2 (zero-cross)."""
        return self.lib.wave_phase_speed(opt)

    def mean_period(self, opt=1):
        return self.lib.wave_mean_period(opt)

    def mean_length(self, opt=1):
        return self.lib.wave_mean_length(opt)

    def water_depth(self):
        return self.lib.wave_water_depth()


if __name__ == "__main__":
    # Minimal smoke test: load, initialize, query one point, clean up.
    import sys
    lib = sys.argv[1] if len(sys.argv) > 1 else \
        "../../builds/linux64/libCFDwavemaker_openmp.so"
    workdir = sys.argv[2] if len(sys.argv) > 2 else None

    wm = CFDwavemaker(lib, workdir=workdir)
    if wm.initialize() != 0 and not wm.is_initialized():
        raise SystemExit("wave_Initialize failed -- is waveinput.dat present?")

    print("water depth :", wm.water_depth())
    print("eta(0,0,0)  :", wm.surface_elevation(0.0, 0.0, 0.0))
    print("u(0,0,-1,0) :", wm.velocity(0.0, 0.0, -1.0, 0.0))
    wm.cleanup()
