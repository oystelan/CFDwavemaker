MAKE_DIR = $(PWD)

.DEFAULT_GOAL := default
.PHONY: default

clean:
	@$(MAKE) -f make_openmp.mk clean
	@$(MAKE) -f make_openmp_vtk.mk clean


# build CFDwavemaker_mpi (base + MPI, distributes lsgrid kinematics across ranks)
mpi:
	@$(MAKE) -f make_mpi.mk clean
	@$(MAKE) -f make_mpi.mk

# build CFDwavemaker_swd (base + swd, openmp) -- the default
default:
	@$(MAKE) -f make_openmp.mk clean
	@$(MAKE) -f make_openmp.mk

# build entire CFDwavemaker + vtk (base + swd + vtk)
vtk:
	@$(MAKE) -f make_openmp_vtk.mk clean
	@$(MAKE) -f make_openmp_vtk.mk
