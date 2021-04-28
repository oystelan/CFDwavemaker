MAKE_DIR = $(PWD)

all:
	@$(MAKE) -f make_complete.mk

.PHONY: clean all

clean:
	@$(MAKE) -f make_openmp.mk clean
	@$(MAKE) -f make_swd.mk clean
	@$(MAKE) -f make_vtk.mk clean

openmp:
	@$(MAKE) -f make_openmp.mk

swd:
	@$(MAKE) -f make_swd.mk

vtk:
	@$(MAKE) -f make_vtk.mk

vtk_cluster:
	@$(MAKE) -f make_vtk_cluster.mk
