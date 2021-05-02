MAKE_DIR = $(PWD)

.PHONY: clean openmp

clean:
	@$(MAKE) -f make_openmp.mk clean
	@$(MAKE) -f make_swd.mk clean

openmp:
	@$(MAKE) -f make_openmp.mk

swd:
	@$(MAKE) -f make_swd.mk
