Build Instructions
==================


Build in Linux/Unix
-------------------

Start by cloning the git repository into a suitable location on you computer

.. code:: none

	git clone git@github.com:oystelan/CFDwavemaker.git

open a terminal and navigate to the location of the CFDwavemaker folder. a Makefile is provided for quick building. The following make commands are available:

Single CPU builds
.................

Make both static (.a) and shared (.so) library

.. code:: none

	make all


Multi-threading OpenMP builds
.............................

Make both static (.a) and shared (.so) library with OpenMP parallel support

.. code:: none

	make openmp

Other commands
..............

As usual,

.. code:: none

	make clean

will remove all previously combiled object files (.o) which are used when assembling the library

.. warning::

	A few words of advice:

	* Be sure to run "make clean" before running either "make openmp" or "make all"
	* Be careful/try to avoid using the OpenMP compiled library in combination with MPI. 

Build in Windows
----------------

