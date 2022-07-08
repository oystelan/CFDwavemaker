Build from source
=================


Build in Linux/Unix
-------------------

Start by cloning the git repository into a suitable location on you computer

.. code:: none

	git clone git@github.com:oystelan/CFDwavemaker.git

open a terminal and navigate to the location of the CFDwavemaker folder. a Makefile is provided for quick building. CFDwavemaker can be built in three different versions:

- **default**: contains core kinematics library + the SWD library extension.
- **basic**: The core kinematics library, excluding SWD. SWD requires a fairly new version of gcc or intel compiler on windows. Therefore it can sometimes be useful to exclude SWD from the build if you dont plan on using it.
- **all**: Complete build, linking core library, swd and the vtk library into one dynamic library. This does however require a bit of preparation, like builting the vtk library upfront, and modifying the make script accordingly to be in line with the version of vtk you have built. Dont attempt this if your not experienced with compiling. 

Right, now you build.

.. code:: none

	make clean
	make basic

If the target is omitted, the default version is built. If everything went according to plan, the resulting static and dynamic link library should now have turned up in the folder ../build/linux64/*
Thats it, in terms of building the library. Now its time to link the library to your software of choice. See.

Dependencies
............

So why did i choose to build the basic version?
The answer is, for the default version (containing the swd extension) and for the all package which links to vtk, there are some dependencies you need which does not by default ship with the compiler. 

For the default version of CFDwavemaker, the only additional dependency is fftw3. This is easily installed on a debian based system with the following command:

.. code:: none

	sudo apt-get install fftw3 libfftw3-dev
		
The **all** version on the other hand requires the installation and compilation of the vtk library. A separate section on how to do this will come.


Building the VTK library for linking with CFDwavemaker
.......................................................

.. note::

	This stepwise procedure is tested on VTK-9.1.1 and VTK-9.2.0. Obviously, its no guarantee that this will work exactly the same in future releases of VTK.

1. Download the latest source code from http://www.vtk.org
2. Extract the source code to a suitable location.
3. install Cmake (typically `sudo apt install cmake cmake-curses-gui` on a debian system)
4. go into where you extracted your source files (hint: it should contain a CMakeLists.txt file), and make a directory to build in. lets call it `./build` for now. now go into the build folder and type `ccmake ../.`. This should open the cmake build gui. first configure by hitting c. You should now have all cmake parameters available.
5. Now, we want to build the static library.  so switch BUILD_SHARED_LIBS=Off, and go into advanced mode by hiting the t-key. Under the flag CMAKE_CXX_FLAGS_RELEASE, add `-fPIC`. 
6. Finally, set an installation directory. This is not strictly required, but makes it much easier when linking later on with CFDwavemaker, since the install function gathers libraries and include files into a very organized folder. The install dir is set under the flag CMAKE_INSTALL_PREFIX, which by default is set to /usr/local. To not mix this build with other system files, lets set it to our local home directory for now, i.e. `CMAKE_INSTALL_PREFIX=~/vtk-9.2.0`. Note: The installation directory can be updated later on without having to rebuild the entire library.
7. Now, hit the g-key to generate, and when finished, exit.
8. You can now build by typing `make`. This step will take some time to complete, so be patient (go have lunch or something).
9. Now, the VTK building part is done. What remains is to set the right paths in the CFDwavemaker make files. Now, remember the path you set in step 5? 



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

Its possible to do, but rarely need to. 










