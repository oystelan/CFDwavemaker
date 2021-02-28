Build from source
=================


Build in Linux/Unix
-------------------

.. warning:: 

  These build instructions are somewhat outdated. Will be updated soon.

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

The CFDwavemaker library may be built and used on windows computers. A Visual studio solution file (.sln) is provided on the git hub repo for this purpose.
To get this up and running, there are at least two ways to go about in order to clone the repo and get solution up and running

Clone using Visual Studio
.........................

	1.	Launch Visual studio (community edition will do. i use v2019 at the time of writing this, but that should not matter) and select "clone a Repository" from the startup screen. This can also be found in the new tab if you already have VS up and running.
	2.	Specify the link to the repo. Generally i find the HTML tag to work best

		.. code:: none

			https://github.com/oystelan/CFDwavemaker.git

	3.	specify user name and password and all that (if required). Once you finish the setup wizard, VS should now fix the rest.

		.. warning::

			VS may require that you install git on your computer if you have not used this before.
	4.	The solution is built by pressing Ctrl+Shift+B. Be sure to select "Release" and "x64" when building.

Clone using git
...............

This is essentially the same method as the one above. except the cloning of the repo is done a bit more manually.

	1.	Start by installing git to your computer.
	2.	find a suitable folder to place the cloned git repo on your computer
	3.	start the git command prompt and navigate to your folder.
	4.	use the following command to clone the git repo

		.. code:: none

			git clone git@github.com:oystelan/CFDwavemaker.git
	5. Once done, you should find a file named "CFDwavemaker.sln". Open this in Visual studio and you should be up and running










