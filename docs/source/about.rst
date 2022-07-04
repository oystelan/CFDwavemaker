Introduction
============

About
-----

CFDwavemaker is a C++ library developed with the soule purpose of providing wave kinematics data for

* initializing a CFD domain with millions of cells in a fast and efficient manner
* work as a wave maker, prescribing kinematics at the boundaries of the CFD domain during runtime

The specific focus of this code is on higher order irregular wave theories for deep and intermediate water depths. A few showcases illustrating what can be achieved using this code in combination with some state of the art CFD codes are shown below.


Showcases
---------

Wave Flume
..........

*to be completed soon*


Short crested wave example
..........................

CFDwavemaker and `Basilisk`_ Navier Stokes solver. Large periodic CFD domain initialized with second order wave kinematics
	
.. _`Basilisk`: http://basilisk.fr

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/1KRlpboGX-A" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


Contribute?
-----------
Constructive feedback or suggestions of improvements are always welcome. If you have codes or not yet implemented wave theory which can fit into this library, `let me know about it <mailto:oystelan@gmail.com>`_


Bug reporting
-------------
Bug in the code you say? `let me know about it <mailto:oystelan@gmail.com>`_. 

Reference to CFDwavemaker
-------------------------

The code is open source and free to use by anyone. If you find the code useful and decide to used it in projects or publications, make sure to reference to

.. code-block:: none

   @inproceedings{landeCFDwavemaker,
      title={CFDWAVEMAKER: An open-source library for efficient generation of higher order wave kinematics},
      author={Lande, Oystein and Helmers, Jens Bloch},
      booktitle={International Conference on Offshore Mechanics and Arctic Engineering},
      volume={57656},
      pages={V03AT02A029},
      year={2022},
      organization={American Society of Mechanical Engineers}
   }
