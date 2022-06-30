WaveDampToIncident
##################

Force incoming and outgoing waves towards the incident wave field by use
of a forcing-zone-approach. A penalty is included in the momentum equation
to penalize deviations from the incoming wave velocity field::

  LHS(U) == RHS - lambda * ramp * rho * alpha * (U - U_incident_wave)

This source term does not force the density (alpha) field, only U.

See, e.g., Perić and Abdel-Maksoud (2016 & 2018) for more information
about the method, references included below. They include recomendations
for damping/forcing zone sizes and the magnitude of the penalty parameter.

Usage
-----

Extra source terms in the momentum equation can be plugged into OpenFOAM
by use of the fvOptions dictionary input file. For speed and to avoid
problems with overlapping zones it is recommended to only include one
dictionary entry. An example is given below::

    TheWaveDampingZone
    {
        type            WaveDampToIncident;

        // Define one or more origins and directions
        // The damping is zero at the origin and grows to full damping
        // in the direction of the corresponding normal
        origins         ((5 0 0));
        directions      ((1 0 0));

        // The ramp is used to define the transition from no damping to
        // full damping. The start is 0 (origin) and the duration is the
        // distance to full damping (the ramps were probably first made
        // for use along the time axis, hence "duration" not "distance")
        ramp
        {
            type        halfCosineRamp;
            start       0;    // distance from origin, typically 0
            duration    5;  // this is the distance from origin to max damping
        }

        // The damping coefficient/penalty factor at maximum damping
        // this is multiplied by the difference between the computed
        // solution and the wanted solution (incident wave velocity).
        // The forcing zone definition (ramps) are also multiplied in to
        // limit the spatial influence of the damping
        penalty 10;
    }

    // This is a different zone since we use a different ramp duration
    // If the duration was the same the origin and directions here could
    // be added to the lists above
    WaveCreationZone
    {
        type            WaveDampToIncident;

        // Define one or more origins and directions
        // The damping is zero at the origin and grows to full damping
        // in the direction of the corresponding normal
        origins         ((2 0 0));
        directions      ((-1 0 0));

        // The ramp is used to define the transition from no damping to
        // full damping. The start is 0 (origin) and the duration is the
        // distance to full damping (the ramps were probably first made
        // for use along the time axis, hence "duration" not "distance")
        ramp
        {
            type        halfCosineRamp;
            start       0;    // distance from origin, typically 0
            duration    2;  // this is the distance from origin to max damping
        }

        // The damping coefficient/penalty factor at maximum damping
        // this is multiplied by the difference between the computed
        // solution and the wanted solution (incident wave velocity).
        // The forcing zone definition (ramps) are also multiplied in to
        // limit the spatial influence of the damping
        penalty 10;
    }


A rant on naming
----------------

This wave-damping method has many names. Among them is the Euler Overlay
Method, EOM, which I (Tormod) dislike since the name is misleading and
the paper introducing the method cites no-one else who has done the exact
same thing before.

We are forcing the solution towards the incident wave solution, which
can be analytical, stream function, HOSM, Stokes, ... and NOT towards a
solution from a simpler Euler solver (where the EOM method got its name).

Copyright
---------

DNV GL 2020 - Tormod Landet

References
----------

@article{peric_analytical_2018,
	title = {Analytical prediction of reflection coefficients for wave absorbing layers in flow simulations of regular free-surface waves},
	volume = {147},
	issn = {0029-8018},
	url = {http://www.sciencedirect.com/science/article/pii/S0029801817306066},
	doi = {10.1016/j.oceaneng.2017.10.009},
	urldate = {2018-02-28},
	journal = {Ocean Engineering},
	author = {Perić, Robinson and Abdel-Maksoud, Moustafa},
	year = {2018},
	pages = {132--147},
}

@article{peric_reliable_2016,
	title = {Reliable damping of free-surface waves in numerical simulations},
	volume = {63},
	issn = {0937-7255},
	url = {https://doi.org/10.1080/09377255.2015.1119921},
	doi = {10.1080/09377255.2015.1119921},
	number = {1},
	urldate = {2018-02-28},
	journal = {Ship Technology Research},
	author = {Perić, Robinson and Abdel-Maksoud, Moustafa},
	year = {2016},
	pages = {1--13},
}
