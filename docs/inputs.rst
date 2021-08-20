.. inputs

Code Inputs
##############

.. _star-list:

***************************
Stellar Observation Lists
***************************

The stellar observation lists are the core source of data for the code to run. They are generated (**how are they generated?**).

Different prescriptions are possible with only minor modifications to the code, however the default assumption is that the files will follow the following prescription:

* One datafile per magnitude bin (:ref:`Nm <Nm>` datafiles in total)
* Filenames will be of the form "``n``.csv", where ``n`` runs from ``0`` to ``Nm-1``. 
* Each line of the datafile corresponds to a single star
* The lines will be comma-separated lists with the following entries::

	element 1: gaia_astrometric_detections
	element 2: predicted_gaia_visitations
	element 3: vistation_time_1
	element 4: visitation_time_2
	....
	element N, final_visitation_time
* The length of the line (i.e. ``N``) can vary from star to star, but ``N-2`` should not exceed :ref:`NumberLargerThanMaxObservations <obs-n>`.
* The lines within the datafiles are randomly shuffled

*********
Needlets
*********

**************************
Needlet-Healpix Mapping
**************************

