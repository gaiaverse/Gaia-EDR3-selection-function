.. inputs

Code Inputs
##############

.. _star-list:

***************************
Stellar Observation Lists
***************************

The stellar observation lists are the core source of data for the code to run. They are generated (**how are they generated?**).

Different prescriptions are possible with only minor modifications to the code, however the default assumption is that the files will follow the following prescription:

* All datafiles within the same directory (no subdirectories) *
* One magnitude bin per datafile (:ref:`Nm <Nm>` datafiles in total)
* Filenames will be of the form "``n``.csv", where ``n`` runs from ``0`` to ``Nm-1``. * 
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

Lines marked with a '*' are assumptions of the :doc:`assignment-creation` protocol, instead of the core Theia code, so are less stringent as long as the protocols are updated.

.. _stellar-directory:

------------------
Stellar Directory
------------------

In order to load in #DataLoadCount stars, it is useful to know a-priori how many stars are within the file. However, parsing such large text files on the fly can take an exceedingly long time, and would be annoying to run every time the system initialised. 

We therefore require that, in the same directory as the datafiles, there exists a file called ``directory.dat``, which contains a precalculated directory of the number of stars in each file, in the format of:

.. code-block:: text

	1234 file1.csv
	2927 file2.csv
	9238 file3.csv
	.
	.
	.
	19282 fileN.csv

Such a file can be easily generated with the following bash command:

.. code-block:: bash
	
	wc -l *.csv > directory.dat

.. _needlet-files:

*********
Needlets
*********

(I don't quite understand enough about how these are generated....would be useful to have help on that!)

**************************
Needlet-Healpix Mapping
**************************

(I don't quite understand enough about how these are generated....would be useful to have help on that!)


.. _gap-list:

**************************
Gap Catalogued 
**************************

Stored in gaps_prior.dat
