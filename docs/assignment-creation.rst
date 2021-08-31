.. assignment-creation

##########################
File-Core Allocation
##########################

The File-Core Allocation protocol is a separate, external script which is called by :ref:`LoadData() <load-data>` function. It is deliberately kept separate to allow for the possibility of differing protocols regarding the generation of the :ref:`star-list`.

The file is kept within a python script, ``starAllocation.py``. Replacing or modifying this code will alter the behaviour of the Theia code without requiring a recompilation.

The purpose of this code is to determine which of the (potentially wildly differing in sizes) of the :ref:`star-list` files should be assigned to which of the :doc:`MPI workers <parallelisation>` running the code. 

Input
--------------

The protocol needs to know the location of the data, and the number of cores over which to distribute the files.


The script should be callable from the command line using the following syntax:

.. code-block:: text

	python3 path/to/starAllocation.py root/to/Data NCores
	
	//for example:
	
	python3 src/DataLoading/starAllocation.py ../Data/ 20


Internal Logic 
--------------------

The assignment is balanced so that the memory and computational balance is roughly equal between the workers, and hence avoid wasted cycles whilst workers wait for their compatriots to catch up. 

This is achieved by using the on-disc file size as a metric for the eventual memory footprint (a poor measure, but better than the alternatives), and then assigning the largest remaining file to the worker with the smallest workload, updating the workload count and interating until all the files have been used. The worker with the resulting smallest workload is assigned to the rootID, as the root must perform other tasks. 

The system also should have some way of assigning each file to a magnitude bin (1 bin per file, but multiple files can share the same bin). In the default assumption this is trivial as the files names are ``N.csv``, where ``N`` is the magnitude bin the stars belong to. 

Output
--------------

The output is a file ``coreAssignments.dat``, which should be placed in the ``ModelInputs`` directory.

This file should have one line per worker (including the root). The first element of each line is the ID of the worker that line belongs to. Each subsequent element should be a list of pairs: a file name, following by the integer representing the magnitude bin of that file. There is no restriction on line length. 

An example of the output of the default script, for 4 cores splitting up 0.csv -> 10.csv:

.. code-block:: text

	0,5.csv,5,3.csv,3
	1,0.csv,0,1.csv,1
	2,7.csv,7,2.csv,2,6.csv,6
	3,8.csv,8,9.csv,9,4.csv,4

