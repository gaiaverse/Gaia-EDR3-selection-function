.. command-arguments

#############################
Command Line Arguments
#############################

In addition to the :doc:`compile-variables`, we make use of the `JSL::Argument class <https://jackstandardlibrary.readthedocs.io/en/latest/argument.html>`_, allowing easy integration of command-line arguments for a number of properties. 

Arguments & Triggers
***************************


.. csv-table:: 
   :file: commandargs.csv
   :widths: 25 15 25 40 25
   :header-rows: 1
   
Usage
**************

Commands can be passed into the code in one of two ways: either through the commandline interface::

	> ./theia -minibatch 16 -gradlim 0.1 -harness-release 12

Or through a configuration file, written as:

.. code-block:: text
	:emphasize-lines: 1
	
	//configuration.txt
	minibatch 16
	gradlim 0.1
	data ../new/data/location
	random-seed 299

And launched as::

	> ./theia -config configuration.txt 

Note that if a configuration file is attached, all other command line arguments will be ignored. 

CommandArgs Holder Class
***************************

Once loaded, the command-line arguments are stored within a central CommandArgs object

.. doxygenclass:: CommandArgs
