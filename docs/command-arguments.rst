.. command-arguments

#############################
Command Line Arguments
#############################

In addition to the :doc:`compile-variables`, we make use of the `JSL::Argument class <https://jackstandardlibrary.readthedocs.io/en/latest/argument.html>`_, allowing easy integration of command-line arguments for a number of properties. 

Arguments & Triggers
***************************

.. code-block:: c++

	Argument<std::string> StartVectorLocation = Argument<std::string>("__null_location__","restart");



CommandArgs Holder Class
***************************

.. doxygenclass:: CommandArgs
