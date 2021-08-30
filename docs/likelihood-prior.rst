.. likelihood-prior

###########################
LogLikelihoodPrior
###########################


.. doxygenclass:: LogLikelihoodPrior
	:private-members:


Prior Functions
******************

Because the priors are often very small, repeatedble sections of code which we might want to swap out trivially, there are several small snippets designed to allow this. 

Holder Struct
---------------------

.. doxygenstruct:: F_dF


Predefined Priors
-----------------------

.. doxygenfunction:: Normal

.. doxygenfunction:: StudentT

.. doxygenfunction:: TemporalBetaPrior
