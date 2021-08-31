.. probability-functions-poisson-exact



#################################################
Probability Functions for Exact Poisson Binomial 
#################################################

These functions are used by the :doc:`ExactPoissonContribution() call of the LogLikelihood class <likelihood>` to calculate a precise version of the Poisson Binomial probability. This version is very slow, so unless significant deviations are expeected, the :doc:`inexact versions should be used <probability-functions-poisson>`.

Forward Loop
***************

.. doxygenfunction:: poisson_binomial_lpmf_forward

Backward Loop
****************

.. doxygenfunction:: poisson_binomial_lpmf_backward


Sub-pmf Loop
****************

.. doxygenfunction:: poisson_binomial_sublpmf
