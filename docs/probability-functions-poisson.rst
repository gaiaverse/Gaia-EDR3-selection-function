.. probability-functions-poisson



###############################################
Probability Functions for Poisson Binomial 
###############################################

These functions are used by the :doc:`PoissonContribution() call of the LogLikelihood class <likelihood>` to calculate a slightly truncated / rounded version of the Poisson Binomial probability. If these functions misbehave (as they are liable to do when a member of {p_i} is close to 0 or 1), then the :doc:`exact versions should be used <probability-functions-poisson-exact>`.

Forward Loop
*********************

.. doxygenfunction:: poisson_binomial_pmf_forward

Backward Loop
*********************

.. doxygenfunction:: poisson_binomial_pmf_backward

Sub-pmf Loop
*********************

.. doxygenfunction:: poisson_binomial_subpmf
