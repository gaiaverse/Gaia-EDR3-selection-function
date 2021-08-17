.. background

Background
############

For the full background and theory underlying this project, please see **(link to papers)**. We provide here a basic summary of the function of the code. 

The purpose of this code is to infer as a function of time :math:`t` and position on the sky :math:`(\alpha,\delta)` the operating efficiency of the `Gaia Space Telescope <https://www.gaiaverse.space/home>`_, where the operating efficiency is defined as the probability that a star at that location on the sky entered into the Astrometric pipeline.

We do this by maximising the Likelihood that the data within the :ref:`Stellar Data <star-list>` would be achieved given a proposed efficiency vector :math:`\vec{x}`, and the appropriate priors:

.. math::
	\mathcal{L}(\vec{x} | ~\text{data}) \propto  p(\text{data} | \vec{x}) \times \text{prior}(\vec{x})

The maximal likelihood is found using a modified `ADAM optimizer <https://arxiv.org/abs/1412.6980>`_, 
