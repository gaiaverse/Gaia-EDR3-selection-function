.. background

Background
############

For the full background and theory underlying this project, please see **(link to papers)**. We provide here a basic summary of the function of the code. 

The end goal of the project is to determine, given a star (real or hypothetical) with magnitude :math:`m` and position on the sky :math:`(\alpha,\delta)`, the probability that such a star would make it into the data pipeline of the `Gaia Space Telescope <https://www.gaiaverse.space/home>`_. In order to do so, we must determine as a function of time, :math:`t`, and position on the sky, :math:`\vec{\theta}`, the operating efficiency of the telescope.

We do this by maximising the Likelihood that the data within the :ref:`Stellar Data <star-list>` would be achieved given a proposed efficiency, :math:`\vec{x}`, and the appropriate priors:

.. math::
	\mathcal{L}(\vec{x} | ~\text{data}) \propto  p(\text{data} | \vec{x}) \times \text{prior}(\vec{x})

The maximal likelihood is found using a modified `ADAM optimizer <https://arxiv.org/abs/1412.6980>`_, which makes use of standard Stochastic Gradient Descent methods. 

Probability Model 
*******************

A proposed efficiency :math:`\vec{x}` contains a list of efficiency parameters, :math:`\{x_{t_i}\}`, for each of the time bins :math:`t_i`, as well as a set of spatial efficiencies, :math:`\{x_{m\ell}\}`, for each magnitude bin, :math:`m`, and `HEALPix <https://healpix.sourceforge.io/>`_ location :math:`\ell`. 

Each star is determined to have been 'visited' by Gaia at a set of time :math:`\tau_j`, which are determined using the *nominal* `Scanning Law <https://www.cosmos.esa.int/web/gaia/scanning-law>`_ and the location of the star on the sky. The probabilty that this visitation leads to a detection within the Gaia pipeline is therefore given as:

.. math::
	p(\text{detect} | \text{visited}, m, \tau_j) = \frac{1}{1 + \exp{(x_{\tau_j})}} \times \exp\left( - 0.5\log(2) \left( e^{x_{m\ell_1(\tau_j)}} + e^{x_{m\ell_2(\tau_j)}} \right) \right)

The terms :math:`\ell_1` and :math:`\ell_2` are both included because Gaia has two fields of view separated by :math:`106.5^\circ`: :math:`\ell_1(\tau_j)` is the HEALPix location being viewed by the first Field of View at the time of the visitation, whilst :math:`\ell_2(\tau_j)` is the viewing location of the second FoV. 

Given the set :math:`\{\tau_j\}` of visitations for a star, we may therefore generate a set :math:`\{p_j\}` of associated probabilities. The observational data associated with each star tells us how many times the star was entered into the astrometric pipelines: :math:`k` successful detections, in contrast to our :math:`n` visitations. 

In general, if there are :math:`n` events each with a probability :math:`p_i` of success, the total number of successes, :math:`K` follows a Poisson Binomial Distribution:

.. math::
	p(K = k | \{p_i\}, n) = \sum_{A \in F_k^n} \prod_{i \in A} p_i \prod_{j \in A^c} (1 - p_j)

Where :math:`F_k^n` is the set of all subsets of :math:`k` integers that can be picked from :math:`\{1,...,n\}` and :math:`A^c` denotes the complement of set :math:`A`.

However, whilst analytically correct, this bears two problems for us: firstly it is expensive to compute, and secondly it is extremely unforgiving: if :math:`n = k`, then the allowed probabilities for :math:`p_i` are extremely constrained - and if :math:`k > n`, the distribution becomes meaningless. This is important for us, as the nominal scanning law is imperfect, and our inferred visitation times may very well diverge strongly from the actual visitations. 

We therefore elect to use the more lenient normal approximation to the  Poisson Binomial, which allows for some additional variance. 
