.. background

Background
############

For the full background and theory underlying this project, please see **(link to papers)**. We provide here a basic summary of the function of the code. 

The end goal of the project is to determine, given a star (real or hypothetical) with magnitude :math:`m` and position on the sky :math:`(\alpha,\delta)`, the probability that such a star would make it into the data pipeline of the `Gaia Space Telescope <https://www.gaiaverse.space/home>`_. In order to do so, we must determine as a function of time, :math:`t`, and position on the sky, :math:`\vec{\theta}`, the operating efficiency of the telescope.

We do this by maximising the Likelihood that the data within the :ref:`Stellar Data <star-list>` would be achieved given a proposed efficiency, :math:`\vec{x}`, and the appropriate priors:

.. math::
	\mathcal{L}(\vec{x} | ~\text{data}) \propto  p(\text{data} | \vec{x}) \times \text{prior}(\vec{x})

The maximal likelihood is found using a modified `ADAM optimizer <https://arxiv.org/abs/1412.6980>`_, which makes use of standard Stochastic Gradient Descent methods. 

.. _model:

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


Variance Model
===================

.. _property-spaces:

Property Spaces
*********************

Within this project, we encode our :doc:`efficiency-vector` on two different spaces, termed the ``Raw`` and ``Transformed`` spaces:

* ``Raw`` space is the version handled by the optimizer, and the space on which the actual optimization occurs. This space is often referred to as :math:`z`-space, and the variables associated with it defined accordingly.
* ``Transform`` space is the more physically/mathematically meaningful space, and the space which the Likelihood function operates within. This space is often referred to as :math:`x`-space, or :math:`p`-space.

The spaces are linked by Forward and Backward transforms. The splitting of the spaces grants us a number of advantages:

* We avoid a complex, interconnected prior by having a simple prior in ``Raw`` space
* We can enforce bounds on our parameters with appropriate transforms (i.e. :math:`x_i > 0` can be enforced by :math:`x_i = e^{x_i}`)



.. _forward-transform:

Forward Transform 
====================

The Forward Transform converts the ``Raw`` vector into the ``Transformed`` vector, such that :math:`\vec{x} = \text{ForwardTransform}(\vec{z})`.

The Forward Transform has 3 components: Temporal, Spatial and Hyper. 


.. _forward-transform-temporal:

Temporal Forward Transform 
-------------------------------

With :math:`Nt` components of both the temporal part of :math:`\vec{x}` and :math:`\vec{z}` (denoted :math:`\vec{x}^t` and :math:`\vec{z}^t` respectively), the transform is given by:

.. math::
	\begin{align}
		q_{Nt-1} & = z^t_{Nt-1}
		\\
		q_i &  = \sqrt{1 - e^{- 2/\ell_t} } \times z^t_i + e^{-1/\ell_t} q_{i+1}
		\\
		x^t_i & = \mu_t + \sigma_t \times q_i
	\end{align}
	
As the prior on :math:`\vec{z}^t` is simply the zero-mean, unit-normal Gaussian, :math:`\mu_t` and :math:`\sigma_t` are the corresponding :ref:`mean <mut>` and :ref:`standard deviations <sigmat>` of the prior on :math:`x_t`. The quantity :math:`\ell_t` is the :ref:`coupling lengthscale <lt>`, which enforces correlation between the temporal components. 

.. _forward-transform-spatial:

Spatial Forward Transform
-------------------------------

We use spherical needlets to decompose the HEALPix-mapped sky into correlated units: our Raw spatial vector, :math:`\vec{z}_{ms}`, contains a needlet-weighting for the :math:`m`-band sky map, whilst the corresponding :math:`\vec{x}_{ml}` contains the efficiency parameter for the :math:`l`-th HEALPix location of the :math:`m`-band sky map. 

They are related to each other by:

.. math::

	x_{ml, p} = \mu_p + \sigma \sum_{j = 0}^{\texttt{needlet_order}} \sum_{k = 0}^{N_j} 
	
**Need to go over this is some more detail**

.. _forward-transform-hyper:

Hyper Forward Transform 
-------------------------------

The hyperparameters associated with the coefficients of the `Variance Model`_ are unconstrained and hence unaltered by the transform:

.. math::

	x_{\text{coef}~i}^h = z_\text{coef}^h

The hyperparameters associated with the population weightings, however, are constrained by the fact that they must be :math:`x_{\text{frac}~i}^h > 0` and :math:`\sum_i x_{\text{frac}~i}^h = 1`. The transform maps the unconstraintd :math:`\vec{z}` such that:

.. math::

	x_{\text{frac}~i}^h = \frac{\exp(z_{\text{frac}~i}^h)}{\sum_i \exp(z_{\text{frac}~i}^h)}

This necessarily removes a degree of freedom, so there is an inherent degeneracy in this transform. 


.. _backward-transform:

Backward Transform
====================

The Backward Transform is **not quite** the inverse of the Forward Transform -- instead of recovering :math:`z` from :math:`x`, we recover the associated *gradients*, such that :math:`\nabla_\vec{z} \mathcal{L} = \text{BackwardTransform}(\nabla_\vec{x} \mathcal{L})`.

.. _backward-transform-temporal:

Backward Forward Transform 
-------------------------------


.. _backward-transform-spatial:

Spatial Backward Transform
-------------------------------

.. _backward-transform-hyper:

Hyper Backward Transform 
-------------------------------
