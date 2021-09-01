.. optimiser


###########
Optimiser
###########


Implementation Details
-----------------------------

.. toctree::
	optimiser-implementation
	optimiser-data
	optimizer-misc
	likelihood-functor
	:maxdepth: 4
	:caption: Members:
	

.. _adabadam-theory:

ADABADAM Theory
-----------------

The ADABADAM (**Ada**\ ptive **B**\ atched, **Ada**\ ptive **M**\ oment) Optimizer is based on the `ADAM Optimizer of Kingma and Ba <https://arxiv.org/abs/1412.6980>`_, and is designed as a batched SGD optimizer. 



.. _minibatching:

Minibatching
+++++++++++++++++++

Minibatching can drastically speed up calculations wherein a function calculates a value given a parameter :math:`\vec{x}` and a set of data (in our case, the :ref:`star-list`). 

The function can also act on a subset of that data - if the data is randomly shuffled, then this value is approximately representative of the full function. By making changes to :math:`\vec{x}` after calculating every small chunk (rather than waiting for all the data to be processed), we can end up with much more rapid convergence. 

A calculation on a small chunk is called a *minibatch*, whilst the set of minibatches which cover the whole dataset is called an *epoch*. 




Optimization Step
+++++++++++++++++++++++

At each minibatch step, we use the usual ADAM optimisation scheme:

.. math::
	
	\vec{G}_t & = \nabla F(\vec{x}_{t -1})
	\\
	\vec{m}_t & = \beta_1 \vec{m}_{t -1} + (1 - \beta_1)\vec{G}_t
	\\
	\vec{v}_t & = \beta_2 \vec{g}_{t-1} + (1 - \beta_2) (\vec{G}_t \otimes \vec{G}_t)
	\\
	\hat{m}_t & = \frac{\vec{m}_t}{1 - \beta_1^t} ~~~~~~~~~~ \hat{v}_t  = \frac{\vec{v}_t}{1 - \beta_2^t}
	\\
	[\vec{x}_{t}]_i & = [\vec{x}_{t - 1}]_i  - \frac{[\hat{m}_t]_i}{\sqrt{[\hat{v}_t]_i} + \epsilon} \times \mathcal{R}

Where :math:`F` is the function being optimised, :math:`\vec{v}\otimes\vec{w}` is the Hadamard pointwise product, :math:`\beta_1` and :math:`\beta_2` are the 'ADAM parameters' (we choose 0.7 and 0.99, based on empirical testing), and :math:`\epsilon` is a small offset to prevent dividing by numbers close to 0. 


Optimizing Rate
+++++++++++++++++++++

The parameter :math:`\mathcal{R}` is the current Optimizing Rate - which we modify from the usual constant of the ADAM optimizer. The optimizing rate is the product of the LearningRate parameter and the Harness parameter.

The Harness is usually set to 1 (and hence has no effect). However, when a discontinuous event occurs within the optimiser - such as a minibatch downstep, or the initialisation of the code - we set the harness to a small number (usually 1/10), thereby slowing the progress of the code and limiting the impact of the discontinuity. This enables the code to then smoothly learn how the surface has changed and adapt without large jumps which take a long time to recover. The Harness returns to 1 over a set number of epochs. 

The LearningRate is nominally set to the step size of the optimizer, however, it can change on two conditions:

1. If, over the course of an epoch, the net change in the function is positive (I.e. minimization did not occur), the LearningRate is halved.
2. If, over the course of an epoch, the net change in the function is positive, the LearningRate increases by 3%. This cannot increase the LearningRate above the original stepsize. 

The rationale for this is that, if the net change is positive across an epoch, the function must be close to the optimum, and therefore needs to take smaller steps to find progressively better estimates. On the other hand, if the function is consistently taking negative (i.e. good) steps, then it can increase its step size to speed up the progress. 
 


