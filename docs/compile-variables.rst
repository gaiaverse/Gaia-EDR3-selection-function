.. compile-variables


########################
Compile-Time Variables
########################

Defined within ``src/Main/GlobalVariables.h``.

The following variables are defined at compile-time within the file and control important user-variables which must be set at compile time for performance reasons. 


Efficiency Vector Parameters
***************************************

These parameters primarily control the number of ordering of the elements within the :doc:`Efficiency Vector <efficiency-vector>`. 

.. doxygenvariable:: Nt
.. doxygenvariable:: healpix_order	
.. doxygenvariable:: needlet_order
.. doxygenvariable:: NVariancePops
.. doxygenvariable:: hyperOrder


Priors & Lengthscales
*************************************************
.. doxygenvariable:: xtPriorNonGap
.. doxygenvariable:: xtPriorInsideGap
.. doxygenvariable:: studentNu
.. doxygenvariable:: sigmat
.. doxygenvariable:: lt_revs

.. doxygenvariable:: xmPrior
.. doxygenvariable:: lm

.. doxygenvariable:: gapPriorAlpha
.. doxygenvariable:: gapPriorPeak
.. doxygenvariable:: gapPriorBeta



Initialisation
****************************
.. doxygenvariable:: initialisationBounds
.. doxygenvariable:: xmInitialised

Data Properties
*****************************

.. doxygenvariable:: TempDirName
.. doxygenvariable:: DataLoadCount
.. doxygenvariable:: magOffset
