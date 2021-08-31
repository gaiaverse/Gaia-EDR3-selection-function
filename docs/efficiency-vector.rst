.. efficiency-vector

#####################
Efficiency Vector
#####################


   
.. doxygenclass:: EfficiencyVector
	:project: theia
	:private-members:


**********
Examples
**********

.. code:: c++
	
	#include "EfficiencyVector.h"
	
	int main()
	{
		//null loading location --> initialises a random number between initialisationBounds (+/- 0.3)
		std::string loadLocation = "__null_location__";
		std::string saveLocation = "Output";
			
		EfficiencyVector V(loadLocation, saveLocation);
		
		
		//Temporal access
		std::cout << "Examine first 5 temporal elements: " << std::endl;
		for (int i = 0; i < 5; ++i)
		{
			std::cout << "z_t[" << i << "] = " << V.Access(V.Raw, V.Temporal, V.Position, i) << std::endl;
		}
		std::cout << "\n\n";
		
		//Spatial access
		std::cout << "Spatial/Magnitude elements have two indices:" << std::endl;
		for (int j = 0; j < 3;++j)
		{
			for (int m = 0; m < 2;++m)
			{
				std::cout << "z_lm[" << j << ", " << m << "] = " << V.Access(V.Raw, V.Spatial, V.Position, j,m) << std::endl;
			}
		}
		std::cout << "\n\n";
		
		//Assignation 
		std::cout << "Assignation uses similar syntax: " << std::endl;
		for (int j = 0; j < 3;++j)
		{
			for (int m = 0; m < 2;++m)
			{
				 V.Assign(V.Raw, V.Spatial, V.Position, j,m,j + m);
				
				std::cout << "z_lm[" << j << ", " << m << "] = " << V.Access(V.Raw, V.Spatial, V.Position, j,m) << std::endl;
			}
		}
	
		
		
	}
	
Gives the following output:

.. code:: text

	Examine first 5 temporal elements: 
	z_t[0] = 0.204113
	z_t[1] = -0.0633702
	z_t[2] = 0.16986
	z_t[3] = 0.179064
	z_t[4] = 0.246988
	
	
	Spatial/Magnitude elements have two indices:
	z_lm[0, 0] = -0.192967
	z_lm[0, 1] = 0.103349
	z_lm[1, 0] = 0.102755
	z_lm[1, 1] = 0.0303147
	z_lm[2, 0] = 0.0325592
	z_lm[2, 1] = -0.0851858
	
	
	Assignation uses similar syntax: 
	z_lm[0, 0] = 0
	z_lm[0, 1] = 1
	z_lm[1, 0] = 1
	z_lm[1, 1] = 2
	z_lm[2, 0] = 2
	z_lm[2, 1] = 3


