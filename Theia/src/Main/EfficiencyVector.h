#pragma once
#include <vector>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../Likelihood/ProbabilityFunctions.h"


/*!
 * 
 * The Efficiency Vector is the central object for the Theia code. It encodes a proposed operating efficiency history for Gaia, and allows for the appropriate transforms between \verbatim embed:rst:inline :ref:`property-spaces` \endverbatim . 
 * 
 */

class EfficiencyVector
{
	public:
		enum VectorMode {Raw, Transformed};
		enum VectorComponent {Temporal, Spatial, Hyper};
		enum VectorType {Position, Gradient};
	
		EfficiencyVector(std::string load_location,std::string save_location);
	
		//Access Functions
		
		double Access(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index);
		void Assign(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient,int index, double newValue);
		void Increment(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index, double value);
		double Access(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index1, int index2);
		void Assign(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index1, int index2, double value);
		void Increment(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index1, int index2, double newValue);
		

		
		
		std::vector<double> RawPosition;
		std::vector<double> TransformedPosition;
	
		std::vector<double> RawGradient;
		std::vector<double> TransformedGradient;

		std::vector<bool> KnownGapList;
		
		/*!
		 * Hello! See \verbatim embed:rst:inline :ref:`forward-transform` \endverbatim 
		*/
		
		void ForwardTransform(const std::vector<double> &z);
		void BackwardTransform();
		void Save(bool finalSave,int saveStep,bool uniqueSave);
	private:
		void ForwardTransform_Hyper();
		void ForwardTransform_Spatial();
		void ForwardTransform_Temporal();
		
		void BackwardTransform_Hyper();
		void BackwardTransform_Spatial();
		void BackwardTransform_Temporal();
		
		//Needlet stuff - has to be public
		int needletN;
		std::vector<int> needlet_u;
    	std::vector<int> needlet_v;
    	std::vector<double> needlet_w;
    	std::vector<double> bVector;
    	
    	//Cholesky stuff
    	Eigen::Matrix<double, Nm, Nm> CholeskyKg;
    	double cholesky_tol = 1e-4;
	    int choleskyN;
	    std::vector<int> cholesky_u;
    	std::vector<int> cholesky_v;
    	std::vector<double> cholesky_w; 
    	
    	
    	
    	void GenerateVector();
    	void LoadVector(std::string load_location);
    	void LoadNeedlets();
    	void LoadCholesky();
    	void Reset();
    	
    	std::string SaveLocation;
    	
};
