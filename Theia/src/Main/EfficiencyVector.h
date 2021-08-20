#pragma once
#include <vector>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../Likelihood/ProbabilityFunctions.h"




class EfficiencyVector
{
	public:
		enum VectorMode {Raw, Transformed};
		enum VectorComponent {Temporal, Spatial, Hyper};
	
	
		EfficiencyVector(std::string load_location);
	
		//Access Functions
		
		double Access(int index, VectorComponent component, VectorMode mode);
		void Assign(int index, VectorComponent, VectorMode mode, double newValue);
		void Increment(int index, VectorComponent, VectorMode mode, double value);
		double Access(int index1, int index2, VectorComponent component, VectorMode mode);
		void Assign(int index1, int index2, VectorComponent, VectorMode mode, double value);
		void Increment(int index1, int index2, VectorComponent, VectorMode mode, double newValue);
		
		
		std::vector<double> RawPosition;
		std::vector<double> TransformedPosition;
	
		std::vector<double> RawGradient;
		std::vector<double> TransformedGradient;

		std::vector<bool> KnownGapList;
		void ForwardTransform(const VectorXd &z);
		void BackwardTransform();
	private:
		//hyperbuffer?
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
};
