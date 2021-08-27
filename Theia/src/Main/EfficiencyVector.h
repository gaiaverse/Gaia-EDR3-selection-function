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
	
		/*!The internal representation of the \verbatim embed:rst:inline :math:`\vec{z}` vector \endverbatim. Elements can be accessed manually, but internal structure is only guaranteed through the Access(), Assign() and Increment() commands. */
		std::vector<double> RawPosition;
		
		/*! The internal representation of the \verbatim embed:rst:inline :math:`\vec{x}` vector \endverbatim. Elements can be accessed manually, but internal structure is only guaranteed through the Access(), Assign() and Increment() commands. */
		std::vector<double> TransformedPosition;
	
		/*! The internal representation of the \verbatim embed:rst:inline :math:`\nabla_z \mathcal{L}` vector \endverbatim. Elements can be accessed manually, but internal structure is only guaranteed through the Access(), Assign() and Increment() commands. */
		std::vector<double> RawGradient;
		
		/*! The internal representation of the \verbatim embed:rst:inline :math:`\nabla_x \mathcal{L}` vector \endverbatim. Elements can be accessed manually, but internal structure is only guaranteed through the Access(), Assign() and Increment() commands. */
		std::vector<double> TransformedGradient;
	
	
	
		/*!
		 * The first component of any Access/Assign/Increment call: determines which side of the \verbatim embed:rst:inline :ref:`Property Space <property-spaces>` \endverbatim is accessed.
		*/ 
		enum VectorMode {Raw, Transformed};
		
		/*!
		 * The second component of any Access/Assign/Increment call: determines which component of the vector is accessed.
		*/ 
		enum VectorComponent {Temporal, Spatial, Hyper};
		
		/*!
		 * The third component of any Access/Assign/Increment call: determines if you are accessing \verbatim embed:rst:inline :math:`\vec{x}` or :math:`\nabla_x \mathcal{L}` \endverbatim
		*/ 
		enum VectorType {Position, Gradient};
	
	
		/*! 
		 * Constructor function.
		 * \param load_location The location of a Raw position vector to load in (via a LoadVector() call). If set to a null value, #GenerateVector() is called instead.
		 * \param save_location The directory into which the vector is saved
		 */
		EfficiencyVector(std::string load_location,std::string save_location);
	
		//Access Functions
		
		/*!
		 * Access the specified component without needing to worry about the internal structure of the vector
		 * \param rawOrTransformed a #VectorMode enum determining which side of  \verbatim embed:rst:inline :ref:`Property Space <property-spaces>` \endverbatim is accessed
		 * \param component a #VectorComponent enum determining which of temporal/spatial/hyper you are accessing
		 * \param positionOrGradient a #VectorType enum determining if access is requested for position or gradient
		 * \param index The single-element index of the requested element
		 * \returns The requested member of the Raw/Transformed Position/Gradient
		*/
		double Access(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index);
		
		/*!
		 * Assign a value to the specified component without needing to worry about the internal structure of the vector
		 * \param rawOrTransformed a #VectorMode enum determining which side of  \verbatim embed:rst:inline :ref:`Property Space <property-spaces>` \endverbatim is accessed
		 * \param component a #VectorComponent enum determining which of temporal/spatial/hyper you are accessing
		 * \param positionOrGradient a #VectorType enum determining if access is requested for position or gradient
		 * \param index The single-element index of the requested element
		 * \param newValue The value to be assigned to the chosen element
		*/
		void Assign(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient,int index, double newValue);
		
		/*!
		 * Increments the value of a specified component without needing to worry about the internal structure of the vector
		 * \param rawOrTransformed a #VectorMode enum determining which side of  \verbatim embed:rst:inline :ref:`Property Space <property-spaces>` \endverbatim is accessed
		 * \param component a #VectorComponent enum determining which of temporal/spatial/hyper you are accessing
		 * \param positionOrGradient a #VectorType enum determining if access is requested for position or gradient
		 * \param index The single-element index of the requested element
		 * \param value The value to be added to the chosen element
		*/
		void Increment(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index, double value);
		
		/*!
		 *  An overload for Access(), but for vector access which requires two indices for easy access. 
		*/
		double Access(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index1, int index2);
		
		/*!
		 *  An overload for Assign(), but for vector access which requires two indices for easy access. 
		*/
		void Assign(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index1, int index2, double value);
		
		/*!
		 *  An overload for Increment(), but for vector access which requires two indices for easy access. 
		*/
		void Increment(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index1, int index2, double newValue);
		

			
		/*!
		 * An implementation of \verbatim embed:rst:inline :ref:`forward-transform` \endverbatim , assigns ``z`` to the #RawPosition and then transforms it into #TransformedPosition
		 * \param z The new value of #RawPosition
		*/
		void ForwardTransform(const std::vector<double> &z);
		
		/*!
		 * An implementation of \verbatim embed:rst:inline :ref:`backward-transform` \endverbatim , transforms #TransformedGradient into #RawGradient
		*/
		void BackwardTransform();
		
		/*!
		 * Saves the current value of #RawPosition and #TransformedPosition to file 
		 * \param finalSave If true, saves the files as "FinalPosition_<>.dat" in a parent directory, else saves them as "TempPosition_<>.dat" in a labelled directory. 
		 * \param saveStep The current epoch-id, if unique-saving is active, saves the Transformed vectors as TempPosition``saveStep``_<>.dat
		 * \param uniqueSave If true, allows multiple copies of the TransformedVector to be saved, instead of only the most recent one  
		*/
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
    	
    	
    	/*!
    	 * Generates a random vector 
    	*/
    	void GenerateVector();
    	void LoadVector(std::string load_location);
    	void LoadNeedlets();
    	void LoadCholesky();
    	void Reset();
    	
    	std::string SaveLocation;
    	
};
