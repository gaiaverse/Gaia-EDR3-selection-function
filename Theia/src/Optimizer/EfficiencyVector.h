#pragma once
#include <vector>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../Likelihood/MiscFunctions.h"
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
		 * Constructor function for the master version of the EfficiencyVector
		 * \param load_location The location of a Raw position vector to load in (via a LoadVector() call). If set to a null value, #GenerateVector() is called instead.
		 * \param save_location The directory into which the vector is saved
		 */
		EfficiencyVector(std::string load_location,std::string save_location);
	
		/*!
		 * Constructor for the basic version of the EfficiencyVector - when initialised in this way it has only the ability to index appropriately into the TransformedVector -- everything else is uninitiliased memory.
		 */
		EfficiencyVector(std::vector<double> newTransformed);
		
		//! Default constructor (does nothing except sit there like a lemon)
		EfficiencyVector();
		
		
		/*!
		 * Access the specified component without needing to worry about the internal structure of the vector
		 * \param rawOrTransformed a #VectorMode enum determining which side of  \verbatim embed:rst:inline :ref:`Property Space <property-spaces>` \endverbatim is accessed
		 * \param component a #VectorComponent enum determining which of temporal/spatial/hyper you are accessing
		 * \param positionOrGradient a #VectorType enum determining if access is requested for position or gradient
		 * \param index The single-element index of the requested element
		 * \returns The requested member of the Raw/Transformed Position/Gradient
		*/
		
		double Access(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index) const;
		
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
		double Access(VectorMode rawOrTransformed, VectorComponent component, VectorType positionOrGradient, int index1, int index2) const;
		
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
	
		/*! Executes the Temporal part of the ForwardTransform(). See \verbatim embed:rst:inline :ref:`Temporal Forward Transform <forward-transform-temporal>` \endverbatim for implementation details. */
		void ForwardTransform_Temporal();
		
		/*! Executes the Spatial part of the ForwardTransform(). See \verbatim embed:rst:inline :ref:`Spatial Forward Transform <forward-transform-spatial>` \endverbatim for implementation details. */
		void ForwardTransform_Spatial();
		
		/*! Executes the Hyperparameter part of the ForwardTransform(). See \verbatim embed:rst:inline :ref:`Hyperparameter Forward Transform <forward-transform-hyper>` \endverbatim for implementation details. */
		void ForwardTransform_Hyper();
		
		/*! Executes the Temporal part of the BackwardTransform(). See \verbatim embed:rst:inline :ref:`Temporal Backward Transform <backward-transform-temporal>` \endverbatim for implementation details. */
		void BackwardTransform_Hyper();
		
		/*! Executes the Spatial part of the BackwardTransform(). See \verbatim embed:rst:inline :ref:`Spatial Backward Transform <backward-transform-spatial>` \endverbatim for implementation details. */
		void BackwardTransform_Spatial();
		
		/*! Executes the Hyperparameter part of the BackwardTransform(). See \verbatim embed:rst:inline :ref:`Hyperparameter Backward Transform <backward-transform-hyper>` \endverbatim for implementation details. */
		void BackwardTransform_Temporal();
		
		//! The number of needlet elements - unknown at run time as it depends on the spatial resolution chosen, and the chosen zero-rounding cutoff
		int needletN;
		
		//! Non-zero i-indices in the needlet matrix 
		std::vector<int> needlet_u;
    	
    	//! Non-zero j-indices in the needlet matrix
    	std::vector<int> needlet_v;
    	
    	/*! The values of P_ij, with ij determined from #needlet_u and #needlet_v */  
    	std::vector<double> needlet_w;
    	
    	//! A temporary holder (to prevent needless memory initialisations) for the spatial transform
    	std::vector<double> bVector;
    	
    	//! The #BackwardTransform() requires the use of a cholesky-decomposed matrix for the magnitude-coupling (encoded by #lm). This required many cycles over elements which were near-zero. We have decomposed this into a list of #choleskyN elements which are significant, to cut the runtime of this expensive operation.
	    int choleskyN;
	    
	    //! Non-zero i-indices of the cholesky-matrix Kg
	    std::vector<int> cholesky_u;
	    
	     //! Non-zero j-indices of the cholesky-matrix Kg
    	std::vector<int> cholesky_v;
    	
    	 //! The value of CholeskyKg_ij, with ij determined from #cholesky_u and #cholesky_v
    	std::vector<double> cholesky_w; 
    	
    	/*!
    	 * Generates a random Raw vector. All elements are initialised to +/- #initialisationBounds. The zeroth-order spatial modes have an additional offset equal to #xmInitialised.
    	*/
    	void GenerateVector();
    	
    	/*!
    	 * Searches for a valid savefile given by the input. If it is of the correct length, loads it in as the #RawPosition, otherwise throws an error
    	 * \param load_location The name of the savefile to attempt to load`  
    	*/
    	void LoadVector(std::string load_location);
    	
    	/*!
    	 * Searches for a valid needlet file and loads the contents into #needlet_u, #needlet_v and #needlet_w. See  \verbatim embed:rst:inline :ref:`Needlet Files <needlet-files>` \endverbatim for the required file properties. Run only at object initialisation.
    	*/
    	void LoadNeedlets();
    	
    	/*!
    	 * Generates a matrix CholeskyKg used for magnitude-correlations, then decomposes it into #choleskyN non-zero elements encodes by #cholesky_w. Run only at object initialisation.
    	 */
    	void LoadCholesky();
    	
    	/*!
    	 * Sets the #RawPosition, #RawGradient, #TransformedPosition and #TransformedGradient to empty vectors -- required because the transforms are additive so need a zero-base to work from (rather than simply overlaying the old values with the new ones) 
    	*/
    	void Reset();
    	
    	//!The directory where the #Save() function places its output
    	std::string SaveLocation;
    	
};
