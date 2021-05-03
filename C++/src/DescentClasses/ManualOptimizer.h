#pragma once
#include <vector>
#include <math.h>

#include <iostream>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"
using Eigen::VectorXd;


struct Conditions
{
	int MaxSteps;
	double StepSize;
	double xConvergence;
	double gConvergence;
	double fConvergence;
};
struct Statuses
{
	int CurrentSteps;
	bool TooManySteps;
	bool ReachedGradConvergence;
	bool ReachedStepConvergence;
	bool ReachedFunctionConvergence;
};

template<class T>
class Optimizer
{
	private:
		T Functor;
		int Dimensions;
	public:
		bool Converged;
		Conditions Condition;
		Statuses Status;
		
		Optimizer<T>(int n) : Functor(n)
		{	
			Dimensions = n;
			SetDefaults();	
		}
		
		
		
		void SetDefaults()
		{
			Converged = false;
			Condition.MaxSteps = 10000;
			Condition.xConvergence = 0;
			Condition.gConvergence = 1e-7;
			Condition.fConvergence = 0;
			Condition.StepSize = 1;
			Status.CurrentSteps = 0;
			Status.TooManySteps = false;
			Status.ReachedGradConvergence = false;
			Status.ReachedStepConvergence = false;
			Status.ReachedFunctionConvergence = false;
		}
		
		void Minimize(VectorXd & x)
		{
			if (x.size() != Dimensions)
			{
				std::cout << "OPTIMIZER ERROR: Initial position vector is not of the provided size." << std::endl;
				exit(2);
			}
			BinaryWolfeDescent(x);
		}
		
		void DumbDescent(VectorXd & x)
		{
			double prevF = 0;
			Functor(x);
			bool minimiseContinues = true;
			double alpha = Condition.StepSize;
			while (minimiseContinues)
			{
				Functor(x);
				
				VectorXd dx = alpha * Functor.Gradient;
				if(dx.hasNaN())
				{
					alpha *= 0.1;
				}
				else
				{
					x -= dx;
				}
	
				++Status.CurrentSteps;				
				double df = Functor.Value - prevF;
				minimiseContinues = CheckContinues(dx,Functor.Gradient,df);
				prevF = Functor.Value;
			}
		}
		
		void BinaryWolfeDescent(VectorXd & x)
		{
			double prevF = 0;
			Functor(x);
			bool minimiseContinues = true;
			double alpha;
			double alphaInit = 1;
			double c1 = 1e-4;
			double c2 = 0.9;

			Functor(x);
			while (minimiseContinues)
			{
				
				double OriginalValue = Functor.Value;
				VectorXd Grad = Functor.Gradient;
				double alpha = alphaInit;
				bool alphaFound = false;
				double armijoValue = Grad.norm();
				armijoValue *= armijoValue;
				VectorXd dx;
				int alphaSteps = 0;
				
				
				while (!alphaFound)
				{
					dx = alpha * Grad;
					VectorXd xHyp = x - dx;
					Functor(xHyp);
					bool armijoSuccess = (Functor.Value <= OriginalValue - alpha*c1*armijoValue);
					bool curvatureSuccess = ( Grad.dot(Functor.Gradient) <= c2* armijoValue);
					bool nanSuccess = ! (std::isnan(Functor.Value) || Functor.Gradient.hasNaN() );
					
					if (armijoSuccess && curvatureSuccess && nanSuccess)
					{
						alphaFound = true;
						x = xHyp;
					}
					else
					{
						//~ std::cout << Grad.transpose() << std::endl;
						++alphaSteps;
						alpha = alpha*0.5;
					}
					
					if (alphaSteps > 100)
					{
						std::cout << "To many alpha steps. Dead in the water" << std::endl;
						exit(2);
					}
				}
				if (alpha == alphaInit)
				{
					alphaInit *= 2;
				}
				else
				{
					alphaInit = alpha;
				}
				++Status.CurrentSteps;				
				double df = Functor.Value - prevF;
				minimiseContinues = CheckContinues(dx,Functor.Gradient,df);
				prevF = Functor.Value;
			}
			
		}
		
		
				
		bool CheckContinues(const VectorXd & dx,const VectorXd & dg, double df)
		{
		
			if (Status.CurrentSteps > Condition.MaxSteps)
			{
				return false;
			}
			if (Condition.xConvergence > 0 && dx.norm() < Condition.xConvergence)
			{
				Converged = true;
				return false;
			}
			if (Condition.gConvergence > 0 && dg.norm() < Condition.gConvergence)
			{
				Converged = true;
				return false;
			}
			if (Condition.fConvergence > 0 && abs(df) < Condition.fConvergence)
			{
				Converged = true;
				return false;
			}

			
			return true;
		} 

};



class TestFunctor
{
	public:
		double Value;
		VectorXd Gradient;
		
		
		TestFunctor(int n);
		
		void Calculate(const VectorXd & x);
		void operator () (const VectorXd & x);
	private:
		int Dimensions; 
};
