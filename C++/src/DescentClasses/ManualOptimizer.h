#pragma once
#include <vector>
#include <math.h>
#include "../GenericFunctions/timeCodes.h"
#include "Star.h"
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
	int SaveSteps;
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
		
		Optimizer<T>(int n,const std::vector<Star> & data, int nTransformParams, int nRawParams,std::string outdir, int nStars) : Functor(n,data, nTransformParams,outdir,nStars)
		{	
			Dimensions = nRawParams;
			SetDefaults();	
		}
		
		
		
		void SetDefaults()
		{
			Converged = false;
			Condition.MaxSteps = 10000;
			Condition.xConvergence = 0;
			Condition.gConvergence = 1e-7;
			Condition.fConvergence = 0;
			Condition.StepSize = 0.01;
			Condition.SaveSteps = 5;
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
			
			double c1_orig = 1e-3;
			double c2_orig = 0.999;
			

			Functor.Calculate(x);
			while (minimiseContinues)
			{
				double c1 = c1_orig;
				double c2 = c2_orig;
			
				double OriginalValue = Functor.Value;
				VectorXd Grad = Functor.Gradient;
				double alpha = alphaInit;
				bool alphaFound = false;
				double gNorm = Grad.norm();
				VectorXd pk = -1*Grad/gNorm; 
				double armijoValue = pk.dot(Grad);
				VectorXd dx;
				int alphaSteps = 0;
				
				
				while (!alphaFound)
				{
					dx = alpha * pk;
					VectorXd xHyp = x + dx;
					Functor.Calculate(xHyp);
					bool armijoSuccess = (Functor.Value <= OriginalValue + alpha*c1*armijoValue);
					bool curvatureSuccess = ( - pk.dot(Functor.Gradient) <= - c2* armijoValue);
					bool nanSuccess = ! (std::isnan(Functor.Value) || Functor.Gradient.hasNaN() );
					
					std::cout << "\t\t\tTrying alpha = " << alpha << " which gives |dx| = " << dx.norm() << " \n\t\t\t\tL = " << Functor.Value << " <=! " << OriginalValue - alpha*c1*armijoValue;
					std::cout << "\n\t\t\t\tGrad: " <<   Grad.dot(Functor.Gradient) << "<=! " << c2*armijoValue << "\n";
					
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
					
					if (alphaSteps > 10)
					{
						alphaSteps = 0;
						std::cout << "Reducing convergence conditions" << std::endl;
						c1 = c1*0.5;
						c2 = c2 * 1.1;
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
				if (Status.CurrentSteps % Condition.SaveSteps == 0)
				{
					Functor.SavePosition(false);
					
				}
				std::cout << "\t\tStep " << Status.CurrentSteps << " Taken, at Calculation Evaluation " << Functor.LoopID << " (L,Gradnorm,df) = (" <<prevF << ", " <<  Functor.Gradient.norm() << ", " << df << ")\n\t\t\t"; printTime();
			}
			
			Functor.SavePosition(true);
		}
		
		
				
		bool CheckContinues(const VectorXd & dx,const VectorXd & dg, double df)
		{
		
			if (Status.CurrentSteps > Condition.MaxSteps)
			{
				//~ std::cout << "STEPS" << std::endl;
				Status.TooManySteps = true;
				return false;
			}
			if (Condition.xConvergence > 0 && dx.norm() < Condition.xConvergence)
			{
				//~ std::cout << "X" << std::endl;
				Status.ReachedStepConvergence = true;
				Converged = true;
				return false;
			}
			if (Condition.gConvergence > 0 && dg.norm() < Condition.gConvergence)
			{
				//~ std::cout << "G" << std::endl;
				Status.ReachedGradConvergence = true;
				Converged = true;
				return false;
			}
			if (Condition.fConvergence > 0 && abs(df) < Condition.fConvergence)
			{
				Status.ReachedFunctionConvergence = true;
				Converged = true;
				return false;
			}

			
			return true;
		} 

		std::string GetStatus()
		{
			std::string s = "";
			s += "Steps Taken: " + std::to_string(Status.CurrentSteps) + " / " + std::to_string(Condition.MaxSteps);
			s += "\nHalt conditions: ";
			std::vector<std::string> titles = {"Too many steps", "Reached Gradient Convergence", "Reached Step Convergence", "Reached Functional Convergence"};
			std::vector<bool> values = {Status.TooManySteps, Status.ReachedGradConvergence, Status.ReachedStepConvergence, Status.ReachedFunctionConvergence};
			for (int i = 0; i < titles.size(); ++i)
			{
				s +=  "\n\t" + titles[i] + ": ";
				std::string answer = "no";
				if (values[i])
				{
					answer = "yes";
				} 
				s += answer;
			}
			s+= "\n";
			return s;
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
