#pragma once
#include <vector>
#include <math.h>
#include <iomanip>
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
		T & Functor;
		int Dimensions;
	public:
		bool Converged;
		Conditions Condition;
		Statuses Status;
		
		Optimizer<T>(int nRawParams, T& functor) : Functor(functor)
		
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
			Condition.StepSize = 0.02;
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
			GradientTester(x);
		}
		
			
		void BinaryWolfeDescent(VectorXd & x)
		{
			double prevF = 0;
			Functor.Calculate(x);
			bool minimiseContinues = true;
			double alpha;
			double alphaInit = 1;
			
			double c1_orig = 1e-2;
			double c2_orig = 0.9;
			

			while (minimiseContinues)
			{
				double c1 = c1_orig;
				double c2 = c2_orig;
			
				double OriginalValue = Functor.Value;
				VectorXd Grad = Functor.Gradient;
				double alpha = alphaInit;
				bool alphaFound = false;
				double gNorm = Grad.norm();
				VectorXd pk = -1*Grad; 
				double armijoValue = pk.dot(Grad);
				VectorXd dx;
				int alphaSteps = 0;
				
				
				while (!alphaFound)
				{
					dx = alpha * pk;
					VectorXd xHyp = x + dx;
					Functor.Calculate(xHyp);
					bool armijoSuccess = (Functor.Value <= OriginalValue);
					//~ bool curvatureSuccess = ( - pk.dot(Functor.Gradient) <= - c2* armijoValue);
					bool nanSuccess = ! (std::isnan(Functor.Value) || Functor.Gradient.hasNaN() );
					
					std::cout << "\t\t\tTrying alpha = " << alpha << " which gives |dx| = " << dx.norm() << " \n\t\t\t\tL = " << Functor.Value << " <=! " << OriginalValue - alpha*c1*armijoValue;
					std::cout << "\n\t\t\t\tGrad: " <<   Grad.dot(Functor.Gradient) << "<=! " << c2*armijoValue << "\n";
					
					if (armijoSuccess && nanSuccess)
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
					alphaInit *= 1.3;
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
		
		void ADAM(VectorXd &x)
		{
			
			
			//initialise ADAM vectors
			VectorXd m = VectorXd::Zero(Dimensions);
			VectorXd v = VectorXd::Zero(Dimensions);
			
			double beta1 = 0.9;
			double beta2 = 0.99;
			double eps = 1e-10;
			bool minimiseContinues = true;
			double prevF = 0;
			VectorXd ones = VectorXd::Constant(Dimensions,1.0);
			while (minimiseContinues)
			{
				//~ std::cout << x.transpose() << std::endl;
				int t = Status.CurrentSteps + 1;
				Functor.Calculate(x);
				double b1Mod = 1.0/(1.0 - pow(beta1,t));
				double b2Mod = 1.0/(1.0 - pow(beta2,t));
				m = (  beta1 *m + (1.0-beta1)*Functor.Gradient);
				
				
				VectorXd gSq = Functor.Gradient.array() * Functor.Gradient.array(); 
				v = ( beta2 * v + (1.0-beta2)* ( gSq) );
			
				VectorXd dx = b1Mod * m * Condition.StepSize;

				for (int i = 0; i < Dimensions; ++i)
				{
					dx[i] /= (sqrt(v[i]*b2Mod) + eps);
				}

				x -= dx;
				
				double df = Functor.Value - prevF;
				minimiseContinues = CheckContinues(dx,Functor.Gradient,df);
				prevF = Functor.Value;
				++Status.CurrentSteps;
				if (Status.CurrentSteps % Condition.SaveSteps == 0)
				{
					Functor.SavePosition(false);
				}
				std::cout << "\t\tStep " << Status.CurrentSteps << " Taken, at Calculation Evaluation " << Functor.LoopID << "\n";
				std::cout << "\t\t\t(L,Gradnorm,df) = (" << std::setprecision(10) << prevF << ", " <<  std::setprecision(10) <<Functor.Gradient.norm() << ", " << std::setprecision(10) <<df << ")\n"; 
				std::cout << "\t\t\t"; printTime();
				
				if (std::isnan(Functor.Value))
				{
					exit(10);
				}
			}
		}
		
		
		void GradientTester(VectorXd &x)
		{
			double dx = 1e-6;
			
			
			std::cout << "Performing a gradient test at the following location: \n\t x = " << x.transpose() << std::endl;
			Functor.Calculate(x);
			double F = Functor.Value;
			VectorXd Grad = Functor.Gradient;
			for (int i = 0; i < Dimensions;++i)
			{
				VectorXd xH = x;
				xH[i] += dx;
				Functor.Calculate(x);
				
				double measureGrad = (Functor.Value - F)/dx;
				double err = (measureGrad - Grad[i])/Grad[i];
				std::cout << "dF/dx_" << i << " = " << measureGrad << "   (theory value = " << Grad[i] << ", err = " << err << ")\n";
				
			}
			
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
		
		TestFunctor();
		TestFunctor(int n);
		int LoopID;
		void Calculate(const VectorXd & x);
		void operator () (const VectorXd & x);
		void SavePosition(bool finalSave);
	private:
		int Dimensions; 
};
