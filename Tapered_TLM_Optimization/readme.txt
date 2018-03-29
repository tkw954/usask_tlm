This code generates the weighting factors used for the Tapered TLM model proposed in ven der Buhs and Wiens (2017). Futher information provided in the help section of each function.

The main m-file is TableGenerate.m:
	- Set range of dissipation number (beta), and taper ratio (lambda)
	- Selects initial guess for optimization
	- Performs optimization for range of beta and lambda
	- Outputs mE, mG, and tau weighting factors for range of beta and lambda

gettlmcoeffs.m
	- Looks up parameters for uniform TLM
	- From Johnston (2014), reproduced by permission

tlm_weights.mat
	- Parameters for unifrom Enhanced TLM
	- From Johnston (2014), reproduced by permission

OptimizationsForTable.m
	- Function that performs optimization for a single lambda and beta pair.
	- Requires an initial guess.
	- Returns mE, mG, and tau for the single case.

TaperedObjectiveFunction.m
	- Objective function to be minimized.
	- Outputs optimal weighting factors.

TaperedTLMFunctions.m
	- Calculates the TLM transfer functions and parameters (E,F,G,Z_c, and T) given the weighting factors

TaperedTLMTransferMatrix.m
	- Calculates the transfer matrix (t11*,t12*,t21*, and t22*) given the TLM transfer functions

ExactSolutionforOpt.m
	- Collects the numerical solution to the system of differential equations.
	- Outputs the numerical exact tranmission matrix.
	- Can be computed in parallel if the user has a Parallel Computing Toolbox license.

t11t21venderbuhsExact.m and t12t22venderbuhsExact.m
	- Functions that calculate the numerical exact solution.
	- Outputs the tranmission matrix terms.

Q_odefun.m, Q_bcfun.m, and PQ_bcfun.m
	- Diferential equation and boundary condition functions required by the bvp4c algotithm.


References:
J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

N Johnston. Simulink Models. http://people.bath.ac.uk/ensdnj/models/newtlm.html, 2014