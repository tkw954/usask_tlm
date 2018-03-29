function [ mE,mG,tau ] = OptimizationsForTable(lambda, beta, params0)
%[ mE, mG, tau ] = OptimizationsForTable(lambda, beta, params0)
% Inputs: lambda  = Taper ratio (dimensionless)
%         beta    = dissipation number (dimensionless)
%         params0 = Initial guess parameter vector of the form [ mE0, mG0, tau0 ]
%
% Outputs: Optimized weighting factors
%
% Note: Optimized parameters are not sensitive to the assumed radii and
% fluid properties since this model is linear. The optimal weighting
% factors only depend on dissipation number (beta) and taper ratio
% (lambda).
%
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

r1 = 20e-3; %Assume some fixed radius value (m)
r2 = r1*lambda; %(m)
nu=100e-6;%(m^2/s) kinematic viscosity
K=1.5e9;%(Pa) bulk modulus
rho=890;%(kg/m^3) density
c=sqrt(K/rho);%(m/s) sonic speed
lambda = min(r1,r2)/max(r1,r2); % Taper ratio (dimensionless)

% Assign length, L (m), based off of given dissipation number
L = (beta*c*max(r1,r2)^2/nu)*((9*lambda^3)/((lambda^2+lambda+1)^2));

T=L/c;% Wave propagation time (s)

k=6; %number of weighting function terms
n=nan(1,k);%weighting function coefficient
n(1)=0.3/(1+3*beta);% Equation (19)
for i=2:k 
    n(i)=n(i-1)*3;% Equation (19)
end
N_per_decade=50;% Number of frequency points per decade
omegaT_min=0.01;% Minimum omega*T for frequency points

%Log space between omegaT_min and n(end)
omegaT=logspace(log10(omegaT_min),log10(n(end)),round(N_per_decade*(log10(n(end)/T)-log10(omegaT_min/T))));

omega=omegaT/T;%(rad/s) frequency

% Calculate exact transmission matrix using boundary value solver.
[ t11, t21, t12, t22 ] = ExactSolutionforOpt( omega, L, r1, r2, nu, rho, K );

fcn_min=@(params) TaperedObjectiveFunction( params, omega, L, r1, r2, nu, rho, K,t11,t21,t12,t22);%objective function to minimize
options = optimset('MaxFunEvals',500000,'MaxIter',500000);

% Ensure parameters are positive or 0 (i.e. set lower bound to 0)
lb = zeros(1,numel(params0));

% Perform optimization
params=fmincon(fcn_min,params0,[],[],[],[],lb,[],[],options);%optimize

mE = params(1:k);
mG = params((k+1):(2*k));
tau = params(end);