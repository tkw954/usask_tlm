  function [ epsilon ] = ...
    TaperedObjectiveFunction(params, omega, L, r1, r2, nu, rho, K, t11, t21, t12, t22)
%[ epsilon, t11_star, t12_star, t21_star, t22_star] = TaperedObjectiveFunctionErrors(params, omega, L, r1, r2, nu, rho, K, t11, t21, t12, t22)
%  Objective function for parameter optimization.
%
% Inputs: params = optmization parameters
%          omega = freqeuncy (rad/s)
%              L = pipeline length (m)
%             r1 = inlet radius (m)
%             r2 = outlet radius (m)
%             nu = kinematic viscosity (m^2/s)
%            rho = density (kg/m^3)
%              K = bulk modulus (Pa)
%t11,t21,t12,t22 = transmission matrix terms
%
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017


%Extract optimization parameters from parameter vector
k = (numel(params)-1)/2;%number of weighting functions
mE = params(1:k);%coefficients for transfer function G_1

mG = params((k+1):(2*k));%coefficients for transfer function G

tau = params(end);

[ E_1, E_2, F_1, F_2, G_1, G_2, Zc_1, Zc_2, T_1, T_2 ] = TaperedTLMFunctions( omega, L, r1, r2, nu, rho, K, mG, mE, tau);

[ t11_star, t12_star, t21_star, t22_star ] = TaperedTLMTransferMatrix( omega, E_1, E_2, F_1, F_2, G_1, G_2, Zc_1, Zc_2, T_1, T_2);

c=sqrt(K/rho);%(m/s) sonic speed
T=L/c;%(s) transmission time

eps_12=sum((abs((t12-t12_star)/Zc_1)).^2./(omega*T));%error in T12
eps_21=sum((abs((t21-t21_star)*Zc_1)).^2./(omega*T));%error in T21
eps_11=sum(abs(t11-t11_star).^2./(omega*T));%error in T12
eps_22=sum(abs(t22-t22_star).^2./(omega*T));%error in T21


eps_E=sum(max(0, mE(3:end)-3*mE(2:(end-1))).^2);%constraint on mE. eq(26)
eps_G=10*max(0,sum(mG)-1).^2;%constraint on mG. eq (25)

epsilon=eps_11+eps_22+eps_12+eps_21+eps_E+eps_G;%total objective function. eq (24)

end

