function [ E_1, E_2, F_1, F_2, G_1, G_2, Zc_1, Zc_2, T_1, T_2 ] = TaperedTLMFunctions( omega, L, r1, r2, nu, rho, K, mG, mE, tau)
%[ E_1, E_2, F_1, F_2, G_1, G_2, Zc_1, Zc_2, T_1, T_2 ] = TaperedTLMFunctions( omega, L, r1, r2, nu, rho, K, mG, mE, tau)
%Calculates the E, F, and G transfer functions for the tapered TLM
%
% Inputs:  omega = freqeuncy (rad/s)
%              L = pipeline length (m)
%             r1 = inlet radius (m)
%             r2 = outlet radius (m)
%             nu = kinematic viscosity (m^2/s)
%            rho = density (kg/m^3)
%              K = bulk modulus (Pa)
%      mG,mE,tau = parameters
%
% Outputs: Transfer functions, characteristic impedance, and wave
%          propagation times.
% NOTE: Outputs are separated to allow investigation into using separate
% weighiting factors, propagation times, and impedance values.
% 
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017


c = sqrt(K/rho);
lambda = min(r1,r2)/max(r1,r2); % Taper ratio (dimensionless). eq(11)
beta = ((nu*L)/(c*max(r1,r2)^2))*((lambda^2+lambda+1)^2/(9*lambda^3)); % Dimensonless dissipation number. eq(13)
Zc = ((3*c*rho)/(pi*max(r1,r2)^2))/(lambda^2+lambda+1); % Characteristic impedance. eq(10)

T=L/c;%(s) Nominal propagation time
T_1=T*tau;%(s) adjusted transmission time forward propagation. eq(18)
T_2=T_1;%(s) adjusted transmission time reverse propagation
Zc_1 = Zc;
Zc_2 = Zc;

k = max([numel(mG), numel(mE)]);

n=nan(1,k);%weighting function coefficient
n(1)=0.3/(1+3*beta);% eq(19)
for i=2:k
    n(i)=n(i-1)*3;% eq(19)
end

%calculate transfer function E. eq(14)
tmpsum=zeros(size(omega));
for i=1:numel(mE)
    tmpsum=tmpsum+mE(i)./(n(i)+1j.*omega*T);
end
E_1=tmpsum*Zc_1;
E_2=E_1;

%calculate transfer function F, eq (15)
tmpsum=0;
for i=1:numel(mE)
    tmpsum=tmpsum+mE(i)/(n(i));
end
b=1-8*beta/tmpsum; % eq (16)
F_1=Zc_2+b*E_1;
F_2=F_1;

%calculate transfer function G, eq (17)
tmpsum=zeros(size(omega));
for i=1:numel(mG)
    tmpsum=tmpsum+mG(i)*1j*omega*T./(n(i)+1j.*omega*T);
end
G_1=1-tmpsum;
G_2=G_1;

end
