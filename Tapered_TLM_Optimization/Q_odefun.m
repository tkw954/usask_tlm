function dydx=Q_odefun(x,y,K,nu,rho,r1,r2,L,s)
%dydx=Q_odefun(x,y,K,nu,rho,r1,r2,L,s)
% ODE for flow in tapered pipe. 
%
% Inputs:    x = X;
%            y = [dQdx;Q]
%            K = bulk modulus (Pa)
%           nu = kinematic viscosity (m^2/s)
%          rho = density (kg/m^3)
%           r1 = inlet radius (m)
%           r2 = outlet radius (m)
%            L = pipeline length (m)
%            s = Laplace frequency, j*omega
%            
% Output: dydx = [d2Qdx2;dQdx]
%
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

dQdx=y(1);
Q=y(2);

r = r1+(r2-r1).*x/L;
N = -besselj(0,1j.*r.*sqrt(s./nu))./besselj(2,1j.*r.*sqrt(s./nu));

A = pi*(r1+(r2-r1)*x/L).^2;
dAdx = 2*pi*(r1+(r2-r1)*x/L)*(r2-r1)/L;
d2Qdx2 = dAdx./A.*dQdx+N.*rho.*s.^2.*Q./K;

dydx=[d2Qdx2;dQdx];