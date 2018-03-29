function res=PQ_bcfun(ya,yb,Q0,P1,r1,r2,L,K,s)
%res=PQ_bcfun(ya,yb,Q0,P1,r1,r2,L,K,s)
% Pressure/Flow boundary conditions for ODE solver
%
% Inputs: ya,yb = States at inlet and outlet
%         Q0,P1 = boundary conditions
%         r1 = inlet radius (m)
%         r2 = outlet radius (m)
%          L = pipeline length (m)
%          K = bulk modulus (Pa)
%          s = Laplace frequency, j*omega
%
% Output: res = residual
%
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

%Outlet Pressure
A1 = pi*r2^2;
dAdx1 = 2*pi*r2*(r2-r1)/L;
%P_1 = -dAdx1.*K.*yb(2)./A1.^2./s-yb(1).*K./s./A1;
P_1 = -yb(1).*K./s./A1; 
    
res=[ya(2)-Q0;
    P_1-P1];