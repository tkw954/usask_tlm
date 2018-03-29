function res=Q_bcfun(ya,yb,Qa,Qb)
%res=Q_bcfun(ya,yb,Qa,Qb)
% Flow boundary conditions for ODE solver
%
% Inputs: ya,yb = States at inlet and outlet
%         Qa,Qb = boundary conditions
%
% Output: res = residual
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

% Calculate residual, NOTE: y(2) is Q and y(1) is dQ/dx.
res=[ya(2)-Qa;
    yb(2)-Qb];