function [ t12, t22 ] = t12t22venderBuhsExact( omega, L, r1, r2, nu, rho, K )
%[ t12, t22 ] = t12t22venderBuhsExact( omega, L, r1, r2, nu, rho, K ) 
%   Solves the boundary value problem for flow in tapered pipe and returns
%   t12 and t22 of the transmission matrix. This solver assumes a blocked
%   outlet with inlet flow excitation.
% Inputs: omega = freqeuncy (rad/s)
%             L = pipeline length (m)
%            r1 = inlet radius (m)
%            r2 = outlet radius (m)
%            nu = kinematic viscosity (m^2/s)
%           rho = density (kg/m^3)
%             K = bulk modulus (Pa)
%
% Outputs = Terms of the tranmission matrix [t12, t22 ] in eq (7).
%
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

%% Boundary conditions
Q0=1;% Flow frequency component at X=0 (Unit Input)
P1=0;% Pressure frequency component at X=1 (Constant 0 Pressure Outlet)

%% Initial guess
N_grid=100;
X0=0;
X1=L;
X=linspace(X0,X1,N_grid);
Q1=0.5;
y_bv0=[(Q1-Q0)/(X1-X0)*[1 1];
    Q0 Q1];%initial guess for y at boundaries

%% set up loop
P0_sol=nan(1,numel(omega));
P1_sol=nan(1,numel(omega));
Q0_sol=nan(1,numel(omega));
Q1_sol=nan(1,numel(omega));
dQ_dX0_sol= nan(1,numel(omega));
dQ_dX1_sol = nan(1,numel(omega));

for i=1:numel(omega)
    s=1j.*omega(i);
    
    %% set up BVP
    odefcn=@(x,y) Q_odefun(x,y,K,nu,rho,r1,r2,L,s);
    bcfun=@(ya,yb) PQ_bcfun(ya,yb,Q0,P1,r1,r2,L,K,s);
    bvpinitfun=@(X) [(X-X0)/(X1-X0)*(y_bv0(1,2)-y_bv0(1,1))+y_bv0(1,1);
        (X-X0)/(X1-X0)*(y_bv0(2,2)-y_bv0(2,1))+y_bv0(2,1)];
    if i==1
    solinit = bvpinit(X,bvpinitfun);
    end
    options = bvpset('RelTol',1e-3);
    %% solve
    sol = bvp4c(odefcn,bcfun,solinit,options);
    
    y=deval(sol,X);
    dQ_dX=y(1,:);
    Q=y(2,:);
    solinit.x = X;
    solinit.y = [dQ_dX;
                 Q];
    y_bv0=[y(:,1) y(:,end)];%update initial guess
   
    
    A = pi*(r1+(r2-r1)*X/L).^2;

    P = -dQ_dX.*K./s./A;

    dQ_dX0_sol(i)=dQ_dX(1);
    dQ_dX1_sol(i)=dQ_dX(end);
    Q0_sol(i)=Q(1);
    Q1_sol(i)=Q(end);
    P0_sol(i)=P(1);
    P1_sol(i)=P(end);
end

t12 = P0_sol./Q1_sol;
t22 = -Q0_sol./Q1_sol;
