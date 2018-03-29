function [ t11, t21 ] = t11t21venderBuhsExact( omega, L, r1, r2, nu, rho, K )
%[ t11, t21 ] = t11t21venderBuhsExact( omega, L, r1, r2, nu, rho, K ) 
%   Solves the boundary value problem for flow in tapered pipe and returns
%   t11 and t21 of the transmission matrix. This solver assumes a blocked
%   outlet with inlet flow excitation.
% Inputs: omega = freqeuncy (rad/s)
%             L = pipeline length (m)
%            r1 = inlet radius (m)
%            r2 = outlet radius (m)
%            nu = kinematic viscosity (m^2/s)
%           rho = density (kg/m^3)
%             K = bulk modulus (Pa)
%
% Outputs = Terms of the tranmission matrix [t11, t21 ] in eq (7).
%
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

%% Boundary conditions
Q0=1;% Flow frequency component at X=0 (Unit Input)
Q1=0;% Flow frequency component at X=1 (Blocked Outlet)

%% Initial guess
N_grid=100; % 100 points over the pipe length
X0=0;
X1=L;
X=linspace(X0,X1,N_grid);

y_bv0=[(Q1-Q0)/(X1-X0)*[1 1];
    Q0 Q1];%initial guess for y at boundaries

%% Set up loop
P0_sol=nan(1,numel(omega));
P1_sol=nan(1,numel(omega));
Q0_sol=nan(1,numel(omega));
Q1_sol=nan(1,numel(omega));

for i=1:numel(omega)
    s=1j.*omega(i);
    
    %% set up BVP
    odefcn=@(x,y) Q_odefun(x,y,K,nu,rho,r1,r2,L,s);
    bcfun=@(ya,yb) Q_bcfun(ya,yb,Q0,Q1);
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
   
    Q0_sol(i)=Q(1);
    Q1_sol(i)=Q(end);
    P0_sol(i)=P(1);
    P1_sol(i)=P(end);
end

t11 = P0_sol./P1_sol;

t21 = Q0_sol./P1_sol;

