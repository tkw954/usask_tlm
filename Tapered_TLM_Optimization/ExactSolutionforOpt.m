function [ t11, t21, t12, t22 ] = ExactSolutionforOpt( omega, L, r1, r2, nu, rho, K )
%[ t11, t21, t12, t22 ] = ExactSolutionforOpt( omega, L, r1, r2, nu, rho, K )
% Calculates the exact transmission matrix for a tapered pipeline using a
% numerical boundary value solver. NOTE: The Parallel Computing Toolbox
% required.
%
% Inputs: omega = freqeuncy (rad/s)
%             L = pipeline length (m)
%            r1 = inlet radius (m)
%            r2 = outlet radius (m)
%            nu = kinematic viscosity (m^2/s)
%           rho = density (kg/m^3)
%             K = bulk modulus (Pa)
%
% Outputs = Four terms of the tranmission matrix [t11, t21, t12, t22] in eq
% (7).
%
% Reference:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017



A = nan(numel(omega),2);
B = nan(numel(omega),2);

% Loop can be computed in parallel for significant time savings.
if license('test', 'distrib_computing_toolbox')
    parfor i=1:2
        if i==1
        % Calculate first 2 terms of the transmission matrix
        [A(:,i), B(:,i) ] = t11t21venderBuhsExact( omega, L, r1, r2, nu, rho, K );

        else
        % Calculate last 2 terms of the transmission matrix
        [ A(:,i), B(:,i) ] = t12t22venderBuhsExact( omega, L, r1, r2, nu, rho, K );
        end
    end
    
else
     for i=1:2
        if i==1
        % Calculate first 2 terms of the transmission matrix
        [A(:,i), B(:,i) ] = t11t21venderBuhsExact( omega, L, r1, r2, nu, rho, K );

        else
        % Calculate last 2 terms of the transmission matrix
        [ A(:,i), B(:,i) ] = t12t22venderBuhsExact( omega, L, r1, r2, nu, rho, K );
        end
     end
end

% Outputs
t11(1,1:numel(omega)) = A(:,1);
t21(1,1:numel(omega)) = B(:,1);
t12(1,1:numel(omega)) = A(:,2);
t22(1,1:numel(omega)) = B(:,2);

end

