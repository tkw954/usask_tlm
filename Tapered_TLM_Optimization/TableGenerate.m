% This code will generate a 3D lookup table of wighting factors for the
% Tapered Transmission Line method as developed in ven der Buhs et. al
% (2017).

% References:

% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017


% Number of weighing factors used
k = 6;

% Initilize counters
opt_failed = 0;
h = 1;

% Range of dissipation number (beta) (dimensionless)
beta = logspace(-4,0,33);

% Range of taper ratio (lambda) (dimensionless)
lambda = linspace(1,0.5,10);

% Preallocate parameters
mE = NaN(numel(lambda),k,numel(beta));
mG = NaN(numel(lambda),k,numel(beta));
tau = NaN(numel(lambda),numel(beta));

% Loop through the range of dissipation number
for i=1:numel(beta)
    % Loop through the range of taper ratio. Always starting at no taper
    % (i.e. lambda = 1) and moving towards larger taper.
    for j=1:numel(lambda)
          if j==1
           % Initial guess of paramters found from reference:
           % N Johnston. Simulink models, http://people.bath.ac.uk/ensdnj/models/newtlm.html, 2014.
           [ni, k2, mgi, mei, deltai] = gettlmcoeffs(beta(i), 6);
           
           % Assign initial parameters to vector "params"
            params = [mei mgi deltai];
          elseif opt_failed==1
           % If optimization failed on previos iteration, select initial
           % paramters again from Johnston et al. solution.
              [ni, k2, mgi, mei, deltai] = gettlmcoeffs(beta(i), 6);
            params = [mei mgi deltai];
          else
           % If optimization was successful on previous iteration, use
           % paramters from previous iteration as initial guess.
            params = [mE(j-1,:,i) mG(j-1,:,i) tau(j-1,i)];
            
          end

          % The numerical solver will not arrive at a solution for small
          % values of dissipation number (beta) and will return an error. 
          % As a result try the optimization, if an error occurs assign NaN
          % to the paramters and try again for the next value.
          try
          [mE(j,:,i), mG(j,:,i), tau(j,i)] = OptimizationsForTable(lambda(j),beta(i),params);
          
          % opt_failed = 0 if the optimization/solution was successful.
          opt_failed = 0;
          
          catch
          mE(j,:,i) = NaN(1,k);
          mG(j,:,i) = NaN(1,k);
          tau(j,i) = NaN;
          warning(['Optimization for beta = ' num2str(beta(i)) ' and lambda = ' num2str(lambda(j)) ' failed since the exact numerical solution could not be found.'])
          % opt_failed = 1 if the optimization/solution failed.
          opt_failed = 1;
          end
          
          % Displays progress to the user (i.e. how many optimizations 
          % completed vs. how many in total)
          disp([num2str(h) ' of ' num2str(numel(beta)*numel(lambda))]);
          
          % Number of optimizations counter
          h = h+1;
    end
end

mE_lookup = permute(mE,[1 3 2]);
mG_lookup = permute(mG,[1 3 2]);
tau_lookup = tau;
